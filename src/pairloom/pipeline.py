"""End-to-end pipeline: streaming, adjustment, merging, statistics, publish.

Processing streams templates in bounded chunks, applies repeat-aware joint
adjustment, overlap analysis, classification, and (in merge mode) consensus
merging, then publishes all outputs atomically. Worker failures propagate
and abort the whole run; no partial output is ever presented as complete.

Determinism: templates are processed in query-name stream order, parallel
chunks are consumed in submission order, and all statistics are integer
aggregations, so results do not depend on scheduling or thread count.
"""

from __future__ import annotations

import contextlib
import os
import sys
import tempfile
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam

from pairloom import __version__
from pairloom.alignment import Placement, build_read_structure
from pairloom.candidates import discover_candidates, normalize_record_cigar, single_read_leftmost
from pairloom.cigarops import compute_md_nm
from pairloom.classify import (
    BOTH_ALIGNED,
    LEFT_ALIGNED_ONLY,
    MERGEABLE,
    NEITHER_ALIGNED,
    RIGHT_ALIGNED_ONLY,
    StatsCounter,
    assign_status_tag,
)
from pairloom.config import STALE_TAGS_REMOVED, PairLoomConfig
from pairloom.exceptions import (
    DataInvariantError,
    OutputError,
    PairLoomError,
    WorkerError,
)
from pairloom.mapq import adjusted_mapq
from pairloom.merge import MergeSuccess, try_merge
from pairloom.overlap import analyze_overlap
from pairloom.pairing import select_pair
from pairloom.parquetout import DisagreementAggregator, write_read_stats_parquet
from pairloom.records import (
    check_pair_flags,
    detect_sort_order,
    iter_templates,
    name_sort_to_temp,
    open_input_bam,
    split_primaries,
)
from pairloom.validate import validate_paths

# Best-effort close during error rollback is intentional (never suppresses
# the in-flight original exception, which is re-raised by the caller).

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAPPED = 0x4
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80


@dataclass(slots=True)
class TemplateResult:
    """Everything one template contributes to the outputs."""

    classification: str
    merge_reason: str | None
    #: (stream name, SAM line) pairs in deterministic template order
    out_lines: list[tuple[str, str]] = field(default_factory=list)
    #: (tid, rpos0, ref_base, observed_base, strand_reverse, pair_disagrees)
    observations: list[tuple[int, int, str, str, bool, bool]] = field(default_factory=list)


@dataclass(slots=True)
class RunSummary:
    classification_counts: dict[str, int]
    merge_reasons: dict[str, int]
    reported_disagreement_positions: int
    outputs: dict[str, Path]
    sorted: bool
    indexed: bool


class _WorkerState:
    """Per-process configuration, header, and reference cache."""

    def __init__(self, cfg: PairLoomConfig, header: Any) -> None:
        from pairloom.reference import ReferenceCache

        self.cfg = cfg
        self.header = header
        self.refcache = ReferenceCache(
            cfg.reference_path,
            max_contigs=cfg.reference_cache_contigs,
            max_bases=cfg.reference_cache_bases,
        )
        self.ref_names: list[str] = [sq["SN"] for sq in (header.to_dict().get("SQ") or [])]

    def contig(self, tid: int) -> str:
        try:
            return self.ref_names[tid]
        except IndexError:
            raise DataInvariantError(f"reference id {tid} missing from BAM header") from None


_STATE: _WorkerState | None = None


def _worker_init(cfg: PairLoomConfig, header_dict: dict) -> None:
    global _STATE
    header = pysam.AlignmentHeader.from_dict(header_dict)
    _STATE = _WorkerState(cfg, header)


def _ensure_state() -> _WorkerState:
    if _STATE is None:
        raise WorkerError("worker state not initialized")
    return _STATE


def _process_chunk_strings(chunk: list[list[str]]) -> list[TemplateResult]:
    """Pool entry point: parse SAM lines and process each template."""
    state = _ensure_state()
    results: list[TemplateResult] = []
    for lines in chunk:
        records = [pysam.AlignedSegment.fromstring(line, state.header) for line in lines]
        results.append(process_template(records, state))
    return results


# ---------------------------------------------------------------------------
# Per-template processing
# ---------------------------------------------------------------------------


def process_template(records: list[Any], state: _WorkerState) -> TemplateResult:
    """Process one template into output lines, classification, and observations."""
    if not records:
        raise DataInvariantError("empty template")
    name = records[0].query_name
    primaries, others = split_primaries(records)
    check_pair_flags(primaries, name)
    # Validate the primary group as a whole BEFORE assigning slots: record
    # order must not decide whether a malformed template is rejected.
    paired_primaries = [r for r in primaries if r.flag & FLAG_PAIRED]
    unpaired_primaries = [r for r in primaries if not (r.flag & FLAG_PAIRED)]
    if paired_primaries and unpaired_primaries:
        raise DataInvariantError(f"template {name!r} mixes paired and unpaired primary records")
    if len(unpaired_primaries) > 1:
        raise DataInvariantError(
            f"template {name!r} has {len(unpaired_primaries)} unpaired primary "
            "records; at most one is allowed"
        )
    r1 = r2 = None
    for record in primaries:
        if record.flag & FLAG_PAIRED:
            if record.flag & FLAG_READ1:
                r1 = record
            else:
                r2 = record
        else:
            r1 = record  # single non-paired primary: treated as R1-like orphan

    classification: str
    merge_reason: str | None = None
    merged_line: str | None = None
    observations: list[tuple[int, int, str, str, bool, bool]] = []

    if r1 is not None and r2 is not None:
        m1 = not r1.flag & FLAG_UNMAPPED
        m2 = not r2.flag & FLAG_UNMAPPED
        if m1 and m2:
            merged_line, classification, merge_reason, observations = _both_mapped(
                r1, r2, others, state
            )
        elif m1:
            classification = LEFT_ALIGNED_ONLY
            _adjust_single(r1, r2, state)
        elif m2:
            classification = RIGHT_ALIGNED_ONLY
            _adjust_single(r2, r1, state)
        else:
            classification = NEITHER_ALIGNED
    elif r1 is None and r2 is None:
        # No primary records at all (e.g. a filtered BAM holding only
        # secondary/supplementary alignments): documented pass-through. No
        # primary alignment evidence exists, so the template counts as
        # neither_aligned; the records below stay untagged and unchanged.
        classification = NEITHER_ALIGNED
    else:
        record = r1 if r1 is not None else r2
        if record.flag & FLAG_UNMAPPED:
            classification = NEITHER_ALIGNED
        else:
            classification = LEFT_ALIGNED_ONLY if record is r1 else RIGHT_ALIGNED_ONLY
            _adjust_single(record, None, state)

    for record in primaries:
        assign_status_tag(record, state.cfg, classification)

    out_lines: list[tuple[str, str]] = []
    if state.cfg.mode == "adjust":
        for record in records:
            out_lines.append(("adjusted", record.to_string()))
    elif merged_line is not None:
        # The consensus record replaces the primary mates; only non-primary
        # records (which by policy prevent merging) could continue, so this
        # template never duplicates into the unmerged output.
        out_lines.append(("merged", merged_line))
        for record in others:
            out_lines.append(("unmerged", record.to_string()))
    else:
        for record in records:
            out_lines.append(("unmerged", record.to_string()))
    return TemplateResult(classification, merge_reason, out_lines, observations)


def _both_mapped(
    r1: Any, r2: Any, others: list[Any], state: _WorkerState
) -> tuple[str | None, str, str | None, list[tuple[int, int, str, str, bool, bool]]]:
    """Jointly adjust a both-mapped pair; returns (merged SAM line, class, reason, obs)."""
    cfg = state.cfg
    tid1, tid2 = r1.reference_id, r2.reference_id
    if tid1 < 0 or tid2 < 0:
        raise DataInvariantError(
            f"template {r1.query_name!r} has mapped records without reference id"
        )
    cands1 = _discover(r1, state)
    cands2 = _discover(r2, state)
    # Eligibility filters (min_mapq/min_base_quality) apply to the ORIGINAL
    # aligner MAPQ, not to the conservative values PairLoom writes after
    # adjustment: capping our own output MAPQ must not disqualify a pair
    # whose input evidence was strong.
    orig_mapq1 = r1.mapping_quality
    orig_mapq2 = r2.mapping_quality
    choice = None
    if tid1 == tid2:
        choice = select_pair(cands1, cands2, r1.is_reverse, r2.is_reverse, cfg)

    if choice is not None:
        _apply_placement(r1, choice.p1, choice.n_equiv1, choice.truncated, state)
        _apply_placement(r2, choice.p2, choice.n_equiv2, choice.truncated, state)
        proper = True
    else:
        # No valid candidate-pair combination: keep the ORIGINAL BAM
        # placements (never slide mates independently).
        _apply_placement(r1, cands1.current, len(cands1.placements), cands1.truncated, state)
        _apply_placement(r2, cands2.current, len(cands2.placements), cands2.truncated, state)
        proper = False
    _sync_mate_fields(r1, r2)
    _set_proper_flag(r1, proper)
    _set_proper_flag(r2, proper)

    observations: list[tuple[int, int, str, str, bool, bool]] = []
    if r1.reference_id != r2.reference_id:
        # Cross-contig pairs keep their adjusted records; no overlap analysis.
        return None, BOTH_ALIGNED, "different_contig", observations

    struct1 = build_read_structure(r1)
    struct2 = build_read_structure(r2)
    # The records now carry the conservative ADJUSTED MAPQs; the consensus
    # record takes its output MAPQ from those, never from the originals.
    merged_mapq = min(struct1.mapq, struct2.mapq)
    # Eligibility filters (min_mapq) continue to apply to the ORIGINAL
    # aligner MAPQ (policy comment above).
    struct1.mapq = orig_mapq1
    struct2.mapq = orig_mapq2
    contig = state.contig(r1.reference_id)
    # Observations exist only where BOTH mates align a base (the aligned-
    # span intersection), so only that reference slice is fetched: distant
    # mates cost a read-sized fetch instead of the whole span between them.
    ov_lo = max(struct1.start0, struct2.start0)
    ov_hi = min(struct1.end0, struct2.end0)
    observations: list[tuple[int, int, str, str, bool, bool]] = []
    if ov_lo < ov_hi:
        ref_seq = state.refcache.fetch(contig, ov_lo, ov_hi)
        overlap = analyze_overlap(struct1, struct2, ref_seq, ov_lo, cfg)
        observations = [
            (o.tid, o.rpos0, o.ref_base, o.base, o.strand_reverse, o.pair_disagrees)
            for o in overlap.observations
        ]

    outcome = try_merge(
        struct1, struct2, cfg, has_non_primary=bool(others), merged_mapq=merged_mapq
    )
    if isinstance(outcome, MergeSuccess):
        merged_record = _build_merged_record(r1.query_name, outcome, contig, state)
        return merged_record.to_string(), MERGEABLE, None, observations
    return None, BOTH_ALIGNED, outcome.reason, observations


def _discover(record: Any, state: _WorkerState) -> Any:
    """Discover candidate placements.

    BAM SEQ/CIGAR are stored in reference orientation for both strands, so
    base-level verification and exact-block sliding need no strand-specific
    transform.
    """
    cfg = state.cfg
    cigar = normalize_record_cigar(record.cigartuples)
    seq = record.query_sequence
    if seq is None:
        raise DataInvariantError(f"record {record.query_name!r} has no SEQ")
    contig = state.contig(record.reference_id)
    return discover_candidates(
        tid=record.reference_id,
        pos0=record.reference_start,
        cigar=cigar,
        query_seq=seq.upper(),
        contig=contig,
        refcache=state.refcache,
        cfg=cfg,
        read_name=record.query_name,
    )


def _apply_placement(
    record: Any,
    placement: Placement,
    n_equiv: int,
    truncated: bool,
    state: _WorkerState,
) -> None:
    """Move a record to ``placement`` and update MAPQ/NM/MD consistently.

    Only exact-match placements are ever moved (discovery guarantees this),
    so the rewritten tags are orientation-independent: NM=0 and MD=<span>.
    Stale aligner scoring tags are removed rather than guessed at.
    """
    cfg = state.cfg
    changed = placement.pos0 != record.reference_start
    if changed:
        if not placement.exact:
            raise DataInvariantError(
                f"internal invariant violated: non-exact placement change for "
                f"read {record.query_name!r}"
            )
        record.reference_start = placement.pos0
        record.set_tag("MD", str(placement.ref_span), value_type="Z")
        record.set_tag("NM", 0, value_type="i")
        for tag in STALE_TAGS_REMOVED:
            if record.has_tag(tag):
                # pysam (0.24) has no remove_tag(); set_tag(tag, None) deletes.
                record.set_tag(tag, None)
    record.mapping_quality = adjusted_mapq(
        orig_mapq=record.mapping_quality,
        n_equiv=n_equiv,
        truncated=truncated,
        exact=placement.exact,
        placement_changed=changed,
        force_mapq60_exact=cfg.force_mapq60_exact,
    )


def _adjust_single(record: Any, mate: Any | None, state: _WorkerState) -> None:
    """Leftmost-normalize a mapped orphan (or mapped mate of an unmapped mate).

    ``mate``, when given, is the unmapped record of the same template.
    """
    if record.flag & FLAG_UNMAPPED or record.reference_id < 0:
        return
    cands = _discover(record, state)
    placement = single_read_leftmost(cands)
    _apply_placement(record, placement, len(cands.placements), cands.truncated, state)
    if mate is None:
        # No mate record in the template: a reciprocal MQ cannot exist.
        record.set_tag("MQ", None)
        return
    # SAM convention: an unmapped mate carries the mapped read's locus.
    if mate.reference_id == -1 or mate.reference_id == record.reference_id:
        mate.reference_id = record.reference_id
        mate.reference_start = record.reference_start
    record.next_reference_id = mate.reference_id
    record.next_reference_start = mate.reference_start
    # The unmapped mate's RNEXT/PNEXT must follow the mapped mate's NEW
    # locus; keeping the input value would point at the pre-adjustment
    # position (SAM spec: an unmapped read's mate fields reference the
    # mapped mate).
    mate.next_reference_id = record.reference_id
    mate.next_reference_start = record.reference_start
    # MQ semantics: the mapped record's mate is unmapped (tag removed); the
    # unmapped record's MQ reflects the mapped mate's ADJUSTED MAPQ.
    record.set_tag("MQ", None)
    mate.set_tag("MQ", record.mapping_quality, value_type="i")
    _set_proper_flag(record, False)
    _set_proper_flag(mate, False)


def _sync_mate_fields(r1: Any, r2: Any) -> None:
    """Refresh mate coordinates, reciprocal MQ, and TLEN after adjustment."""
    r1.next_reference_id = r2.reference_id
    r1.next_reference_start = r2.reference_start
    r2.next_reference_id = r1.reference_id
    r2.next_reference_start = r1.reference_start
    # Reciprocal MQ (mate MAPQ) follows the ADJUSTED values — also when only
    # the MAPQ changed (repeat-ambiguous placements that did not move).
    r1.set_tag("MQ", r2.mapping_quality, value_type="i")
    r2.set_tag("MQ", r1.mapping_quality, value_type="i")
    if r1.reference_id == r2.reference_id and r1.reference_id >= 0:
        start = min(r1.reference_start, r2.reference_start)
        end = max(r1.reference_start + _ref_len(r1), r2.reference_start + _ref_len(r2))
        tlen = end - start
        if r1.reference_start < r2.reference_start:
            r1.template_length, r2.template_length = tlen, -tlen
        elif r1.reference_start > r2.reference_start:
            r1.template_length, r2.template_length = -tlen, tlen
        else:
            # SAM TLEN: signs must be opposite within a template; with
            # identical coordinates the first segment (read1) is positive.
            r1.template_length, r2.template_length = tlen, -tlen
    else:
        r1.template_length = 0
        r2.template_length = 0


def _ref_len(record: Any) -> int:
    from pairloom.cigarops import cigar_ref_len

    cigar = record.cigartuples
    if not cigar:
        return 0
    return cigar_ref_len(tuple(cigar))


def _set_proper_flag(record: Any, proper: bool) -> None:
    if proper:
        record.flag = record.flag | FLAG_PROPER
    else:
        record.flag = record.flag & ~FLAG_PROPER


def _build_merged_record(
    qname: str,
    merged: Any,
    contig: str,
    state: _WorkerState,
) -> Any:
    """Create the single consensus BAM record for a merged pair."""
    record = pysam.AlignedSegment(state.header)
    record.query_name = qname
    record.flag = 0
    record.reference_id = state.ref_names.index(contig)
    record.reference_start = merged.start0
    record.cigartuples = list(merged.cigar)
    record.query_sequence = merged.seq
    record.query_qualities = list(merged.quals)
    record.mapping_quality = merged.mapq
    record.next_reference_id = -1
    record.next_reference_start = -1
    record.template_length = 0
    record.set_tag(state.cfg.status_tag, MERGEABLE, value_type="Z")
    if merged.read_group is not None:
        record.set_tag("RG", merged.read_group, value_type="Z")
    ref_end = merged.start0 + _merged_ref_len(merged.cigar)
    ref_seq = state.refcache.fetch(contig, merged.start0, ref_end)
    md, nm = compute_md_nm(
        merged.start0,
        merged.cigar,
        merged.seq,
        ref_seq,
        f"merged read {qname} at {contig}:{merged.start0 + 1}",
    )
    record.set_tag("MD", md, value_type="Z")
    record.set_tag("NM", nm, value_type="i")
    return record


def _merged_ref_len(cigar: tuple[tuple[int, int], ...]) -> int:
    from pairloom.cigarops import cigar_ref_len

    return cigar_ref_len(cigar)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


class _BamWriter:
    """BAM writer to a temp file on the destination filesystem."""

    def __init__(self, final_path: Path, header: Any) -> None:
        self.final_path = final_path
        self.tmp_path = Path(str(final_path) + f".pairloom-tmp-{os.getpid()}")
        try:
            self._writer = pysam.AlignmentFile(str(self.tmp_path), "wb", header=header)
        except (OSError, ValueError) as exc:
            raise OutputError(f"cannot create output BAM {self.tmp_path}: {exc}") from exc

    def write(self, record: Any) -> None:
        try:
            self._writer.write(record)
        except (OSError, ValueError) as exc:
            raise OutputError(f"failed writing BAM record to {self.tmp_path}: {exc}") from exc

    def close(self) -> None:
        self._writer.close()

    def cleanup(self) -> None:
        with contextlib.suppress(OSError, ValueError):
            self._writer.close()
        self.tmp_path.unlink(missing_ok=True)


def run_pipeline(cfg: PairLoomConfig) -> RunSummary:
    """Run the full pipeline and publish outputs atomically."""
    validate_paths(cfg)
    stats = StatsCounter()
    reasons: Counter[str] = Counter()
    published: list[Path] = []
    backups: list[tuple[Path, Path]] = []
    writers: dict[str, _BamWriter] = {}
    tmp_files: list[Path] = []
    committed = False

    try:
        with tempfile.TemporaryDirectory(prefix="pairloom-") as tmpdir:
            input_file = open_input_bam(str(cfg.input_path))
            try:
                header_dict = input_file.header.to_dict()
                sort_order = detect_sort_order(input_file)
                if sort_order != "queryname":
                    input_file.close()
                    sorted_input = name_sort_to_temp(str(cfg.input_path), tmpdir)
                    input_file = open_input_bam(sorted_input)
                out_header_dict = _output_header(header_dict, cfg, sort_output=cfg.sort_output)
                header_obj = pysam.AlignmentHeader.from_dict(out_header_dict)
                ref_names = [sq["SN"] for sq in (header_dict.get("SQ") or [])]

                main_state = _WorkerState(cfg, header_obj)
                try:
                    for contig in ref_names:
                        main_state.refcache.length(contig)
                finally:
                    main_state.refcache.close()

                outputs = cfg.mode_outputs()
                # Writer keys are stream names ("adjusted"/"merged"/"unmerged"),
                # matching TemplateResult.out_lines roles.
                for stream, role in (
                    ("adjusted", "adjusted_bam"),
                    ("merged", "merged_bam"),
                    ("unmerged", "unmerged_bam"),
                ):
                    if role in outputs:
                        writers[stream] = _BamWriter(outputs[role], header_obj)

                aggregator = DisagreementAggregator(cfg)
                aggregator.set_spill_dir(tmpdir)

                def handle(result: TemplateResult) -> None:
                    stats.add(result.classification)
                    if result.merge_reason:
                        reasons[result.merge_reason] += 1
                    for stream, line in result.out_lines:
                        writers[stream].write(pysam.AlignedSegment.fromstring(line, header_obj))
                    for obs in result.observations:
                        aggregator.add_observation(*obs)

                chunks = _iter_chunks(input_file, cfg.chunk_size)
                if cfg.threads == 1:
                    state = _WorkerState(cfg, header_obj)
                    try:
                        for chunk in chunks:
                            for records in chunk:
                                handle(process_template(records, state))
                    finally:
                        state.refcache.close()
                else:
                    _run_parallel(chunks, cfg, header_dict, handle)

                for writer in writers.values():
                    writer.close()
                    tmp_files.append(writer.tmp_path)

                disag_tmp = Path(str(cfg.disagreement_parquet) + f".pairloom-tmp-{os.getpid()}")
                stats_tmp = Path(str(cfg.read_stats_parquet) + f".pairloom-tmp-{os.getpid()}")
                # Register temps before writing so failures clean them up.
                tmp_files.extend([disag_tmp, stats_tmp])
                aggregator.finalize(disag_tmp, ref_names)
                write_read_stats_parquet(stats_tmp, stats.counts(), cfg.parquet_compression)

                # Optional coordinate sort (in place on the temp files).
                if cfg.sort_output:
                    for writer in writers.values():
                        sorted_tmp = Path(str(writer.tmp_path) + ".sorted")
                        try:
                            pysam.sort(
                                "-@",
                                str(min(cfg.threads, 4)),
                                "-T",
                                os.path.join(tmpdir, "sortout"),
                                "-o",
                                str(sorted_tmp),
                                str(writer.tmp_path),
                            )
                        except (OSError, ValueError) as exc:
                            raise OutputError(
                                f"coordinate sorting failed for {writer.tmp_path}: {exc}"
                            ) from exc
                        writer.tmp_path.unlink()
                        writer.tmp_path = sorted_tmp
                        tmp_files.append(sorted_tmp)

                index_files: list[tuple[Path, Path]] = []
                if cfg.index_output:
                    for writer in writers.values():
                        idx_tmp = Path(str(writer.tmp_path) + ".bai")
                        try:
                            pysam.index("-o", str(idx_tmp), str(writer.tmp_path))
                        except (OSError, ValueError) as exc:
                            raise OutputError(
                                f"indexing failed for {writer.tmp_path}: {exc}"
                            ) from exc
                        final_idx = Path(str(writer.final_path) + ".bai")
                        index_files.append((idx_tmp, final_idx))
                        tmp_files.append(idx_tmp)

                # Atomic publish: BAMs, then indexes, then Parquet tables.
                # Pre-existing outputs (possible only under --force) are
                # moved aside first and restored on failure, so a failed
                # run never destroys the previous results.
                pairs: list[tuple[Path, Path]] = []
                for writer in writers.values():
                    pairs.append((writer.tmp_path, writer.final_path))
                for idx_tmp, final_idx in index_files:
                    pairs.append((idx_tmp, final_idx))
                pairs.append((disag_tmp, cfg.disagreement_parquet))
                pairs.append((stats_tmp, cfg.read_stats_parquet))
                for tmp, final in pairs:
                    if final.exists():
                        backup = Path(str(final) + f".pairloom-bak-{os.getpid()}")
                        os.replace(final, backup)
                        backups.append((backup, final))
                    os.replace(tmp, final)
                    published.append(final)
                # Commit point: every requested output is published and
                # complete. From here rollback is FORBIDDEN — destroying
                # committed outputs (or their backups) would lose data
                # instead of protecting it, whatever fails later.
                committed = True
                # Backup removal is post-commit housekeeping: a failure
                # keeps the outputs and reports the leftovers loudly.
                leftover: list[Path] = []
                for backup, _final in backups:
                    try:
                        backup.unlink()
                    except OSError:
                        leftover.append(backup)
                if leftover:
                    print(
                        "pairmerge: warning: could not remove backup file(s); "
                        "delete them manually: " + ", ".join(str(p) for p in leftover),
                        file=sys.stderr,
                    )
            finally:
                input_file.close()
    except PairLoomError:
        if not committed:
            _rollback(writers, tmp_files, published, backups)
        else:
            _warn_committed()
        raise
    except OSError as exc:
        if not committed:
            _rollback(writers, tmp_files, published, backups)
            raise OutputError(f"I/O failure: {exc}") from exc
        _warn_committed()
        raise OutputError(
            f"post-commit failure after publication ({exc}); published outputs were kept"
        ) from exc
    except Exception:
        if not committed:
            _rollback(writers, tmp_files, published, backups)
        else:
            _warn_committed()
        raise

    print(
        "pairmerge: done. "
        + ", ".join(f"{k}={v}" for k, v in stats.counts().items())
        + (f"; merge-reasons={dict(reasons)}" if reasons else "")
        + f"; disagreement-positions={aggregator.reported_positions}",
        file=sys.stderr,
    )
    return RunSummary(
        classification_counts=stats.counts(),
        merge_reasons=dict(reasons),
        reported_disagreement_positions=aggregator.reported_positions,
        outputs=cfg.mode_outputs(),
        sorted=cfg.sort_output,
        indexed=cfg.index_output,
    )


def _iter_chunks(input_file: Any, chunk_size: int) -> Any:
    """Yield lists of templates (each a list of records) of bounded size."""
    chunk: list[list[Any]] = []
    for template in iter_templates(input_file):
        chunk.append(template)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _run_parallel(chunks: Any, cfg: PairLoomConfig, header_dict: dict, handle: Any) -> None:
    """Process chunks in a process pool with bounded, ordered submission."""
    import multiprocessing

    ctx = multiprocessing.get_context("spawn")
    max_inflight = cfg.threads * 2
    with ProcessPoolExecutor(
        max_workers=cfg.threads,
        mp_context=ctx,
        initializer=_worker_init,
        initargs=(cfg, header_dict),
    ) as pool:
        pending: deque = deque()

        def submit_next(chunk: list[list[Any]]) -> None:
            payload = [[record.to_string() for record in template] for template in chunk]
            pending.append(pool.submit(_process_chunk_strings, payload))

        try:
            for chunk in chunks:
                submit_next(chunk)
                while len(pending) >= max_inflight:
                    _drain_one(pending, handle)
            while pending:
                _drain_one(pending, handle)
        except BaseException:
            for future in pending:
                future.cancel()
            raise


def _drain_one(pending: deque, handle: Any) -> None:
    future = pending.popleft()
    try:
        results = future.result()
    except PairLoomError:
        raise
    except Exception as exc:
        raise WorkerError(f"worker chunk failed: {exc}") from exc
    for result in results:
        handle(result)


def _warn_committed() -> None:
    """Post-commit failure notice: outputs are complete and stay in place."""
    print(
        "pairmerge: warning: failure after publication; all outputs were "
        "already published and are kept (no rollback).",
        file=sys.stderr,
    )


def _rollback(
    writers: dict[str, _BamWriter],
    tmp_files: list[Path],
    published: list[Path],
    backups: list[tuple[Path, Path]] = (),
) -> None:
    """Best-effort cleanup: remove temps and already-published outputs,
    then restore the backed-up previous outputs (present only when the run
    overwrote existing files under --force)."""
    for writer in writers.values():
        writer.cleanup()
    for path in tmp_files:
        path.unlink(missing_ok=True)
    for path in published:
        path.unlink(missing_ok=True)
    for backup, final in backups:
        os.replace(backup, final)


def _output_header(header_dict: dict, cfg: PairLoomConfig, sort_output: bool) -> dict:
    """Copy the input header, set SO, and append a @PG record."""
    out = dict(header_dict)
    hd = dict(out.get("HD") or {})
    hd["SO"] = "coordinate" if sort_output else "queryname"
    out["HD"] = hd
    pg_list = list(out.get("PG") or [])
    pg_id = "pairmerge"
    existing_ids = {pg.get("ID") for pg in pg_list}
    suffix = 0
    while pg_id in existing_ids:
        suffix += 1
        pg_id = f"pairmerge.{suffix}"
    pg_list.append(
        {
            "ID": pg_id,
            "PN": "pairmerge",
            "VN": __version__,
            "CL": (
                f"pairmerge mode={cfg.mode} input={cfg.input_path} reference={cfg.reference_path}"
            ),
        }
    )
    out["PG"] = pg_list
    return out
