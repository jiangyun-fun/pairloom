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
import gc
import os
import subprocess
import sys
import tempfile
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pairloom_core as _core
import pysam

from pairloom import __version__
from pairloom.alignment import (
    Placement,
    expand_reference_matches,
    revcomp,
)
from pairloom.candidates import discover_candidates, normalize_record_cigar, single_read_leftmost
from pairloom.cigarops import CIGAR_H, CIGAR_I, CIGAR_M, CIGAR_N, CIGAR_S, reverse_cigar
from pairloom.classify import (
    BOTH_ALIGNED,
    LEFT_ALIGNED_ONLY,
    MERGEABLE,
    NEITHER_ALIGNED,
    RIGHT_ALIGNED_ONLY,
    StatsCounter,
    assign_rescue_tag,
    assign_status_tag,
)
from pairloom.config import RESCUE_TAG, RESCUE_VALUE, STALE_TAGS_REMOVED, PairLoomConfig
from pairloom.exceptions import (
    DataInvariantError,
    OutputError,
    PairLoomError,
    WorkerError,
)
from pairloom.mapq import adjusted_mapq
from pairloom.pairing import select_pair
from pairloom.parquetout import (
    FastDisagreementAccumulator,
    MergeStatsCollector,
    write_merge_stats_parquet,
    write_read_stats_parquet,
)
from pairloom.records import (
    check_pair_flags,
    detect_sort_order,
    iter_templates,
    name_sort_to_temp,
    open_input_bam,
    split_primaries,
)
from pairloom.rescue import STATUS_RESCUED, RescuePlacement, rescue_eligible, rescue_pair
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
    #: batched overlap aggregates as Rust-core blobs:
    #: (keys, disg, ref, altkeys, altbases, altstrands), or None
    observations: tuple[bytes, bytes, bytes, bytes, bytes, bytes] | None = None
    #: aligned reference overlap of a merged pair (merged templates only)
    overlap_len: int | None = None
    #: True when a MERGED template carried secondary/supplementary
    #: evidence (consensus carries the evidence tag; those records dropped)
    non_primary_evidence: bool = False


@dataclass(slots=True)
class RunSummary:
    classification_counts: dict[str, int]
    merge_reasons: dict[str, int]
    merged_non_primary: int
    reported_disagreement_positions: int
    outputs: dict[str, Path]
    sorted: bool
    indexed: bool


class _WorkerState:
    """Per-process configuration, header, and reference cache."""

    def __init__(self, cfg: PairLoomConfig, header: Any | None) -> None:
        from pairloom.reference import ReferenceCache

        self.cfg = cfg
        self.header = header
        ref_names = [sq["SN"] for sq in (header.to_dict().get("SQ") or [])] if header else []
        self.refcache = ReferenceCache(
            cfg.reference_path,
            max_contigs=cfg.reference_cache_contigs,
            max_bases=cfg.reference_cache_bases,
        )
        self.ref_names: list[str] = ref_names
        #: Reverse mapping for rescue's global search: reference FASTA names
        #: that the BAM header can express as output placements.
        self.tid_by_name: dict[str, int] = {name: tid for tid, name in enumerate(ref_names)}
        #: Discovery memo: identical (placement, read sequence) reads repeat
        #: massively on oligo panels; discovery is a pure function of these
        #: inputs, so results are shared read-only. Bounded; cleared on cap.
        self.discover_cache: dict[tuple, Any] = {}
        #: Rescue scan memo: (block, reverse, span) -> (occurrences, truncated).
        self.occ_cache: dict[tuple, tuple] = {}
        #: Contig-key memo for rescue's global search.
        self.contig_keys_cache: dict[str, int] | None = None
        #: Whole-contig string cache for rescue scans, built lazily once per
        #: worker when the reference fits the cache base budget; None means
        #: not yet attempted, False means too large (fall back to refcache).
        self.contig_seqs: dict[str, str] | None | bool = None
        #: Lazily opened BAM handle for chunk-range streaming (parallel path).
        self.bam: Any = None

    def contig(self, tid: int) -> str:
        try:
            return self.ref_names[tid]
        except IndexError:
            raise DataInvariantError(f"reference id {tid} missing from BAM header") from None

    def bind_header(self, header: Any) -> None:
        """Adopt the streamed BAM's header object (workers, lazily).

        Sharing the file's exact header object makes record construction and
        writer emission take pysam's fast same-header path; the output @PG
        header is stamped once at publish time by ``samtools reheader``.
        """
        self.header = header
        self.ref_names = [sq["SN"] for sq in (header.to_dict().get("SQ") or [])]
        self.tid_by_name = {name: tid for tid, name in enumerate(self.ref_names)}


_STATE: _WorkerState | None = None


def _worker_init(cfg: PairLoomConfig, out_header_dict: dict) -> None:
    global _STATE
    state = _WorkerState(cfg, None)
    # The output header (input header + PairLoom @PG records) is used by
    # chunk 0's writers: `samtools cat` keeps the FIRST input's header, so
    # the concatenated stream already carries the output header and no
    # post-hoc reheader pass is needed.
    state.out_header = pysam.AlignmentHeader.from_dict(out_header_dict)
    _STATE = state


def _ensure_state() -> _WorkerState:
    if _STATE is None:
        raise WorkerError("worker state not initialized")
    return _STATE


_BLOB_DTYPES = (
    np.uint64,
    np.uint64,
    np.uint8,
    np.uint8,
    np.uint64,
    np.uint8,
    np.uint8,
    np.uint64,
)


def _chunk_observation_arrays(
    obs: list[tuple[bytes, bytes, bytes, bytes, bytes, bytes]],
) -> tuple[np.ndarray, ...]:
    """Merge per-template overlap-aggregate blobs into chunk arrays (worker side).

    The Rust core merges the blobs; numpy views them without copying.
    """
    blobs = _core.merge_aggregates(obs)
    return tuple(np.frombuffer(b, dtype=dt) for b, dt in zip(blobs, _BLOB_DTYPES, strict=True))


def _process_chunk_range(
    bam_path: str,
    start_vo: int,
    end_vo: int,
    chunk_index: int,
    tmp_dir: str,
    streams: tuple[str, ...],
) -> tuple[list[tuple[str, str | None, int | None]], dict[str, str], tuple]:
    """Pool entry point: stream one template-group-aligned BAM range.

    Workers read the (name-sorted) BAM themselves at a virtual offset and
    write their output records into per-chunk BAM files, so records never
    cross the process boundary in either direction. Returns per-template
    metadata (classification, merge reason), the chunk's output file paths,
    and the chunk's aggregated observation arrays.
    """
    state = _ensure_state()
    # One BAM handle per worker: opening parses the full header (tens of ms
    # on large panels), which per-chunk opens would repeat hundreds of times.
    if state.bam is None:
        state.bam = pysam.AlignmentFile(bam_path, "rb")
        state.bind_header(state.bam.header)
    bam = state.bam
    writers: dict[str, Any] = {}
    paths: dict[str, str] = {}
    try:
        for stream in streams:
            path = os.path.join(tmp_dir, f"{stream}_{chunk_index:06d}.bam")
            header = state.out_header if chunk_index == 0 else state.header
            writers[stream] = pysam.AlignmentFile(path, "wb", header=header)
            paths[stream] = path
        bam.seek(start_vo)
        meta: list[tuple[str, str | None, int | None, bool]] = []
        obs: list[tuple[bytes, bytes, bytes, bytes, bytes, bytes]] = []
        group: list[Any] = []
        current_name: str | None = None
        # Virtual offsets only need checking at template-group boundaries:
        # records within a group always belong to the chunk, and ``end_vo``
        # sits on a boundary by construction.
        while True:
            if bam.tell() >= end_vo:
                break
            record = next(bam)
            name = record.query_name
            if name is None:
                raise DataInvariantError("BAM record without QNAME")
            if group and name != current_name:
                result = process_template(group, state, emit_writers=writers)
                meta.append(
                    (
                        result.classification,
                        result.merge_reason,
                        result.overlap_len,
                        result.non_primary_evidence,
                    )
                )
                if result.observations is not None:
                    obs.append(result.observations)
                group = []
            group.append(record)
            current_name = name
        if group:
            result = process_template(group, state, emit_writers=writers)
            meta.append(
                (
                    result.classification,
                    result.merge_reason,
                    result.overlap_len,
                    result.non_primary_evidence,
                )
            )
            if result.observations is not None:
                obs.append(result.observations)
        return meta, paths, _chunk_observation_arrays(obs)
    finally:
        for writer in writers.values():
            with contextlib.suppress(OSError, ValueError):
                writer.close()
        # pysam writer objects participate in reference cycles; without an
        # explicit collection each chunk's writers (and their internal
        # header copies) linger in the worker for the whole run.
        writers.clear()
        gc.collect()


def _iter_chunk_ranges(input_file: Any, chunk_size: int) -> Any:
    """Yield (start, end) virtual-offset ranges covering template-aligned chunks.

    ``end`` is the virtual offset just after the last record of the chunk's
    final template group, which is exactly where the next chunk begins.
    """
    prev_vo = input_file.tell()
    start_vo = prev_vo
    last_name: str | None = None
    groups = 0
    for record in input_file.fetch(until_eof=True):
        name = record.query_name
        if name is None:
            raise DataInvariantError("BAM record without QNAME")
        if name != last_name:
            if last_name is not None:
                groups += 1
                if groups >= chunk_size:
                    yield start_vo, prev_vo
                    start_vo = prev_vo
                    groups = 0
            last_name = name
        prev_vo = input_file.tell()
    if last_name is not None:
        yield start_vo, prev_vo


# ---------------------------------------------------------------------------
# Per-template processing
# ---------------------------------------------------------------------------


def process_template(
    records: list[Any], state: _WorkerState, emit_writers: dict[str, Any] | None = None
) -> TemplateResult:
    """Process one template into output lines, classification, and observations.

    With ``emit_writers`` (stream name -> object with ``write(record)``),
    output records are written directly and ``out_lines`` stays empty —
    the caller owns serialization-free emission.
    """
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
    obs_payload: tuple[bytes, bytes, bytes, bytes, bytes, bytes] | None = None
    overlap_len: int | None = None

    if r1 is not None and r2 is not None:
        m1 = not r1.flag & FLAG_UNMAPPED
        m2 = not r2.flag & FLAG_UNMAPPED
        if m1 and m2:
            (
                merged_line,
                classification,
                merge_reason,
                obs_payload,
                overlap_len,
            ) = _both_mapped(r1, r2, others, state)
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
        if merge_reason is not None:
            record.set_tag(state.cfg.reason_tag, merge_reason, value_type="Z")

    out_lines: list[tuple[str, str]] = []
    if state.cfg.mode == "adjust":
        for record in records:
            if emit_writers is None:
                out_lines.append(("adjusted", record.to_string()))
            else:
                emit_writers["adjusted"].write(record)
    elif merged_line is not None:
        # The consensus record replaces the primary mates, and the
        # template's non-primary records are dropped by policy: their
        # evidence survives in the consensus record's preserved outer
        # soft clips plus the evidence tag. Refused templates keep every
        # record (the branch below), so nothing is silently discarded.
        if emit_writers is None:
            out_lines.append(("merged", merged_line.to_string()))
        else:
            emit_writers["merged"].write(merged_line)
    else:
        for record in records:
            if emit_writers is None:
                out_lines.append(("unmerged", record.to_string()))
            else:
                emit_writers["unmerged"].write(record)
    return TemplateResult(
        classification,
        merge_reason,
        out_lines,
        obs_payload,
        overlap_len,
        non_primary_evidence=(
            merged_line is not None and bool(others) and state.cfg.merge_non_primary
        ),
    )


def _build_core_structure(record: Any) -> Any:
    """Build the Rust-core read structure for a mapped record."""
    cigar = normalize_record_cigar(record.cigartuples)
    seq = record.query_sequence
    if seq is None:
        raise DataInvariantError(f"record {record.query_name!r} has no SEQ")
    quals = record.query_qualities
    rg = record.get_tag("RG") if record.has_tag("RG") else None
    return _core.build_structure(
        record.reference_id,
        record.reference_start,
        bool(record.is_reverse),
        record.mapping_quality,
        record.query_name,
        cigar,
        seq,
        list(quals) if quals is not None else None,
        rg,
    )


def _core_merge(
    s1: Any,
    s2: Any,
    cfg: PairLoomConfig,
    merged_mapq: int | None,
    rescue_bypass_mapq: bool,
    has_non_primary: bool,
) -> Any:
    """Call the Rust-core merge with the configuration's gate parameters."""
    return _core.merge_pair(
        s1,
        s2,
        cfg.min_mapq,
        cfg.min_overlap,
        cfg.min_fragment_length,
        cfg.max_fragment_length,
        cfg.min_quality_difference,
        cfg.max_consensus_quality,
        cfg.conflict_policy,
        cfg.library_orientation,
        merged_mapq,
        rescue_bypass_mapq,
        has_non_primary,
    )


#: Minimum outer-clip length that counts as read-through evidence.
_CHIMERIC_MIN_CLIP = 10
#: Memoization cap for the clip-containment test (sequence-pair keyed).
_CHIMERIC_CACHE_CAP = 65536


def _outer_clips(record: Any) -> tuple[str, str]:
    """A record's leading and trailing soft-clip bases, stored orientation."""
    cigar = record.cigartuples or ()
    seq = record.query_sequence or ""
    lead = trail = ""
    qi = 0
    for op, length in cigar:
        if op == CIGAR_S:
            if qi == 0:
                lead = seq[:length]
            else:
                trail = seq[qi : qi + length]
        if op in (CIGAR_M, CIGAR_I, CIGAR_S, 7, 8):
            qi += length
    return lead, trail


def _is_chimeric_read(r1: Any, r2: Any, state: _WorkerState) -> bool:
    """True when cross-contig mates show read-through evidence.

    A concatemer or chimera makes each mate read past its oligo into the
    partner's sequence, so a mate's outer clip occurs as a substring of the
    other mate's stored sequence. Pure placement-level check; memoized on
    the sequence triple (concatemers repeat massively on oligo pools).
    """
    key = (r1.query_name, r1.query_sequence, r2.query_sequence)
    cached = getattr(state, "chimeric_cache", None)
    if cached is not None and key in cached:
        return cached[key]
    seq1 = r1.query_sequence or ""
    seq2 = r2.query_sequence or ""
    verdict = False
    if seq1 and seq2:
        for clip in _outer_clips(r1):
            if len(clip) >= _CHIMERIC_MIN_CLIP and clip in seq2:
                verdict = True
                break
        if not verdict:
            for clip in _outer_clips(r2):
                if len(clip) >= _CHIMERIC_MIN_CLIP and clip in seq1:
                    verdict = True
                    break
    if cached is None:
        cached = {}
        state.chimeric_cache = cached
    if len(cached) >= _CHIMERIC_CACHE_CAP:
        cached.clear()
    cached[key] = verdict
    return verdict


def _both_mapped(
    r1: Any, r2: Any, others: list[Any], state: _WorkerState
) -> tuple[
    Any | None,
    str,
    str | None,
    tuple[bytes, bytes, bytes, bytes, bytes, bytes] | None,
    int | None,
]:
    """Jointly adjust a both-mapped pair; returns (merged record, class, reason, obs, overlap)."""
    cfg = state.cfg
    tid1, tid2 = r1.reference_id, r2.reference_id
    if tid1 < 0 or tid2 < 0:
        raise DataInvariantError(
            f"template {r1.query_name!r} has mapped records without reference id"
        )
    # SAM '=' SEQ shorthand (base identical to the reference) is expanded
    # before any scoring, statistics, or merging sees the record.
    expand_reference_matches(r1, state.contig(tid1), state.refcache)
    expand_reference_matches(r2, state.contig(tid2), state.refcache)
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

    obs_payload = None
    if r1.reference_id != r2.reference_id:
        # Cross-contig pairs keep their adjusted records; no overlap analysis.
        if rescue_eligible(
            "different_contig",
            orig_mapq1,
            orig_mapq2,
            cfg,
            shape_blocked=_rescue_shape_blocked(r1, r2, others, cfg),
        ):
            return _attempt_rescue(r1, r2, cands1, cands2, others, obs_payload, state)
        if _is_chimeric_read(r1, r2, state):
            return None, BOTH_ALIGNED, "chimeric_read", obs_payload, None
        return None, BOTH_ALIGNED, "different_contig", obs_payload, None

    s1 = _build_core_structure(r1)
    s2 = _build_core_structure(r2)
    # The records now carry the conservative ADJUSTED MAPQs; the consensus
    # record takes its output MAPQ from those, never from the originals.
    # Templates merged despite non-primary evidence emit MAPQ 0 instead
    # (rescue-consistent honesty policy for contested placements).
    evidence = bool(others) and cfg.merge_non_primary
    gate_non_primary = bool(others) and not cfg.merge_non_primary
    merged_mapq = 0 if evidence else min(s1.mapq, s2.mapq)
    # Eligibility filters (min_mapq) continue to apply to the ORIGINAL
    # aligner MAPQ (policy comment above).
    s1.mapq = orig_mapq1
    s2.mapq = orig_mapq2
    contig = state.contig(r1.reference_id)
    # Observations exist only where BOTH mates align a base (the aligned-
    # span intersection), so only that reference slice is fetched: distant
    # mates cost a read-sized fetch instead of the whole span between them.
    ov_lo = max(s1.start0, s2.start0)
    ov_hi = min(s1.end0, s2.end0)
    obs_payload = None
    if ov_lo < ov_hi:
        ref_seq = state.refcache.fetch(contig, ov_lo, ov_hi)
        blobs = _core.aggregate_overlap(
            s1, s2, r1.reference_id, ref_seq, ov_lo, cfg.min_mapq, cfg.min_base_quality
        )
        if blobs[0]:
            obs_payload = blobs

    outcome = _core_merge(s1, s2, cfg, merged_mapq, False, gate_non_primary)
    overlap_len = outcome.overlap_len
    if not outcome.ok:
        if rescue_eligible(
            outcome.reason,
            orig_mapq1,
            orig_mapq2,
            cfg,
            shape_blocked=_rescue_shape_blocked(r1, r2, others, cfg),
        ):
            return _attempt_rescue(r1, r2, cands1, cands2, others, obs_payload, state)
        return None, BOTH_ALIGNED, outcome.reason, obs_payload, None
    merged_record = _build_merged_record(
        r1.query_name, outcome, contig, state, non_primary=evidence
    )
    return merged_record, MERGEABLE, None, obs_payload, overlap_len


def _rescue_shape_blocked(r1: Any, r2: Any, others: list[Any], cfg: PairLoomConfig) -> bool:
    """Placement-independent merge exclusions, checked before any rescue search.

    These are the gates ``try_merge`` applies before its placement-level
    reasons on the same-contig path: non-primary records (only when
    ``--no-merge-non-primary`` restores the refusal), hard clips,
    missing qualities, reference skips, and conflicting read groups. The
    cross-contig path reaches rescue without running any merge gate, so
    the exclusions are enforced here uniformly — a template failing one
    keeps its original placements on every path.
    """
    if others and not cfg.merge_non_primary:
        return True
    for record in (r1, r2):
        cigar = record.cigartuples or ()
        if any(op in (CIGAR_H, CIGAR_N) for op, _length in cigar):
            return True
        if record.query_qualities is None:
            return True
    rg1 = r1.get_tag("RG") if r1.has_tag("RG") else None
    rg2 = r2.get_tag("RG") if r2.has_tag("RG") else None
    return rg1 is not None and rg2 is not None and rg1 != rg2


def _attempt_rescue(
    r1: Any,
    r2: Any,
    cands1: Any,
    cands2: Any,
    others: list[Any],
    obs_payload: tuple[bytes, bytes, bytes, bytes, bytes, bytes] | None,
    state: _WorkerState,
) -> tuple[
    Any | None,
    str,
    str | None,
    tuple[bytes, bytes, bytes, bytes, bytes, bytes] | None,
    int | None,
]:
    """Attempt joint low-MAPQ rescue for an eligible template.

    Rescue relocates the mates to the unique compatible joint placement
    (when one exists) and re-runs every merge gate except the original-MAPQ
    eligibility gate. Returns the template's full result tuple; the reason
    is a rescue outcome status when the pair was not relocated or a
    post-rescue merge-rejection reason when it was.
    """
    outcome = rescue_pair(r1, r2, cands1, cands2, state)
    if outcome.status != STATUS_RESCUED:
        return None, BOTH_ALIGNED, outcome.status, obs_payload, None
    _relocate_record(r1, outcome.p1, state)
    _relocate_record(r2, outcome.p2, state)
    # A relocation may flip a mate's strand; the mate-reverse flags (0x20)
    # must then follow the mates' actual strands in the output. Ordinary
    # (unrescued) pairs never change strand, so their flags pass through.
    for record, mate in ((r1, r2), (r2, r1)):
        if mate.is_reverse:
            record.flag = record.flag | 0x20
        else:
            record.flag = record.flag & ~0x20
    _sync_mate_fields(r1, r2)
    _set_proper_flag(r1, True)
    _set_proper_flag(r2, True)
    s1 = _build_core_structure(r1)
    s2 = _build_core_structure(r2)
    # Rescued pairs contribute no disagreement observations by policy: the
    # statistics keep testing the ORIGINAL aligner MAPQ, and an eligible
    # rescue always has a sub-threshold mate. Placements are not
    # re-observed at the rescued locus (documented in the README).
    merged = _core_merge(
        s1,
        s2,
        state.cfg,
        0,
        True,
        bool(others) and not state.cfg.merge_non_primary,
    )
    if not merged.ok:
        return None, BOTH_ALIGNED, merged.reason, None, None
    merged_record = _build_merged_record(
        r1.query_name,
        merged,
        state.contig(r1.reference_id),
        state,
        rescued=True,
        non_primary=bool(others) and state.cfg.merge_non_primary,
    )
    return merged_record, MERGEABLE, None, None, merged.overlap_len


_DISCOVER_CACHE_CAP = 262_144


def _discover(record: Any, state: _WorkerState) -> Any:
    """Discover candidate placements.

    BAM SEQ/CIGAR are stored in reference orientation for both strands, so
    base-level verification and exact-block sliding need no strand-specific
    transform. Results are memoized per worker on
    (tid, pos0, cigar, sequence): identical reads repeat massively on oligo
    panels, and discovery is a pure function of these inputs.
    """
    cigar = normalize_record_cigar(record.cigartuples)
    seq = record.query_sequence
    if seq is None:
        raise DataInvariantError(f"record {record.query_name!r} has no SEQ")
    key = (record.reference_id, record.reference_start, cigar, seq.upper())
    cached = state.discover_cache.get(key)
    if cached is not None:
        return cached
    if len(state.discover_cache) >= _DISCOVER_CACHE_CAP:
        state.discover_cache.clear()
    result = discover_candidates(
        tid=record.reference_id,
        pos0=record.reference_start,
        cigar=cigar,
        query_seq=seq.upper(),
        contig=state.contig(record.reference_id),
        refcache=state.refcache,
        cfg=state.cfg,
        read_name=record.query_name,
    )
    state.discover_cache[key] = result
    return result


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
    expand_reference_matches(record, state.contig(record.reference_id), state.refcache)
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


def _relocate_record(record: Any, placement: RescuePlacement, state: _WorkerState) -> None:
    """Move a record to a rescued global placement (possibly new contig/strand).

    A strand change rewrites the stored sequence into the new alignment
    orientation — SAM stores SEQ/QUAL in reference orientation for both
    strands — by reverse-complementing SEQ, reversing QUAL, and reversing
    the CIGAR op order (clips swap ends). The rescued placement is exact by
    construction, so NM/MD are orientation-independent; the output MAPQ is
    0 (a unique joint placement is not a calibrated confidence); stale
    aligner scoring tags are removed rather than guessed at.
    """
    strand_changed = placement.is_reverse != bool(record.is_reverse)
    if strand_changed:
        seq = record.query_sequence
        if seq is None:  # pragma: no cover - discovery fails earlier on missing SEQ
            raise DataInvariantError(f"record {record.query_name!r} has no SEQ to relocate")
        quals = record.query_qualities
        record.query_sequence = revcomp(seq, f"rescue of {record.query_name!r}")
        # pysam's query_sequence setter drops query_qualities; restore them
        # in the new orientation (their order follows the stored sequence).
        if quals is not None:
            record.query_qualities = quals[::-1]
        record.cigartuples = list(reverse_cigar(tuple(record.cigartuples or ())))
    record.reference_id = placement.tid
    record.reference_start = placement.pos0
    if placement.is_reverse:
        record.flag = record.flag | 0x10
    else:
        record.flag = record.flag & ~0x10
    record.mapping_quality = 0
    record.set_tag("MD", str(placement.span), value_type="Z")
    record.set_tag("NM", 0, value_type="i")
    for tag in STALE_TAGS_REMOVED:
        if record.has_tag(tag):
            # pysam (0.24) has no remove_tag(); set_tag(tag, None) deletes.
            record.set_tag(tag, None)
    assign_rescue_tag(record)


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
    rescued: bool = False,
    non_primary: bool = False,
) -> Any:
    """Create the single consensus BAM record for a merged pair."""
    record = pysam.AlignedSegment(state.header)
    record.query_name = qname
    record.flag = 0
    record.reference_id = state.tid_by_name[contig]
    record.reference_start = merged.start0
    record.cigartuples = list(merged.cigar)
    record.query_sequence = merged.seq
    record.query_qualities = list(merged.quals)
    record.mapping_quality = merged.mapq
    record.next_reference_id = -1
    record.next_reference_start = -1
    record.template_length = 0
    record.set_tag(state.cfg.status_tag, MERGEABLE, value_type="Z")
    if rescued:
        # Jointly rescued pairs: output MAPQ 0 (passed by the caller) plus
        # the rescue marker so downstream tools can tell them apart.
        record.set_tag(RESCUE_TAG, RESCUE_VALUE, value_type="i")
    if non_primary:
        # Merged despite secondary/supplementary evidence: the consensus
        # carries output MAPQ 0 (passed by the caller) and the evidence
        # marker naming what was dropped (the clips stay in the CIGAR).
        record.set_tag(state.cfg.evidence_tag, 1, value_type="i")
    if merged.read_group is not None:
        record.set_tag("RG", merged.read_group, value_type="Z")
    ref_end = merged.start0 + _merged_ref_len(merged.cigar)
    ref_seq = state.refcache.fetch(contig, merged.start0, ref_end)
    md, nm = _core.compute_md_nm(merged.start0, list(merged.cigar), merged.seq, ref_seq)
    record.set_tag("MD", md, value_type="Z")
    record.set_tag("NM", nm, value_type="i")
    return record


def _merged_ref_len(cigar: tuple[tuple[int, int], ...]) -> int:
    from pairloom.cigarops import cigar_ref_len

    return cigar_ref_len(cigar)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


_CAT_FANIN = 64


def _concat_bams(dest: Path, files: list[str], tmp_dir: str) -> None:
    """Concatenate BAM chunks into ``dest`` with bounded fan-in.

    ``samtools cat`` keeps every input open simultaneously and parses each
    input's header, so a single call over hundreds of chunk files would
    scale memory with the chunk count (megabytes per contig-heavy header).
    Cascading 64-way batches keeps the open set bounded; intermediate
    batches are consumed as they are produced.
    """
    current = list(files)
    generation = 0
    try:
        while len(current) > _CAT_FANIN:
            reduced: list[str] = []
            for start in range(0, len(current), _CAT_FANIN):
                group = current[start : start + _CAT_FANIN]
                if len(group) == 1:
                    reduced.append(group[0])
                    continue
                merged = os.path.join(tmp_dir, f"cat_{dest.stem}_{generation}_{start:06d}.bam")
                pysam.samtools.cat("-o", merged, *group)
                for path in group:
                    os.unlink(path)
                reduced.append(merged)
            current = reduced
            generation += 1
        pysam.samtools.cat("-o", str(dest), *current)
    except (OSError, ValueError) as exc:
        raise OutputError(f"output concatenation failed for {dest}: {exc}") from exc
    finally:
        for path in current:
            with contextlib.suppress(OSError):
                if path != str(dest):
                    os.unlink(path)


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
    merged_non_primary = 0
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
                stream_path = str(cfg.input_path)
                if sort_order != "queryname":
                    input_file.close()
                    sorted_input = name_sort_to_temp(str(cfg.input_path), tmpdir)
                    input_file = open_input_bam(sorted_input)
                    stream_path = sorted_input
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
                        # Writers share the INPUT header object so records
                        # stream through pysam's fast same-header path; the
                        # output @PG header is stamped at publish time.
                        writers[stream] = _BamWriter(outputs[role], input_file.header)

                aggregator = FastDisagreementAccumulator(cfg)
                aggregator.set_spill_dir(tmpdir)

                merge_stats = MergeStatsCollector()

                def on_meta(
                    classification: str,
                    merge_reason: str | None,
                    overlap_len: int | None,
                    non_primary_evidence: bool = False,
                ) -> None:
                    nonlocal merged_non_primary
                    stats.add(classification)
                    if merge_reason:
                        reasons[merge_reason] += 1
                    if non_primary_evidence:
                        merged_non_primary += 1
                    merge_stats.add(classification, merge_reason, overlap_len)

                chunks = _iter_chunks(input_file, cfg.chunk_size)
                chunk_files: dict[str, list[str]] = {}
                if cfg.threads == 1:
                    state = _WorkerState(cfg, input_file.header)
                    try:
                        for chunk in chunks:
                            obs = []
                            for records in chunk:
                                result = process_template(records, state, emit_writers=writers)
                                on_meta(
                                    result.classification,
                                    result.merge_reason,
                                    result.overlap_len,
                                    result.non_primary_evidence,
                                )
                                if result.observations is not None:
                                    obs.append(result.observations)
                            aggregator.add_chunk(*_chunk_observation_arrays(obs))
                    finally:
                        state.refcache.close()
                else:
                    chunk_files = _run_parallel(
                        input_file,
                        stream_path,
                        cfg,
                        out_header_dict,
                        on_meta,
                        aggregator,
                        tmpdir,
                        tuple(writers.keys()),
                    )

                for writer in writers.values():
                    writer.close()
                    tmp_files.append(writer.tmp_path)
                needs_reheader = cfg.threads == 1 or not any(
                    chunk_files.get(stream) for stream in writers
                )
                if writers and needs_reheader and str(header_obj) != str(input_file.header):
                    header_sam = os.path.join(tmpdir, "pairloom_out_header.sam")
                    with open(header_sam, "w") as fh:
                        fh.write(str(header_obj))
                    for writer in writers.values():
                        reheadered = Path(str(writer.tmp_path) + ".rh")
                        try:
                            with open(reheadered, "wb") as fh:
                                subprocess.run(
                                    ["samtools", "reheader", header_sam, str(writer.tmp_path)],
                                    stdout=fh,
                                    check=True,
                                )
                        except (OSError, subprocess.CalledProcessError) as exc:
                            raise OutputError(
                                f"output reheader failed for {writer.tmp_path}: {exc}"
                            ) from exc
                        os.replace(reheadered, writer.tmp_path)
                if cfg.threads != 1:
                    # Parallel path: concatenate the workers' per-chunk BAMs
                    # (identical headers; chunk order == template order) into
                    # each stream's temp file, replacing the empty header-only
                    # stream the writer created.
                    for stream, writer in writers.items():
                        files = chunk_files.get(stream) or []
                        if not files:
                            continue
                        _concat_bams(writer.tmp_path, files, tmpdir)

                disag_tmp = Path(str(cfg.disagreement_parquet) + f".pairloom-tmp-{os.getpid()}")
                stats_tmp = Path(str(cfg.read_stats_parquet) + f".pairloom-tmp-{os.getpid()}")
                mstats_tmp: Path | None = None
                if cfg.merge_stats_parquet is not None:
                    mstats_tmp = Path(str(cfg.merge_stats_parquet) + f".pairloom-tmp-{os.getpid()}")
                    write_merge_stats_parquet(
                        mstats_tmp, merge_stats.rows(), cfg.parquet_compression
                    )
                # Register temps before writing so failures clean them up.
                tmp_files.extend([disag_tmp, stats_tmp] + ([mstats_tmp] if mstats_tmp else []))
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
                if mstats_tmp is not None:
                    assert cfg.merge_stats_parquet is not None  # noqa: S101
                    pairs.append((mstats_tmp, cfg.merge_stats_parquet))
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
                        "pairloom: warning: could not remove backup file(s); "
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
        "pairloom: done. "
        + ", ".join(f"{k}={v}" for k, v in stats.counts().items())
        + (f"; merge-reasons={dict(reasons)}" if reasons else "")
        + f"; merged_non_primary={merged_non_primary}"
        + f"; disagreement-positions={aggregator.reported_positions}",
        file=sys.stderr,
    )
    return RunSummary(
        classification_counts=stats.counts(),
        merge_reasons=dict(reasons),
        merged_non_primary=merged_non_primary,
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


def _run_parallel(
    input_file: Any,
    stream_path: str,
    cfg: PairLoomConfig,
    header_dict: dict,
    on_meta: Any,
    aggregator: Any,
    tmp_dir: str,
    streams: tuple[str, ...],
) -> dict[str, list[str]]:
    """Process chunks in a process pool with bounded, ordered submission.

    Chunks are virtual-offset ranges; each worker streams its range from
    ``stream_path`` and writes its output records into per-chunk BAM files,
    so records never serialize across the IPC boundary in either direction.
    Returns each stream's chunk file list in chunk order for the final
    concatenation. The pool forks when the platform allows it (warm
    imports; no re-exec of the interpreter per worker).
    """
    import multiprocessing

    methods = multiprocessing.get_all_start_methods()
    ctx = multiprocessing.get_context("fork" if "fork" in methods else "spawn")
    max_inflight = cfg.threads * 2
    chunk_files: dict[str, list[str]] = {stream: [] for stream in streams}
    with ProcessPoolExecutor(
        max_workers=cfg.threads,
        mp_context=ctx,
        initializer=_worker_init,
        initargs=(cfg, header_dict),
    ) as pool:
        pending: deque = deque()
        chunk_index = 0

        def submit_next(start_vo: int, end_vo: int) -> None:
            nonlocal chunk_index
            pending.append(
                pool.submit(
                    _process_chunk_range,
                    stream_path,
                    start_vo,
                    end_vo,
                    chunk_index,
                    tmp_dir,
                    streams,
                )
            )
            chunk_index += 1

        try:
            for start_vo, end_vo in _iter_chunk_ranges(input_file, cfg.chunk_size):
                submit_next(start_vo, end_vo)
                while len(pending) >= max_inflight:
                    _drain_one(pending, on_meta, aggregator, chunk_files)
            while pending:
                _drain_one(pending, on_meta, aggregator, chunk_files)
        except BaseException:
            for future in pending:
                future.cancel()
            raise
    return chunk_files


def _drain_one(pending: deque, on_meta: Any, aggregator: Any, chunk_files: dict) -> None:
    future = pending.popleft()
    try:
        meta, paths, arrays = future.result()
    except PairLoomError:
        raise
    except Exception as exc:
        raise WorkerError(f"worker chunk failed: {exc}") from exc
    for classification, merge_reason, overlap_len, non_primary_evidence in meta:
        on_meta(classification, merge_reason, overlap_len, non_primary_evidence)
    for stream, path in paths.items():
        chunk_files[stream].append(path)
    aggregator.add_chunk(*arrays)


def _warn_committed() -> None:
    """Post-commit failure notice: outputs are complete and stay in place."""
    print(
        "pairloom: warning: failure after publication; all outputs were "
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
    pg_id = "pairloom"
    existing_ids = {pg.get("ID") for pg in pg_list}
    suffix = 0
    while pg_id in existing_ids:
        suffix += 1
        pg_id = f"pairloom.{suffix}"
    pg_list.append(
        {
            "ID": pg_id,
            "PN": "pairloom",
            "VN": __version__,
            "CL": (
                f"pairloom mode={cfg.mode} input={cfg.input_path} reference={cfg.reference_path}"
            ),
        }
    )
    out["PG"] = pg_list
    return out
