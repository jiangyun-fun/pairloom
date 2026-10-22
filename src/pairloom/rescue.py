"""Joint low-MAPQ pair rescue: global cross-reference placement search.

When a both-mapped template fails merging for a placement-level reason and
at least one mate's ORIGINAL MAPQ is below ``--min-mapq``, the aligner's
primary placements carry no usable joint evidence. Rescue re-places the
template jointly:

- both mates' exact-match ungapped aligned blocks are searched over the
  WHOLE reference (every BAM-header contig, both strands), not just the
  original contig and window;
- placements are combined into candidate pairs under the same constraints
  as ordinary selection (same reference, library orientation model,
  fragment-length bounds);
- the template is rescued only when EXACTLY ONE compatible pair exists and
  the search was complete (neither mate's occurrence list truncated by
  ``--max-candidates-per-read``, pair enumeration not capped by
  ``--max-pair-evaluations``). A deterministic tie-break never counts as
  evidence of uniqueness: two equally supported placements mean the pair
  stays ambiguous and unmerged.

Uniqueness claims therefore come from a global search by construction —
the ``--realign-window`` of ordinary sliding never bounds them.

Search cost: one C-speed substring scan of the reference per rescued mate
block (chunked at ``_SCAN_CHUNK`` fetches so memory stays bounded on large
contigs; contigs at or below the chunk size are fetched whole once). The
scan aborts as soon as the occurrence cap is exceeded — the outcome
(truncated) is already decided. Only templates that failed the placement
checks pay this cost; ``--candidate-source bam`` and ``--no-rescue`` avoid
it entirely.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pairloom.alignment import Placement, aligned_block, revcomp, slidable_cigar
from pairloom.exceptions import DataInvariantError
from pairloom.merge import (
    REASON_CONTIG,
    REASON_FRAGMENT,
    REASON_MAPQ,
    REASON_ORIENTATION,
)
from pairloom.pairing import fragment_ok, leftmost_is_read1, orientation_ok

#: Rescue outcomes double as documented merge-rejection reasons.
STATUS_RESCUED = "rescue_applied"
STATUS_AMBIGUOUS = "rescue_ambiguous"
STATUS_NO_PAIR = "rescue_no_compatible_pair"
STATUS_TRUNCATED = "rescue_search_truncated"
STATUS_UNSUPPORTED = "rescue_unsupported_alignment"

#: Merge-failure reasons that joint re-placement can address. Record-shape
#: failures (non-primary, hard clips, missing qualities, reference skips,
#: read-group conflicts) and evidence conflicts are placement-independent
#: and never trigger rescue.
TRIGGER_REASONS: frozenset[str] = frozenset(
    {REASON_MAPQ, REASON_CONTIG, REASON_ORIENTATION, REASON_FRAGMENT}
)

#: Reference fetch granularity for occurrence scans on large contigs.
_SCAN_CHUNK = 1_000_000


@dataclass(frozen=True, slots=True)
class RescuePlacement:
    """One global placement of a mate's aligned block."""

    tid: int
    pos0: int
    is_reverse: bool
    span: int

    @property
    def end0(self) -> int:
        return self.pos0 + self.span


@dataclass(frozen=True, slots=True)
class RescueOutcome:
    status: str
    p1: RescuePlacement | None = None
    p2: RescuePlacement | None = None


def rescue_eligible(failure_reason: str, orig_mapq1: int, orig_mapq2: int, cfg: Any) -> bool:
    """True when the failed template should attempt joint rescue.

    Requires a placement-level failure, at least one mate below the
    original-MAPQ gate (high-MAPQ discordant pairs keep their placements —
    their discordance is evidence, not ambiguity), reference-source
    discovery, and rescue not disabled.
    """
    return (
        cfg.rescue_enabled
        and cfg.candidate_source == "reference"
        and failure_reason in TRIGGER_REASONS
        and (orig_mapq1 < cfg.min_mapq or orig_mapq2 < cfg.min_mapq)
    )


def rescue_pair(r1: Any, r2: Any, cands1: Any, cands2: Any, state: Any) -> RescueOutcome:
    """Search globally for a unique compatible joint placement."""
    blocks: list[str] = []
    spans: list[int] = []
    strands: list[bool] = []
    for read, cands in ((r1, cands1), (r2, cands2)):
        current: Placement | None = cands.current
        if current is None or not current.exact or not slidable_cigar(current.cigar):
            return RescueOutcome(STATUS_UNSUPPORTED)
        seq = read.query_sequence
        if seq is None:  # pragma: no cover - discovery fails earlier on missing SEQ
            return RescueOutcome(STATUS_UNSUPPORTED)
        block = aligned_block(current.cigar, seq.upper())
        if len(block) != current.ref_span:
            # Slidable cigars are a single M block, so this cannot happen;
            # guarded because a wrong span would corrupt every coordinate.
            raise DataInvariantError(
                f"rescue block length {len(block)} != reference span {current.ref_span} "
                f"for read {read.query_name!r}"
            )
        blocks.append(block)
        spans.append(current.ref_span)
        strands.append(bool(read.is_reverse))

    occ1, truncated = _global_occurrences(blocks[0], strands[0], spans[0], state)
    if truncated:
        return RescueOutcome(STATUS_TRUNCATED)
    occ2, truncated = _global_occurrences(blocks[1], strands[1], spans[1], state)
    if truncated:
        return RescueOutcome(STATUS_TRUNCATED)

    pairs, complete = _compatible_pairs(occ1, occ2, spans, state.cfg)
    if not complete:
        return RescueOutcome(STATUS_TRUNCATED)
    if not pairs:
        return RescueOutcome(STATUS_NO_PAIR)
    if len(pairs) > 1:
        return RescueOutcome(STATUS_AMBIGUOUS)
    return RescueOutcome(STATUS_RESCUED, pairs[0][0], pairs[0][1])


def _global_occurrences(
    block: str, record_reverse: bool, span: int, state: Any
) -> tuple[list[RescuePlacement], bool]:
    """All placements of one mate's aligned block (both strands), globally.

    Strand semantics: BAM stores the aligned block in reference orientation,
    so ``block`` equals the reference bases at the CURRENT placement — for a
    forward-strand record that locus is a forward placement, for a
    reverse-strand record a reverse placement. The raw read is
    ``revcomp(block)`` for reverse records; its forward placements are the
    occurrences of the raw sequence and its reverse placements the
    occurrences of ``revcomp(raw)``. Labelling the probes the other way
    would place reads on strands their sequence cannot support.

    Returns ``(occurrences, truncated)``; ``truncated`` means the count
    exceeded ``--max-candidates-per-read`` (the scan stopped early — the
    rescue outcome is already decided and uniqueness can never be claimed).
    """
    cap = state.cfg.max_candidates_per_read
    raw = revcomp(block, "rescue block") if record_reverse else block
    probes = ((raw, False), (revcomp(raw, "rescue block"), True))
    occurrences: list[RescuePlacement] = []
    for name in state.refcache.contig_names:
        tid = state.tid_by_name.get(name)
        if tid is None:
            # Reference contigs absent from the BAM header cannot hold
            # expressible output placements; the pipeline ignores them
            # everywhere else, and so does rescue.
            continue
        contig_len = state.refcache.length(name)
        for probe, is_reverse in probes:
            for pos in _scan_contig(state.refcache, name, contig_len, probe):
                occurrences.append(RescuePlacement(tid, pos, is_reverse, span))
            if len(occurrences) > cap:
                return occurrences, True
    return occurrences, False


def _scan_contig(refcache: Any, name: str, contig_len: int, probe: str) -> list[int]:
    """Ascending start positions of ``probe`` in one contig, chunked fetches.

    Chunks overlap by ``len(probe) - 1`` so occurrences spanning a chunk
    boundary are found exactly once (positions at or beyond the next
    chunk's base are deferred to it).
    """
    span = len(probe)
    if span > contig_len:
        return []
    hits: list[int] = []
    if contig_len <= _SCAN_CHUNK:
        seq = refcache.fetch(name, 0, contig_len)
        idx = seq.find(probe)
        while idx != -1:
            hits.append(idx)
            idx = seq.find(probe, idx + 1)
        return hits
    base = 0
    while base + span <= contig_len:
        window_end = min(base + _SCAN_CHUNK, contig_len)
        chunk = refcache.fetch(name, base, min(window_end + span - 1, contig_len))
        idx = chunk.find(probe)
        while idx != -1:
            pos = base + idx
            if pos >= window_end:
                break
            hits.append(pos)
            idx = chunk.find(probe, idx + 1)
        if window_end >= contig_len:
            break
        base = window_end
    return hits


def _compatible_pairs(
    occ1: list[RescuePlacement],
    occ2: list[RescuePlacement],
    spans: list[int],
    cfg: Any,
) -> tuple[list[tuple[RescuePlacement, RescuePlacement]], bool]:
    """Constraint-satisfying joint placements plus search completeness.

    Enumeration order is deterministic: candidates sorted by
    (tid, pos0, strand). The evaluation budget mirrors ``select_pair``:
    same-contig combinations count, and ending exactly at the cap counts
    as truncated (unexplored combinations could hold another valid pair).
    """
    list1 = sorted(occ1, key=lambda p: (p.tid, p.pos0, p.is_reverse))
    list2 = sorted(occ2, key=lambda p: (p.tid, p.pos0, p.is_reverse))
    budget = cfg.max_pair_evaluations
    evaluated = 0
    pairs: list[tuple[RescuePlacement, RescuePlacement]] = []
    for p1 in list1:
        for p2 in list2:
            if p1.tid != p2.tid:
                continue
            if evaluated >= budget:
                break
            evaluated += 1
            left_is_read1 = leftmost_is_read1(
                p1.pos0, p2.pos0, p1.is_reverse, p2.is_reverse, cfg.library_orientation
            )
            if not orientation_ok(
                p1.is_reverse, p2.is_reverse, left_is_read1, cfg.library_orientation
            ):
                continue
            if not fragment_ok(_as_placement(p1, spans[0]), _as_placement(p2, spans[1]), cfg):
                continue
            pairs.append((p1, p2))
        if evaluated >= budget:
            break
    return pairs, evaluated < budget


def _as_placement(rp: RescuePlacement, span: int) -> Placement:
    """Adapter so the shared ``fragment_ok`` bounds check applies unchanged."""
    return Placement(tid=rp.tid, pos0=rp.pos0, cigar=((0, span),), score=span, exact=True)
