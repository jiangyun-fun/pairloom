"""Repeat-aware candidate placement discovery.

A BAM carries one primary placement per read even when the read matches
repeatedly. Candidate discovery therefore distinguishes:

- placements present in the BAM (``--candidate-source bam``), and
- additional equivalent placements found by sliding an exact-match ungapped
  block along the reference inside ``--realign-window`` (default source).

Discovery is bounded and never claims global uniqueness: if the candidate cap
is hit, the candidate set is flagged truncated and the MAPQ policy applies a
conservative cap.
"""

from __future__ import annotations

from dataclasses import dataclass

from pairloom.alignment import (
    Placement,
    aligned_block,
    score_placement,
    slidable_cigar,
)
from pairloom.cigarops import CIGAR_EQUAL, CIGAR_X, cigar_ref_len, normalize_cigar
from pairloom.config import PairLoomConfig
from pairloom.exceptions import DataInvariantError
from pairloom.reference import ReferenceCache


@dataclass(frozen=True, slots=True)
class CandidateSet:
    """Candidate placements for one read plus discovery status."""

    placements: tuple[Placement, ...]
    truncated: bool


def normalize_record_cigar(cigar: list[tuple[int, int]] | None) -> tuple[tuple[int, int], ...]:
    """Normalize a record CIGAR: map =/X onto M, merge adjacent equal ops."""
    if not cigar:
        raise DataInvariantError("record has no CIGAR")
    mapped = tuple((0 if op in (CIGAR_EQUAL, CIGAR_X) else op, length) for op, length in cigar)
    return normalize_cigar(mapped)


def discover_candidates(
    tid: int,
    pos0: int,
    cigar: tuple[tuple[int, int], ...],
    query_seq: str,
    contig: str,
    refcache: ReferenceCache,
    cfg: PairLoomConfig,
    read_name: str,
) -> CandidateSet:
    """Return all deterministic candidate placements for one read.

    Only exact-match ungapped blocks are slid: for such reads every candidate
    is equivalent by construction, so sliding cannot override alignment
    evidence. Mismatch- or indel-bearing alignments keep only their current
    placement (documented limitation).
    """
    context = f"read {read_name} at {contig}:{pos0 + 1}"
    contig_len = refcache.length(contig)
    span = cigar_ref_len(cigar)
    if pos0 < 0 or pos0 + span > contig_len:
        raise DataInvariantError(
            f"alignment [{pos0}, {pos0 + span}) out of bounds for {contig} "
            f"(length {contig_len}) in {context}"
        )
    ref_seq = refcache.fetch(contig, pos0, pos0 + span)
    score, exact = score_placement(pos0, cigar, query_seq, ref_seq, context)
    current = Placement(tid=tid, pos0=pos0, cigar=cigar, score=score, exact=exact)

    if cfg.candidate_source == "bam":
        return CandidateSet((current,), False)
    if not exact or not slidable_cigar(cigar):
        return CandidateSet((current,), False)

    block = aligned_block(cigar, query_seq)
    if len(block) != span:  # pragma: no cover - slidable cigars guarantee this
        raise DataInvariantError(f"aligned block length mismatch in {context}")

    lo = max(0, pos0 - cfg.realign_window)
    hi = min(contig_len - span, pos0 + cfg.realign_window)
    if hi < lo:
        return CandidateSet((current,), False)
    window_seq = refcache.fetch(contig, lo, hi + span)
    occurrences: list[int] = []
    idx = window_seq.find(block)
    while idx != -1:
        occurrences.append(lo + idx)
        idx = window_seq.find(block, idx + 1)
    if not occurrences:  # pragma: no cover - current placement always matches
        raise DataInvariantError(f"exact block not found at own position in {context}")
    truncated = len(occurrences) > cfg.max_candidates_per_read
    kept = sorted(occurrences)[: cfg.max_candidates_per_read]
    placements = tuple(
        Placement(tid=tid, pos0=p, cigar=cigar, score=len(block), exact=True) for p in kept
    )
    return CandidateSet(placements, truncated)


def single_read_leftmost(cands: CandidateSet) -> Placement:
    """Deterministic leftmost choice for orphan/single-read normalization."""
    return min(cands.placements, key=lambda p: p.pos0)
