"""Consensus merging of overlapping mate pairs.

Merging operates on reference coordinates via :class:`ReadStructure`
objects, never by naive string concatenation. The merged record represents
the fragment in reference orientation:

- outer soft clips of the leftmost/rightmost reads are preserved;
- internal soft clips are dropped (the other mate's aligned bases cover
  those reference positions);
- overlapping aligned bases appear exactly once, resolved by the consensus
  policy;
- deletions observed by one or both mates survive; an indel interpretation
  that conflicts with the other mate's aligned bases makes the pair
  unmergeable with an explicit reason;
- no sequence is synthesized across uncovered fragment gaps.

The merged CIGAR is reconstructed from the per-reference-position walk, not
by concatenating mate CIGARs.
"""

from __future__ import annotations

from dataclasses import dataclass

from pairloom.alignment import ACGT, ReadStructure, covers
from pairloom.cigarops import CIGAR_D, CIGAR_I, CIGAR_M, CIGAR_S, normalize_cigar
from pairloom.config import PairLoomConfig
from pairloom.pairing import orientation_ok

#: Documented quality of an unresolved conflicting base (policy ``quality``).
UNRESOLVED_QUALITY = 2


@dataclass(frozen=True, slots=True)
class MergeSuccess:
    start0: int
    cigar: tuple[tuple[int, int], ...]
    seq: str
    quals: tuple[int, ...]
    mapq: int
    overlap_len: int
    read_group: str | None


@dataclass(frozen=True, slots=True)
class MergeFailure:
    reason: str


MergeOutcome = MergeSuccess | MergeFailure

# Documented merge-rejection reasons (expected biological outcomes, not errors).
REASON_NON_PRIMARY = "non_primary_records"
REASON_HARD_CLIP = "hard_clipped_mate"
REASON_NO_QUALITIES = "qualities_missing"
REASON_MAPQ = "mapq_below_min"
REASON_SKIP = "reference_skip_unsupported"
REASON_CONTIG = "different_contig"
REASON_ORIENTATION = "orientation_incompatible"
REASON_FRAGMENT = "fragment_length_out_of_bounds"
REASON_OVERLAP = "insufficient_overlap"
REASON_RG = "read_group_conflict"
REASON_GAP = "fragment_gap"
REASON_INDEL = "indel_conflict"
REASON_INSERTION = "insertion_conflict"
REASON_STRICT = "strict_conflict"
REASON_EMPTY = "empty_consensus"


class _StrictSignal:
    """Sentinel: strict policy rejects an unresolved base conflict."""


STRICT = _StrictSignal()


def try_merge(
    struct1: ReadStructure,
    struct2: ReadStructure,
    cfg: PairLoomConfig,
    has_non_primary: bool = False,
) -> MergeOutcome:
    """Attempt to merge two adjusted mates into a consensus fragment."""
    if has_non_primary:
        return MergeFailure(REASON_NON_PRIMARY)
    if struct1.hard_clipped or struct2.hard_clipped:
        return MergeFailure(REASON_HARD_CLIP)
    if not struct1.has_qualities or not struct2.has_qualities:
        return MergeFailure(REASON_NO_QUALITIES)
    if struct1.mapq < cfg.min_mapq or struct2.mapq < cfg.min_mapq:
        return MergeFailure(REASON_MAPQ)
    if struct1.skipped or struct2.skipped:
        return MergeFailure(REASON_SKIP)
    if struct1.tid != struct2.tid:
        return MergeFailure(REASON_CONTIG)
    left_is_read1 = struct1.start0 <= struct2.start0
    if not orientation_ok(
        struct1.is_reverse, struct2.is_reverse, left_is_read1, cfg.library_orientation
    ):
        return MergeFailure(REASON_ORIENTATION)
    start = min(struct1.start0, struct2.start0)
    end = max(struct1.end0, struct2.end0)
    span = end - start
    if not cfg.min_fragment_length <= span <= cfg.max_fragment_length:
        return MergeFailure(REASON_FRAGMENT)
    overlap = len(struct1.bases.keys() & struct2.bases.keys())
    if overlap < cfg.min_overlap:
        return MergeFailure(REASON_OVERLAP)
    if (
        struct1.read_group is not None
        and struct2.read_group is not None
        and struct1.read_group != struct2.read_group
    ):
        return MergeFailure(REASON_RG)

    seq_parts: list[str] = []
    qual_parts: list[int] = []
    ops: list[tuple[int, int]] = []

    left = struct1 if left_is_read1 else struct2
    right = struct2 if left_is_read1 else struct1

    def emit(op: int, length: int, bases: str, quals: tuple[int, ...] | list[int]) -> None:
        if length <= 0:
            return
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + length)
        else:
            ops.append((op, length))
        seq_parts.append(bases)
        qual_parts.extend(quals)

    # Outer soft clips: the leftmost read's leading clip and the rightmost
    # read's trailing clip, both already in reference orientation.
    emit(CIGAR_S, len(left.lead_clip[0]), left.lead_clip[0], left.lead_clip[1])

    for rpos in range(start, end):
        ins = _merge_insertion(rpos, struct1, struct2)
        if isinstance(ins, MergeFailure):
            return ins
        if ins is not None:
            emit(CIGAR_I, len(ins[0]), ins[0], ins[1])
        in1 = struct1.bases.get(rpos)
        in2 = struct2.bases.get(rpos)
        d1 = rpos in struct1.deleted
        d2 = rpos in struct2.deleted
        if in1 is None and in2 is None:
            # Neither mate aligns a base: a deletion observed by either mate
            # spans this position (the other mate cannot contradict it here).
            if d1 or d2:
                emit(CIGAR_D, 1, "", [])
                continue
            return MergeFailure(REASON_GAP)
        if (in1 is not None and d2) or (in2 is not None and d1):
            # One mate aligns a base where the other deletes it: the indel
            # interpretations are irreconcilably inconsistent.
            return MergeFailure(REASON_INDEL)
        resolved = _consensus_base(in1, in2, cfg)
        if resolved is STRICT:
            return MergeFailure(REASON_STRICT)
        base, qual = resolved
        emit(CIGAR_M, 1, base, [qual])

    ins = _merge_insertion(end, struct1, struct2)
    if isinstance(ins, MergeFailure):
        return ins
    if ins is not None:
        emit(CIGAR_I, len(ins[0]), ins[0], ins[1])

    emit(CIGAR_S, len(right.trail_clip[0]), right.trail_clip[0], right.trail_clip[1])
    if not any(op in (CIGAR_M, CIGAR_I, CIGAR_D) for op, _ in ops):
        return MergeFailure(REASON_EMPTY)
    return MergeSuccess(
        start0=start,
        cigar=normalize_cigar(tuple(ops)),
        seq="".join(seq_parts),
        quals=tuple(qual_parts),
        mapq=min(struct1.mapq, struct2.mapq),
        overlap_len=overlap,
        read_group=struct1.read_group or struct2.read_group,
    )


def _consensus_base(
    in1: tuple[str, int] | None,
    in2: tuple[str, int] | None,
    cfg: PairLoomConfig,
) -> tuple[str, int] | None | _StrictSignal:
    """Resolve one reference position; None when no mate aligns a base here.

    Quality model (independence of sequencing errors assumed, documented in
    the README): agreeing bases combine as ``min(q1 + q2, cap)``; a
    disagreement is resolved toward the higher-quality base when the gap
    reaches ``--min-quality-difference``; unresolved disagreements emit ``N``
    at low quality under policy ``quality`` and reject the pair under
    policy ``strict``.
    """
    if in1 is None and in2 is None:
        return None
    if in2 is None:
        return in1
    if in1 is None:
        return in2
    base1, qual1 = in1
    base2, qual2 = in2
    if base1 == base2:
        if base1 in ACGT:
            return base1, min(qual1 + qual2, cfg.max_consensus_quality, 93)
        return base1, min(qual1, qual2)
    if base1 == "N":
        return in2
    if base2 == "N":
        return in1
    hi, lo = (in1, in2) if qual1 >= qual2 else (in2, in1)
    if hi[1] - lo[1] >= cfg.min_quality_difference:
        return hi[0], min(hi[1], cfg.max_consensus_quality)
    if cfg.conflict_policy == "strict":
        return STRICT
    return "N", UNRESOLVED_QUALITY


def _merge_insertion(
    boundary: int,
    struct1: ReadStructure,
    struct2: ReadStructure,
) -> tuple[str, tuple[int, ...]] | MergeFailure | None:
    """Resolve the insertions both mates attach to ``boundary``.

    Both mates inserting: identical sequences merge once (per-base quality
    maximum), anything else is a conflict. A single-mate insertion is
    authoritative when the other mate does not span the boundary's interior;
    when the other mate aligns across it, the interpretations conflict.
    """
    ins1 = struct1.insertions.get(boundary)
    ins2 = struct2.insertions.get(boundary)
    if ins1 is None and ins2 is None:
        return None
    if ins1 is not None and ins2 is not None:
        if ins1[0] == ins2[0]:
            quals = tuple(max(a, b) for a, b in zip(ins1[1], ins2[1], strict=True))
            return ins1[0], quals
        return MergeFailure(REASON_INSERTION)
    owner, other, ins = (struct1, struct2, ins1) if ins1 is not None else (struct2, struct1, ins2)
    if covers(other, boundary) and other.start0 < boundary < other.end0:
        return MergeFailure(REASON_INSERTION)
    return ins
