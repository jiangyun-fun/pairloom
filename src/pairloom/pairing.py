"""Joint R1/R2 placement optimization.

Candidate pairs are evaluated against the documented preference hierarchy:

1. both mates on the same reference (candidates always are; cross-contig
   pairs therefore have no valid combination and keep their BAM placements);
2. library orientation and fragment-length plausibility;
3. best combined alignment score: candidate pairs scoring below
   ``best - --alignment-score-tolerance`` are discarded as clearly inferior;
4. largest aligned reference overlap (the tolerance defines which pairs
   count as near-equivalent, so overlap never overrides clearly superior
   evidence);
5. deterministic tie-breaks: leftmost fragment start, rightmost fragment
   end, then R1/R2 start coordinates.

``pair-aware`` (default) applies the score-tolerance filter first and then
prefers overlap inside the surviving set. ``max-overlap`` skips the filter:
among all constraint-satisfying pairs it prefers overlap first, which is
more aggressive about maximizing overlap. Both strategies are deterministic.
"""

from __future__ import annotations

from dataclasses import dataclass

from pairloom.alignment import Placement
from pairloom.candidates import CandidateSet
from pairloom.config import PairLoomConfig


@dataclass(frozen=True, slots=True)
class PairChoice:
    p1: Placement
    p2: Placement
    n_equiv1: int
    n_equiv2: int
    truncated: bool


def orientation_ok(rev1: bool, rev2: bool, left_is_read1: bool, model: str) -> bool:
    """Check mate orientation against the library model.

    ``left_is_read1`` refers to reference-coordinate order, not read number.
    """
    left_rev = rev1 if left_is_read1 else rev2
    right_rev = rev2 if left_is_read1 else rev1
    if model == "fr":
        return (not left_rev) and right_rev
    if model == "rf":
        return left_rev and (not right_rev)
    if model == "ff":
        return (not rev1) and (not rev2)
    if model == "rr":
        return rev1 and rev2
    if model == "any":
        return True
    raise ValueError(f"unknown library orientation model {model!r}")


def overlap_len(p1: Placement, p2: Placement) -> int:
    return max(0, min(p1.end0, p2.end0) - max(p1.pos0, p2.pos0))


def fragment_ok(p1: Placement, p2: Placement, cfg: PairLoomConfig) -> bool:
    span = max(p1.end0, p2.end0) - min(p1.pos0, p2.pos0)
    return cfg.min_fragment_length <= span <= cfg.max_fragment_length


def select_pair(
    cands1: CandidateSet,
    cands2: CandidateSet,
    rev1: bool,
    rev2: bool,
    cfg: PairLoomConfig,
) -> PairChoice | None:
    """Choose the joint placement, or None when no combination satisfies constraints."""
    list1 = sorted(cands1.placements, key=lambda p: (-p.score, p.pos0, p.cigar))
    list2 = sorted(cands2.placements, key=lambda p: (-p.score, p.pos0, p.cigar))
    evaluated = 0
    valid: list[tuple[Placement, Placement, int, int, int, int]] = []
    for p1 in list1:
        for p2 in list2:
            if evaluated >= cfg.max_pair_evaluations:
                break
            if p1.tid != p2.tid:
                continue
            evaluated += 1
            left_is_read1 = p1.pos0 <= p2.pos0
            if not orientation_ok(rev1, rev2, left_is_read1, cfg.library_orientation):
                continue
            if not fragment_ok(p1, p2, cfg):
                continue
            frag_start = min(p1.pos0, p2.pos0)
            frag_end = max(p1.end0, p2.end0)
            valid.append((p1, p2, p1.score + p2.score, overlap_len(p1, p2), frag_start, frag_end))
        if evaluated >= cfg.max_pair_evaluations:
            break
    if not valid:
        return None
    best_score = max(entry[2] for entry in valid)
    tol = cfg.alignment_score_tolerance
    near = [entry for entry in valid if entry[2] >= best_score - tol]
    # Both strategies rank overlap first, then residual score, then the
    # documented coordinate tie-breaks; they differ only in the candidate
    # pool (score-tolerance filtered vs. all constraint-satisfying pairs).
    pool = near if cfg.pair_strategy == "pair-aware" else valid

    def key(e: tuple[Placement, Placement, int, int, int, int]) -> tuple:
        return (-e[3], -e[2], e[4], -e[5], e[0].pos0, e[1].pos0)

    chosen = min(pool, key=key)
    n_equiv1 = len({e[0].pos0 for e in near})
    n_equiv2 = len({e[1].pos0 for e in near})
    return PairChoice(
        p1=chosen[0],
        p2=chosen[1],
        n_equiv1=n_equiv1,
        n_equiv2=n_equiv2,
        truncated=cands1.truncated or cands2.truncated,
    )


def single_read_ambiguity(cands: CandidateSet) -> int:
    """Number of equivalent placements for an independently assessed read."""
    return len(cands.placements)
