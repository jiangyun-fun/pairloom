"""Reference-coordinate overlap extraction for mate-pair base comparison.

Two bases overlap when they align to the same reference coordinate (CIGAR-
aware mapping; soft clips, insertions, deletions, and skips never produce
observations). Only A/C/G/T observations passing the configured base- and
mapping-quality filters are eligible. Strand is taken from the actual
alignment orientation, not the read number. Compared bases are expressed in
reference orientation, so a reverse-strand mate's stored base is
complemented before comparison.
"""

from __future__ import annotations

from dataclasses import dataclass

from pairloom.alignment import ACGT, ReadStructure
from pairloom.config import PairLoomConfig


@dataclass(frozen=True, slots=True)
class OverlapObservation:
    """One eligible read-base observation inside a mate overlap."""

    tid: int
    rpos0: int
    ref_base: str
    base: str
    strand_reverse: bool
    pair_disagrees: bool


@dataclass(frozen=True, slots=True)
class OverlapResult:
    observations: tuple[OverlapObservation, ...]
    #: number of reference positions covered by both mates' aligned blocks
    #: (quality-independent; used for the --min-overlap merge criterion)
    aligned_overlap_len: int


def analyze_overlap(
    struct1: ReadStructure,
    struct2: ReadStructure,
    ref_seq: str,
    ref_offset: int,
    cfg: PairLoomConfig,
) -> OverlapResult:
    """Join two mates' observations by reference position.

    ``ref_seq`` is the slice ``ref_seq[i] == reference[ref_offset + i]``.
    """
    common = struct1.bases.keys() & struct2.bases.keys()
    overlap_len = len(common)
    if not common:
        return OverlapResult((), 0)
    if struct1.mapq < cfg.min_mapq or struct2.mapq < cfg.min_mapq:
        return OverlapResult((), overlap_len)
    min_bq = cfg.min_base_quality
    obs: list[OverlapObservation] = []
    tid = struct1.tid
    for rpos in sorted(common):
        base1, qual1 = struct1.bases[rpos]
        base2, qual2 = struct2.bases[rpos]
        ok1 = base1 in ACGT and qual1 >= min_bq
        ok2 = base2 in ACGT and qual2 >= min_bq
        if not (ok1 and ok2):
            continue
        ref_base = ref_seq[rpos - ref_offset]
        disagrees = base1 != base2
        obs.append(OverlapObservation(tid, rpos, ref_base, base1, struct1.is_reverse, disagrees))
        obs.append(OverlapObservation(tid, rpos, ref_base, base2, struct2.is_reverse, disagrees))
    return OverlapResult(tuple(obs), overlap_len)


def aggregate_overlap(
    struct1: ReadStructure,
    struct2: ReadStructure,
    ref_seq: str,
    ref_offset: int,
    cfg: PairLoomConfig,
) -> dict[int, tuple[str, int, bool, dict[str, list[int]]]]:
    """Aggregate eligible overlap observations per position (batched path).

    Produces exactly the aggregates the per-observation path
    (``analyze_overlap`` feeding ``DisagreementAggregator.add_observation``)
    would produce: ``{rpos0: (ref_base, dp, has_disag, alts)}`` where ``dp``
    counts eligible observations (two per eligible base pair), ``alts`` maps
    each non-reference base to ``[forward, reverse]`` counts, and
    ``has_disag`` marks a base-pair disagreement. Gating (mapq, base
    quality, ACGT) is identical to ``analyze_overlap``.
    """
    common = struct1.bases.keys() & struct2.bases.keys()
    if not common:
        return {}
    if struct1.mapq < cfg.min_mapq or struct2.mapq < cfg.min_mapq:
        return {}
    min_bq = cfg.min_base_quality
    rev1 = struct1.is_reverse
    rev2 = struct2.is_reverse
    out: dict[int, tuple[str, int, bool, dict[str, list[int]]]] = {}
    for rpos in sorted(common):
        base1, qual1 = struct1.bases[rpos]
        base2, qual2 = struct2.bases[rpos]
        ok1 = base1 in ACGT and qual1 >= min_bq
        ok2 = base2 in ACGT and qual2 >= min_bq
        if not (ok1 and ok2):
            continue
        ref_base = ref_seq[rpos - ref_offset]
        disagrees = base1 != base2
        alts: dict[str, list[int]] = {}
        if base1 != ref_base:
            alts[base1] = [0, 1] if rev1 else [1, 0]
        if base2 != ref_base:
            counts = alts.get(base2)
            if counts is None:
                alts[base2] = [0, 1] if rev2 else [1, 0]
            else:
                counts[1 if rev2 else 0] += 1
        out[rpos] = (ref_base, 2, disagrees, alts)
    return out
