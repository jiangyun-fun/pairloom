"""Unit tests for joint pair optimization (scenarios 2, 3).

select_pair is exercised directly with hand-built placements so that score,
overlap, orientation, fragment, and tie-break behaviour are pinned down
independently of discovery.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pairloom.alignment import Placement
from pairloom.candidates import CandidateSet
from pairloom.config import PairLoomConfig
from pairloom.exceptions import ConfigurationError
from pairloom.pairing import orientation_ok, overlap_len, select_pair


def placement(pos: int, score: int, exact: bool = True, span: int = 20) -> Placement:
    return Placement(tid=0, pos0=pos, cigar=((0, span),), score=score, exact=exact)


def cfg(**overrides) -> PairLoomConfig:
    base = dict(
        input_path=Path("in.bam"),
        reference_path=Path("ref.fa"),
        mode="adjust",
        adjusted_bam=Path("adj.bam"),
        disagreement_parquet=Path("d.parquet"),
        read_stats_parquet=Path("s.parquet"),
    )
    base.update(overrides)
    return PairLoomConfig(**base)


class TestOrientation:
    @pytest.mark.parametrize(
        ("rev1", "rev2", "left_r1", "model", "expected"),
        [
            # FR: leftmost forward, rightmost reverse
            (False, True, True, "fr", True),
            (True, False, True, "fr", False),
            (False, False, True, "fr", False),
            (True, True, True, "fr", False),
            (True, False, False, "fr", True),  # R2 leftmost: R2 fwd, R1 rev
            (False, True, False, "fr", False),
            # RF: leftmost reverse, rightmost forward
            (True, False, True, "rf", True),
            # FF / RR
            (False, False, True, "ff", True),
            (True, True, False, "rr", True),
            # any
            (True, True, True, "any", True),
        ],
    )
    def test_models(self, rev1, rev2, left_r1, model, expected) -> None:
        assert orientation_ok(rev1, rev2, left_r1, model) is expected


class TestSelectPair:
    def test_score_beats_overlap_within_tolerance_zero(self) -> None:
        """Scenario 3: alignment evidence outranks overlap preference.

        (80,115) scores 200 with zero overlap; (100,115) scores 190 with 5
        overlapping bases. Tolerance 0 discards the weaker pair regardless
        of its larger overlap.
        """
        c1 = CandidateSet((placement(80, 100), placement(100, 90)), False)
        c2 = CandidateSet((placement(115, 100),), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg())
        assert choice is not None
        assert choice.p1.pos0 == 80
        assert choice.p2.pos0 == 115

    def test_near_equivalent_set_prefers_overlap(self) -> None:
        """Scenario 2: inside the tolerance window, larger overlap wins."""
        c1 = CandidateSet((placement(80, 100), placement(100, 90)), False)
        c2 = CandidateSet((placement(115, 100),), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg(alignment_score_tolerance=10))
        assert choice is not None
        assert choice.p1.pos0 == 100  # overlap outranks the 10-point score gap

    def test_max_overlap_strategy_ignores_tolerance_filter(self) -> None:
        """max-overlap prefers overlap among all constraint-satisfying pairs."""
        c1 = CandidateSet((placement(80, 100), placement(100, 90)), False)
        c2 = CandidateSet((placement(115, 100),), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg(pair_strategy="max-overlap"))
        assert choice is not None
        assert choice.p1.pos0 == 100

    def test_orientation_filters_invalid_combinations(self) -> None:
        c1 = CandidateSet((placement(100, 100),), False)
        c2 = CandidateSet((placement(80, 100),), False)
        # R2 leftmost and reverse violates FR -> no valid pair.
        assert select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg()) is None

    def test_fragment_bounds_filter(self) -> None:
        c1 = CandidateSet((placement(100, 100), placement(400, 100)), False)
        c2 = CandidateSet((placement(500, 100),), False)
        # (100,500) span 420; (400,500) span 120.
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg(min_fragment_length=250))
        assert choice is not None
        assert choice.p1.pos0 == 100  # only the long fragment passes

    def test_tiebreak_leftmost_start_then_rightmost_end(self) -> None:
        c1 = CandidateSet((placement(100, 50), placement(200, 50)), False)
        c2 = CandidateSet((placement(110, 50), placement(210, 50)), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg())
        assert choice is not None
        assert (choice.p1.pos0, choice.p2.pos0) == (100, 110)

    def test_n_equiv_counts_distinct_near_placements(self) -> None:
        c1 = CandidateSet((placement(100, 50), placement(104, 50)), False)
        c2 = CandidateSet((placement(110, 50),), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg())
        assert choice is not None
        assert choice.n_equiv1 == 2
        assert choice.n_equiv2 == 1

    def test_max_pair_evaluations_caps_work(self) -> None:
        c1 = CandidateSet(tuple(placement(100 + 4 * i, 50) for i in range(50)), False)
        c2 = CandidateSet(tuple(placement(400 + 4 * i, 50) for i in range(50)), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg(max_pair_evaluations=10))
        # Fragment bounds default 0..1000 keep most of the first 10 valid.
        assert choice is not None

    def test_eval_budget_truncation_flags_choice(self) -> None:
        """Hitting --max-pair-evaluations must mark the search truncated.

        Unexplored combinations could hold equivalent placements, so the
        MAPQ policy must not claim uniqueness (min(orig, 10) applies).
        """
        c1 = CandidateSet((placement(100, 50), placement(104, 50)), False)
        c2 = CandidateSet((placement(110, 50),), False)
        full = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg())
        assert full is not None
        assert full.truncated is False
        capped = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg(max_pair_evaluations=1))
        assert capped is not None
        assert capped.truncated is True

    def test_discovery_truncation_propagates(self) -> None:
        c1 = CandidateSet((placement(100, 50),), True)
        c2 = CandidateSet((placement(110, 50),), False)
        choice = select_pair(c1, c2, rev1=False, rev2=True, cfg=cfg())
        assert choice is not None
        assert choice.truncated is True

    def test_overlap_len(self) -> None:
        a = placement(100, 50, span=20)
        b = placement(110, 50, span=20)
        assert overlap_len(a, b) == 10
        c = placement(200, 50, span=20)
        assert overlap_len(a, c) == 0


class TestConfigValidation:
    def test_bad_tag_rejected(self) -> None:
        with pytest.raises(ConfigurationError, match="status_tag"):
            cfg(status_tag="ZC1")

    def test_bad_orientation_rejected(self) -> None:
        with pytest.raises(ConfigurationError, match="library_orientation"):
            cfg(library_orientation="xx")  # type: ignore[arg-type]

    def test_fragment_bounds_consistency(self) -> None:
        with pytest.raises(ConfigurationError, match="fragment length bounds"):
            cfg(min_fragment_length=100, max_fragment_length=50)

    def test_mode_output_wiring(self) -> None:
        with pytest.raises(ConfigurationError, match="--mode merge requires"):
            cfg(mode="merge", adjusted_bam=None)  # no merged/unmerged outputs
        with pytest.raises(ConfigurationError, match="only valid with --mode adjust"):
            cfg(
                mode="merge",
                merged_bam=Path("m.bam"),
                unmerged_bam=Path("u.bam"),
                adjusted_bam=Path("a.bam"),
            )
