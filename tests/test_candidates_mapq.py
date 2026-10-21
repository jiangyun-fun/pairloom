"""Repeat-aware adjustment and MAPQ policy (scenarios 1, 4, 5, 21).

The long tandem repeat at chr1:400-463 (period 4) gives every read fully
inside it multiple exact placements; the unique region at chr1:500 gives
single placements. Expected placements are derived in-test from the fixture
reference strings.
"""

from __future__ import annotations

from conftest import LONG_REPEAT, REPEAT2, UNIQUE
from pairloom.mapq import adjusted_mapq


def alt_base(ref_base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[ref_base]


class TestMapqPolicyUnit:
    def test_ambiguous_capped_at_3(self) -> None:
        assert (
            adjusted_mapq(60, n_equiv=3, truncated=False, exact=True, placement_changed=True) == 3
        )

    def test_truncated_capped_at_10(self) -> None:
        assert (
            adjusted_mapq(60, n_equiv=1, truncated=True, exact=True, placement_changed=True) == 10
        )

    def test_unchanged_preserved(self) -> None:
        assert (
            adjusted_mapq(42, n_equiv=1, truncated=False, exact=True, placement_changed=False) == 42
        )

    def test_exact_corroborated_change_gets_60(self) -> None:
        assert (
            adjusted_mapq(60, n_equiv=1, truncated=False, exact=True, placement_changed=True) == 60
        )

    def test_exact_without_corroboration_capped(self) -> None:
        assert (
            adjusted_mapq(20, n_equiv=1, truncated=False, exact=True, placement_changed=True) == 20
        )
        assert (
            adjusted_mapq(60, n_equiv=1, truncated=False, exact=False, placement_changed=True) == 25
        )

    def test_force_mapq60_exact(self) -> None:
        assert (
            adjusted_mapq(
                5,
                n_equiv=9,
                truncated=True,
                exact=True,
                placement_changed=False,
                force_mapq60_exact=True,
            )
            == 60
        )


class TestSingleReadLeftmost:
    def test_orphan_in_repeat_moves_leftmost(self, fx) -> None:
        """Scenario 1: equivalent placements prefer the leftmost position."""
        ref = fx.ref_seqs["chr1"]
        block = "AACCAACCAACC"
        read = fx.make_record("orphan1", 0x1 | 0x40, REPEAT2, "12M", block, mapq=60)
        fx.write_bam_and_set([read])
        code, out = fx.run(mode="adjust")
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        assert len(records) == 1
        rec = records[0]
        expected_leftmost = ref.find(block)
        assert 0 <= expected_leftmost < REPEAT2  # the first repeat copy
        assert rec.reference_start == expected_leftmost
        # Multiple equivalent placements: conservative MAPQ, never 60.
        assert rec.mapping_quality <= 3
        tags = dict(rec.get_tags())
        assert tags["NM"] == 0
        assert tags["MD"] == "12"
        assert tags["ZC"] == "left_aligned_only"

    def test_deterministic_across_runs(self, fx) -> None:
        block = "AACCAACCAACC"
        results = []
        for _ in range(3):
            read = fx.make_record("orphan1", 0x1 | 0x40, REPEAT2, "12M", block)
            fx.write_bam_and_set([read])
            code, out = fx.run(mode="adjust")
            assert code == 0
            rec = fx.read_bam(out["adjusted_bam"])[0]
            results.append((rec.reference_start, rec.mapping_quality))
        assert len(set(results)) == 1


class TestAmbiguousRepeatPair:
    def test_pair_in_long_repeat_conservative_mapq(self, fx) -> None:
        """Scenario 4: ambiguous exact repeats never earn MAPQ 60."""
        pair = fx.pair("rep", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        for rec in fx.read_bam(out["adjusted_bam"]):
            assert rec.mapping_quality <= 3
            assert rec.mapping_quality != 60

    def test_capped_pair_search_caps_mapq(self, fx) -> None:
        """A truncated pair search must not imply uniqueness: no false 60.

        With the default budget this repeat pair earns MAPQ 3 (many
        equivalent placements); capping the search after one evaluation
        hides those placements, so the truncation cap (min(orig, 10))
        must apply instead.
        """
        pair = fx.pair("rep", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16)
        fx.write_bam_and_set(pair)
        code, out = fx.run(
            mode="adjust",
            extra=["--max-pair-evaluations", "1", "--max-candidates-per-read", "100"],
        )
        assert code == 0
        for rec in fx.read_bam(out["adjusted_bam"]):
            assert rec.mapping_quality <= 10


class TestMateMapqTag:
    def test_mq_refreshed_when_only_mapq_changes(self, fx) -> None:
        """Repeat-ambiguous unchanged placements: MQ follows the adjusted MAPQ.

        The placement does not move, so the old code never entered the
        stale-tag removal branch and left MQ:i:60 next to MAPQ 3.
        """
        pair = fx.pair("mq", LONG_REPEAT, LONG_REPEAT, read_len=16, tags={"MQ": 60})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        for rec in fx.read_bam(out["adjusted_bam"]):
            assert rec.reference_start == LONG_REPEAT  # unchanged placement
            assert rec.mapping_quality == 3
            assert dict(rec.get_tags())["MQ"] == 3

    def test_mq_follows_movement(self, fx) -> None:
        """After joint movement, MQ reflects the mate's adjusted MAPQ."""
        pair = fx.pair("mq2", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16, tags={"MQ": 60})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        for rec in fx.read_bam(out["adjusted_bam"]):
            assert rec.reference_start == LONG_REPEAT
            assert rec.mapping_quality == 3
            assert dict(rec.get_tags())["MQ"] == 3


class TestUniqueExact:
    def test_unique_read_preserves_mapq(self, fx) -> None:
        """Unique placement, unchanged: original MAPQ survives."""
        ref = fx.ref_seqs["chr1"]
        read = fx.make_record("uniq", 0x1 | 0x40, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        fx.write_bam_and_set([read])
        code, out = fx.run(mode="adjust")
        assert code == 0
        rec = fx.read_bam(out["adjusted_bam"])[0]
        assert rec.reference_start == UNIQUE
        assert rec.mapping_quality == 60

    def test_corroborated_move_earns_mapq60(self, fx) -> None:
        """Scenario 5: a moved exact read with unique pair evidence earns 60.

        Mini reference: R1 block "AAAAGGGG" occurs at 4 and 34; R2 block
        "CCCCGGGG" occurs only at 44. With --max-fragment-length 30 the only
        valid FR combination is (34, 44), so R1 moves 4 -> 34 with exactly
        one equivalent placement and full corroboration.
        """
        ref_seqs = {
            "mini": "TTTT"
            + "AAAAGGGG"
            + "TTTTTT"
            + "CCGTACGT"
            + "ACGTACGT"
            + "AAAAGGGG"
            + "TT"
            + "CCCCGGGG"
            + "TTTT"
        }
        fx.ref_seqs = ref_seqs
        fx.__post_init__()  # rewrite reference files
        r1 = fx.make_record(
            "cor",
            0x1 | 0x40 | 0x2,
            4,
            "8M",
            "AAAAGGGG",
            contig="mini",
            mate_contig="mini",
            mate_pos=44,
        )
        r2 = fx.make_record(
            "cor",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            44,
            "8M",
            "CCCCGGGG",
            contig="mini",
            mate_contig="mini",
            mate_pos=4,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="adjust", extra=["--max-fragment-length", "30"])
        assert code == 0
        by_flag = {r.flag: r for r in fx.read_bam(out["adjusted_bam"])}
        r1_out = by_flag[0x1 | 0x40 | 0x2]
        r2_out = by_flag[0x1 | 0x80 | 0x10 | 0x20 | 0x2]
        assert r1_out.reference_start == 34  # moved by pair evidence
        assert r1_out.mapping_quality == 60  # exact + unique + corroborated
        assert r2_out.reference_start == 44  # single candidate, unchanged
        assert r2_out.mapping_quality == 60
