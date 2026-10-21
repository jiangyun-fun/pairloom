"""Classification, tagging, and record policy (scenarios 16-19, 22)."""

from __future__ import annotations

from conftest import UNIQUE


class TestClassification:
    def test_all_categories_mutually_exclusive_and_summing(self, fx) -> None:
        """Scenario 22: one template per category, total equals the sum."""
        records: list = []
        # mergeable
        records += fx.pair("a_merge", UNIQUE, UNIQUE + 5, read_len=20)
        # both_aligned: mates do not overlap
        records += fx.pair("b_both", UNIQUE, UNIQUE + 200, read_len=20)
        # left_aligned_only: R1 mapped, R2 unmapped
        r1 = fx.make_record(
            "c_left",
            0x1 | 0x40 | 0x8 | 0x2,
            UNIQUE,
            "20M",
            fx.ref_seqs["chr1"][UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        r2 = fx.make_record("c_left", 0x1 | 0x80 | 0x4, UNIQUE, "*", "ACGTACGTACGTACGTACGT")
        records += [r1, r2]
        # right_aligned_only: R1 unmapped, R2 mapped
        r1 = fx.make_record(
            "d_right",
            0x1 | 0x40 | 0x4,
            UNIQUE,
            "*",
            "ACGTACGTACGTACGTACGT",
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        r2 = fx.make_record(
            "d_right",
            0x1 | 0x80 | 0x8 | 0x2,
            UNIQUE,
            "20M",
            fx.ref_seqs["chr1"][UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        records += [r1, r2]
        # neither_aligned
        records += [
            fx.make_record("e_none", 0x1 | 0x40 | 0x4, -1, "*", "ACGTACGTAC"),
            fx.make_record("e_none", 0x1 | 0x80 | 0x4, -1, "*", "ACGTACGTAC"),
        ]
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type == {
            "total": 5,
            "mergeable": 1,
            "left_aligned_only": 1,
            "right_aligned_only": 1,
            "both_aligned": 1,
            "neither_aligned": 1,
        }
        parts = by_type["total"]
        assert parts == sum(
            by_type[k]
            for k in (
                "mergeable",
                "left_aligned_only",
                "right_aligned_only",
                "both_aligned",
                "neither_aligned",
            )
        )

    def test_tags_on_both_mates_and_merged_record(self, fx) -> None:
        records: list = []
        records += fx.pair("m1", UNIQUE, UNIQUE + 5, read_len=20)
        records += fx.pair("m2", UNIQUE, UNIQUE + 200, read_len=20)
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [dict(r.get_tags())["ZC"] for r in merged] == ["mergeable"]
        m2_tags = [dict(r.get_tags())["ZC"] for r in unmerged if r.query_name == "m2"]
        assert m2_tags == ["both_aligned", "both_aligned"]

    def test_adjust_mode_tags_without_merging(self, fx) -> None:
        pair = fx.pair("am", UNIQUE, UNIQUE + 5, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        assert len(records) == 2  # separate mates retained
        assert all(dict(r.get_tags())["ZC"] == "mergeable" for r in records)


class TestTagCollision:
    def test_unrelated_existing_tag_fails_loudly(self, fx) -> None:
        pair = fx.pair("col", UNIQUE, UNIQUE + 5, read_len=20, tags={"ZC": "unrelated-value"})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 1  # loud failure, not a silent overwrite

    def test_rerun_idempotent_classification_tag(self, fx) -> None:
        pair = fx.pair("re", UNIQUE, UNIQUE + 5, read_len=20, tags={"ZC": "mergeable"})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        assert all(dict(r.get_tags())["ZC"] == "mergeable" for r in records)

    def test_custom_status_tag(self, fx) -> None:
        pair = fx.pair("ct", UNIQUE, UNIQUE + 200, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust", extra=["--status-tag", "XY"])
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        assert all(dict(r.get_tags())["XY"] == "both_aligned" for r in records)


class TestRecordPolicy:
    def test_secondary_and_supplementary_passthrough(self, fx) -> None:
        """Scenario 19: non-primary records pass through untouched."""
        ref = fx.ref_seqs["chr1"]
        pair = fx.pair("sp", UNIQUE, UNIQUE + 5, read_len=20)
        secondary = fx.make_record(
            "sp", 0x1 | 0x40 | 0x100, UNIQUE + 2, "20M", ref[UNIQUE + 2 : UNIQUE + 22]
        )
        supplementary = fx.make_record(
            "sp", 0x1 | 0x40 | 0x800, UNIQUE + 2, "15M", ref[UNIQUE + 2 : UNIQUE + 17]
        )
        fx.write_bam_and_set([*pair, secondary, supplementary])
        code, out = fx.run(mode="adjust")
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        assert len(records) == 4
        non_primary = [r for r in records if r.flag & 0x900]
        assert len(non_primary) == 2
        for rec in non_primary:
            assert not rec.has_tag("ZC")  # never tagged
            assert rec.reference_start == UNIQUE + 2  # never adjusted
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        # Non-primary records block merging by policy.
        assert by_type["both_aligned"] == 1

    def test_non_primary_blocks_merging_in_merge_mode(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        pair = fx.pair("sp2", UNIQUE, UNIQUE + 5, read_len=20)
        supplementary = fx.make_record(
            "sp2", 0x1 | 0x40 | 0x800, UNIQUE + 2, "15M", ref[UNIQUE + 2 : UNIQUE + 17]
        )
        fx.write_bam_and_set([*pair, supplementary])
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        names = [r.query_name for r in fx.read_bam(out["unmerged_bam"])]
        assert names == ["sp2", "sp2", "sp2"]  # nothing silently discarded

    def test_mates_on_different_contigs(self, fx) -> None:
        """Scenario 16: cross-contig pair -> both_aligned, no overlap rows."""
        ref = fx.ref_seqs
        r1 = fx.make_record(
            "cc",
            0x1 | 0x40 | 0x2,
            UNIQUE,
            "20M",
            ref["chr1"][UNIQUE : UNIQUE + 20],
            contig="chr1",
            mate_contig="chr2",
            mate_pos=50,
        )
        r2 = fx.make_record(
            "cc",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            50,
            "20M",
            ref["chr2"][50:70],
            contig="chr2",
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["both_aligned"] == 1
        assert fx.read_parquet(out["disagreement_parquet"])["pos"] == []
        assert fx.read_bam(out["merged_bam"]) == []
        # Both adjusted records present with proper-pair flag cleared.
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert len(unmerged) == 2
        assert all(not (r.flag & 0x2) for r in unmerged)

    def test_more_than_two_primaries_fails_loudly(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        trio = fx.pair("t3", UNIQUE, UNIQUE + 5, read_len=20)
        extra = fx.make_record("t3", 0x1 | 0x40, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        fx.write_bam_and_set([*trio, extra])
        code, out = fx.run(mode="adjust")
        assert code == 1

    def test_orphan_single_record(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        orphan = fx.make_record("orphan", 0x1 | 0x40, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        fx.write_bam_and_set([orphan])
        code, out = fx.run(mode="adjust")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["left_aligned_only"] == 1
        assert by_type["total"] == 1

    def test_unpaired_record_treated_as_r1(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        unpaired = fx.make_record("single", 0, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        fx.write_bam_and_set([unpaired])
        code, out = fx.run(mode="adjust")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["left_aligned_only"] == 1
