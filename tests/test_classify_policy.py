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

    def test_secondary_only_template_passes_through(self, fx) -> None:
        """A template with no primary records passes through unchanged.

        Possible in filtered BAMs; the documented non-primary pass-through
        policy applies — no adjustment, no tagging, classified
        neither_aligned (no primary alignment evidence exists).
        """
        ref = fx.ref_seqs["chr1"]
        s1 = fx.make_record(
            "sonly",
            0x1 | 0x40 | 0x100,
            UNIQUE,
            "20M",
            ref[UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        s2 = fx.make_record(
            "sonly",
            0x1 | 0x80 | 0x100,
            UNIQUE,
            "20M",
            ref[UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([s1, s2])
        code, out = fx.run(mode="adjust")
        assert code == 0
        recs = fx.read_bam(out["adjusted_bam"])
        assert len(recs) == 2
        for rec in recs:
            assert rec.reference_start == UNIQUE  # never adjusted
            assert not rec.has_tag("ZC")  # never tagged
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["neither_aligned"] == 1
        assert by_type["total"] == 1
        # Merge mode: nothing merged; everything retained in the unmerged stream.
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert len(unmerged) == 2
        assert all(not r.has_tag("ZC") for r in unmerged)

    def test_orphan_stale_mq_removed(self, fx) -> None:
        """Orphan MQ semantics: removed on the mapped record, refreshed on
        the unmapped mate."""
        import pysam

        r1 = fx.make_record(
            "omq",
            0x1 | 0x40,
            UNIQUE,
            "20M",
            fx.ref_seqs["chr1"][UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE,
            tags={"MQ": 60},
        )
        hdr = pysam.AlignmentHeader.from_dict(fx.header_dict())
        r2 = pysam.AlignedSegment(hdr)
        r2.query_name = "omq"
        r2.flag = 0x1 | 0x80 | 0x4
        r2.reference_id = 0
        r2.reference_start = UNIQUE
        r2.next_reference_id = 0
        r2.next_reference_start = UNIQUE
        r2.mapping_quality = 0
        r2.query_sequence = "A" * 20
        r2.query_qualities = pysam.qualitystring_to_array("I" * 20)
        r2.set_tag("MQ", 60, value_type="i")
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="adjust")
        assert code == 0
        recs = {r.flag & 0xC0: r for r in fx.read_bam(out["adjusted_bam"])}
        assert "MQ" not in dict(recs[0x40].get_tags())  # mate is unmapped
        assert dict(recs[0x80].get_tags())["MQ"] == 60  # mapped mate's MAPQ

    def test_mixed_paired_unpaired_rejected_in_any_order(self, fx) -> None:
        """Record order must not decide whether the template is rejected."""
        ref = fx.ref_seqs["chr1"]
        unpaired = fx.make_record("mix", 0, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        paired_r1 = fx.make_record(
            "mix",
            0x1 | 0x40,
            UNIQUE + 100,
            "20M",
            ref[UNIQUE + 100 : UNIQUE + 120],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        for order in ([unpaired, paired_r1], [paired_r1, unpaired]):
            fx.write_bam_and_set(list(order))
            code, out = fx.run(mode="adjust")
            assert code == 1, order
            assert not out["adjusted_bam"].exists()

    def test_two_unpaired_primaries_rejected(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        a = fx.make_record("two", 0, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 20])
        b = fx.make_record("two", 0, UNIQUE + 100, "20M", ref[UNIQUE + 100 : UNIQUE + 120])
        fx.write_bam_and_set([a, b])
        code, out = fx.run(mode="adjust")
        assert code == 1
        assert not out["adjusted_bam"].exists()

    def test_duplicate_primary_read_numbers_fail_loudly(self, fx) -> None:
        """Two primary R1 records must not silently overwrite each other."""
        ref = fx.ref_seqs["chr1"]
        ra = fx.make_record(
            "dup",
            0x1 | 0x40,
            UNIQUE,
            "20M",
            ref[UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE + 100,
        )
        rb = fx.make_record(
            "dup",
            0x1 | 0x40,
            UNIQUE + 100,
            "20M",
            ref[UNIQUE + 100 : UNIQUE + 120],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([ra, rb])
        code, out = fx.run(mode="adjust")
        assert code == 1
        assert not out["adjusted_bam"].exists()

    def test_unmapped_mate_rnext_follows_movement(self, fx) -> None:
        """After the mapped mate slides, RNEXT/PNEXT update on BOTH records.

        The unmapped mate sits at the mapped read's locus; after the mapped
        mate moves 404 -> 400 (leftmost placement in the tandem repeat), the
        unmapped record's own next_reference_start must follow to 400 too —
        not stay at the stale pre-adjustment position.
        """
        import pysam

        from conftest import LONG_REPEAT

        ref = fx.ref_seqs["chr1"]
        r1 = fx.make_record(
            "half",
            0x1 | 0x40,
            LONG_REPEAT + 4,
            "16M",
            ref[LONG_REPEAT + 4 : LONG_REPEAT + 20],
            mate_contig="chr1",
            mate_pos=LONG_REPEAT + 4,
        )
        hdr = pysam.AlignmentHeader.from_dict(fx.header_dict())
        r2 = pysam.AlignedSegment(hdr)
        r2.query_name = "half"
        r2.flag = 0x1 | 0x80 | 0x4
        r2.reference_id = 0
        r2.reference_start = LONG_REPEAT + 4  # unmapped mate at mapped locus
        r2.next_reference_id = 0
        r2.next_reference_start = LONG_REPEAT + 4  # input points at OLD pos
        r2.mapping_quality = 0
        r2.query_sequence = "A" * 16
        r2.query_qualities = pysam.qualitystring_to_array("I" * 16)
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="adjust")
        assert code == 0
        recs = {r.flag & 0xC0: r for r in fx.read_bam(out["adjusted_bam"])}
        mapped, unmapped = recs[0x40], recs[0x80]
        assert mapped.reference_start == LONG_REPEAT  # slid to leftmost
        assert unmapped.reference_start == LONG_REPEAT  # SAM convention locus
        for rec in (mapped, unmapped):
            assert rec.next_reference_id == 0
            assert rec.next_reference_start == LONG_REPEAT
