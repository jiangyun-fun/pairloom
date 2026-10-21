"""Consensus merging (scenarios 2, 8.x, 10, 12, 13, 14, 15)."""

from __future__ import annotations

from conftest import LONG_REPEAT, UNIQUE


def alt_bases(ref_base: str, n: int) -> list[str]:
    choices = [b for b in "ACGT" if b != ref_base]
    return choices[:n]


class TestMergeOnce:
    def test_spec_example_overlap_appears_once(self, fx) -> None:
        """Spec 8.1: R1 AACCGGTT + R2 GGTTAA -> AACCGGTTAA."""
        fx.ref_seqs = {"mini": "AACCGGTTAA" + "CCGGTTACGAC" * 3}
        fx.__post_init__()
        r1 = fx.make_record(
            "ex",
            0x1 | 0x40 | 0x2,
            0,
            "8M",
            "AACCGGTT",
            contig="mini",
            mate_contig="mini",
            mate_pos=4,
        )
        r2 = fx.make_record(
            "ex",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            4,
            "6M",
            "GGTTAA",
            contig="mini",
            mate_contig="mini",
            mate_pos=0,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge", extra=["--min-overlap", "4"])
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        assert len(merged) == 1
        rec = merged[0]
        assert rec.query_sequence == "AACCGGTTAA"
        assert rec.cigarstring == "10M"
        assert rec.reference_start == 0

    def test_mates_extending_different_distances(self, fx) -> None:
        """Scenario 12: R1 extends 10bp left, R2 extends 30bp right."""
        ref = fx.ref_seqs["chr1"]
        # R1 [500,520) forward; R2 [510,540) reverse (30bp): overlap 10.
        r1 = fx.make_record(
            "ext",
            0x1 | 0x40 | 0x2,
            UNIQUE,
            "20M",
            ref[UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE + 10,
        )
        r2 = fx.make_record(
            "ext",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE + 10,
            "30M",
            ref[UNIQUE + 10 : UNIQUE + 40],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        assert len(merged) == 1
        rec = merged[0]
        assert rec.reference_start == UNIQUE
        assert rec.cigarstring == "40M"
        assert rec.query_sequence == ref[UNIQUE : UNIQUE + 40]
        tags = dict(rec.get_tags())
        assert tags["NM"] == 0
        assert tags["MD"] == "40"


class TestOuterSoftClips:
    def test_outer_clips_preserved_internal_dropped(self, fx) -> None:
        """Scenario 13: keep left read's lead clip and right read's trail clip.

        R1 also carries an internal (trailing) clip, which must be dropped.
        """
        ref = fx.ref_seqs["chr1"]
        pair = fx.pair(
            "cl", UNIQUE, UNIQUE + 8, read_len=20, r1_lead_clip=4, r1_trail_clip=3, r2_trail_clip=5
        )
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        assert len(merged) == 1
        rec = merged[0]
        # Aligned union [500, 528) = 28M; only outer clips survive.
        assert rec.cigarstring == "4S28M5S"
        expected_seq = "NNNN" + ref[UNIQUE : UNIQUE + 28] + "NNNNN"
        assert rec.query_sequence == expected_seq
        # Qualities: clips and single-mate flanks stay 40; the 12 agreeing
        # overlap bases combine to min(40+40, cap 60) = 60.
        quals = list(rec.query_qualities)
        assert quals == [40] * 12 + [60] * 12 + [40] * 13

    def test_reverse_mate_outer_clip_mirrored(self, fx) -> None:
        """Reverse-strand R2's outer (ref-orientation trailing) clip survives.

        Stored CIGAR is in reference orientation for both strands, so R2's
        ref-orientation trailing 6S is stored verbatim. Its internal 3S
        (ref-orientation leading) clip is dropped per policy.
        """
        ref = fx.ref_seqs["chr1"]
        pair = fx.pair("rv", UNIQUE, UNIQUE + 10, read_len=20, r2_lead_clip=3, r2_trail_clip=6)
        fx.write_bam_and_set(pair)
        # Sanity: R2 stored cigar is the reference-orientation cigar.
        stored = fx.read_bam(fx.last_bam)
        assert any(r.cigarstring == "3S20M6S" for r in stored)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.cigarstring == "30M6S"
        assert rec.query_sequence == ref[UNIQUE : UNIQUE + 30] + "NNNNNN"

    def test_container_owns_trailing_outer_clip(self, fx) -> None:
        """A contained mate's trailing clip is internal and must be dropped.

        R1 (40M3S at UNIQUE) fully contains R2's aligned span (10M2S at
        UNIQUE+10). R2 starts later but ends earlier: the fragment's outer
        trailing clip is R1's 3S, never R2's internal 2S.
        """
        ref = fx.ref_seqs["chr1"]
        r1 = fx.make_record(
            "cont",
            0x1 | 0x40 | 0x2,
            UNIQUE,
            "40M3S",
            ref[UNIQUE : UNIQUE + 40] + "NNN",
            mate_contig="chr1",
            mate_pos=UNIQUE + 10,
        )
        r2 = fx.make_record(
            "cont",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE + 10,
            "10M2S",
            ref[UNIQUE + 10 : UNIQUE + 20] + "NN",
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.cigarstring == "40M3S"
        assert rec.query_sequence == ref[UNIQUE : UNIQUE + 40] + "NNN"


class TestMergedMapqPolicy:
    def test_merged_mapq_is_conservative_adjusted(self, fx) -> None:
        """The consensus MAPQ comes from the ADJUSTED values, not the input.

        Both mates of this repeat pair are ambiguous (adjusted MAPQ 3), so
        the merged record must carry 3 — never the original 60. Eligibility
        still uses the original MAPQ (60 >= 20), so the pair merges at all.
        """
        pair = fx.pair("rm", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.reference_start == LONG_REPEAT
        assert rec.mapping_quality <= 3

    def test_unique_pair_merged_mapq_still_60(self, fx) -> None:
        """Unique unchanged placements keep their original MAPQ in the consensus."""
        pair = fx.pair("um", UNIQUE, UNIQUE + 5, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.mapping_quality == 60


class TestQualityModel:
    def test_unequal_quality_conflict_resolved(self, fx) -> None:
        """Scenario 10: the higher-quality base wins."""
        ref = fx.ref_seqs["chr1"]
        p = UNIQUE + 12
        ref_base = ref[p]
        alt1 = alt_bases(ref_base, 1)[0]
        pair = fx.pair(
            "qq", UNIQUE, UNIQUE + 5, read_len=20, r1_mutations={p: alt1}, r1_qual=40, r2_qual=20
        )
        # R2 keeps the reference base with quality 20; R1's alt has 40.
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        expected = ref[UNIQUE:p] + alt1 + ref[p + 1 : UNIQUE + 25]
        assert rec.query_sequence == expected
        # The conflict resolves toward R1's alt, which differs from the
        # reference: one mismatch in the merged record.
        tags = dict(rec.get_tags())
        assert tags["NM"] == 1
        assert tags["MD"] == f"{p - UNIQUE}{ref_base}{UNIQUE + 25 - p - 1}"

    def test_equal_quality_conflict_emits_n(self, fx) -> None:
        """Scenario 10: equal qualities cannot resolve -> N at low quality."""
        ref = fx.ref_seqs["chr1"]
        p = UNIQUE + 12
        alt1, alt2 = alt_bases(ref[p], 2)
        pair = fx.pair(
            "nn", UNIQUE, UNIQUE + 5, read_len=20, r1_mutations={p: alt1}, r2_mutations={p: alt2}
        )
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.query_sequence[p - UNIQUE] == "N"
        assert rec.query_qualities[p - UNIQUE] == 2

    def test_strict_policy_rejects_unresolved_conflict(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        p = UNIQUE + 12
        alt1, alt2 = alt_bases(ref[p], 2)
        pair = fx.pair(
            "st", UNIQUE, UNIQUE + 5, read_len=20, r1_mutations={p: alt1}, r2_mutations={p: alt2}
        )
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge", extra=["--conflict-policy", "strict"])
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["mergeable"] == 0
        assert by_type["both_aligned"] == 1
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [r.flag & 0xC0 for r in unmerged] == [0x40, 0x80]

    def test_agreeing_qualities_combine_capped(self, fx) -> None:
        fx.ref_seqs["chr1"]
        pair = fx.pair("cap", UNIQUE, UNIQUE + 5, read_len=20, r1_qual=40, r2_qual=40)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        # Overlap [505,520): q = min(40+40, cap 60) = 60.
        assert list(rec.query_qualities)[5:20] == [60] * 15
        # Non-overlap flanks keep single-mate quality 40.
        assert list(rec.query_qualities)[:5] == [40] * 5


class TestIndels:
    def _del_records(self, fx, r2_cigar: str, r2_seq: str):
        ref = fx.ref_seqs["chr1"]
        r1_seq = ref[UNIQUE : UNIQUE + 20] + ref[UNIQUE + 23 : UNIQUE + 40]
        r1 = fx.make_record(
            "ind",
            0x1 | 0x40 | 0x2,
            UNIQUE,
            "20M3D17M",
            r1_seq,
            mate_contig="chr1",
            mate_pos=UNIQUE + 10,
        )
        r2 = fx.make_record(
            "ind",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE + 10,
            r2_cigar,
            r2_seq,
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        return [r1, r2]

    def test_shared_deletion_merged(self, fx) -> None:
        """Scenario 14: identical deletions in both mates survive once."""
        ref = fx.ref_seqs["chr1"]
        r2_seq = ref[UNIQUE + 10 : UNIQUE + 20] + ref[UNIQUE + 23 : UNIQUE + 50]
        records = self._del_records(fx, "10M3D27M", r2_seq)
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        assert len(merged) == 1
        rec = merged[0]
        assert rec.reference_start == UNIQUE
        assert rec.cigarstring == "20M3D27M"
        assert rec.query_sequence == (ref[UNIQUE : UNIQUE + 20] + ref[UNIQUE + 23 : UNIQUE + 50])
        tags = dict(rec.get_tags())
        assert tags["NM"] == 3
        assert tags["MD"] == f"20^{ref[UNIQUE + 20 : UNIQUE + 23]}27"

    def test_conflicting_deletion_unmergeable(self, fx) -> None:
        """R2 aligns across R1's deletion -> indel_conflict, pair retained."""
        ref = fx.ref_seqs["chr1"]
        r2_seq = ref[UNIQUE + 10 : UNIQUE + 45]
        records = self._del_records(fx, "35M", r2_seq)
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert len(unmerged) == 2
        stats = fx.read_parquet(out["read_stats_parquet"])
        assert dict(zip(stats["type"], stats["count"], strict=False))["both_aligned"] == 1

    def test_shared_insertion_merged_once(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        ins = "GGA"
        r1_seq = ref[UNIQUE : UNIQUE + 10] + ins + ref[UNIQUE + 10 : UNIQUE + 20]
        r1 = fx.make_record(
            "ins", 0x1 | 0x40 | 0x2, UNIQUE, "10M3I10M", r1_seq, mate_contig="chr1", mate_pos=UNIQUE
        )
        r2_seq = ref[UNIQUE : UNIQUE + 10] + ins + ref[UNIQUE + 10 : UNIQUE + 20]
        r2 = fx.make_record(
            "ins",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE,
            "10M3I10M",
            r2_seq,
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.cigarstring == "10M3I10M"
        assert rec.query_sequence == r1_seq
        tags = dict(rec.get_tags())
        assert tags["NM"] == 3
        assert tags["MD"] == "20"

    def test_conflicting_insertion_unmergeable(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        r1_seq = ref[UNIQUE : UNIQUE + 10] + "GGA" + ref[UNIQUE + 10 : UNIQUE + 20]
        r2_seq = ref[UNIQUE : UNIQUE + 10] + "CCT" + ref[UNIQUE + 10 : UNIQUE + 20]
        r1 = fx.make_record(
            "ic", 0x1 | 0x40 | 0x2, UNIQUE, "10M3I10M", r1_seq, mate_contig="chr1", mate_pos=UNIQUE
        )
        r2 = fx.make_record(
            "ic",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE,
            "10M3I10M",
            r2_seq,
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert len(fx.read_bam(out["unmerged_bam"])) == 2


class TestNonOverlap:
    def test_non_overlapping_mates_kept_separate(self, fx) -> None:
        """Scenario 15: no overlap -> both_aligned, no exception."""
        pair = fx.pair("far", UNIQUE, UNIQUE + 100, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert len(fx.read_bam(out["unmerged_bam"])) == 2
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["both_aligned"] == 1

    def test_fragment_gap_not_synthesized(self, fx) -> None:
        """Mates separated by an uncovered middle gap stay unmerged."""
        ref = fx.ref_seqs["chr1"]
        r1 = fx.make_record(
            "gap",
            0x1 | 0x40 | 0x2,
            UNIQUE,
            "20M",
            ref[UNIQUE : UNIQUE + 20],
            mate_contig="chr1",
            mate_pos=UNIQUE + 60,
        )
        r2 = fx.make_record(
            "gap",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            UNIQUE + 60,
            "20M",
            ref[UNIQUE + 60 : UNIQUE + 80],
            mate_contig="chr1",
            mate_pos=UNIQUE,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []


class TestMergedRecordSemantics:
    def test_flags_pos_tags_and_mate_removal(self, fx) -> None:
        fx.ref_seqs["chr1"]
        pair = fx.pair("sem", UNIQUE, UNIQUE + 5, read_len=20, tags={"RG": "rg1"})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.flag == 0
        assert rec.query_name == "sem"
        assert rec.next_reference_start == -1
        assert rec.next_reference_id == -1
        assert rec.template_length == 0
        assert rec.mapping_quality == 60
        tags = dict(rec.get_tags())
        assert tags["RG"] == "rg1"
        assert tags["ZC"] == "mergeable"

    def test_read_group_conflict_unmergeable(self, fx) -> None:
        fx.ref_seqs["chr1"]
        r1 = fx.pair("rgc", UNIQUE, UNIQUE + 5, read_len=20, tags={"RG": "rg1"})[0]
        r2 = fx.pair("rgc", UNIQUE, UNIQUE + 5, read_len=20, tags={"RG": "rg2"})[1]
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert len(fx.read_bam(out["unmerged_bam"])) == 2

    def test_merged_not_duplicated_in_unmerged(self, fx) -> None:
        pair = fx.pair("nd", UNIQUE, UNIQUE + 5, read_len=20)
        r1u = fx.make_record("un", 0x1 | 0x40 | 0x4, -1, "*", "ACGTACGTAC")
        r2u = fx.make_record("un", 0x1 | 0x80 | 0x4, -1, "*", "ACGTACGTAC")
        fx.write_bam_and_set([*pair, r1u, r2u])
        code, out = fx.run(mode="merge")
        assert code == 0
        merged_names = [r.query_name for r in fx.read_bam(out["merged_bam"])]
        unmerged_names = [r.query_name for r in fx.read_bam(out["unmerged_bam"])]
        assert merged_names == ["nd"]
        assert "nd" not in unmerged_names
        assert sorted(unmerged_names) == ["un", "un"]
