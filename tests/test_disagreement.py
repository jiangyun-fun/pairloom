"""Overlapping-base disagreement statistics (scenarios 6-9, 11, 24, 25, 5.4).

All expectations are derived from the fixture reference strings. Reference
positions below are 0-based; Parquet ``pos`` is 1-based.
"""

from __future__ import annotations

import pyarrow.parquet as pq

from conftest import UNIQUE

# Overlap window used throughout: R1 [500,520) forward, R2 [510,530) reverse.
# Overlapping reference positions: [510,520) -> 10 bases >= --min-overlap.
P = 515  # the disagreement position inside the overlap


def alt_bases(ref_base: str, n: int) -> list[str]:
    choices = [b for b in "ACGT" if b != ref_base]
    return choices[:n]


def read_disagreement(fx, out) -> list[dict]:
    table = pq.read_table(str(out["disagreement_parquet"]))
    return table.to_pylist()


class TestSchema:
    def test_exact_schema_and_order(self, fx) -> None:
        pair = fx.pair("p", UNIQUE, UNIQUE + 10, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        schema = pq.read_schema(str(out["disagreement_parquet"]))
        assert [f.name for f in schema] == [
            "seq_name",
            "pos",
            "ref",
            "alt",
            "dp",
            "adf",
            "adr",
            "af",
        ]
        types = [str(f.type) for f in schema]
        assert types == ["string", "int64", "string", "string", "int64", "int64", "int64", "double"]

    def test_stats_schema(self, fx) -> None:
        pair = fx.pair("p", UNIQUE, UNIQUE + 10, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        schema = pq.read_schema(str(out["read_stats_parquet"]))
        assert [f.name for f in schema] == ["type", "count"]
        assert [str(f.type) for f in schema] == ["string", "int64"]

    def test_empty_input_empty_tables_with_schema(self, fx) -> None:
        """Scenario 24: empty Parquet outputs keep the correct schemas."""
        fx.write_bam_and_set([])
        code, out = fx.run(mode="merge")
        assert code == 0
        dis = pq.read_table(str(out["disagreement_parquet"]))
        assert dis.num_rows == 0
        assert [f.name for f in dis.schema] == [
            "seq_name",
            "pos",
            "ref",
            "alt",
            "dp",
            "adf",
            "adr",
            "af",
        ]
        stats = fx.read_parquet(out["read_stats_parquet"])
        assert stats["type"] == [
            "total",
            "mergeable",
            "left_aligned_only",
            "right_aligned_only",
            "both_aligned",
            "neither_aligned",
        ]
        assert stats["count"] == [0, 0, 0, 0, 0, 0]


class TestDisagreementCounting:
    def test_perfect_agreement_reports_nothing(self, fx) -> None:
        """Scenario 6: agreeing overlaps never appear in the table."""
        pair = fx.pair("ok", UNIQUE, UNIQUE + 10, read_len=20)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert read_disagreement(fx, out) == []

    def test_single_disagreement(self, fx) -> None:
        """Scenario 7: one conflicting base produces exactly one row."""
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1 = alt_bases(ref_base, 1)[0]
        pair = fx.pair("d1", UNIQUE, UNIQUE + 10, read_len=20, r1_mutations={P: alt1})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        assert rows == [
            {
                "seq_name": "chr1",
                "pos": P + 1,
                "ref": ref_base,
                "alt": alt1,
                "dp": 2,
                "adf": 1,
                "adr": 0,
                "af": 0.5,
            }
        ]

    def test_reverse_strand_r2_counts_adr(self, fx) -> None:
        """Scenario 11: strand comes from alignment orientation (R2 reverse)."""
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt2 = alt_bases(ref_base, 1)[0]
        # R2 carries the non-reference base; R1 agrees with the reference.
        pair = fx.pair("d2", UNIQUE, UNIQUE + 10, read_len=20, r2_mutations={P: alt2})
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        assert rows == [
            {
                "seq_name": "chr1",
                "pos": P + 1,
                "ref": ref_base,
                "alt": alt2,
                "dp": 2,
                "adf": 0,
                "adr": 1,
                "af": 0.5,
            }
        ]

    def test_both_differ_from_reference_two_rows(self, fx) -> None:
        """Scenario 8 / spec 5.4: R1 and R2 disagree with each other AND ref."""
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1, alt2 = alt_bases(ref_base, 2)
        pair = fx.pair(
            "d3", UNIQUE, UNIQUE + 10, read_len=20, r1_mutations={P: alt1}, r2_mutations={P: alt2}
        )
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        # Sorted by alt allele.
        expected = [
            {
                "seq_name": "chr1",
                "pos": P + 1,
                "ref": ref_base,
                "alt": alt1,
                "dp": 2,
                "adf": 1,
                "adr": 0,
                "af": 0.5,
            },
            {
                "seq_name": "chr1",
                "pos": P + 1,
                "ref": ref_base,
                "alt": alt2,
                "dp": 2,
                "adf": 0,
                "adr": 1,
                "af": 0.5,
            },
        ]
        expected.sort(key=lambda r: r["alt"])
        assert rows == expected

    def test_multiple_alternatives_and_agreeing_pairs(self, fx) -> None:
        """Scenario 9: dp includes observations from agreeing overlapping pairs."""
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1, alt2, alt3 = alt_bases(ref_base, 3)
        records = []
        records += fx.pair(
            "p1", UNIQUE, UNIQUE + 10, read_len=20, r1_mutations={P: alt1}, r2_mutations={P: alt2}
        )
        # p2: both mates agree on a third non-reference allele.
        records += fx.pair(
            "p2",
            UNIQUE + 5,
            UNIQUE + 15,
            read_len=20,
            r1_mutations={P: alt3},
            r2_mutations={P: alt3},
        )
        # p3: both mates agree with the reference (adds depth only).
        records += fx.pair("p3", UNIQUE + 3, UNIQUE + 13, read_len=20)
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        by_alt = {row["alt"]: row for row in rows}
        assert set(by_alt) == {alt1, alt2, alt3}
        assert by_alt[alt1]["dp"] == 6 and by_alt[alt1]["adf"] == 1
        assert by_alt[alt1]["adr"] == 0 and by_alt[alt1]["af"] == 1 / 6
        assert by_alt[alt2]["adr"] == 1 and by_alt[alt2]["af"] == 1 / 6
        assert by_alt[alt3]["adf"] == 1 and by_alt[alt3]["adr"] == 1
        assert by_alt[alt3]["af"] == 2 / 6

    def test_ordering_by_header_pos_alt(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        # Two disagreeing positions on chr1 (alts sorted), one on chr2.
        ref_base_a, ref_base_b = ref[P], ref[P + 2]
        alt_a1, alt_a2 = alt_bases(ref_base_a, 2)
        alt_b = alt_bases(ref_base_b, 1)[0]
        chr2 = fx.ref_seqs["chr2"]
        alt_c = alt_bases(chr2[50], 1)[0]
        records = []
        records += fx.pair(
            "o1",
            UNIQUE,
            UNIQUE + 10,
            read_len=20,
            r1_mutations={P: alt_a1},
            r2_mutations={P: alt_a2},
        )
        records += fx.pair(
            "o2", UNIQUE + 10, UNIQUE + 20, read_len=20, r1_mutations={P + 12: alt_b}
        )
        chr2_pair = fx.pair("o3", 40, 50, read_len=20, contig="chr2", r1_mutations={50: alt_c})
        records += chr2_pair
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        keys = [(row["seq_name"], row["pos"], row["alt"]) for row in rows]
        assert keys == sorted(keys)
        assert keys[0][0] == "chr1"
        assert keys[-1][0] == "chr2"


class TestFilters:
    def test_low_base_quality_excluded(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1 = alt_bases(ref_base, 1)[0]
        # R1's whole read at quality 15 < min 20: its bases are ineligible.
        pair = fx.pair("q1", UNIQUE, UNIQUE + 10, read_len=20, r1_mutations={P: alt1}, r1_qual=15)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert read_disagreement(fx, out) == []

    def test_low_mapq_excluded(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1 = alt_bases(ref_base, 1)[0]
        pair = fx.pair("q2", UNIQUE, UNIQUE + 10, read_len=20, r1_mutations={P: alt1}, mapq=10)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert read_disagreement(fx, out) == []

    def test_disagreements_computed_from_adjusted_pre_merge_mates(self, fx) -> None:
        """Scenario: statistics reflect the two original mates, never a consensus."""
        ref = fx.ref_seqs["chr1"]
        ref_base = ref[P]
        alt1, alt2 = alt_bases(ref_base, 2)
        pair = fx.pair(
            "pm",
            UNIQUE,
            UNIQUE + 10,
            read_len=20,
            r1_mutations={P: alt1},
            r2_mutations={P: alt2},
            r1_qual=40,
            r2_qual=20,
        )  # R1 wins the conflict
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = read_disagreement(fx, out)
        assert len(rows) == 2  # both original alleles reported pre-consensus
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["mergeable"] == 1

    def test_soft_clips_and_indels_never_counted(self, fx) -> None:
        """Clipped and indel bases do not create overlapping observations."""
        fx.ref_seqs["chr1"]
        pair = fx.pair("sc", UNIQUE, UNIQUE + 10, read_len=20, r1_lead_clip=4, r2_trail_clip=4)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert read_disagreement(fx, out) == []
        merged = fx.read_bam(out["merged_bam"])
        assert merged[0].cigarstring == "4S30M4S"
