"""Default-on merging of templates carrying non-primary records."""

from __future__ import annotations

from conftest import UNIQUE


def _supplementary(fx, qname: str):
    ref = fx.ref_seqs["chr1"]
    return fx.make_record(
        qname, 0x1 | 0x40 | 0x800, UNIQUE + 2, "15M", ref[UNIQUE + 2 : UNIQUE + 17]
    )


class TestDefaultOn:
    def test_supplementary_template_merges_with_evidence_tag(self, fx) -> None:
        """Concordant primaries + supplementary evidence merge by default."""
        pair = fx.pair("np", UNIQUE, UNIQUE + 5, read_len=20)
        fx.write_bam_and_set([*pair, _supplementary(fx, "np")])
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = fx.read_bam(out["merged_bam"])
        assert len(merged) == 1
        rec = merged[0]
        assert rec.query_name == "np"
        assert rec.get_tag("ZC") == "mergeable"
        assert rec.get_tag("ZN") == 1  # evidence marker
        assert rec.mapping_quality == 0  # honesty policy, like rescued pairs
        # Nothing for this template reaches the unmerged stream (drop policy)
        assert [r.query_name for r in fx.read_bam(out["unmerged_bam"])] == []
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["mergeable"] == 1

    def test_refused_case_keeps_all_records_with_real_reason(self, fx) -> None:
        """A supplementary template failing a real gate keeps everything."""
        # Mates 400 bp apart: same contig, FR, MAPQ 60, but an uncovered
        # fragment gap -> fragment_gap (not rescue-triggering, not shape).
        pair = fx.pair("gap", UNIQUE, UNIQUE + 400, read_len=20)
        fx.write_bam_and_set([*pair, _supplementary(fx, "gap")])
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert len(unmerged) == 3  # primaries + supplementary, nothing dropped
        prim = [r for r in unmerged if r.has_tag("ZM")]
        assert {r.get_tag("ZM") for r in prim} == {"fragment_gap"}
        assert all(not r.has_tag("ZN") for r in unmerged)


class TestOptOut:
    def test_no_merge_non_primary_restores_refusal(self, fx) -> None:
        """--no-merge-non-primary reproduces the old non_primary_records gate."""
        pair = fx.pair("opt", UNIQUE, UNIQUE + 5, read_len=20)
        fx.write_bam_and_set([*pair, _supplementary(fx, "opt")])
        code, out = fx.run(mode="merge", extra=["--no-merge-non-primary"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [r.query_name for r in unmerged] == ["opt", "opt", "opt"]
        prim = [r for r in unmerged if r.has_tag("ZM")]
        assert {r.get_tag("ZM") for r in prim} == {"non_primary_records"}
