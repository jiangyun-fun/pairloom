"""Per-read merge-failure classification: ZM tag, chimeric_read, merge-stats."""

from __future__ import annotations

import pysam


def _header(fx) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict(fx.header_dict())


def _rec(header, name, flag, ref, pos, cigar, seq, quals="I" * 200, mapq=60):
    rec = pysam.AlignedSegment(header)
    rec.query_name = name
    rec.flag = flag
    tid = header.get_tid(ref) if ref else -1
    rec.reference_id = tid
    rec.next_reference_id = tid
    rec.reference_start = pos
    rec.cigarstring = cigar
    rec.query_sequence = seq
    rec.query_qualities = pysam.qualitystring_to_array(quals[: len(seq)])
    rec.mapping_quality = mapq
    return rec


class TestChimericDetection:
    def test_chimeric_read_split_from_different_contig(self, fx) -> None:
        """Cross-contig mates with read-through clips classify chimeric_read."""
        # R1: 20 bp on ctg1 then 20 bp clip that also appears in R2's sequence.
        r1_seq = fx.ref_seqs["chr1"][:20] + fx.ref_seqs["chr2"][:20]
        r2_seq = fx.ref_seqs["chr2"][:20] + fx.ref_seqs["chr1"][:20]
        r1 = _rec(_header(fx), "chim", 0x1 | 0x40, "chr1", 0, "20M20S", r1_seq)
        r2 = _rec(_header(fx), "chim", 0x1 | 0x80 | 0x10, "chr2", 0, "20M20S", r2_seq)
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        recs = list(fx.read_bam(out["unmerged_bam"]))
        assert all(r.get_tag("ZM") == "chimeric_read" for r in recs)

    def test_plain_cross_contig_stays_different_contig(self, fx) -> None:
        """Cross-contig mates without read-through evidence keep the old reason."""
        r1_seq = fx.ref_seqs["chr1"][:40]
        r2_seq = fx.ref_seqs["chr2"][10:50]
        r1 = _rec(_header(fx), "plain", 0x1 | 0x40, "chr1", 0, "40M", r1_seq)
        r2 = _rec(_header(fx), "plain", 0x1 | 0x80 | 0x10, "chr2", 10, "40M", r2_seq)
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge")
        assert code == 0
        recs = list(fx.read_bam(out["unmerged_bam"]))
        assert all(r.get_tag("ZM") == "different_contig" for r in recs)


class TestReasonTag:
    def test_merged_and_unaligned_templates_carry_no_zm(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        seq = ref[:60]
        r1 = _rec(_header(fx), "ok", 0x1 | 0x40, "chr1", 0, "60M", seq)
        r2 = _rec(_header(fx), "ok", 0x1 | 0x80 | 0x10, "chr1", 20, "60M", ref[20:80])
        u1 = _rec(_header(fx), "un", 0x1 | 0x40, "chr1", 0, "40M", ref[:40])
        u2 = _rec(_header(fx), "un", 0x1 | 0x80 | 0x4, None, -1, None, ref[:40], mapq=0)
        fx.write_bam_and_set([r1, r2, u1, u2])
        code, out = fx.run(mode="merge")
        assert code == 0
        merged = list(fx.read_bam(out["merged_bam"]))
        assert merged and all(not r.has_tag("ZM") for r in merged)
        unmerged = list(fx.read_bam(out["unmerged_bam"]))
        by_name: dict[str, list] = {}
        for r in unmerged:
            by_name.setdefault(r.query_name, []).append(r)
        assert all(not r.has_tag("ZM") for r in by_name["un"])
        # unaligned-mate templates classify left_aligned_only: no reason tag

    def test_reason_tag_name_configurable(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        r1 = _rec(_header(fx), "x", 0x1 | 0x40, "chr1", 0, "40M", ref[:40])
        r2 = _rec(_header(fx), "x", 0x1 | 0x80 | 0x10, "chr2", 0, "40M", fx.ref_seqs["chr2"][:40])
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge", extra=["--reason-tag", "XY"])
        assert code == 0
        recs = list(fx.read_bam(out["unmerged_bam"]))
        assert all(r.has_tag("XY") for r in recs)


class TestMergeStatsParquet:
    def test_tables_reconcile_with_done_line(self, fx) -> None:
        ref = fx.ref_seqs["chr1"]
        r1 = _rec(_header(fx), "ok", 0x1 | 0x40, "chr1", 0, "60M", ref[:60])
        r2 = _rec(_header(fx), "ok", 0x1 | 0x80 | 0x10, "chr1", 20, "60M", ref[20:80])
        c1 = _rec(_header(fx), "c", 0x1 | 0x40, "chr1", 0, "40M", ref[:40])
        c2 = _rec(_header(fx), "c", 0x1 | 0x80 | 0x10, "chr2", 0, "40M", fx.ref_seqs["chr2"][:40])
        fx.write_bam_and_set([r1, r2, c1, c2])
        code, out = fx.run(
            mode="merge",
            extra=["--merge-stats-parquet", str(fx.tmp_path / "m.parquet")],
        )
        assert code == 0
        t = fx.read_parquet(str(fx.tmp_path / "m.parquet"))
        rows = {(ty, k): c for ty, k, c in zip(t["type"], t["key"], t["count"], strict=True)}
        assert rows[("classification", "mergeable")] == 1
        assert rows[("reason", "different_contig")] == 1
        assert rows[("family", "chimeric_split")] + rows[("family", "merged")] == 2
        hist = [c for (ty, _k), c in rows.items() if ty == "overlap_hist"]
        assert sum(hist) == 1  # one merged template
        assert rows[("reason", "mergeable")] == 1  # successes keyed by classification
