"""Joint low-MAPQ pair rescue and adjacent joining.

Spec cases: raw reads R1=AA and R2=GG from an FR pair (R2 contributes CC
in reference orientation). Against decoy references (``seq2=AATT``,
``seq3=CCTT``, a duplicate ``seq4=AACC``) the joint placement on
``seq1=AACC`` is uniquely supported even though every individual read has
several equally good alignments and both original MAPQs are 0.
"""

from __future__ import annotations

import ast
import shutil
import subprocess

import pysam
import pytest

from pairloom.config import PairLoomConfig
from pairloom.exceptions import ConfigurationError
from pairloom.reference import ReferenceCache
from pairloom.rescue import _scan_contig

CASE_A = {"seq1": "AACC", "seq2": "AATT"}
CASE_B = {"seq1": "AACC", "seq2": "AATT", "seq3": "CCTT"}
CASE_DUP = {"seq1": "AACC", "seq2": "AATT", "seq3": "CCTT", "seq4": "AACC"}

#: (contig, pos0, reverse) primary placements for R1/R2. Every entry is a
#: valid exact placement of the raw reads AA (R1) / GG (R2): a forward
#: placement needs ref == raw, a reverse placement ref == revcomp(raw).
R1_FWD_SITES = [("seq1", 0, False), ("seq2", 0, False)]
R1_REV_SITES = [("seq2", 2, True), ("seq3", 2, True)]  # ref == TT
R2_REV_SITES = [("seq1", 2, True), ("seq3", 0, True)]  # ref == CC


def build_pair(
    fx,
    ref_seqs: dict[str, str],
    qname: str,
    r1_site: tuple[str, int, bool],
    r2_site: tuple[str, int, bool],
    mapq: int = 0,
    r1_quals: tuple[int, ...] = (30, 31),
    r2_quals: tuple[int, ...] = (21, 22),
) -> list[pysam.AlignedSegment]:
    """Build one FR template at the given exact primary placements.

    ``fx.pair`` cannot express cross-contig or reverse R1 geometry, so the
    records are built directly. Stored SEQ is the reference-orientation
    string (exact alignments), quals follow the stored orientation.
    """
    fx.ref_seqs = ref_seqs
    fx.__post_init__()
    records = []
    for site, mate, read_flag, mate_quals in (
        (r1_site, r2_site, 0x1 | 0x40, r1_quals),
        (r2_site, r1_site, 0x1 | 0x80, r2_quals),
    ):
        contig, pos, rev = site
        mcontig, mpos, mrev = mate
        span = 2
        flag = read_flag | 0x2 | (0x10 if rev else 0) | (0x20 if mrev else 0)
        rec = fx.make_record(
            qname,
            flag,
            pos,
            "2M",
            fx.ref_seqs[contig][pos : pos + span],
            mapq=mapq,
            contig=contig,
            mate_contig=mcontig,
            mate_pos=mpos,
        )
        rec.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + q) for q in mate_quals))
        records.append(rec)
    fx.write_bam_and_set(records)
    return records


def merged_single(fx, out) -> pysam.AlignedSegment:
    merged = fx.read_bam(out["merged_bam"])
    assert len(merged) == 1
    return merged[0]


def reason_summary(capsys) -> dict[str, int]:
    """Parse the stderr merge-reason summary of the last run."""
    err = capsys.readouterr().err
    line = [ln for ln in err.splitlines() if "merge-reasons=" in ln]
    assert line, f"no merge-reason summary in stderr: {err!r}"
    body = line[0].split("merge-reasons=", 1)[1].split(";", 1)[0].strip()
    return ast.literal_eval(body)


class TestRescueCases:
    def test_case_a_rescued_from_either_r1_primary(self, fx, capsys) -> None:
        """Required case A: R1's MAPQ-0 primary on seq1 or on seq2."""
        for r1_site in R1_FWD_SITES:
            build_pair(fx, CASE_A, "a", r1_site, ("seq1", 2, True))
            code, out = fx.run(mode="merge")
            assert code == 0
            rec = merged_single(fx, out)
            assert rec.reference_name == "seq1"
            assert rec.reference_start == 0  # SAM POS 1
            assert rec.query_sequence == "AACC"
            assert rec.cigarstring == "4M"

    def test_case_b_rescued_across_ambiguous_primaries(self, fx) -> None:
        """Required case B: both MAPQs 0, all valid primary combinations.

        Includes reverse-strand primaries (R1 on TT loci, R2 on CC loci)
        and the RR decoy pair on seq3 that the FR model must reject.
        """
        sites = [
            (r1, r2)
            for r1 in R1_FWD_SITES + R1_REV_SITES
            for r2 in R2_REV_SITES
            if r1[0] != r2[0] or r1[1] != r2[1] or r1[2] != r2[2]
        ]
        assert len(sites) >= 6
        for r1_site, r2_site in sites:
            build_pair(fx, CASE_B, "b", r1_site, r2_site)
            code, out = fx.run(mode="merge")
            assert code == 0
            rec = merged_single(fx, out)
            assert (rec.reference_name, rec.reference_start, rec.query_sequence) == (
                "seq1",
                0,
                "AACC",
            )

    def test_merged_record_exact_contract(self, fx) -> None:
        """Required default-parameter output: seq1, POS 1, AACC, 4M, NM 0, MD 4."""
        build_pair(fx, CASE_A, "x", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = merged_single(fx, out)
        tags = dict(rec.get_tags())
        assert rec.flag == 0
        assert rec.cigarstring == "4M"
        assert tags["NM"] == 0
        assert tags["MD"] == "4"
        # Rescued records never claim calibrated joint confidence.
        assert rec.mapping_quality == 0
        assert tags["ZC"] == "mergeable"
        assert tags["ZR"] == 1
        assert rec.next_reference_id == -1
        assert rec.next_reference_start == -1
        assert rec.template_length == 0

    def test_reference_oriented_qualities_concatenated(self, fx) -> None:
        """Qualities follow the reference orientation, reversing when needed.

        R1 AA carries (30, 31); R2's stored CC carries (21, 22) in reference
        order, so the consensus is [30, 31, 21, 22]. When R1's own primary
        sits on a TT locus (reverse, stored TT with quals (30, 31)), the
        flip to the forward rescue placement reverses them to (31, 30).
        """
        build_pair(fx, CASE_A, "q", ("seq1", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge")
        assert code == 0
        assert list(merged_single(fx, out).query_qualities) == [30, 31, 21, 22]

        build_pair(fx, CASE_B, "q", ("seq3", 2, True), ("seq1", 2, True))
        code, out = fx.run(mode="merge")
        assert code == 0
        assert list(merged_single(fx, out).query_qualities) == [31, 30, 21, 22]


class TestHonestAmbiguity:
    def test_duplicate_reference_stays_unmerged(self, fx, capsys) -> None:
        """Required negative: a second AACC copy makes the pair ambiguous.

        No unique rescue, seq1 is never chosen as confidently resolved, the
        records keep their original placements, and no rescue marker is set.
        """
        build_pair(fx, CASE_DUP, "d", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [(r.reference_name, r.reference_start) for r in unmerged] == [
            ("seq2", 0),
            ("seq1", 2),
        ]
        for rec in unmerged:
            assert not rec.has_tag("ZR")
        assert reason_summary(capsys) == {"rescue_ambiguous": 1}

    def test_candidate_cap_blocks_uniqueness_claim(self, fx, capsys) -> None:
        """Truncated discovery must prevent a uniqueness claim.

        R1 has four global occurrences in case B (AA fwd on seq1/seq2, TT
        rev on seq2/seq3); a cap of 3 truncates the search, so rescue must
        refuse even though a unique compatible pair exists.
        """
        build_pair(fx, CASE_B, "c", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--max-candidates-per-read", "3"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert reason_summary(capsys) == {"rescue_search_truncated": 1}

    def test_pair_evaluation_budget_blocks_uniqueness_claim(self, fx, capsys) -> None:
        """Exhausting --max-pair-evaluations must prevent a uniqueness claim."""
        build_pair(fx, CASE_B, "e", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--max-pair-evaluations", "1"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert reason_summary(capsys) == {"rescue_search_truncated": 1}


class TestAdjacentJoining:
    def test_uncovered_base_prevents_joining(self, fx, capsys) -> None:
        """One uncovered reference base rejects the join (fragment_gap).

        seq1=AATCC: R1 AA at [0,2) and R2 (CC, reverse) at [3,5) leave
        position 2 uncovered. The joint placement is unique, so the records
        are still rescued (ZR set) — but no base is synthesized.
        """
        build_pair(fx, {"seq1": "AATCC"}, "g", ("seq1", 0, False), ("seq1", 3, True))
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [(r.reference_start, r.query_sequence) for r in unmerged] == [
            (0, "AA"),
            (3, "CC"),
        ]
        assert all(r.has_tag("ZR") for r in unmerged)
        assert reason_summary(capsys) == {"fragment_gap": 1}

    def test_explicit_min_overlap_1_rejects_adjacency(self, fx, capsys) -> None:
        """Explicit --min-overlap 1 must reject the zero-overlap examples."""
        build_pair(fx, CASE_A, "o", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--min-overlap", "1"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert reason_summary(capsys) == {"insufficient_overlap": 1}
        # The rescue itself succeeded: the relocation stands, only the merge
        # is rejected by the user's explicit overlap requirement.
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert all(r.reference_name == "seq1" and r.has_tag("ZR") for r in unmerged)

    def test_min_overlap_zero_is_valid_and_negative_rejected(self) -> None:
        base = dict(
            input_path="in.bam",
            reference_path="ref.fa",
            mode="merge",
            merged_bam="m.bam",
            unmerged_bam="u.bam",
            disagreement_parquet="d.parquet",
            read_stats_parquet="s.parquet",
        )
        assert PairLoomConfig(**base).min_overlap == 0
        with pytest.raises(ConfigurationError, match="min_overlap"):
            PairLoomConfig(**base, min_overlap=-1)


class TestOptOuts:
    def test_bam_source_preserves_placements(self, fx, capsys) -> None:
        """BAM-only discovery keeps original placements: no rescue search."""
        build_pair(fx, CASE_A, "s", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--candidate-source", "bam"])
        assert code == 0
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [(r.reference_name, r.reference_start) for r in unmerged] == [
            ("seq2", 0),
            ("seq1", 2),
        ]
        assert all(not r.has_tag("ZR") for r in unmerged)
        assert reason_summary(capsys) == {"different_contig": 1}

    def test_bam_source_same_contig_reports_mapq_gate(self, fx, capsys) -> None:
        build_pair(fx, CASE_A, "s2", ("seq1", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--candidate-source", "bam"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert reason_summary(capsys) == {"mapq_below_min": 1}

    def test_no_rescue_flag_disables(self, fx, capsys) -> None:
        build_pair(fx, CASE_A, "n", ("seq2", 0, False), ("seq1", 2, True))
        code, out = fx.run(mode="merge", extra=["--no-rescue"])
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        assert reason_summary(capsys) == {"different_contig": 1}

    def test_high_mapq_discordant_pair_not_rescued(self, fx, capsys) -> None:
        """High-MAPQ discordance is evidence, not ambiguity: never relocated."""
        build_pair(fx, CASE_A, "h", ("seq2", 0, False), ("seq1", 2, True), mapq=60)
        code, out = fx.run(mode="merge")
        assert code == 0
        assert fx.read_bam(out["merged_bam"]) == []
        unmerged = fx.read_bam(out["unmerged_bam"])
        assert [(r.reference_name, r.reference_start) for r in unmerged] == [
            ("seq2", 0),
            ("seq1", 2),
        ]
        assert reason_summary(capsys) == {"different_contig": 1}


class TestMateExchange:
    def test_read_number_exchange_preserves_result(self, fx) -> None:
        """Swapping read numbers yields the same physical fragment.

        With R1=GG and R2=AA the FR geometry still puts AA leftmost-forward
        on seq1 and GG rightmost-reverse: the same AACC consensus.
        """
        # R2 (AA read) forward on seq1:0; R1 (GG read) reverse on seq1:2.
        build_pair(fx, CASE_B, "x", ("seq1", 2, True), ("seq1", 0, False))
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = merged_single(fx, out)
        assert rec.query_sequence == "AACC"
        assert (rec.reference_name, rec.reference_start) == ("seq1", 0)
        assert rec.mapping_quality == 0
        assert dict(rec.get_tags())["ZR"] == 1


class TestStrandFlipAndClips:
    def test_rf_pair_flip_carries_clips_through(self, fx) -> None:
        """rf library: R2 flips strand on rescue and its clip survives.

        s1=AACCGG, s2=AAT. R1 (block AA, reverse, 2S2M stored NNAA) starts
        on s2:0; R2 (block GG, reverse, 2S2M stored NNGG) starts on s1:4 —
        a cross-contig pair. The unique rf joint placement is R1 reverse at
        [0,2) and R2 forward at [2,4): R2's stored sequence flips to CCNN
        (cigar 2M2S), so its clip becomes the fragment's outer trailing
        clip next to R1's outer leading clip.
        """
        fx.ref_seqs = {"s1": "AACCGG", "s2": "AAT"}
        fx.__post_init__()
        r1 = fx.make_record(
            "rf",
            0x1 | 0x40 | 0x2 | 0x10 | 0x20,
            0,
            "2S2M",
            "NNAA",
            mapq=0,
            contig="s2",
            mate_contig="s1",
            mate_pos=4,
        )
        r2 = fx.make_record(
            "rf",
            0x1 | 0x80 | 0x2 | 0x10 | 0x20,
            4,
            "2S2M",
            "NNGG",
            mapq=0,
            contig="s1",
            mate_contig="s2",
            mate_pos=0,
        )
        fx.write_bam_and_set([r1, r2])
        code, out = fx.run(mode="merge", extra=["--library-orientation", "rf"])
        assert code == 0
        rec = merged_single(fx, out)
        assert rec.cigarstring == "2S4M2S"
        assert rec.query_sequence == "NNAACCNN"
        assert (rec.reference_name, rec.reference_start) == ("s1", 0)
        tags = dict(rec.get_tags())
        assert (tags["NM"], tags["MD"], tags["ZR"]) == (0, "4", 1)
        assert rec.mapping_quality == 0


class TestAdjustModeRelocation:
    def test_rescued_records_relocated_in_adjust_mode(self, fx) -> None:
        """Adjust mode emits the rescued placements as two records."""
        build_pair(fx, CASE_B, "adj", ("seq3", 2, True), ("seq3", 0, True))
        code, out = fx.run(mode="adjust")
        assert code == 0
        records = fx.read_bam(out["adjusted_bam"])
        by_flag = {r.flag & 0xC0: r for r in records}
        r1_out, r2_out = by_flag[0x40], by_flag[0x80]
        # R1: seq3:2 reverse (stored TT) -> seq1:0 forward (stored AA).
        assert (r1_out.reference_name, r1_out.reference_start) == ("seq1", 0)
        assert r1_out.query_sequence == "AA"
        assert not r1_out.is_reverse
        assert r1_out.is_proper_pair
        # R2 stays reverse on seq1 at 2 (stored CC).
        assert (r2_out.reference_name, r2_out.reference_start) == ("seq1", 2)
        assert r2_out.is_reverse
        assert r2_out.is_proper_pair
        # Mate fields, TLEN signs, and rescue markers all follow.
        assert r1_out.next_reference_start == 2 and r2_out.next_reference_start == 0
        assert (r1_out.template_length, r2_out.template_length) == (4, -4)
        for rec in records:
            assert rec.mapping_quality == 0
            tags = dict(rec.get_tags())
            assert (tags["NM"], tags["MD"], tags["ZR"]) == (0, "2", 1)
            assert tags["ZC"] == "mergeable"
        stats = fx.read_parquet(out["read_stats_parquet"])
        assert dict(zip(stats["type"], stats["count"], strict=False))["mergeable"] == 1


class TestStatisticsExclusion:
    def test_rescued_pairs_do_not_add_disagreement_observations(self, fx) -> None:
        """Rescued pairs contribute no observations — eligibility unchanged.

        seq1=ATCG holds two templates overlapping position 1: a normal
        MAPQ-60 pair disagreeing there (T->C on the forward mate) and a
        rescued MAPQ-0 pair (R1 AT forward at 0, R2 reverse at 1 with
        stored TC) whose bases agree with the reference at the same
        position. With rescue excluded from statistics the row reports
        dp=2 from the normal pair only; wrongly counting the rescued pair
        would give 4.
        """
        # Rescued pair: unique FR joint placement (0 fwd, 1 rev) — the
        # dinucleotides here are deliberately non-self-complementary so
        # each block supports only one strand.
        rescued = build_pair(fx, {"seq1": "ATCG"}, "res", ("seq1", 0, False), ("seq1", 1, True))
        # Normal pair with an SNV at position 1 on the forward mate.
        normal = fx.pair("norm", 0, 1, read_len=2, r1_mutations={1: "C"}, mapq=60, contig="seq1")
        fx.write_bam_and_set([*rescued, *normal])
        code, out = fx.run(mode="merge")
        assert code == 0
        rows = fx.read_parquet(out["disagreement_parquet"])
        assert len(rows["pos"]) == 1
        assert rows["pos"] == [2]  # 1-based
        assert rows["alt"] == ["C"]
        assert rows["dp"] == [2]
        assert rows["adf"] == [1]
        assert rows["adr"] == [0]
        assert rows["af"] == [0.5]
        # Both templates merged; the rescued one carries the marker.
        merged = {r.query_name: r for r in fx.read_bam(out["merged_bam"])}
        assert set(merged) == {"res", "norm"}
        assert dict(merged["res"].get_tags())["ZR"] == 1
        assert not merged["norm"].has_tag("ZR")


class TestDeterminism:
    def test_identical_results_across_order_and_threads(self, fx) -> None:
        """All templates share one reference (case B); four outcomes mix:
        two rescued-and-merged, one rescue-unsupported (mutated mate), one
        ordinary high-MAPQ adjacent merge."""
        rescued1 = build_pair(fx, CASE_B, "b1", ("seq2", 0, False), ("seq1", 2, True))
        rescued2 = build_pair(fx, CASE_B, "a1", ("seq3", 2, True), ("seq3", 0, True))
        normal = build_pair(fx, CASE_B, "h1", ("seq1", 0, False), ("seq1", 2, True), mapq=60)
        unsupported = build_pair(fx, CASE_B, "m1", ("seq2", 0, False), ("seq1", 2, True))
        # Mismatch at the second base: rescue-unsupported (the pysam SEQ
        # setter drops qualities — restore them).
        unsupported[0].query_sequence = "AT"
        unsupported[0].query_qualities = pysam.qualitystring_to_array('!"')
        pairs = [*rescued1, *rescued2, *normal, *unsupported]
        signatures = []
        for order in (pairs, pairs[::-1]):
            fx.write_bam_and_set(order)
            for extra in (
                ["--threads", "1", "--chunk-size", "1"],
                ["--threads", "4", "--chunk-size", "1"],
            ):
                code, out = fx.run(mode="merge", extra=extra)
                assert code == 0
                sig = sorted(
                    (
                        stream,
                        r.query_name,
                        r.reference_name,
                        r.reference_start,
                        r.flag,
                        r.cigarstring,
                        r.query_sequence,
                        r.mapping_quality,
                        sorted(r.get_tags()),
                    )
                    for stream, path in (
                        ("merged", out["merged_bam"]),
                        ("unmerged", out["unmerged_bam"]),
                    )
                    for r in fx.read_bam(path)
                )
                signatures.append(sig)
        assert signatures[0] == signatures[1] == signatures[2] == signatures[3]
        assert any(r[0] == "merged" and r[1] == "b1" for r in signatures[0])
        assert any(r[0] == "unmerged" and r[1] == "m1" for r in signatures[0])


class TestScanPrimitives:
    def test_chunked_contig_scan_finds_boundary_hits_once(self, fx) -> None:
        """Chunked occurrence scan: cross-boundary hits appear exactly once."""
        span = 4
        # All planted loci are at least `span` apart, and 1_000_002 crosses
        # the first 1 Mbp chunk boundary (fetch overlap must not duplicate
        # or lose it).
        planted = [999_997, 1_000_002, 1_999_999, 2_400_001]
        seq = "A" * 2_400_005
        for pos in planted:
            seq = seq[:pos] + "TGCA" + seq[pos + span :]
        fx.ref_seqs = {"big": seq}
        fx.__post_init__()
        refcache = ReferenceCache(fx.ref_path, max_contigs=1, max_bases=1)
        try:
            hits = _scan_contig(refcache, "big", len(seq), "TGCA")
        finally:
            refcache.close()
        assert hits == planted


class TestIndependentValidation:
    @pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools")
    def test_samtools_quickcheck_rescued_outputs(self, fx) -> None:
        build_pair(fx, CASE_B, "qc", ("seq3", 2, True), ("seq3", 0, True))
        code, out = fx.run(mode="merge", extra=["--sort-output"])
        assert code == 0
        for key in ("merged_bam", "unmerged_bam"):
            result = subprocess.run(
                ["samtools", "quickcheck", "-v", str(out[key])],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, result.stderr
