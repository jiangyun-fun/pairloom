"""The classic gallery BAM: every documented outcome, checked end to end.

Builds examples/build_classic_bam.py's deterministic gallery and asserts
the exact per-template results, so the gallery (and OUTCOMES.md) cannot
drift from the implementation without a failing test.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]

EXPECTED_MERGED = {
    "ok_overlap",
    "ok_adjacent",
    "ok_equal_qual_N",
    "ok_rescued",
    "ok_nonprimary",
}

EXPECTED_ZM = {
    "fail_cross_contig": "different_contig",
    "fail_chimeric": "chimeric_read",
    "fail_gap": "fragment_gap",
    "fail_frag_len": "fragment_length_out_of_bounds",
    "fail_orientation": "orientation_incompatible",
    "fail_indel": "indel_conflict",
    "fail_insertion": "insertion_conflict",
    "fail_hardclip": "hard_clipped_mate",
    "fail_noqual": "qualities_missing",
    "fail_refskip": "reference_skip_unsupported",
    "fail_rg": "read_group_conflict",
    "fail_rescue_unsupported": "rescue_unsupported_alignment",
    "fail_rescue_ambiguous": "rescue_ambiguous",
    "fail_rescue_no_pair": "rescue_no_compatible_pair",
    "fail_rescue_truncated": "rescue_search_truncated",
}

EXPECTED_ZC = {
    "cls_left_only": "left_aligned_only",
    "cls_right_only": "right_aligned_only",
    "cls_neither": "neither_aligned",
}


@pytest.fixture(scope="module")
def gallery(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("gallery")
    subprocess.run(
        [sys.executable, str(REPO / "examples" / "build_classic_bam.py"), "--out-dir", str(out)],
        check=True,
    )
    subprocess.run(["samtools", "faidx", str(out / "classic.fa")], check=True)
    return out


def _run(gallery: Path, dest: Path, extra: list[str] | None = None) -> dict[str, dict]:
    from pairloom.cli import main

    dest.mkdir(exist_ok=True)
    code = main(
        [
            "--input",
            str(gallery / "classic.bam"),
            "--reference",
            str(gallery / "classic.fa"),
            "--mode",
            "merge",
            "--merged-bam",
            str(dest / "merged.bam"),
            "--unmerged-bam",
            str(dest / "unmerged.bam"),
            "--disagreement-parquet",
            str(dest / "disg.parquet"),
            "--read-stats-parquet",
            str(dest / "stats.parquet"),
            "--reason-tag",
            "ZM",
            "--threads",
            "1",
            "--force",
            *(extra or []),
        ]
    )
    assert code == 0
    out: dict[str, dict] = {}
    with pysam.AlignmentFile(str(dest / "merged.bam")) as fh:
        for r in fh:
            out[r.query_name] = {
                "stream": "merged",
                "ZC": r.get_tag("ZC"),
                "ZR": r.get_tag("ZR") if r.has_tag("ZR") else None,
                "ZN": r.get_tag("ZN") if r.has_tag("ZN") else None,
                "mapq": r.mapping_quality,
                "seq": r.query_sequence,
            }
    with pysam.AlignmentFile(str(dest / "unmerged.bam")) as fh:
        for r in fh:
            if r.has_tag("ZM"):
                out[r.query_name] = {
                    "stream": "unmerged",
                    "ZC": r.get_tag("ZC"),
                    "ZM": r.get_tag("ZM"),
                }
            elif r.has_tag("ZC") and r.query_name not in out:
                out[r.query_name] = {"stream": "unmerged", "ZC": r.get_tag("ZC"), "ZM": None}
    return out


class TestGalleryDefaultRun:
    def test_every_template_accounted(self, gallery, tmp_path) -> None:
        out = _run(gallery, tmp_path / "default")
        assert set(out) == EXPECTED_MERGED | set(EXPECTED_ZM) | set(EXPECTED_ZC)
        assert {q for q, r in out.items() if r["stream"] == "merged"} == EXPECTED_MERGED

    def test_success_paths_and_tags(self, gallery, tmp_path) -> None:
        out = _run(gallery, tmp_path / "default")
        assert out["ok_rescued"]["ZR"] == 1 and out["ok_rescued"]["mapq"] == 0
        assert out["ok_nonprimary"]["ZN"] == 1 and out["ok_nonprimary"]["mapq"] == 0
        assert "N" in out["ok_equal_qual_N"]["seq"]  # equal-quality conflict -> N
        for q in ("ok_overlap", "ok_adjacent"):
            assert out[q]["mapq"] == 60 and out[q]["ZC"] == "mergeable"

    def test_every_failure_reason(self, gallery, tmp_path) -> None:
        out = _run(gallery, tmp_path / "default")
        for q, reason in EXPECTED_ZM.items():
            assert out[q]["ZM"] == reason, f"{q}: expected {reason}, got {out[q]['ZM']}"
            assert out[q]["ZC"] == "both_aligned"

    def test_classifications(self, gallery, tmp_path) -> None:
        out = _run(gallery, tmp_path / "default")
        for q, cls in EXPECTED_ZC.items():
            assert out[q]["ZC"] == cls
        # only-non-primary-records template: its single secondary record
        # passes through to the unmerged stream unchanged and untagged
        with pysam.AlignmentFile(str(tmp_path / "default" / "unmerged.bam")) as fh:
            recs = [r for r in fh.fetch(until_eof=True) if r.query_name == "cls_only_nonprimary"]
        assert len(recs) == 1
        assert recs[0].flag & 0x100 and not recs[0].has_tag("ZC")


class TestGalleryFlagRuns:
    @pytest.mark.parametrize(
        ("extra", "qname", "reason"),
        [
            (["--no-merge-non-primary"], "ok_nonprimary", "non_primary_records"),
            (["--no-rescue"], "ok_rescued", "mapq_below_min"),
            (["--conflict-policy", "strict"], "ok_equal_qual_N", "strict_conflict"),
            (["--min-overlap", "30"], "ok_overlap", "insufficient_overlap"),
        ],
    )
    def test_policy_gated_reasons(self, gallery, tmp_path_factory, extra, qname, reason) -> None:
        out = _run(gallery, tmp_path_factory.mktemp("flag"), extra)
        assert out[qname]["stream"] == "unmerged"
        assert out[qname]["ZM"] == reason
