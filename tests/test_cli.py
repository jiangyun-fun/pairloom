"""CLI validation, output handling, sorting, indexing (scenarios 27, 28)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from conftest import UNIQUE


def base_args(fx, out_dir: Path, mode: str = "adjust") -> list[str]:
    args = [
        "--input",
        str(fx.last_bam),
        "--reference",
        str(fx.ref_path),
        "--mode",
        mode,
        "--disagreement-parquet",
        str(out_dir / "d.parquet"),
        "--read-stats-parquet",
        str(out_dir / "s.parquet"),
    ]
    if mode == "adjust":
        args += ["--adjusted-bam", str(out_dir / "a.bam")]
    else:
        args += ["--merged-bam", str(out_dir / "m.bam"), "--unmerged-bam", str(out_dir / "u.bam")]
    return args


@pytest.fixture()
def one_pair_bam(fx):
    fx.write_bam_and_set(fx.pair("cli", UNIQUE, UNIQUE + 5, read_len=20))
    return fx


class TestPathValidation:
    def test_missing_required_output_for_mode(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o1"
        out.mkdir()
        args = [
            "--input",
            str(fx.last_bam),
            "--reference",
            str(fx.ref_path),
            "--mode",
            "adjust",
            "--disagreement-parquet",
            str(out / "d.parquet"),
            "--read-stats-parquet",
            str(out / "s.parquet"),
            # --adjusted-bam missing
        ]
        assert main(args) == 2

    def test_wrong_mode_output_rejected(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o2"
        out.mkdir()
        args = base_args(fx, out, mode="merge")
        args += ["--adjusted-bam", str(out / "a.bam")]  # adjust-only output
        assert main(args) == 2

    def test_existing_output_without_force(self, one_pair_bam, fx, tmp_path) -> None:
        """Scenario 28: existing outputs refuse to be overwritten."""
        from pairloom.cli import main

        out = tmp_path / "o3"
        out.mkdir()
        existing = out / "a.bam"
        existing.write_bytes(b"precious")
        assert main(base_args(fx, out)) == 2
        assert existing.read_bytes() == b"precious"
        # --force overwrites.
        assert main(base_args(fx, out) + ["--force"]) == 0
        assert existing.read_bytes() != b"precious"

    def test_duplicate_output_paths_rejected(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o4"
        out.mkdir()
        args = base_args(fx, out)
        args[args.index("--read-stats-parquet") + 1] = str(out / "d.parquet")
        assert main(args) == 2

    def test_output_over_input_rejected(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o5"
        out.mkdir()
        args = base_args(fx, out)
        args[args.index("--adjusted-bam") + 1] = str(fx.last_bam)
        assert main(args) == 2

    def test_unwritable_directory(self, one_pair_bam, fx, tmp_path) -> None:
        import os

        from pairloom.cli import main

        if os.geteuid() == 0:
            pytest.skip("root ignores directory permissions")
        locked = tmp_path / "locked"
        locked.mkdir()
        locked.chmod(0o500)
        try:
            args = base_args(fx, locked)
            assert main(args) == 2
        finally:
            locked.chmod(0o700)

    def test_missing_input(self, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o6"
        out.mkdir()
        args = [
            "--input",
            str(out / "nope.bam"),
            "--reference",
            str(fx.ref_path),
            "--mode",
            "adjust",
            "--adjusted-bam",
            str(out / "a.bam"),
            "--disagreement-parquet",
            str(out / "d.parquet"),
            "--read-stats-parquet",
            str(out / "s.parquet"),
        ]
        assert main(args) == 2

    def test_missing_reference_index(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        ref = tmp_path / "ref2.fa"
        shutil.copy(fx.ref_path, ref)
        out = tmp_path / "o7"
        out.mkdir()
        args = base_args(fx, out)
        args[args.index("--reference") + 1] = str(ref)
        assert main(args) == 2

    def test_invalid_status_tag(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o8"
        out.mkdir()
        assert main(base_args(fx, out) + ["--status-tag", "TOOLONG"]) == 2

    def test_bad_compression_choice(self, one_pair_bam, fx, tmp_path) -> None:
        from pairloom.cli import main

        out = tmp_path / "o9"
        out.mkdir()
        args = base_args(fx, out)
        # choices= enforced by argparse -> SystemExit(2)
        with pytest.raises(SystemExit) as excinfo:
            main(args + ["--parquet-compression", "brotli"])
        assert excinfo.value.code == 2

    def test_existing_index_without_force_refused(self, one_pair_bam, fx, tmp_path) -> None:
        """A stale .bai is an output like any other: overwriting needs --force."""
        from pairloom.cli import main

        out = tmp_path / "oi"
        out.mkdir()
        args = base_args(fx, out) + ["--index-output"]
        assert main(args) == 0
        idx = out / "a.bam.bai"
        assert idx.read_bytes()[:4] == b"BAI\x01"
        # Keep only the stale index; declared outputs are gone.
        (out / "a.bam").unlink()
        (out / "d.parquet").unlink()
        (out / "s.parquet").unlink()
        idx.write_bytes(b"stale")
        assert main(args) == 2
        assert idx.read_bytes() == b"stale"
        assert main(args + ["--force"]) == 0
        assert idx.read_bytes()[:4] == b"BAI\x01"

    def test_index_colliding_with_output_rejected(self, one_pair_bam, fx, tmp_path) -> None:
        """An index path equal to another requested output fails validation."""
        from pairloom.cli import main

        out = tmp_path / "oc"
        out.mkdir()
        args = base_args(fx, out, mode="merge")
        # The unmerged BAM sits exactly where the merged BAM's index goes.
        args[args.index("--unmerged-bam") + 1] = str(out / "m.bam.bai")
        args += ["--index-output"]
        assert main(args) == 2
        assert not (out / "m.bam").exists()


class TestSortAndIndex:
    def test_default_output_is_template_order_unindexed(self, one_pair_bam, fx) -> None:
        code, out = fx.run(mode="adjust")
        assert code == 0
        assert not Path(str(out["adjusted_bam"]) + ".bai").exists()
        header = pysam.AlignmentFile(str(out["adjusted_bam"])).header.to_dict()
        assert header["HD"]["SO"] == "queryname"

    def test_sort_output_coordinate_order(self, one_pair_bam, fx) -> None:
        code, out = fx.run(mode="adjust", extra=["--sort-output"])
        assert code == 0
        header = pysam.AlignmentFile(str(out["adjusted_bam"])).header.to_dict()
        assert header["HD"]["SO"] == "coordinate"
        positions = [r.reference_start for r in fx.read_bam(out["adjusted_bam"])]
        assert positions == sorted(positions)

    def test_index_output_implies_sort_and_creates_bai(self, one_pair_bam, fx) -> None:
        code, out = fx.run(mode="adjust", extra=["--index-output"])
        assert code == 0
        idx = Path(str(out["adjusted_bam"]) + ".bai")
        assert idx.exists()
        # pysam can query the index.
        with pysam.AlignmentFile(str(out["adjusted_bam"])) as f:
            assert f.has_index()
            assert len(list(f.fetch("chr1", 0, 1000))) >= 2


class TestSubprocess:
    def test_console_script_entry_point(self, one_pair_bam, fx, tmp_path) -> None:
        """The installed `pairmerge` entry point wires exit codes correctly."""
        exe = shutil.which("pairmerge")
        if exe is None:
            pytest.skip("pairmerge not on PATH")
        out = tmp_path / "o10"
        out.mkdir()
        proc = subprocess.run([exe, *base_args(fx, out)], capture_output=True, text=True)
        assert proc.returncode == 0
        assert (out / "a.bam").exists()

    def test_version_flag(self) -> None:
        exe = shutil.which("pairmerge")
        if exe is None:
            pytest.skip("pairmerge not on PATH")
        proc = subprocess.run([exe, "--version"], capture_output=True, text=True)
        assert proc.returncode == 0
        assert "pairmerge" in proc.stdout
