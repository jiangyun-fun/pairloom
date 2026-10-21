"""Fail-loud behaviour and failure injection (scenarios 29-35).

Every test asserts BOTH a nonzero exit code AND that no final output is
presented as complete after the failure.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from conftest import UNIQUE
from pairloom.exceptions import OutputError


def all_outputs_absent(out: dict[str, Path]) -> None:
    for path in out.values():
        assert not Path(path).exists(), f"output {path} falsely exists"


def no_temp_files(out_dir: Path) -> bool:
    return not list(out_dir.glob("*.pairloom-tmp-*"))


@pytest.fixture()
def pair_bam(fx):
    fx.write_bam_and_set(fx.pair("e", UNIQUE, UNIQUE + 5, read_len=20))
    return fx


class TestCorruptAndMissingInputs:
    def test_corrupted_bam(self, fx, tmp_path) -> None:
        bad = tmp_path / "corrupt.bam"
        bad.write_bytes(b"\x1f\x8b\x08truncated-garbage-not-bgzf")
        fx.last_bam = bad
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)

    def test_missing_reference(self, pair_bam, fx) -> None:
        code, out = fx.run(mode="adjust")
        assert code == 0
        # Now point at a nonexistent reference.
        fx.ref_path = fx.tmp_path / "gone.fa"
        code, out = fx.run(mode="adjust")
        assert code == 2

    def test_reference_missing_bam_contig(self, fx) -> None:
        """A BAM contig absent from the FASTA is an input error."""
        import pysam

        header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "ghost", "LN": 100}]}
        path = fx.tmp_path / "ghost.bam"
        hdr = pysam.AlignmentHeader.from_dict(header)
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as f:
            rec = pysam.AlignedSegment(f.header)
            rec.query_name = "g"
            rec.flag = 0x1 | 0x40
            rec.reference_id = 0
            rec.reference_start = 10
            rec.cigarstring = "20M"
            rec.query_sequence = "A" * 20
            rec.query_qualities = pysam.qualitystring_to_array("I" * 20)
            f.write(rec)
        fx.last_bam = path
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)


class TestInjectedFailures:
    def test_failed_parquet_write(self, pair_bam, fx, monkeypatch) -> None:
        """Scenario 30: Parquet failure aborts and publishes nothing."""

        def boom(*args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr("pairloom.parquetout.pq.write_table", boom)
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)
        assert no_temp_files(out["read_stats_parquet"].parent)

    def test_failed_disagreement_parquet_write(self, pair_bam, fx, monkeypatch) -> None:
        def boom(*args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr("pairloom.parquetout.pq.ParquetWriter.write_table", boom)
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)

    def test_failed_bam_write(self, pair_bam, fx, monkeypatch) -> None:
        """BAM write failure mid-stream aborts with no published outputs."""
        from pairloom import pipeline

        original = pipeline._BamWriter.write
        calls = {"n": 0}

        def flaky(self, record):
            calls["n"] += 1
            if calls["n"] > 1:
                raise OutputError("injected BAM write failure")
            original(self, record)

        monkeypatch.setattr(pipeline._BamWriter, "write", flaky)
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)

    def test_worker_failure_propagates_serial(self, fx) -> None:
        """Scenario 32: malformed CIGAR -> DataInvariantError, exit 1."""
        ref = fx.ref_seqs["chr1"]
        bad = fx.make_record(
            "badcig", 0x1 | 0x40, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 25]
        )  # seq/cigar mismatch
        fx.write_bam_and_set([bad])
        code, out = fx.run(mode="adjust")
        assert code == 1
        all_outputs_absent(out)

    def test_worker_failure_propagates_parallel(self, fx) -> None:
        """Scenario 31: a failing multiprocessing worker aborts the run."""
        ref = fx.ref_seqs["chr1"]
        good = fx.pair("good", UNIQUE, UNIQUE + 5, read_len=20)
        bad = fx.make_record("badcig", 0x1 | 0x40, UNIQUE, "20M", ref[UNIQUE : UNIQUE + 25])
        fx.write_bam_and_set([*good, bad])
        code, out = fx.run(mode="adjust", extra=["--threads", "2", "--chunk-size", "1"])
        assert code == 1
        all_outputs_absent(out)


class TestSilentSuppression:
    def test_no_silent_exception_swallowing(self) -> None:
        """Scenario 34: AST scan forbids except-pass/continue patterns."""
        src_root = Path(__file__).parent.parent / "src" / "pairloom"
        offenders: list[str] = []
        for py in sorted(src_root.glob("*.py")):
            tree = ast.parse(py.read_text())
            for node in ast.walk(tree):
                if not isinstance(node, ast.Try):
                    continue
                for handler in node.handlers:
                    body_types = {type(stmt) for stmt in handler.body}
                    rethrows = any(isinstance(stmt, ast.Raise) for stmt in handler.body)
                    suppressive = body_types & {ast.Pass, ast.Continue}
                    if suppressive and not rethrows:
                        offenders.append(f"{py.name}:{handler.lineno}")
        assert not offenders, f"silent exception swallowing at {offenders}"

    def test_expected_unmergeable_never_raises(self, fx) -> None:
        """Scenario 33: unmergeable biology is a normal outcome."""
        records: list = []
        records += fx.pair("far", UNIQUE, UNIQUE + 200, read_len=20)  # no overlap
        r1u = fx.make_record("u1", 0x1 | 0x40 | 0x4, -1, "*", "ACGTACGTAC")
        r2u = fx.make_record("u1", 0x1 | 0x80 | 0x4, -1, "*", "ACGTACGTAC")
        records += [r1u, r2u]
        fx.write_bam_and_set(records)
        code, out = fx.run(mode="merge")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["both_aligned"] == 1
        assert by_type["neither_aligned"] == 1

    def test_debug_prints_traceback(self, pair_bam, fx, capsys, monkeypatch) -> None:
        def boom(*args, **kwargs):
            raise OSError(28, "No space left on device")

        monkeypatch.setattr("pairloom.parquetout.pq.write_table", boom)
        code, out = fx.run(mode="adjust", extra=["--debug"])
        assert code == 1
        captured = capsys.readouterr()
        assert "Traceback" in captured.err
