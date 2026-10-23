"""Compare two pairloom output sets for exact equivalence.

Compares merged/unmerged BAMs record by record (SAM-string equality in
stream order), disagreement and read-stats parquet frames, and the done
line. Exit 0 only when everything matches.

Usage:
    python benchmarks/equivalence_check.py --a-dir BASE_OUT --b-dir CAND_OUT
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import polars as pl
import pysam


def sam_lines(path: Path) -> list[str]:
    out = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            out.append(rec.to_string())
    return out


def compare_bam(a: Path, b: Path, label: str) -> bool:
    la, lb = sam_lines(a), sam_lines(b)
    if la == lb:
        print(f"BAM {label}: identical ({len(la)} records)")
        return True
    print(f"BAM {label}: MISMATCH ({len(la)} vs {len(lb)} records)")
    for i, (x, y) in enumerate(zip(la, lb, strict=False)):
        if x != y:
            print(f"  first diff at record {i}:\n    A: {x}\n    B: {y}")
            break
    return False


def compare_parquet(a: Path, b: Path, label: str, sort_by: list[str] | None = None) -> bool:
    fa, fb = pl.read_parquet(a), pl.read_parquet(b)
    if sort_by:
        fa = fa.sort(sort_by)
        fb = fb.sort(sort_by)
    if fa.equals(fb):
        print(f"parquet {label}: identical ({fa.height} rows)")
        return True
    print(f"parquet {label}: MISMATCH ({fa.height} vs {fb.height} rows)")
    print(f"  schemas equal: {fa.schema == fb.schema}")
    return False


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a-dir", required=True)
    p.add_argument("--b-dir", required=True)
    args = p.parse_args()
    a, b = Path(args.a_dir), Path(args.b_dir)
    ok = True
    for name in ("merged.bam", "unmerged.bam"):
        pa, pb = a / name, b / name
        if pa.exists() or pb.exists():
            ok &= compare_bam(pa, pb, name)
    ok &= compare_parquet(
        a / "B.disagreements.parquet", b / "B.disagreements.parquet", "disagreements"
    )
    ok &= compare_parquet(a / "B.stats.parquet", b / "B.stats.parquet", "read-stats")
    ta = (a / "done.txt").read_text() if (a / "done.txt").exists() else ""
    tb = (b / "done.txt").read_text() if (b / "done.txt").exists() else ""
    da = re.sub(r"^\d+(\.\d+)?s.*$", "", ta)
    db = re.sub(r"^\d+(\.\d+)?s.*$", "", tb)
    if da == db and da:
        print("done line: identical")
    else:
        print(f"done line: MISMATCH\n  A: {da!r}\n  B: {db!r}")
        ok = False
    print("EQUIVALENCE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
