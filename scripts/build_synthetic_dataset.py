#!/usr/bin/env python
"""Build the small synthetic example dataset under data/synthetic/.

The dataset exercises every classification: a mergeable overlapping pair
(with one disagreement), a repeat-region pair whose placement is jointly
adjusted, a non-overlapping pair, a one-mate-unmapped pair, and a
both-unmapped pair. Expected outputs are documented in EXPECTED.md and are
regression-covered by the test suite.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import pysam


def main() -> int:
    root = Path(__file__).resolve().parent.parent / "data" / "synthetic"
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261021)

    def filler(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    # chr1: [120,132) and [236,248) hold "AACCAACCAACC"; [400,464) is the
    # long tandem repeat "AACCAACC"*8; [500,700) is the unique region.
    chr1 = (
        filler(120)
        + "AACCAACCAACC"
        + filler(104)
        + "AACCAACCAACC"
        + filler(152)
        + "AACCAACC" * 8
        + filler(36)
        + filler(200)
        + filler(300)
    )
    chr2 = filler(300)
    ref_path = root / "ref.fa"
    with open(ref_path, "w") as handle:
        for name, seq in (("chr1", chr1), ("chr2", chr2)):
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(ref_path))

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "queryname"},
            "SQ": [{"SN": "chr1", "LN": len(chr1)}, {"SN": "chr2", "LN": len(chr2)}],
            "RG": [{"ID": "rg1", "SM": "sample1"}],
        }
    )
    bam_path = root / "input.bam"
    records: list[pysam.AlignedSegment] = []

    def rec(qname, flag, pos, cigar, seq, mpos, mapq=60, tags=None):
        r = pysam.AlignedSegment(header)
        r.query_name = qname
        r.flag = flag
        r.reference_id = 0 if pos >= 0 else -1
        r.reference_start = pos
        r.next_reference_id = 0 if pos >= 0 else -1
        r.next_reference_start = mpos if mpos >= 0 else -1
        if cigar != "*":
            r.cigarstring = cigar
        r.query_sequence = seq
        r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        r.mapping_quality = mapq
        for tag, value in (tags or {}).items():
            r.set_tag(tag, value)
        return r

    unique = chr1[500:700]
    # 1. mergeable overlap with one disagreement at chr1:516.
    p = 515
    alt = "C" if unique[p - 500] != "C" else "G"
    r1 = list(unique[0:20])
    r1[p - 500] = alt
    records.append(rec("pair_mergeable", 0x1 | 0x40 | 0x2, 500, "20M", "".join(r1), 510))
    records.append(
        rec("pair_mergeable", 0x1 | 0x80 | 0x10 | 0x20 | 0x2, 510, "20M", unique[10:30], 500)
    )
    # 2. repeat pair, jointly adjusted to maximize overlap.
    records.append(rec("pair_repeat", 0x1 | 0x40 | 0x2, 404, "16M", chr1[404:420], 416))
    records.append(rec("pair_repeat", 0x1 | 0x80 | 0x10 | 0x20 | 0x2, 416, "16M", chr1[416:432], 404))
    # 3. non-overlapping pair.
    records.append(rec("pair_far", 0x1 | 0x40 | 0x2, 600, "20M", chr1[600:620], 650))
    records.append(
        rec("pair_far", 0x1 | 0x80 | 0x10 | 0x20 | 0x2, 650, "20M", chr1[650:670], 600)
    )
    # 4. one mate unmapped.
    records.append(rec("pair_half", 0x1 | 0x40 | 0x8, 620, "20M", chr1[620:640], 620))
    records.append(rec("pair_half", 0x1 | 0x80 | 0x4, -1, "*", "TTTTTTTTTT", -1))
    # 5. both mates unmapped.
    records.append(rec("pair_none", 0x1 | 0x40 | 0x4, -1, "*", "TTTTTTTTTT", -1))
    records.append(rec("pair_none", 0x1 | 0x80 | 0x4, -1, "*", "TTTTTTTTTT", -1))

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for record in records:
            out.write(record)

    print(f"wrote {ref_path}, {ref_path}.fai, {bam_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
