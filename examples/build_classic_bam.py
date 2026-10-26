#!/usr/bin/env python
"""Build the classic gallery BAM: one template per documented outcome.

The gallery pairs a hand-built reference (classic.fa) with a name-sorted
alignment BAM (classic.bam) in which every template demonstrates exactly
one PairLoom outcome: each success path, each classification, and each
merge-failure reason that default settings can surface. Policy-gated
reasons (non_primary_records, mapq_below_min, strict_conflict,
insufficient_overlap) appear when the same BAM runs with their flags —
see examples/gallery/OUTCOMES.md.

Sequences are generated deterministically from --seed.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pysam

READ = 20  # standard gallery read length


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True, help="output directory (created)")
    p.add_argument("--seed", type=int, default=20261026, help="reference seed")
    return p.parse_args()


def rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=n))


def parse_cigar(cigar: str) -> list[tuple[int, int]]:
    import re

    ops = []
    for num, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        ops.append(("MIDNSHP=X".index(op), int(num)))
    return ops


class Builder:
    def __init__(self, contigs: dict[str, str]) -> None:
        self.contigs = contigs
        self.header = {
            "HD": {"VN": "1.6", "SO": "queryname"},
            "SQ": [{"SN": name, "LN": len(seq)} for name, seq in contigs.items()],
            "RG": [{"ID": "rgA"}, {"ID": "rgB"}],
        }
        self.records: list[pysam.AlignedSegment] = []

    def rec(
        self,
        name: str,
        flag: int,
        contig: str | None,
        pos: int,
        cigar: str,
        seq: str,
        mapq: int = 60,
        rg: str | None = None,
        quals: bool = True,
    ) -> pysam.AlignedSegment:
        r = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(self.header))
        r.query_name = name
        r.flag = flag
        if contig is None:
            r.reference_id = -1
            r.next_reference_id = -1
        else:
            r.reference_id = list(self.contigs).index(contig)
            r.next_reference_id = r.reference_id
        r.reference_start = pos
        r.next_reference_start = pos
        r.cigartuples = parse_cigar(cigar) if cigar != "*" else None
        r.query_sequence = seq
        if quals:
            r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        if rg:
            r.set_tag("RG", rg, value_type="Z")
        r.mapping_quality = mapq
        return r

    def add(self, records: list[pysam.AlignedSegment]) -> None:
        self.records.extend(records)

    def write(self, bam_path: Path) -> None:
        with pysam.AlignmentFile(str(bam_path), "wb", header=self.header) as out:
            for r in self.records:
                out.write(r)


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # ---- reference ----------------------------------------------------
    ctga = rand_seq(2000, rng)
    ctgb = rand_seq(2000, rng)
    # jointly-unique rescue blocks: B1 on ctgA AND ctgB, B2 only on ctgA
    b1 = rand_seq(40, rng)
    b2 = rand_seq(40, rng)
    ctga = ctga[:700] + b1 + b2 + ctga[780:]
    ctgb = ctgb[:1300] + b1 + ctgb[1340:]
    # duplicate contig pair for rescue_ambiguous: S1/S2 embedded once each
    s1 = rand_seq(40, rng)
    s2 = rand_seq(40, rng)
    ctgdup = rand_seq(60, rng) + s1 + rand_seq(100, rng) + s2 + rand_seq(80, rng)
    # repeat contig for rescue_search_truncated: 11 exact copies of a unit
    unit = rand_seq(40, rng)
    ctgrep = unit * 11
    contigs = {
        "ctgA": ctga,
        "ctgB": ctgb,
        "ctgDup1": ctgdup,
        "ctgDup2": ctgdup,
        "ctgRep": ctgrep,
    }
    fa = out / "classic.fa"
    with open(fa, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n{seq}\n")

    b = Builder(contigs)
    R = READ

    def pair(name, pos1, pos2, contig="ctgA", mapq=60):
        # Stored SEQ is reference-oriented for BOTH strands (PairLoom/bwa
        # convention); strand lives in the flags. R2 carries 0x10.
        recs = [
            b.rec(name, 0x1 | 0x40 | 0x2, contig, pos1, f"{R}M", contigs[contig][pos1 : pos1 + R], mapq=mapq),
            b.rec(name, 0x1 | 0x80 | 0x10 | 0x2, contig, pos2, f"{R}M", contigs[contig][pos2 : pos2 + R], mapq=mapq),
        ]
        b.add(recs)
        return recs

    # ---- successes ----------------------------------------------------
    pair("ok_overlap", 100, 105)                       # 15 bp overlap -> merge
    pair("ok_adjacent", 200, 220)                      # zero-overlap adjacent join
    # equal-quality conflict in overlap -> N consensus base (quality policy)
    nq = pair("ok_equal_qual_N", 300, 305)
    stored = list(nq[1].query_sequence)                # reference-oriented storage
    stored[10] = "A" if stored[10] != "A" else "C"     # conflicts R1 at equal qual
    nq[1].query_sequence = "".join(stored)
    nq[1].query_qualities = pysam.qualitystring_to_array("I" * len(stored))
    # rescued: individually ambiguous (B1 twice), jointly unique
    b.add([
        b.rec("ok_rescued", 0x1 | 0x40 | 0x2, "ctgA", 700, "40M", b1, mapq=0),
        b.rec("ok_rescued", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 740, "40M", b2, mapq=0),
    ])
    # non-primary evidence: concordant pair + supplementary -> ZN merge
    pair("ok_nonprimary", 400, 405)
    b.add([b.rec("ok_nonprimary", 0x1 | 0x40 | 0x800, "ctgB", 500, f"{R}M", ctgb[500 : 500 + R])])

    # ---- classifications ----------------------------------------------
    b.add([
        b.rec("cls_left_only", 0x1 | 0x40, "ctgA", 600, f"{R}M", ctga[600 : 620]),
        b.rec("cls_left_only", 0x1 | 0x80 | 0x4, None, -1, "*", ctga[600 : 620], mapq=0),
    ])
    b.add([
        b.rec("cls_right_only", 0x1 | 0x40 | 0x4, None, -1, "*", ctga[650 : 670], mapq=0),
        b.rec("cls_right_only", 0x1 | 0x80, "ctgA", 650, f"{R}M", ctga[650 : 670]),
    ])
    b.add([
        b.rec("cls_neither", 0x1 | 0x40 | 0x4, None, -1, "*", "ACGTACGTACGTACGTACGT", mapq=0),
        b.rec("cls_neither", 0x1 | 0x80 | 0x4, None, -1, "*", "ACGTACGTACGTACGTACGT", mapq=0),
    ])
    b.add([  # only non-primary records -> neither_aligned pass-through
        b.rec("cls_only_nonprimary", 0x1 | 0x40 | 0x100, "ctgA", 700, f"{R}M", ctga[700 : 720], mapq=0),
    ])

    # ---- failures under default settings -------------------------------
    b.add([
        b.rec("fail_cross_contig", 0x1 | 0x40 | 0x2, "ctgA", 900, f"{R}M", ctga[900 : 920]),
        b.rec("fail_cross_contig", 0x1 | 0x80 | 0x10 | 0x2, "ctgB", 900, f"{R}M", ctgb[900 : 920]),
    ])
    # chimeric: head-to-tail dimer shape (primary 20M10S, clips = other contig prefix)
    b.add([
        b.rec("fail_chimeric", 0x1 | 0x40, "ctgA", 0, f"{R}M10S", ctga[:R] + ctgb[:10]),
        b.rec("fail_chimeric", 0x1 | 0x80, "ctgB", 0, f"{R}M10S", ctgb[:R] + ctga[:10]),
    ])
    pair("fail_gap", 300, 480)                         # 160 bp gap, span <= 1000 -> fragment_gap
    pair("fail_frag_len", 100, 1600)                   # > max_fragment_length
    b.add([  # same-strand pair -> orientation_incompatible
        b.rec("fail_orientation", 0x1 | 0x40 | 0x2, "ctgA", 1200, f"{R}M", ctga[1200 : 1220]),
        b.rec("fail_orientation", 0x1 | 0x80 | 0x2, "ctgA", 1205, f"{R}M", ctga[1205 : 1225]),
    ])
    # indel conflict: R1 carries a 2 bp deletion across the overlap, R2 exact
    b.add([
        b.rec("fail_indel", 0x1 | 0x40 | 0x2, "ctgA", 100, "10M2D10M", ctga[100:110] + ctga[112:122]),
        b.rec("fail_indel", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 100, "20M", ctga[100:120]),
    ])
    # insertion conflict: both mates insert different bases at the same locus
    ins1 = ctga[100:110] + "GG" + ctga[110:120]
    ins2 = ctga[100:110] + "CC" + ctga[110:120]
    b.add([
        b.rec("fail_insertion", 0x1 | 0x40 | 0x2, "ctgA", 100, "10M2I10M", ins1),
        b.rec("fail_insertion", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 100, "10M2I10M", ins2),
    ])
    # hard-clipped primary
    b.add([
        b.rec("fail_hardclip", 0x1 | 0x40 | 0x2, "ctgA", 100, f"5H{R}M", ctga[100 : 100 + R]),
        b.rec("fail_hardclip", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 105, f"{R}M", ctga[105 : 125]),
    ])
    # missing qualities
    nq1 = b.rec("fail_noqual", 0x1 | 0x40 | 0x2, "ctgA", 100, f"{R}M", ctga[100 : 120], quals=False)
    nq2 = b.rec("fail_noqual", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 105, f"{R}M", ctga[105 : 125])
    nq2.query_qualities = None
    b.add([nq1, nq2])
    # reference skip (N operation)
    b.add([
        b.rec("fail_refskip", 0x1 | 0x40 | 0x2, "ctgA", 100, "10M5N10M", ctga[100:110] + ctga[115:125]),
        b.rec("fail_refskip", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 100, f"{R}M", ctga[100:120]),
    ])
    # read-group conflict
    b.add([
        b.rec("fail_rg", 0x1 | 0x40 | 0x2, "ctgA", 100, f"{R}M", ctga[100 : 120], rg="rgA"),
        b.rec("fail_rg", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 105, f"{R}M", ctga[105 : 125], rg="rgB"),
    ])
    # rescue failures (low original MAPQ)
    b.add([  # gapped alignment -> rescue_unsupported_alignment
        b.rec("fail_rescue_unsupported", 0x1 | 0x40 | 0x2, "ctgA", 100, "10M2D8M", ctga[100:110] + ctga[112:120], mapq=19),
        b.rec("fail_rescue_unsupported", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 100, "18M", ctga[100:118], mapq=19),
    ])
    b.add([  # blocks on both duplicate contigs -> rescue_ambiguous
        b.rec("fail_rescue_ambiguous", 0x1 | 0x40 | 0x2, "ctgDup1", 60, "40M", s1, mapq=0),
        b.rec("fail_rescue_ambiguous", 0x1 | 0x80 | 0x10 | 0x2, "ctgDup1", 200, "40M", s2, mapq=0),
    ])
    b.add([  # unique blocks, no compatible joint placement -> rescue_no_compatible_pair
        b.rec("fail_rescue_no_pair", 0x1 | 0x40 | 0x2, "ctgA", 0, f"{R}M", ctga[:R], mapq=0),
        b.rec("fail_rescue_no_pair", 0x1 | 0x80 | 0x10 | 0x2, "ctgB", 0, f"{R}M", ctgb[:R], mapq=0),
    ])
    b.add([  # 11 occurrences > cap 10 -> rescue_search_truncated
        b.rec("fail_rescue_truncated", 0x1 | 0x40 | 0x2, "ctgRep", 0, "40M", unit, mapq=0),
        b.rec("fail_rescue_truncated", 0x1 | 0x80 | 0x10 | 0x2, "ctgA", 1500, f"{R}M", ctga[1500:1520], mapq=0),
    ])

    b.write(out / "classic.bam")
    print(f"wrote {out/'classic.fa'} ({len(contigs)} contigs) and {out/'classic.bam'} ({len(b.records)} records, {len({r.query_name for r in b.records})} templates)")


if __name__ == "__main__":
    main()
