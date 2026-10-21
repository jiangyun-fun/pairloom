#!/usr/bin/env python
"""Reproducible PairLoom benchmark.

Generates a synthetic paired-end dataset (reference with repeat islands,
BAM with overlapping pairs, disagreements, and unmapped templates), then
measures end-to-end wall time, per-stage profile, and peak RSS for several
worker counts. Results are written as Parquet plus a Markdown table.

Usage:
    python benchmarks/benchmark.py --pairs 20000 --threads 1,2,4 \
        --output-dir benchmarks/results
"""

from __future__ import annotations

import argparse
import cProfile
import pstats
import random
import resource
import sys
import time
from pathlib import Path

import pysam

COMP = str.maketrans("ACGT", "TGCA")


def build_reference(path: Path, length: int = 2_000_000, seed: int = 42) -> str:
    """Reference with a 1 kb period-4 repeat island every 50 kb."""
    rng = random.Random(seed)
    parts: list[str] = []
    pos = 0
    unit = "AACCGGTT"
    while pos < length:
        block = min(50_000, length - pos)
        chunk = [rng.choice("ACGT") for _ in range(block)]
        parts.append("".join(chunk))
        pos += block
        if pos < length:
            island = min(1000, length - pos)
            parts.append(unit * (island // len(unit)))
            pos += island
    seq = "".join(parts)[:length]
    with open(path, "w") as handle:
        handle.write(">chr1\n")
        for i in range(0, len(seq), 60):
            handle.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return seq


def build_bam(bam_path: Path, ref: str, pairs: int, seed: int = 7) -> None:
    """Query-name-sorted BAM: 20% repeat-island pairs, 2% unmapped, ~0.2% SNVs."""
    rng = random.Random(seed)
    length = len(ref)
    islands = []
    pos = 50_000
    while pos < length:
        islands.append((pos, pos + 1000))
        pos += 51_000
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": length}],
    }
    read_len = 100
    hdr = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(str(bam_path), "wb", header=hdr) as out:
        for i in range(pairs):
            name = f"read{i:08d}"
            if rng.random() < 0.02:
                for flag in (0x1 | 0x40 | 0x4, 0x1 | 0x80 | 0x4):
                    rec = pysam.AlignedSegment(hdr)
                    rec.query_name = name
                    rec.flag = flag
                    rec.reference_id = -1
                    rec.query_sequence = "".join(rng.choice("ACGT") for _ in range(read_len))
                    rec.query_qualities = pysam.qualitystring_to_array("I" * read_len)
                    out.write(rec)
                continue
            if rng.random() < 0.20 and islands:
                start, end = rng.choice(islands)
                pos1 = rng.randint(start, max(start, end - 350))
            else:
                pos1 = rng.randint(0, length - 400)
            insert = rng.randint(120, 280)
            pos2 = pos1 + insert - read_len
            if pos2 + read_len > length:
                continue
            r1_seq = list(ref[pos1 : pos1 + read_len])
            r2_seq = list(ref[pos2 : pos2 + read_len])
            for seq, _offset in ((r1_seq, pos1), (r2_seq, pos2)):
                for k in range(read_len):
                    if rng.random() < 0.002:
                        seq[k] = rng.choice([b for b in "ACGT" if b != seq[k]])
            r1 = _record(hdr, name, 0x1 | 0x40 | 0x2, pos1, "".join(r1_seq), pos2)
            r2 = _record(hdr, name, 0x1 | 0x80 | 0x10 | 0x20 | 0x2, pos2, "".join(r2_seq), pos1)
            out.write(r1)
            out.write(r2)


def _record(
    hdr: pysam.AlignmentHeader, name: str, flag: int, pos: int, seq: str, mpos: int
) -> pysam.AlignedSegment:
    rec = pysam.AlignedSegment(hdr)
    rec.query_name = name
    rec.flag = flag
    rec.reference_id = 0
    rec.reference_start = pos
    rec.next_reference_id = 0
    rec.next_reference_start = mpos
    rec.cigarstring = f"{len(seq)}M"
    rec.query_sequence = seq
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    rec.mapping_quality = 60
    return rec


def peak_rss_mb() -> float:
    self_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return max(self_kb, child_kb) / 1024.0


def run_once(
    bam: Path, ref: Path, out_dir: Path, threads: int, chunk: int
) -> tuple[float, float, int]:
    from pairloom.cli import main

    out = out_dir / f"t{threads}"
    out.mkdir(parents=True, exist_ok=True)
    argv = [
        "--input",
        str(bam),
        "--reference",
        str(ref),
        "--mode",
        "merge",
        "--merged-bam",
        str(out / "merged.bam"),
        "--unmerged-bam",
        str(out / "unmerged.bam"),
        "--disagreement-parquet",
        str(out / "disagreements.parquet"),
        "--read-stats-parquet",
        str(out / "stats.parquet"),
        "--threads",
        str(threads),
        "--chunk-size",
        str(chunk),
        "--force",
    ]
    start = time.perf_counter()
    code = main(argv)
    elapsed = time.perf_counter() - start
    if code != 0:
        raise SystemExit(f"benchmark run failed with exit code {code} (threads={threads})")
    import pyarrow.parquet as pq

    total = pq.read_table(str(out / "stats.parquet")).column("count")[0].as_py()
    return elapsed, peak_rss_mb(), total


_STAGE_FUNCTIONS = {
    "bam_read_and_grouping": ("iter_templates", "write"),
    "candidate_discovery": ("discover_candidates",),
    "pair_optimization": ("select_pair",),
    "overlap_analysis": ("analyze_overlap", "build_read_structure"),
    "consensus_merge": ("try_merge",),
    "parquet_aggregation": ("add_observation", "finalize"),
    "record_processing": ("process_template", "fromstring", "to_string"),
}


def stage_profile(bam: Path, ref: Path, out_dir: Path) -> dict[str, float]:
    """Single-thread cProfile run bucketed into pipeline stages (seconds)."""
    from pairloom.cli import main

    out = out_dir / "profile"
    out.mkdir(parents=True, exist_ok=True)
    argv = [
        "--input",
        str(bam),
        "--reference",
        str(ref),
        "--mode",
        "merge",
        "--merged-bam",
        str(out / "merged.bam"),
        "--unmerged-bam",
        str(out / "unmerged.bam"),
        "--disagreement-parquet",
        str(out / "disagreements.parquet"),
        "--read-stats-parquet",
        str(out / "stats.parquet"),
        "--force",
    ]
    profiler = cProfile.Profile()
    profiler.enable()
    code = main(argv)
    profiler.disable()
    if code != 0:
        raise SystemExit("profiled run failed")
    stats = pstats.Stats(profiler)
    # stats.stats maps (file, line, funcname) -> (cc, nc, total, cumulative, callers)
    cumulative_by_func: dict[str, float] = {}
    for (_file, _line, funcname), (_cc, _nc, _tt, cumulative, _callers) in stats.stats.items():
        cumulative_by_func[funcname] = max(cumulative_by_func.get(funcname, 0.0), cumulative)
    stage_seconds: dict[str, float] = {}
    for stage, funcs in _STAGE_FUNCTIONS.items():
        total = max((cumulative_by_func.get(f, 0.0) for f in funcs), default=0.0)
        stage_seconds[stage] = round(total, 3)
    return stage_seconds


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=20000)
    parser.add_argument("--threads", default="1,2,4")
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--reference-length", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("benchmarks/results"))
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = args.output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    ref_path = data_dir / "bench_ref.fa"
    bam_path = data_dir / "bench_input.bam"

    print(f"[bench] building reference ({args.reference_length:,} bp) ...")
    t0 = time.perf_counter()
    ref = build_reference(ref_path, args.reference_length, args.seed)
    print(f"[bench] reference built in {time.perf_counter() - t0:.1f}s")

    print(f"[bench] building BAM ({args.pairs:,} pairs) ...")
    t0 = time.perf_counter()
    build_bam(bam_path, ref, args.pairs, args.seed + 1)
    print(f"[bench] BAM built in {time.perf_counter() - t0:.1f}s")

    thread_list = [int(t) for t in args.threads.split(",")]
    rows: list[dict] = []
    for threads in thread_list:
        elapsed, rss_mb, total = run_once(
            bam_path, ref_path, args.output_dir, threads, args.chunk_size
        )
        pairs_per_s = total / elapsed if elapsed else 0.0
        print(
            f"[bench] threads={threads}: {elapsed:.1f}s, "
            f"{pairs_per_s:,.0f} pairs/s, peak RSS {rss_mb:.0f} MiB, "
            f"templates={total:,}"
        )
        rows.append(
            {
                "threads": threads,
                "seconds": round(elapsed, 2),
                "pairs_per_second": round(pairs_per_s, 1),
                "peak_rss_mib": round(rss_mb, 1),
                "templates": total,
            }
        )

    print("[bench] profiling stages (single thread) ...")
    stages = stage_profile(bam_path, ref_path, args.output_dir)
    for stage, seconds in stages.items():
        print(f"[bench]   {stage}: {seconds:.2f}s")

    import pyarrow as pa
    import pyarrow.parquet as pq

    table = pa.table(
        {
            "threads": pa.array([r["threads"] for r in rows], pa.int64()),
            "seconds": pa.array([r["seconds"] for r in rows], pa.float64()),
            "pairs_per_second": pa.array([r["pairs_per_second"] for r in rows], pa.float64()),
            "peak_rss_mib": pa.array([r["peak_rss_mib"] for r in rows], pa.float64()),
            "templates": pa.array([r["templates"] for r in rows], pa.int64()),
        }
    )
    out_parquet = args.output_dir / "benchmark.parquet"
    pq.write_table(table, str(out_parquet), compression="zstd")

    lines = [
        "# PairLoom benchmark results",
        "",
        f"- pairs: {args.pairs:,}, reference: {args.reference_length:,} bp, "
        f"chunk size: {args.chunk_size}",
        "- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs",
        "",
        "## Scaling",
        "",
        "| threads | seconds | pairs/s | peak RSS (MiB) |",
        "| --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            f"| {r['threads']} | {r['seconds']:.1f} | {r['pairs_per_second']:,.0f} | "
            f"{r['peak_rss_mib']:.0f} |"
        )
    lines += [
        "",
        "## Stage profile (single thread, cumulative seconds)",
        "",
        "| stage | seconds |",
        "| --- | --- |",
    ]
    for stage, seconds in stages.items():
        lines.append(f"| {stage} | {seconds:.2f} |")
    lines.append("")
    (args.output_dir / "RESULTS.md").write_text("\n".join(lines))
    print(f"[bench] wrote {out_parquet} and {args.output_dir / 'RESULTS.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
