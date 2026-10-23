"""``pairloom`` command-line interface.

Exit codes: 0 success, 2 usage/validation errors, 1 processing failures.
Errors are concise by default; ``--debug`` prints full tracebacks.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from pairloom import __version__
from pairloom.config import (
    CANDIDATE_SOURCES,
    CONFLICT_POLICIES,
    LIBRARY_ORIENTATIONS,
    MODES,
    PAIR_STRATEGIES,
    PARQUET_COMPRESSIONS,
    PairLoomConfig,
)
from pairloom.exceptions import ConfigurationError, PairLoomError
from pairloom.pipeline import run_pipeline
from pairloom.validate import validate_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pairloom",
        description=(
            "Pair-aware BAM realignment, overlap analysis, and paired-end read "
            "merging. Every output file is set through its own explicit path "
            "argument; no filenames are inferred."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input", required=True, help="input paired-end BAM")
    parser.add_argument(
        "--reference", required=True, help="indexed reference FASTA (.fai required)"
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=MODES,
        help="adjust: adjust/tag alignments only; merge: also emit consensus reads",
    )
    out = parser.add_argument_group("explicit outputs")
    out.add_argument("--adjusted-bam", help="adjusted unmerged BAM (mode adjust)")
    out.add_argument("--merged-bam", help="consensus reads BAM (mode merge)")
    out.add_argument("--unmerged-bam", help="unmergeable pairs BAM (mode merge)")
    out.add_argument("--disagreement-parquet", required=True, help="overlap disagreement table")
    out.add_argument("--read-stats-parquet", required=True, help="classification counts table")

    pair = parser.add_argument_group("pair optimization")
    pair.add_argument(
        "--candidate-source",
        choices=CANDIDATE_SOURCES,
        default="reference",
        help="placements from BAM only, or BAM plus windowed reference sliding (default)",
    )
    pair.add_argument(
        "--max-candidates-per-read",
        type=int,
        default=10,
        help="candidate placement cap per read (default 10)",
    )
    pair.add_argument(
        "--realign-window",
        type=int,
        default=500,
        help=(
            "max reference window for candidate sliding per side (default 500 bp); "
            "never bounds joint-rescue uniqueness, which searches the whole reference"
        ),
    )
    pair.add_argument(
        "--no-rescue",
        action="store_true",
        help=(
            "disable joint low-MAPQ pair rescue (default: enabled; see the README's rescue section)"
        ),
    )
    pair.add_argument(
        "--max-pair-evaluations",
        type=int,
        default=10000,
        help="max candidate-pair evaluations per template (default 10000)",
    )
    pair.add_argument(
        "--pair-strategy",
        choices=PAIR_STRATEGIES,
        default="pair-aware",
        help="evidence-first with overlap tie-break (default) or aggressive overlap",
    )
    pair.add_argument(
        "--alignment-score-tolerance",
        type=int,
        default=0,
        help="score window treating candidate pairs as near-equivalent (default 0)",
    )
    pair.add_argument(
        "--library-orientation",
        choices=LIBRARY_ORIENTATIONS,
        default="fr",
        help="expected mate orientation model (default fr)",
    )
    pair.add_argument("--min-fragment-length", type=int, default=0, help="default 0")
    pair.add_argument("--max-fragment-length", type=int, default=1000, help="default 1000")

    qual = parser.add_argument_group("quality and merging")
    qual.add_argument("--min-base-quality", type=int, default=20, help="default 20")
    qual.add_argument("--min-mapq", type=int, default=20, help="default 20")
    qual.add_argument(
        "--min-overlap",
        type=int,
        default=0,
        help=(
            "min aligned overlap in reference bases (default 0; 0 explicitly "
            "allows adjacent joining of mates covering a continuous fragment)"
        ),
    )
    qual.add_argument(
        "--conflict-policy",
        choices=CONFLICT_POLICIES,
        default="quality",
        help="resolve conflicts by quality (default) or reject pair (strict)",
    )
    qual.add_argument(
        "--min-quality-difference",
        type=int,
        default=5,
        help="phred gap that resolves a conflicting base (default 5)",
    )
    qual.add_argument(
        "--max-consensus-quality",
        type=int,
        default=60,
        help="consensus quality cap (default 60)",
    )
    qual.add_argument(
        "--force-mapq60-exact",
        action="store_true",
        help=(
            "compatibility mode: MAPQ 60 for every exact match. WARNING: this "
            "overstates mapping confidence for repeated reads."
        ),
    )
    qual.add_argument(
        "--status-tag", default="ZC", help="two-character classification tag (default ZC)"
    )

    io_grp = parser.add_argument_group("output control")
    io_grp.add_argument("--threads", type=int, default=1, help="worker processes (default 1)")
    io_grp.add_argument(
        "--parquet-compression",
        choices=PARQUET_COMPRESSIONS,
        default="zstd",
        help="statistical output compression (default zstd)",
    )
    io_grp.add_argument(
        "--parquet-batch-size",
        type=int,
        default=65536,
        help="rows per Parquet batch (default 65536)",
    )
    io_grp.add_argument(
        "--max-aggregation-records",
        type=int,
        default=1_000_000,
        help="in-memory disagreement positions before spilling (default 1e6)",
    )
    io_grp.add_argument(
        "--chunk-size",
        type=int,
        default=8192,
        help="templates per work unit in parallel runs (default 8192)",
    )
    io_grp.add_argument(
        "--sort-output",
        action="store_true",
        help="coordinate-sort BAM outputs (default: template order)",
    )
    io_grp.add_argument(
        "--index-output",
        action="store_true",
        help="index BAM outputs (implies --sort-output)",
    )
    io_grp.add_argument("--force", action="store_true", help="overwrite existing output files")
    io_grp.add_argument("--debug", action="store_true", help="full tracebacks on error")
    return parser


def build_config(args: argparse.Namespace) -> PairLoomConfig:
    sort_output = args.sort_output or args.index_output
    if args.index_output and not args.sort_output:
        print(
            "pairloom: note: --index-output implies coordinate sorting of BAM outputs",
            file=sys.stderr,
        )
    return PairLoomConfig(
        input_path=Path(args.input),
        reference_path=Path(args.reference),
        mode=args.mode,
        adjusted_bam=Path(args.adjusted_bam) if args.adjusted_bam else None,
        merged_bam=Path(args.merged_bam) if args.merged_bam else None,
        unmerged_bam=Path(args.unmerged_bam) if args.unmerged_bam else None,
        disagreement_parquet=Path(args.disagreement_parquet),
        read_stats_parquet=Path(args.read_stats_parquet),
        threads=args.threads,
        status_tag=args.status_tag,
        min_base_quality=args.min_base_quality,
        min_mapq=args.min_mapq,
        min_overlap=args.min_overlap,
        conflict_policy=args.conflict_policy,
        min_quality_difference=args.min_quality_difference,
        max_consensus_quality=args.max_consensus_quality,
        library_orientation=args.library_orientation,
        min_fragment_length=args.min_fragment_length,
        max_fragment_length=args.max_fragment_length,
        candidate_source=args.candidate_source,
        max_candidates_per_read=args.max_candidates_per_read,
        realign_window=args.realign_window,
        max_pair_evaluations=args.max_pair_evaluations,
        pair_strategy=args.pair_strategy,
        alignment_score_tolerance=args.alignment_score_tolerance,
        parquet_compression=args.parquet_compression,
        parquet_batch_size=args.parquet_batch_size,
        sort_output=sort_output,
        index_output=args.index_output,
        force_mapq60_exact=args.force_mapq60_exact,
        rescue_enabled=not args.no_rescue,
        force=args.force,
        debug=args.debug,
        max_aggregation_records=args.max_aggregation_records,
        chunk_size=args.chunk_size,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg = build_config(args)
        validate_paths(cfg)
    except ConfigurationError as exc:
        print(f"pairloom: error: {exc}", file=sys.stderr)
        return 2
    if cfg.force_mapq60_exact:
        print(
            "pairloom: WARNING: --force-mapq60-exact assigns MAPQ 60 to every "
            "exact match; this overstates mapping confidence for reads in "
            "repeated regions.",
            file=sys.stderr,
        )
    try:
        run_pipeline(cfg)
    except PairLoomError as exc:
        if cfg.debug:
            traceback.print_exc()
        else:
            print(f"pairloom: error: {exc}", file=sys.stderr)
            print("re-run with --debug for a full traceback", file=sys.stderr)
        return 1
    except Exception as exc:
        # Unexpected failures are never swallowed: concise message, nonzero
        # exit, full traceback under --debug.
        if cfg.debug:
            traceback.print_exc()
        else:
            print(
                f"pairloom: internal error ({type(exc).__name__}): {exc}",
                file=sys.stderr,
            )
            print("re-run with --debug for a full traceback", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
