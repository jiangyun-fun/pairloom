"""Validated run configuration shared by the CLI, pipeline, and workers."""

from __future__ import annotations

import re
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Literal

from pairloom.exceptions import ConfigurationError

Mode = Literal["adjust", "merge"]
PairStrategy = Literal["pair-aware", "max-overlap"]
CandidateSource = Literal["bam", "reference"]
ConflictPolicy = Literal["quality", "strict"]
LibraryOrientation = Literal["fr", "rf", "ff", "rr", "any"]

MODES: tuple[str, ...] = ("adjust", "merge")
PAIR_STRATEGIES: tuple[str, ...] = ("pair-aware", "max-overlap")
CANDIDATE_SOURCES: tuple[str, ...] = ("bam", "reference")
CONFLICT_POLICIES: tuple[str, ...] = ("quality", "strict")
LIBRARY_ORIENTATIONS: tuple[str, ...] = ("fr", "rf", "ff", "rr", "any")

# Tags that are recomputed or invalidated when PairLoom rewrites a record.
STALE_TAGS_REMOVED: tuple[str, ...] = ("AS", "XS", "SA", "MC", "MQ", "XA")

#: Fixed diagnostic tag marking records whose placement came from joint
#: low-MAPQ pair rescue (``ZR:i:1``). Lives here so configuration
#: validation can reserve it against ``--status-tag`` without an import
#: cycle (classify imports config, never the reverse).
RESCUE_TAG = "ZR"
RESCUE_VALUE = 1

#: Two-character SAM tag: first character a letter, second alphanumeric.
_TAG_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]$")

#: Parquet compressions supported for statistical outputs.
PARQUET_COMPRESSIONS: tuple[str, ...] = ("zstd", "snappy", "gzip", "lz4", "none")


@dataclass(frozen=True, slots=True)
class PairLoomConfig:
    """All runtime parameters, fully validated at construction time."""

    input_path: Path
    reference_path: Path
    mode: Mode
    disagreement_parquet: Path
    read_stats_parquet: Path
    adjusted_bam: Path | None = None
    merged_bam: Path | None = None
    unmerged_bam: Path | None = None
    threads: int = 1
    status_tag: str = "ZC"
    min_base_quality: int = 20
    min_mapq: int = 20
    #: 0 explicitly allows adjacent (zero-overlap) joining; an uncovered
    #: fragment base still rejects the merge (``fragment_gap``).
    min_overlap: int = 0
    conflict_policy: ConflictPolicy = "quality"
    min_quality_difference: int = 5
    max_consensus_quality: int = 60
    library_orientation: LibraryOrientation = "fr"
    min_fragment_length: int = 0
    max_fragment_length: int = 1000
    candidate_source: CandidateSource = "reference"
    max_candidates_per_read: int = 10
    realign_window: int = 500
    max_pair_evaluations: int = 10000
    pair_strategy: PairStrategy = "pair-aware"
    alignment_score_tolerance: int = 0
    parquet_compression: str = "zstd"
    parquet_batch_size: int = 65536
    sort_output: bool = False
    index_output: bool = False
    force_mapq60_exact: bool = False
    #: Joint low-MAPQ pair rescue across references (default on; the CLI
    #: exposes it as the ``--no-rescue`` opt-out).
    rescue_enabled: bool = True
    force: bool = False
    debug: bool = False
    # Internal tuning knobs (documented; not CLI-facing).
    max_aggregation_records: int = 1_000_000
    chunk_size: int = 256
    reference_cache_contigs: int = 4
    reference_cache_bases: int = 4_000_000

    def __post_init__(self) -> None:
        for name in ("mode",):
            value = getattr(self, name)
            if value not in MODES:
                raise ConfigurationError(f"invalid {name} {value!r}; expected one of {MODES}")
        if self.pair_strategy not in PAIR_STRATEGIES:
            raise ConfigurationError(
                f"invalid pair_strategy {self.pair_strategy!r}; expected one of {PAIR_STRATEGIES}"
            )
        if self.candidate_source not in CANDIDATE_SOURCES:
            raise ConfigurationError(
                f"invalid candidate_source {self.candidate_source!r}; "
                f"expected one of {CANDIDATE_SOURCES}"
            )
        if self.conflict_policy not in CONFLICT_POLICIES:
            raise ConfigurationError(
                f"invalid conflict_policy {self.conflict_policy!r}; "
                f"expected one of {CONFLICT_POLICIES}"
            )
        if self.library_orientation not in LIBRARY_ORIENTATIONS:
            raise ConfigurationError(
                f"invalid library_orientation {self.library_orientation!r}; "
                f"expected one of {LIBRARY_ORIENTATIONS}"
            )
        if self.parquet_compression not in PARQUET_COMPRESSIONS:
            raise ConfigurationError(
                f"invalid parquet_compression {self.parquet_compression!r}; "
                f"expected one of {PARQUET_COMPRESSIONS}"
            )
        if not _TAG_RE.match(self.status_tag):
            raise ConfigurationError(
                f"invalid status_tag {self.status_tag!r}: must be exactly two SAM tag "
                "characters ([A-Za-z][A-Za-z0-9])"
            )
        if self.status_tag == RESCUE_TAG:
            raise ConfigurationError(
                f"status_tag {RESCUE_TAG!r} is reserved for the joint-rescue marker "
                "(ZR:i:1); choose another --status-tag"
            )
        _require_positive(self.threads, "threads")
        _require_non_negative(self.min_overlap, "min_overlap")
        _require_positive(self.max_candidates_per_read, "max_candidates_per_read")
        _require_positive(self.realign_window, "realign_window")
        _require_positive(self.max_pair_evaluations, "max_pair_evaluations")
        _require_positive(self.parquet_batch_size, "parquet_batch_size")
        _require_positive(self.chunk_size, "chunk_size")
        _require_non_negative(self.min_base_quality, "min_base_quality")
        _require_non_negative(self.min_mapq, "min_mapq")
        _require_non_negative(self.min_quality_difference, "min_quality_difference")
        _require_non_negative(self.alignment_score_tolerance, "alignment_score_tolerance")
        if self.max_consensus_quality <= 0:
            raise ConfigurationError("max_consensus_quality must be > 0")
        if not 0 <= self.min_base_quality <= 93:
            raise ConfigurationError("min_base_quality must be within [0, 93]")
        if not 0 <= self.max_consensus_quality <= 93:
            raise ConfigurationError("max_consensus_quality must be within [0, 93]")
        if self.min_fragment_length < 0 or self.max_fragment_length < self.min_fragment_length:
            raise ConfigurationError(
                "fragment length bounds must satisfy 0 <= min <= max "
                f"(got min={self.min_fragment_length}, max={self.max_fragment_length})"
            )
        if self.alignment_score_tolerance < 0:
            raise ConfigurationError("alignment_score_tolerance must be >= 0")
        # Mode-specific output wiring.
        if self.mode == "adjust":
            if self.merged_bam is not None or self.unmerged_bam is not None:
                raise ConfigurationError(
                    "--merged-bam/--unmerged-bam are only valid with --mode merge"
                )
            if self.adjusted_bam is None:
                raise ConfigurationError("--mode adjust requires --adjusted-bam")
        else:
            if self.adjusted_bam is not None:
                raise ConfigurationError("--adjusted-bam is only valid with --mode adjust")
            if self.merged_bam is None or self.unmerged_bam is None:
                raise ConfigurationError("--mode merge requires --merged-bam and --unmerged-bam")

    def mode_outputs(self) -> dict[str, Path]:
        """Explicit output paths required by the selected mode, keyed by role."""
        outputs: dict[str, Path] = {
            "disagreement_parquet": self.disagreement_parquet,
            "read_stats_parquet": self.read_stats_parquet,
        }
        if self.mode == "adjust":
            assert self.adjusted_bam is not None  # noqa: S101 - validated above
            outputs["adjusted_bam"] = self.adjusted_bam
        else:
            assert self.merged_bam is not None and self.unmerged_bam is not None  # noqa: S101
            outputs["merged_bam"] = self.merged_bam
            outputs["unmerged_bam"] = self.unmerged_bam
        return outputs

    def debug_dict(self) -> dict[str, object]:
        """Plain-dict form for log lines and tests."""
        return {f.name: _plain(getattr(self, f.name)) for f in fields(self)}


def _require_positive(value: int, name: str) -> None:
    if value <= 0:
        raise ConfigurationError(f"{name} must be > 0 (got {value})")


def _require_non_negative(value: int, name: str) -> None:
    if value < 0:
        raise ConfigurationError(f"{name} must be >= 0 (got {value})")


def _plain(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    return value
