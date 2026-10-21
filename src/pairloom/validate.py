"""CLI-level validation of input and output paths.

Failures are loud, actionable, and happen before any processing starts:
missing inputs, missing FASTA indexes, unwritable directories, output paths
that collide with each other or with the inputs, and existing outputs
without ``--force``.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from pairloom.config import PairLoomConfig
from pairloom.exceptions import ConfigurationError


def validate_paths(cfg: PairLoomConfig) -> None:
    """Validate all input/output paths; raise ConfigurationError on any issue."""
    _check_readable_file(cfg.input_path, "--input")
    _check_reference(cfg.reference_path)
    outputs = cfg.mode_outputs()
    _check_distinct_outputs(cfg, outputs)
    for role, path in outputs.items():
        _check_output_path(path, role, cfg.force)


def _check_readable_file(path: Path, flag: str) -> None:
    if not path.exists():
        raise ConfigurationError(f"{flag} file does not exist: {path}")
    if not path.is_file():
        raise ConfigurationError(f"{flag} is not a regular file: {path}")
    if not os.access(path, os.R_OK):
        raise ConfigurationError(f"{flag} file is not readable: {path}")


def _check_reference(path: Path) -> None:
    _check_readable_file(path, "--reference")
    fai = Path(str(path) + ".fai")
    if not fai.exists():
        raise ConfigurationError(
            f"reference index {fai} is missing; run 'samtools faidx {path}' first"
        )


def _check_output_path(path: Path, role: str, force: bool) -> None:
    if str(path) in ("", "."):
        raise ConfigurationError(f"{role} path is empty")
    parent = path.parent if str(path.parent) else Path(".")
    try:
        parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ConfigurationError(
            f"cannot create output directory {parent} for {role}: {exc}"
        ) from exc
    try:
        probe = tempfile.NamedTemporaryFile(dir=str(parent), prefix=".pairloom_probe_")  # noqa: SIM115
        probe.close()
    except OSError as exc:
        raise ConfigurationError(
            f"output directory for {role} is not writable: {parent} ({exc})"
        ) from exc
    if path.exists() and not force:
        raise ConfigurationError(f"{role} output already exists: {path}. Use --force to overwrite.")


def _check_distinct_outputs(cfg: PairLoomConfig, outputs: dict[str, Path]) -> None:
    """Reject duplicate output paths and outputs that collide with inputs."""
    resolved: dict[str, str] = {}
    existing = [cfg.input_path, cfg.reference_path, Path(str(cfg.reference_path) + ".fai")]
    for role, path in outputs.items():
        norm = os.path.normpath(os.path.abspath(path))
        for other_role, other_norm in resolved.items():
            if norm == other_norm:
                raise ConfigurationError(
                    f"{role} and {other_role} resolve to the same file: {norm}. "
                    "Each output must have its own explicit path."
                )
        resolved[role] = norm
        for src in existing:
            src_norm = os.path.normpath(os.path.abspath(src))
            if norm == src_norm:
                raise ConfigurationError(f"{role} output {path} would overwrite input {src}")
    # samefile() also catches hardlinks and symlinks among existing paths.
    all_paths = list(outputs.values())
    for i, a in enumerate(all_paths):
        for b in all_paths[i + 1 :]:
            if a.exists() and b.exists() and os.path.samefile(a, b):
                raise ConfigurationError(f"output paths {a} and {b} refer to the same file")
        for src in existing:
            if a.exists() and os.path.samefile(a, src):
                raise ConfigurationError(f"output path {a} refers to input file {src}")
