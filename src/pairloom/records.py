"""Template grouping and record-policy helpers.

Records are grouped by query name from a query-name-sorted stream. Primary
alignments drive pair processing; secondary, supplementary, and unmapped
records pass through unchanged (never silently discarded). A template with
more than two primary records, or with duplicate primary read1/read2
records, violates a documented invariant and fails loudly.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import pysam

from pairloom.exceptions import DataInvariantError, InputError

FLAG_PAIRED = 0x1
FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800


def open_input_bam(path_str: str) -> pysam.AlignmentFile:
    """Open a BAM for streaming, converting failures into InputError."""
    try:
        return pysam.AlignmentFile(path_str, "rb")
    except (OSError, ValueError) as exc:
        raise InputError(f"cannot open input BAM {path_str}: {exc}") from exc


def detect_sort_order(bam: pysam.AlignmentFile) -> str | None:
    header = bam.header.to_dict()
    hd = header.get("HD") or {}
    return hd.get("SO")


def name_sort_to_temp(input_path: str, temp_dir: str) -> str:
    """Query-name-sort a BAM into ``temp_dir``; return the temp path."""
    import os

    out = os.path.join(temp_dir, "qname_sorted.bam")
    try:
        pysam.sort("-n", "-@", "2", "-T", os.path.join(temp_dir, "sorttmp"), "-o", out, input_path)
    except (OSError, ValueError) as exc:
        raise InputError(f"query-name sorting failed for {input_path}: {exc}") from exc
    return out


def iter_templates(bam: pysam.AlignmentFile) -> Iterator[list[Any]]:
    """Yield lists of records sharing a query name, in stream order."""
    current_name: str | None = None
    group: list[Any] = []
    primary_flags: list[int] = []
    for record in bam.fetch(until_eof=True):
        name = record.query_name
        if name is None:
            raise DataInvariantError("BAM record without QNAME")
        if group and name != current_name:
            yield group
            group = []
            primary_flags = []
        group.append(record)
        if not (record.flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)):
            primary_flags.append(record.flag)
            if len(primary_flags) > 2:
                raise DataInvariantError(
                    f"template {name!r} has more than two primary records "
                    f"({len(primary_flags)} found); input violates pairing invariants"
                )
        current_name = name
    if group:
        yield group


def split_primaries(records: list[Any]) -> tuple[list[Any], list[Any]]:
    """Split a template into primary and non-primary (secondary/supplementary) records."""
    primaries = [r for r in records if not (r.flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY))]
    others = [r for r in records if r.flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY)]
    return primaries, others


def check_pair_flags(primaries: list[Any], name: str) -> None:
    """Validate the read1/read2 flag structure of primary records."""
    seen_r1 = 0
    seen_r2 = 0
    for record in primaries:
        if not (record.flag & FLAG_PAIRED):
            continue
        is_r1 = bool(record.flag & 0x40)
        is_r2 = bool(record.flag & 0x80)
        if is_r1 == is_r2:
            raise DataInvariantError(
                f"primary record {name!r} sets both or neither of read1/read2 flags "
                f"(flag {record.flag})"
            )
        seen_r1 += is_r1
        seen_r2 += is_r2
    if seen_r1 > 1 or seen_r2 > 1:
        # Two primaries claiming the same read slot would silently overwrite
        # each other in the processing loop below; reject loudly instead.
        raise DataInvariantError(
            f"template {name!r} has duplicate primary read1/read2 records "
            f"(read1={seen_r1}, read2={seen_r2}); input violates pairing invariants"
        )
