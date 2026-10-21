"""CIGAR and reference-coordinate arithmetic.

All positions are 0-based half-open reference coordinates unless a name says
otherwise (BAM ``pos`` is 0-based; SAM text POS is 1-based and only appears at
the CLI/serialization boundary).

Nothing here trusts ``M`` as evidence of a match: verification functions
compare aligned query bases against the reference base by base.
"""

from __future__ import annotations

import numpy as np

from pairloom.exceptions import DataInvariantError

# BAM CIGAR operation codes (SAM spec v1.7).
CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQUAL = 7
CIGAR_X = 8

#: Operations that consume query bases.
QUERY_OPS = frozenset({CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQUAL, CIGAR_X})
#: Operations that consume reference bases.
REF_OPS = frozenset({CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQUAL, CIGAR_X})
#: Operations that align query bases to reference bases.
ALIGNED_OPS = frozenset({CIGAR_M, CIGAR_EQUAL, CIGAR_X})

_CIGAR_LETTERS = "MIDNSHP=X"


def validate_cigar(
    cigar: list[tuple[int, int]] | tuple[tuple[int, int], ...],
    seq_len: int | None,
    context: str,
) -> None:
    """Fail loudly on malformed CIGARs: bad ops, zero lengths, or seq mismatch."""
    if not cigar:
        raise DataInvariantError(f"empty CIGAR for {context}")
    total_q = 0
    for op, length in cigar:
        if not 0 <= op <= 8 or op == CIGAR_P:
            raise DataInvariantError(f"unsupported CIGAR op {op} (len {length}) for {context}")
        if length <= 0:
            raise DataInvariantError(f"non-positive CIGAR length {length} (op {op}) for {context}")
        if op in QUERY_OPS:
            total_q += length
    if seq_len is not None and total_q != seq_len:
        raise DataInvariantError(
            f"CIGAR consumes {total_q} query bases but sequence has {seq_len} for {context}"
        )


def cigar_ref_len(cigar: tuple[tuple[int, int], ...]) -> int:
    return sum(length for op, length in cigar if op in REF_OPS)


def cigar_query_len(cigar: tuple[tuple[int, int], ...]) -> int:
    return sum(length for op, length in cigar if op in QUERY_OPS)


def ref_interval(pos0: int, cigar: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    """Return the aligned reference interval [start0, end0) covered by ``cigar``."""
    return pos0, pos0 + cigar_ref_len(cigar)


def aligned_pairs(
    pos0: int,
    cigar: tuple[tuple[int, int], ...],
    seq_len: int,
    context: str,
) -> np.ndarray:
    """Return an (n, 2) int64 array of (query_pos, ref_pos) for aligned bases.

    Only M/=/X operations produce pairs; insertions, deletions, skips, and
    clips are excluded by construction.
    """
    validate_cigar(cigar, seq_len, context)
    qblocks: list[np.ndarray] = []
    rblocks: list[np.ndarray] = []
    qpos = 0
    rpos = pos0
    for op, length in cigar:
        if op in ALIGNED_OPS:
            qblocks.append(np.arange(qpos, qpos + length, dtype=np.int64))
            rblocks.append(np.arange(rpos, rpos + length, dtype=np.int64))
        if op in QUERY_OPS:
            qpos += length
        if op in REF_OPS:
            rpos += length
    if not qblocks:
        return np.empty((0, 2), dtype=np.int64)
    return np.stack([np.concatenate(qblocks), np.concatenate(rblocks)], axis=1)


def compute_md_nm(
    pos0: int,
    cigar: tuple[tuple[int, int], ...],
    query_seq: str,
    ref_seq: str,
    context: str,
) -> tuple[str, int]:
    """Recompute the MD tag and NM count by comparing bases to the reference.

    ``query_seq`` must be in BAM storage orientation; ``ref_seq`` must cover
    exactly [pos0, pos0 + cigar_ref_len). Reference skips (N) are emitted as
    deletions in MD, matching samtools behaviour; they do not add to NM.
    """
    need = cigar_ref_len(cigar)
    if len(ref_seq) < need:
        raise DataInvariantError(
            f"reference slice too short for {context}: need {need} bases from "
            f"pos {pos0}, got {len(ref_seq)}"
        )
    md_parts: list[str] = []
    run = 0
    nm = 0
    qi = 0
    ri = 0
    for op, length in cigar:
        if op in ALIGNED_OPS:
            for k in range(length):
                qbase = query_seq[qi + k]
                rbase = ref_seq[ri + k]
                if qbase == rbase:
                    run += 1
                else:
                    md_parts.append(str(run))
                    md_parts.append(rbase)
                    run = 0
                    nm += 1
            qi += length
            ri += length
        elif op == CIGAR_I:
            nm += length
            qi += length
        elif op == CIGAR_D:
            md_parts.append(str(run))
            run = 0
            md_parts.append("^" + ref_seq[ri : ri + length])
            nm += length
            ri += length
        elif op == CIGAR_N:
            md_parts.append(str(run))
            run = 0
            md_parts.append("^" + ref_seq[ri : ri + length])
            ri += length
        elif op == CIGAR_S:
            qi += length
        elif op == CIGAR_H:
            continue
        else:  # pragma: no cover - validate_cigar rejects other ops first
            raise DataInvariantError(f"unsupported CIGAR op {op} in MD recompute for {context}")
    md_parts.append(str(run))
    return "".join(md_parts), nm


def reverse_cigar(cigar: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Reverse CIGAR op order, preserving op types and lengths.

    Reversing the op list yields the valid CIGAR of the reverse-complemented
    read aligned to the same reference interval; leading and trailing clips
    swap ends, which is exactly the reference-orientation transform.
    """
    return tuple(reversed(cigar))


def normalize_cigar(cigar: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Merge adjacent identical operations; keep everything else intact."""
    merged: list[tuple[int, int]] = []
    for op, length in cigar:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return tuple(merged)


def cigar_to_string(cigar: tuple[tuple[int, int], ...]) -> str:
    return "".join(f"{length}{_CIGAR_LETTERS[op]}" for op, length in cigar)
