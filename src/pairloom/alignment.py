"""Alignment evidence and reference-orientation read structures.

Two central abstractions:

- :class:`Placement` - a candidate genomic placement with an alignment score
  that was verified base by base against the reference (never inferred from
  CIGAR ``M`` ops alone).
- :class:`ReadStructure` - a read rewritten in reference orientation: aligned
  bases keyed by reference position, deletions, insertions keyed by boundary,
  and outer soft clips. This makes reverse-strand mates directly comparable
  to forward-strand mates for overlap analysis and merging.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pairloom.cigarops import (
    ALIGNED_OPS,
    CIGAR_D,
    CIGAR_H,
    CIGAR_I,
    CIGAR_N,
    CIGAR_S,
    QUERY_OPS,
    cigar_ref_len,
    reverse_cigar,
    validate_cigar,
)
from pairloom.exceptions import DataInvariantError

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)

ACGT = frozenset("ACGT")


def revcomp(seq: str, context: str) -> str:
    """Reverse-complement an IUPAC nucleotide string; unknown letters fail loudly."""
    unknown = set(seq) - frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
    if unknown:
        raise DataInvariantError(
            f"non-IUPAC sequence character(s) {sorted(unknown)!r} in {context}"
        )
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True, slots=True)
class Placement:
    """A candidate placement with verified alignment evidence."""

    tid: int
    pos0: int
    cigar: tuple[tuple[int, int], ...]
    score: int
    exact: bool
    ref_span: int = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "ref_span", cigar_ref_len(self.cigar))

    @property
    def end0(self) -> int:
        return self.pos0 + self.ref_span


#: Score model (documented in the README): matches earn +1, mismatches cost 3,
#: each inserted or deleted base costs 2.
MATCH_SCORE = 1
MISMATCH_PENALTY = 3
GAP_BASE_PENALTY = 2


def score_placement(
    pos0: int,
    cigar: tuple[tuple[int, int], ...],
    query_seq: str,
    ref_seq: str,
    context: str,
) -> tuple[int, bool]:
    """Verify a placement against the reference and return (score, exact).

    ``exact`` is True only when every aligned query base equals the reference
    base and the alignment has no insertions or deletions. Soft clips never
    count toward identity.
    """
    need = cigar_ref_len(cigar)
    if len(ref_seq) < need:
        raise DataInvariantError(
            f"reference slice too short for {context}: need {need} bases from pos {pos0}, "
            f"got {len(ref_seq)}"
        )
    validate_cigar(cigar, len(query_seq), context)
    qi = 0
    ri = 0
    matches = 0
    mismatches = 0
    ins_bases = 0
    del_bases = 0
    for op, length in cigar:
        if op in ALIGNED_OPS:
            for k in range(length):
                if query_seq[qi + k] == ref_seq[ri + k]:
                    matches += 1
                else:
                    mismatches += 1
            qi += length
            ri += length
        elif op == CIGAR_I:
            ins_bases += length
            qi += length
        elif op == CIGAR_D:
            del_bases += length
            ri += length
        elif op == CIGAR_N:
            ri += length
        elif op == CIGAR_S:
            qi += length
        elif op == CIGAR_H:
            continue
    score = (
        MATCH_SCORE * matches
        - MISMATCH_PENALTY * mismatches
        - GAP_BASE_PENALTY * (ins_bases + del_bases)
    )
    exact = mismatches == 0 and ins_bases == 0 and del_bases == 0
    return score, exact


def slidable_cigar(cigar: tuple[tuple[int, int], ...]) -> bool:
    """True when the CIGAR is a single ungapped M block with optional soft clips."""
    seen_m = False
    for op, _length in cigar:
        if op == 0:
            if seen_m:
                return False
            seen_m = True
        elif op in (CIGAR_S, CIGAR_H):
            continue
        else:
            return False
    return seen_m


def aligned_block(
    cigar: tuple[tuple[int, int], ...],
    query_seq: str,
) -> str:
    """Return the query bases covered by M/=/X operations, in stored orientation."""
    parts: list[str] = []
    qi = 0
    for op, length in cigar:
        if op in ALIGNED_OPS:
            parts.append(query_seq[qi : qi + length])
        if op in QUERY_OPS:
            qi += length
    return "".join(parts)


@dataclass(slots=True)
class ReadStructure:
    """A read expressed in reference orientation."""

    tid: int
    start0: int
    end0: int
    is_reverse: bool
    mapq: int
    query_name: str
    #: reference position -> (base, phred quality) for aligned M bases
    bases: dict[int, tuple[str, int]]
    #: reference positions consumed by D operations
    deleted: frozenset[int]
    #: reference positions consumed by N skips
    skipped: frozenset[int]
    #: insertion boundary -> (inserted bases, phred qualities); an insertion at
    #: boundary b sits between reference positions b-1 and b
    insertions: dict[int, tuple[str, tuple[int, ...]]]
    lead_clip: tuple[str, tuple[int, ...]]
    trail_clip: tuple[str, tuple[int, ...]]
    hard_clipped: bool
    has_qualities: bool
    read_group: str | None
    cigar_ref_orientation: tuple[tuple[int, int], ...]
    aligned_len: int


def build_read_structure(record: object) -> ReadStructure:
    """Build the reference-orientation structure of a mapped pysam record.

    For reverse-strand reads the stored sequence is the reverse complement of
    the reference strand, so we reverse-complement the sequence, reverse the
    quality array, and reverse the CIGAR op order; afterwards every value is
    expressed in reference orientation for both strands.
    """
    seq = record.query_sequence
    if seq is None:
        raise DataInvariantError(
            f"record {record.query_name} has no SEQ; PairLoom requires sequences"
        )
    seq = seq.upper()
    quals = record.query_qualities
    has_qualities = quals is not None
    qual_list: list[int] = list(quals) if has_qualities else [0] * len(seq)
    cigar = tuple(record.cigartuples or ())
    context = f"read {record.query_name}"
    validate_cigar(cigar, len(seq), context)

    is_reverse = record.is_reverse
    if is_reverse:
        seq = revcomp(seq, context)
        qual_list.reverse()
        cigar = reverse_cigar(cigar)

    start0 = record.reference_start
    bases: dict[int, tuple[str, int]] = {}
    deleted: set[int] = set()
    skipped: set[int] = set()
    insertions: dict[int, tuple[str, tuple[int, ...]]] = {}
    lead_bases: list[str] = []
    lead_quals: list[int] = []
    trail_bases: list[str] = []
    trail_quals: list[int] = []
    in_trail = False

    qi = 0
    ri = start0
    aligned_len = 0
    hard_clipped = False
    for op, length in cigar:
        if op in ALIGNED_OPS:
            for k in range(length):
                bases[ri + k] = (seq[qi + k], qual_list[qi + k])
            aligned_len += length
            qi += length
            ri += length
            in_trail = True
        elif op == CIGAR_I:
            insertions[ri] = (seq[qi : qi + length], tuple(qual_list[qi : qi + length]))
            qi += length
        elif op == CIGAR_D:
            deleted.update(range(ri, ri + length))
            ri += length
        elif op == CIGAR_N:
            skipped.update(range(ri, ri + length))
            ri += length
        elif op == CIGAR_S:
            clip_bases = seq[qi : qi + length]
            clip_quals = tuple(qual_list[qi : qi + length])
            if in_trail:
                trail_bases.append(clip_bases)
                trail_quals.append(clip_quals)
            else:
                lead_bases.append(clip_bases)
                lead_quals.append(clip_quals)
            qi += length
        elif op == CIGAR_H:
            hard_clipped = True
    rg = record.get_tag("RG") if record.has_tag("RG") else None
    return ReadStructure(
        tid=record.reference_id,
        start0=start0,
        end0=ri,
        is_reverse=is_reverse,
        mapq=record.mapping_quality,
        query_name=record.query_name,
        bases=bases,
        deleted=frozenset(deleted),
        skipped=frozenset(skipped),
        insertions=insertions,
        lead_clip=("".join(lead_bases), tuple(q for chunk in lead_quals for q in chunk)),
        trail_clip=("".join(trail_bases), tuple(q for chunk in trail_quals for q in chunk)),
        hard_clipped=hard_clipped,
        has_qualities=has_qualities,
        read_group=rg,
        cigar_ref_orientation=cigar,
        aligned_len=aligned_len,
    )


def covers(struct: ReadStructure, boundary: int) -> bool:
    """True when ``boundary`` lies within the read's aligned span [start0, end0]."""
    return struct.start0 <= boundary <= struct.end0
