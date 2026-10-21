"""Cached FASTA reference access.

The cache stores whole contigs up to a configurable base budget (LRU across
contigs); larger contigs are fetched directly through HTSlib, which performs
its own buffered region reads. Bases are always returned uppercase.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import pysam

from pairloom.exceptions import DataInvariantError, InputError


class ReferenceCache:
    """Bounded, thread/process-local cached view of an indexed FASTA."""

    def __init__(
        self,
        path: Path,
        max_contigs: int = 4,
        max_bases: int = 4_000_000,
    ) -> None:
        self._path = Path(path)
        self._max_contigs = max_contigs
        self._max_bases = max_bases
        self._cache: OrderedDict[str, str] = OrderedDict()
        self._cached_bases = 0
        try:
            self._fasta = pysam.FastaFile(str(self._path))
        except (OSError, ValueError) as exc:
            raise InputError(
                f"cannot open reference FASTA {self._path}: {exc}. "
                "If the file exists, an accompanying .fai index is required "
                "(samtools faidx)."
            ) from exc
        names = list(self._fasta.references)
        lengths = list(self._fasta.lengths)
        self._lengths: dict[str, int] = dict(zip(names, lengths, strict=True))

    @property
    def contig_names(self) -> list[str]:
        return list(self._lengths)

    def has_contig(self, name: str) -> bool:
        return name in self._lengths

    def length(self, name: str) -> int:
        try:
            return self._lengths[name]
        except KeyError:
            raise InputError(
                f"reference {self._path} lacks sequence {name!r}; the BAM header "
                "and the FASTA must agree"
            ) from None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        """Return uppercase reference bases for [start0, end0)."""
        contig_len = self.length(name)
        if start0 < 0 or end0 > contig_len or start0 > end0:
            raise DataInvariantError(
                f"reference fetch [{start0}, {end0}) out of bounds for {name} (length {contig_len})"
            )
        if end0 == start0:
            return ""
        seq = self._cache.get(name)
        if seq is None:
            seq = self._fasta.fetch(name, start0, end0).upper()
            if end0 - start0 == contig_len and self._contig_cacheable(contig_len):
                self._store(name, seq)
            return seq
        return seq[start0:end0]

    def close(self) -> None:
        self._fasta.close()
        self._cache.clear()
        self._cached_bases = 0

    def _contig_cacheable(self, contig_len: int) -> bool:
        return contig_len <= self._max_bases

    def _store(self, name: str, seq: str) -> None:
        self._cache[name] = seq
        self._cached_bases += len(seq)
        while len(self._cache) > self._max_contigs or self._cached_bases > self._max_bases:
            _, evicted = self._cache.popitem(last=False)
            self._cached_bases -= len(evicted)
