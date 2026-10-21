"""Unit tests for CIGAR arithmetic and NM/MD recomputation (scenario 21, 32)."""

from __future__ import annotations

import numpy as np
import pytest

from pairloom.alignment import revcomp
from pairloom.cigarops import (
    aligned_pairs,
    cigar_query_len,
    cigar_ref_len,
    cigar_to_string,
    compute_md_nm,
    normalize_cigar,
    reverse_cigar,
    validate_cigar,
)
from pairloom.exceptions import DataInvariantError

M, INS, D, S_ = 0, 1, 2, 4


class TestCigarLengths:
    def test_ref_and_query_len(self) -> None:
        cigar = ((S_, 4), (M, 20), (INS, 2), (M, 5), (D, 3), (M, 10), (S_, 6))
        assert cigar_ref_len(cigar) == 38
        assert cigar_query_len(cigar) == 47

    def test_malformed_cigar_fails_loudly(self) -> None:
        with pytest.raises(DataInvariantError, match="consumes"):
            validate_cigar(((M, 20),), seq_len=25, context="read X")
        with pytest.raises(DataInvariantError, match="non-positive"):
            validate_cigar(((M, 0),), seq_len=0, context="read X")
        with pytest.raises(DataInvariantError, match="unsupported CIGAR op"):
            validate_cigar(((6, 3), (M, 5)), seq_len=5, context="read X")


class TestAlignedPairs:
    def test_pairs_exclude_clips_and_indels(self) -> None:
        cigar = ((S_, 2), (M, 3), (INS, 2), (M, 2), (D, 1), (M, 2))
        pairs = aligned_pairs(100, cigar, 11, "ctx")
        assert pairs.tolist() == [
            [2, 100],
            [3, 101],
            [4, 102],
            [7, 103],
            [8, 104],
            [9, 106],
            [10, 107],
        ]

    def test_reverse_cigar_mirrors(self) -> None:
        cigar = ((S_, 2), (M, 3), (INS, 1), (M, 2))
        assert reverse_cigar(cigar) == ((M, 2), (INS, 1), (M, 3), (S_, 2))
        assert cigar_to_string(normalize_cigar(((M, 2), (M, 3)))) == "5M"


class TestMdNm:
    REF = "AACCGGTTAACCGGTTAACCGGTT"

    def test_exact_match(self) -> None:
        md, nm = compute_md_nm(0, ((M, 10),), "AACCGGTTAA", self.REF, "t")
        assert (md, nm) == ("10", 0)

    def test_mismatches(self) -> None:
        md, nm = compute_md_nm(0, ((M, 10),), "AACCGGTXAAC", self.REF, "t")
        assert md == "7T2"
        assert nm == 1

    def test_deletion(self) -> None:
        # ref AACCGG|TT|AACCGGTT; read skips the middle TT
        seq = "AACCGGAACCGGT"
        md, nm = compute_md_nm(0, ((M, 6), (D, 2), (M, 7)), seq, self.REF, "t")
        assert md == "6^TT7"
        assert nm == 2

    def test_insertion_counts_nm_not_md(self) -> None:
        seq = "AACCGGTTAA" + "GGG" + "CCGGTT"[:6]
        md, nm = compute_md_nm(0, ((M, 10), (INS, 3), (M, 6)), seq, self.REF, "t")
        assert md == "16"
        assert nm == 3

    def test_short_reference_fails_loudly(self) -> None:
        with pytest.raises(DataInvariantError, match="too short"):
            compute_md_nm(0, ((M, 10),), "AACCGGTTAA", "AACCGG", "t")

    def test_revcomp(self) -> None:
        assert revcomp("AACCGGTTN", "t") == "NAACCGGTT"
        with pytest.raises(DataInvariantError):
            revcomp("ACGX", "t")


class TestNumpyPairsEmpty:
    def test_no_aligned_ops(self) -> None:
        pairs = aligned_pairs(10, ((S_, 5),), 5, "ctx")
        assert isinstance(pairs, np.ndarray)
        assert pairs.shape == (0, 2)
