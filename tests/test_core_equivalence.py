"""Rust core vs Python reference parity.

The Rust core (pairloom_core) is a performance port of the Python reference
implementations. These tests hold the two bit-identical on a randomized
corpus: structures, merge outcomes (success fields and every failure
reason), MD/NM computation, and overlap aggregation.
"""

from __future__ import annotations

import random
from dataclasses import replace

import pairloom_core as core
import pysam
import pytest

from pairloom.alignment import build_read_structure
from pairloom.cigarops import compute_md_nm
from pairloom.config import PairLoomConfig
from pairloom.merge import MergeFailure, try_merge
from pairloom.overlap import analyze_overlap
from pairloom.parquetout import FastDisagreementAccumulator
from pairloom.pipeline import _chunk_observation_arrays

HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 200}]}


def _cfg() -> PairLoomConfig:
    return PairLoomConfig(
        input_path="in.bam",
        reference_path="ref.fa",
        mode="merge",
        merged_bam="m.bam",
        unmerged_bam="u.bam",
        disagreement_parquet="d.parquet",
        read_stats_parquet="s.parquet",
    )


def _record(
    name: str, flag: int, pos: int, cigar: list[tuple[int, int]], seq: str, quals: list[int]
):
    rec = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(HEADER))
    rec.query_name = name
    rec.flag = flag
    rec.reference_id = 0
    rec.reference_start = pos
    rec.cigartuples = cigar
    rec.query_sequence = seq
    rec.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + q) for q in quals))
    rec.mapping_quality = 60
    return rec


COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def _rand_pair(rng: random.Random, n: int):
    ref = "".join(rng.choice("ACGT") for _ in range(200))
    frag_start = rng.randrange(5, 60)
    insert = rng.choice([50, 70, 90, 110])
    rl = 40
    r1_fwd = rng.random() < 0.5
    frag = ref[frag_start : frag_start + insert]
    r1_seq = frag[:rl]
    r2_seq = "".join(COMP[b] for b in reversed(frag[insert - rl :]))
    # inject sequencing errors and N calls
    r1_seq = "".join(rng.choice("ACGTN") if rng.random() < 0.05 else b for b in r1_seq)
    r2_seq = "".join(rng.choice("ACGTN") if rng.random() < 0.05 else b for b in r2_seq)
    r1_quals = [rng.randrange(2, 60) for _ in r1_seq]
    r2_quals = [rng.randrange(2, 60) for _ in r2_seq]
    if r1_fwd:
        flag1, flag2 = 0x1 | 0x40, 0x1 | 0x80 | 0x10
    else:
        flag1, flag2 = 0x1 | 0x40 | 0x10, 0x1 | 0x80
        # stored SEQ is reference orientation: reverse record's stored seq is revcomp
        r1_seq = "".join(COMP[b] for b in reversed(r1_seq))
        r1_quals = r1_quals[::-1]
        r2_seq = "".join(COMP[b] for b in reversed(r2_seq))
        r2_quals = r2_quals[::-1]
    # occasional indel cigar on r1 (mid-fragment insertion)
    if rng.random() < 0.3 and insert >= 90:
        r1_seq = r1_seq[:20] + "AC" + r1_seq[20:]
        r1_quals = r1_quals[:20] + [30, 30] + r1_quals[20:]
        cigar1 = [(0, 20), (1, 2), (0, rl - 20)]
    else:
        cigar1 = [(0, rl)]
    cigar2 = [(0, rl)]
    pos1 = frag_start if r1_fwd else frag_start  # reference starts (simplified FR)
    pos2 = frag_start + insert - rl
    if not r1_fwd:
        pos1, pos2 = pos2, pos1
    r1 = _record(f"t{n}", flag1, pos1, cigar1, r1_seq, r1_quals)
    r2 = _record(f"t{n}", flag2, pos2, cigar2, r2_seq, r2_quals)
    return ref, r1, r2


@pytest.mark.parametrize("seed", range(40))
def test_merge_pair_matches_python(seed: int) -> None:
    rng = random.Random(seed)
    ref, r1, r2 = _rand_pair(rng, seed)
    cfg = _cfg()
    py1 = build_read_structure(r1)
    py2 = build_read_structure(r2)
    from pairloom.candidates import normalize_record_cigar

    def core_struct(rec):
        return core.build_structure(
            rec.reference_id,
            rec.reference_start,
            bool(rec.is_reverse),
            rec.mapping_quality,
            rec.query_name,
            normalize_record_cigar(rec.cigartuples),
            rec.query_sequence,
            list(rec.query_qualities),
            rec.get_tag("RG") if rec.has_tag("RG") else None,
        )

    c1, c2 = core_struct(r1), core_struct(r2)
    assert (c1.start0, c1.end0, c1.n_positions) == (py1.start0, py1.end0, len(py1.bases))
    assert (c2.start0, c2.end0, c2.n_positions) == (py2.start0, py2.end0, len(py2.bases))

    for policy in ("quality", "strict"):
        cfg = replace(cfg, conflict_policy=policy)
        py = try_merge(py1, py2, cfg)
        rust = core.merge_pair(
            c1,
            c2,
            cfg.min_mapq,
            cfg.min_overlap,
            cfg.min_fragment_length,
            cfg.max_fragment_length,
            cfg.min_quality_difference,
            cfg.max_consensus_quality,
            cfg.conflict_policy,
            cfg.library_orientation,
            None,
            False,
            False,
        )
        if isinstance(py, MergeFailure):
            assert not rust.ok
            assert rust.reason == py.reason
            continue
        assert rust.ok
        assert rust.start0 == py.start0
        assert [tuple(x) for x in rust.cigar] == list(py.cigar)
        assert rust.seq == py.seq
        assert list(rust.quals) == list(py.quals)
        assert rust.overlap_len == py.overlap_len
        if rust.cigar:
            ref_seq = ref[rust.start0 : rust.start0 + 200]
            assert core.compute_md_nm(rust.start0, list(rust.cigar), rust.seq, ref_seq) == (
                compute_md_nm(rust.start0, tuple(rust.cigar), rust.seq, ref_seq, "ctx")
            )


@pytest.mark.parametrize("seed", range(20))
def test_aggregate_matches_per_observation(seed: int) -> None:
    rng = random.Random(1000 + seed)
    ref, r1, r2 = _rand_pair(rng, seed)
    cfg = _cfg()
    py1 = build_read_structure(r1)
    py2 = build_read_structure(r2)
    ov_lo = max(py1.start0, py2.start0)
    ov_hi = min(py1.end0, py2.end0)
    if ov_lo >= ov_hi:
        pytest.skip("no overlap")
    ref_seq = ref[ov_lo:ov_hi]
    overlap = analyze_overlap(py1, py2, ref_seq, ov_lo, cfg)
    entries: dict[int, list] = {}
    for o in overlap.observations:
        e = entries.setdefault(o.rpos0, [o.ref_base, 0, False, {}])
        e[1] += 1
        e[2] = e[2] or o.pair_disagrees
        if o.base != o.ref_base:
            counts = e[3].setdefault(o.base, [0, 0])
            counts[1 if o.strand_reverse else 0] += 1
    # rust path
    from pairloom.candidates import normalize_record_cigar

    def core_struct(rec):
        return core.build_structure(
            rec.reference_id,
            rec.reference_start,
            bool(rec.is_reverse),
            rec.mapping_quality,
            rec.query_name,
            normalize_record_cigar(rec.cigartuples),
            rec.query_sequence,
            list(rec.query_qualities),
            None,
        )

    blobs = core.aggregate_overlap(
        core_struct(r1), core_struct(r2), 0, ref_seq, ov_lo, cfg.min_mapq, cfg.min_base_quality
    )
    arrays = _chunk_observation_arrays([blobs])
    fast = FastDisagreementAccumulator(cfg)
    fast.add_chunk(*arrays)
    fast._reduce(force=True)
    keys = fast._r_keys
    py_total = sum(e[1] for e in entries.values())
    rust_total = int(fast._r_dp.sum())
    assert rust_total == py_total
    for i in range(keys.size):
        rpos = int(keys[i]) & 0xFFFFFFFF
        e = entries[rpos]
        assert int(fast._r_dp[i]) == e[1]
        assert bool(fast._r_disg[i]) == e[2]
        assert chr(int(fast._r_ref[i])) == e[0]
