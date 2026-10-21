"""End-to-end integration tests (scenarios 2, 20, 21, 23, 26, 12 quickcheck)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from conftest import LONG_REPEAT, UNIQUE


def rich_records(fx) -> list:
    """A mixed fixture: repeats, unique overlaps, disagreements, unmapped."""
    ref = fx.ref_seqs["chr1"]

    def alt(base):
        return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]

    records: list = []
    # Pair-aware optimization target inside the tandem repeat.
    records += fx.pair("rep1", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16)
    # Overlapping unique pairs with disagreements.
    p = UNIQUE + 12
    records += fx.pair("dis1", UNIQUE, UNIQUE + 5, read_len=20, r1_mutations={p: alt(ref[p])})
    records += fx.pair("dis2", UNIQUE + 5, UNIQUE + 15, read_len=20, r1_mutations={p: alt(ref[p])})
    # Clean overlap, soft clips, indels.
    records += fx.pair("clean", UNIQUE, UNIQUE + 10, read_len=20, r1_lead_clip=3, r2_trail_clip=4)
    r1_seq = ref[UNIQUE : UNIQUE + 10] + "GGA" + ref[UNIQUE + 10 : UNIQUE + 20]
    r1 = fx.make_record(
        "ins", 0x1 | 0x40 | 0x2, UNIQUE, "10M3I10M", r1_seq, mate_contig="chr1", mate_pos=UNIQUE
    )
    r2 = fx.make_record(
        "ins",
        0x1 | 0x80 | 0x10 | 0x20 | 0x2,
        UNIQUE,
        "10M3I10M",
        r1_seq,
        mate_contig="chr1",
        mate_pos=UNIQUE,
    )
    records += [r1, r2]
    # Non-overlapping and unmapped.
    records += fx.pair("far", UNIQUE, UNIQUE + 200, read_len=20)
    records += [
        fx.make_record("unm", 0x1 | 0x40 | 0x4, -1, "*", "ACGTACGTAC"),
        fx.make_record("unm", 0x1 | 0x80 | 0x4, -1, "*", "ACGTACGTAC"),
    ]
    return records


class TestPairAwareOptimization:
    def test_joint_overlap_maximization(self, fx) -> None:
        """Scenario 2: joint adjustment favours the larger supported overlap.

        Both mates sit in the period-4 tandem repeat; the original placement
        (404, 416) overlaps on 4 bases. Diagonal candidate pairs (s, s)
        overlap fully (16 bases); the leftmost of them, (400, 400), wins.
        """
        pair = fx.pair("rep", LONG_REPEAT + 4, LONG_REPEAT + 16, read_len=16)
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust", extra=["--min-overlap", "10"])
        assert code == 0
        by_flag = {r.flag & 0xC0: r for r in fx.read_bam(out["adjusted_bam"])}
        r1 = by_flag[0x40]
        r2 = by_flag[0x80]
        assert r1.reference_start == LONG_REPEAT
        assert r2.reference_start == LONG_REPEAT
        # Mate fields follow the adjustment.
        assert r1.next_reference_start == LONG_REPEAT
        assert r2.next_reference_start == LONG_REPEAT
        assert r1.template_length == 16
        assert r2.template_length == 16  # equal starts: both positive
        # Ambiguous repeats: conservative MAPQ on both mates.
        assert r1.mapping_quality <= 3
        assert r2.mapping_quality <= 3
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["mergeable"] == 1  # 16 >= min-overlap 10

    def test_moved_read_with_stale_aligner_tags(self, fx) -> None:
        """Stale aligner tags on moved reads are removed, not crashed on.

        Real aligners (BWA) tag AS/XS/MQ on every primary; before the fix,
        removing them on movement called a nonexistent pysam method and
        aborted the run.
        """
        pair = fx.pair(
            "tagged",
            LONG_REPEAT + 4,
            LONG_REPEAT + 16,
            read_len=16,
            tags={"AS": 16, "XS": 12, "MQ": 40},
        )
        fx.write_bam_and_set(pair)
        code, out = fx.run(mode="adjust")
        assert code == 0
        for rec in fx.read_bam(out["adjusted_bam"]):
            assert rec.reference_start == LONG_REPEAT  # moved to leftmost
            tags = dict(rec.get_tags())
            assert "AS" not in tags and "XS" not in tags and "MQ" not in tags
            assert tags["ZC"] == "mergeable"
            assert tags["NM"] == 0


class TestBoundedReferenceFetch:
    def test_distant_mates_fetch_only_the_overlap(self, fx, monkeypatch) -> None:
        """No reference fetch may span the gap between distant mates.

        Old behaviour fetched the entire inter-mate span (100,002 bases
        here) before the fragment-length rejection; every fetch must stay
        bounded by the discovery window and read size.
        """
        import random

        from pairloom.reference import ReferenceCache

        rng = random.Random(11)
        fx.ref_seqs = {"big": "".join(rng.choice("ACGT") for _ in range(200_000))}
        fx.__post_init__()
        big = fx.ref_seqs["big"]
        r1 = fx.make_record(
            "gap2",
            0x1 | 0x40 | 0x2,
            50_000,
            "20M",
            big[50_000:50_020],
            contig="big",
            mate_contig="big",
            mate_pos=150_000,
        )
        r2 = fx.make_record(
            "gap2",
            0x1 | 0x80 | 0x10 | 0x20 | 0x2,
            150_000,
            "20M",
            big[150_000:150_020],
            contig="big",
            mate_contig="big",
            mate_pos=50_000,
        )
        fx.write_bam_and_set([r1, r2])
        lengths: list[int] = []
        orig = ReferenceCache.fetch

        def spy(self, name, start0, end0):
            lengths.append(end0 - start0)
            return orig(self, name, start0, end0)

        monkeypatch.setattr(ReferenceCache, "fetch", spy)
        code, out = fx.run(mode="adjust")
        assert code == 0
        stats = fx.read_parquet(out["read_stats_parquet"])
        by_type = dict(zip(stats["type"], stats["count"], strict=False))
        assert by_type["both_aligned"] == 1  # kept, no overlap
        assert lengths  # discovery still fetches, but bounded
        assert max(lengths) < 2_000


class TestSortOrderIndependence:
    def test_coordinate_and_qname_sorted_inputs_agree(self, fx) -> None:
        """Scenario 20: both sort orders produce identical results."""
        records = rich_records(fx)
        qname_bam = fx.write_bam(records, so="queryname")
        coord_bam = fx.write_bam(records, so="coordinate")
        fx.last_bam = qname_bam
        code_q, out_q = fx.run(mode="adjust")
        fx.last_bam = coord_bam
        code_c, out_c = fx.run(mode="adjust")
        assert code_q == 0 and code_c == 0

        def signature(path):
            return [
                (
                    r.query_name,
                    r.flag,
                    r.reference_start,
                    r.cigarstring,
                    r.mapping_quality,
                    sorted(r.get_tags()),
                )
                for r in fx.read_bam(path)
            ]

        assert sorted(signature(out_q["adjusted_bam"])) == sorted(signature(out_c["adjusted_bam"]))
        assert fx.read_parquet(out_q["disagreement_parquet"]) == fx.read_parquet(
            out_c["disagreement_parquet"]
        )


class TestDeterminism:
    def test_identical_results_across_thread_counts(self, fx) -> None:
        """Scenario 23: threads 1 vs 4 give byte-identical statistics."""
        fx.write_bam_and_set(rich_records(fx))
        code1, out1 = fx.run(mode="merge", extra=["--threads", "1", "--chunk-size", "4"])
        code4, out4 = fx.run(mode="merge", extra=["--threads", "4", "--chunk-size", "4"])
        assert code1 == 0 and code4 == 0
        for key in ("disagreement_parquet", "read_stats_parquet"):
            assert fx.read_parquet(out1[key]) == fx.read_parquet(out4[key])
        for key in ("merged_bam", "unmerged_bam"):
            sig1 = [
                (
                    r.query_name,
                    r.flag,
                    r.reference_start,
                    r.cigarstring,
                    r.query_sequence,
                    r.mapping_quality,
                    sorted(r.get_tags()),
                )
                for r in fx.read_bam(out1[key])
            ]
            sig4 = [
                (
                    r.query_name,
                    r.flag,
                    r.reference_start,
                    r.cigarstring,
                    r.query_sequence,
                    r.mapping_quality,
                    sorted(r.get_tags()),
                )
                for r in fx.read_bam(out4[key])
            ]
            assert sig1 == sig4


class TestBoundedMemory:
    def test_spilling_aggregator_matches_direct_result(self, fx) -> None:
        """Scenario 26: tiny spill threshold still yields exact tables."""
        ref = fx.ref_seqs["chr1"]
        records: list = []
        # 12 well-separated pairs, each with a distinct disagreement position
        # (30bp spacing so no position is covered by two pairs).
        for i in range(12):
            p = UNIQUE + 20 + 30 * i
            base = {"A": "C", "C": "G", "G": "T", "T": "A"}[ref[p]]
            records += fx.pair(
                f"sp{i}",
                UNIQUE + 10 + 30 * i,
                UNIQUE + 20 + 30 * i,
                read_len=20,
                r1_mutations={p: base},
            )
        fx.write_bam_and_set(records)
        code, out = fx.run(
            mode="merge",
            extra=["--max-aggregation-records", "3", "--parquet-batch-size", "2"],
        )
        assert code == 0
        table = fx.read_parquet(out["disagreement_parquet"])
        assert len(table["pos"]) == 12
        assert all(dp == 2 for dp in table["dp"])
        keys = list(zip(table["seq_name"], table["pos"], table["alt"], strict=False))
        assert keys == sorted(keys)

    def test_aggregator_unit_spill_merge(self, tmp_path) -> None:
        from pairloom.config import PairLoomConfig
        from pairloom.parquetout import DisagreementAggregator

        cfg = PairLoomConfig(
            input_path=Path("i.bam"),
            reference_path=Path("r.fa"),
            mode="adjust",
            adjusted_bam=Path("a.bam"),
            disagreement_parquet=Path("d.parquet"),
            read_stats_parquet=Path("s.parquet"),
            max_aggregation_records=2,
        )
        agg = DisagreementAggregator(cfg)
        agg.set_spill_dir(str(tmp_path))
        # Three pairs disagreeing at the same position; observations spread
        # across many positions to force multiple spills.
        for i in range(10):
            agg.add_observation(0, i, "A", "C", False, True)
            agg.add_observation(0, i, "A", "G", True, True)
            agg.add_observation(0, i, "A", "A", False, False)
        out = tmp_path / "final.parquet"
        reported = agg.finalize(out, ["chr1"])
        assert reported == 10
        import pyarrow.parquet as pq

        table = pq.read_table(str(out)).to_pylist()
        assert len(table) == 20  # two alt rows per position
        assert all(row["dp"] == 3 for row in table)
        for row in table[:2]:
            assert row["af"] == pytest.approx(1 / 3)

    def test_spill_finalize_memory_stays_bounded(self, tmp_path) -> None:
        """Lazy, bounded-fan-in spill merge: finalize must not load the dataset.

        100k positions with a 1000-entry threshold produce ~100 spills. The
        previous eager-list finalize peaked around 145 bytes/row (~15 MB
        here); the lazy cascade stays orders below that. The generous 10 MB
        bound separates the two decisively and holds across machines.
        """
        import tracemalloc

        import pyarrow.parquet as pq

        from pairloom.config import PairLoomConfig
        from pairloom.parquetout import DisagreementAggregator, _spill_stream

        cfg = PairLoomConfig(
            input_path=Path("i.bam"),
            reference_path=Path("r.fa"),
            mode="adjust",
            adjusted_bam=Path("a.bam"),
            disagreement_parquet=Path("d.parquet"),
            read_stats_parquet=Path("s.parquet"),
            max_aggregation_records=1000,
            parquet_batch_size=8192,
        )
        agg = DisagreementAggregator(cfg)
        agg.set_spill_dir(str(tmp_path))
        for i in range(100_000):
            disagrees = i % 500 == 0
            agg.add_observation(0, i, "A", "C" if disagrees else "A", False, disagrees)
        assert len(agg._spill_paths) >= 50
        # Contract: spill files stream as iterators, never materialized lists.
        pf = pq.ParquetFile(agg._spill_paths[0])
        stream = _spill_stream(pf, 8192)
        assert iter(stream) is stream
        assert not isinstance(stream, list)
        pf.close()
        tracemalloc.start()
        reported = agg.finalize(tmp_path / "final.parquet", ["chr1"])
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        assert reported == 200  # every 500th position disagrees
        assert peak < 10_000_000, f"finalize peak {peak / 1e6:.1f} MB is not bounded"


class TestIndependentValidation:
    @pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools")
    def test_samtools_quickcheck_all_outputs(self, fx) -> None:
        import subprocess

        fx.write_bam_and_set(rich_records(fx))
        code, out = fx.run(mode="merge", extra=["--sort-output"])
        assert code == 0
        for key in ("merged_bam", "unmerged_bam"):
            result = subprocess.run(
                ["samtools", "quickcheck", "-v", str(out[key])],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, result.stderr

    def test_merged_records_reference_consistency(self, fx) -> None:
        """Merged seq/quals/cigar/POS all consistent with the reference."""
        ref = fx.ref_seqs["chr1"]
        fx.write_bam_and_set(fx.pair("mc", UNIQUE, UNIQUE + 5, read_len=20))
        code, out = fx.run(mode="merge")
        assert code == 0
        rec = fx.read_bam(out["merged_bam"])[0]
        assert rec.reference_start == UNIQUE
        assert rec.query_sequence == ref[UNIQUE : UNIQUE + 25]
        assert rec.cigarstring == "25M"
        assert rec.is_unmapped is False
        assert rec.is_paired is False
        assert rec.is_reverse is False
        tags = dict(rec.get_tags())
        assert tags == {"ZC": "mergeable", "MD": "25", "NM": 0}
