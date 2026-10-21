"""Shared fixtures: canonical reference, BAM builders, and a run helper.

All expected values in tests are derived from the same reference strings the
fixtures use, so tests never hard-code random filler content.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path

import pysam
import pytest

COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


# Canonical reference layout (0-based intervals):
#   chr1 [0,120)    unique filler
#   chr1 [120,132)  "AACCAACCAACC"  repeat copy 1
#   chr1 [132,236)  unique filler
#   chr1 [236,248)  "AACCAACCAACC"  repeat copy 2
#   chr1 [248,400)  unique filler
#   chr1 [400,464)  "AACCAACC" * 8  long tandem repeat (period 4)
#   chr1 [464,500)  unique filler
#   chr1 [500,700)  unique region for overlap/merge tests
#   chr1 [700,1000) unique filler
#   chr2 [0,300)    unique filler
REPEAT1 = 120
REPEAT2 = 236
LONG_REPEAT = 400
LONG_REPEAT_LEN = 64
UNIQUE = 500
UNIQUE_LEN = 200

_RNG_SEED = 20261021


def canonical_reference() -> dict[str, str]:
    rng = random.Random(_RNG_SEED)

    def filler(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    chr1 = (
        filler(120)
        + "AACCAACCAACC"
        + filler(104)
        + "AACCAACCAACC"
        + filler(LONG_REPEAT - 248)
        + "AACCAACC" * 8
        + filler(UNIQUE - (LONG_REPEAT + LONG_REPEAT_LEN))
        + filler(UNIQUE_LEN)
        + filler(300)
    )
    chr2 = filler(300)
    assert len(chr1) == 1000
    assert len(chr2) == 300
    return {"chr1": chr1, "chr2": chr2}


@dataclass
class Fixture:
    """Everything a test needs to build inputs and check outputs."""

    tmp_path: Path
    ref_seqs: dict[str, str] = field(default_factory=canonical_reference)
    ref_path: Path | None = None
    _bam_count: int = 0
    _run_count: int = 0

    def __post_init__(self) -> None:
        self.ref_path = self.tmp_path / "ref.fa"
        with open(self.ref_path, "w") as handle:
            for name, seq in self.ref_seqs.items():
                handle.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    handle.write(seq[i : i + 60] + "\n")
        pysam.faidx(str(self.ref_path))

    # -- BAM construction ---------------------------------------------------

    def header_dict(self, so: str = "queryname") -> dict:
        return {
            "HD": {"VN": "1.6", "SO": so},
            "SQ": [{"SN": name, "LN": len(seq)} for name, seq in self.ref_seqs.items()],
            "RG": [{"ID": "rg1", "SM": "sample1"}, {"ID": "rg2", "SM": "sample2"}],
        }

    def write_bam(
        self,
        records: list[pysam.AlignedSegment],
        so: str = "queryname",
        name: str | None = None,
    ) -> Path:
        """Write records to a BAM; optionally coordinate-sort afterwards."""
        self._bam_count += 1
        path = self.tmp_path / (name or f"input_{self._bam_count}.bam")
        header = pysam.AlignmentHeader.from_dict(self.header_dict(so=so))
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for record in records:
                out.write(record)
        if so == "coordinate":
            sorted_path = str(path) + ".cs"
            pysam.sort("-o", sorted_path, str(path))
            path.unlink()
            Path(sorted_path).rename(path)
        return path

    def make_record(
        self,
        qname: str,
        flag: int,
        pos: int,
        cigar: str,
        seq: str,
        qual: int = 40,
        mapq: int = 60,
        contig: str = "chr1",
        mate_contig: str | None = None,
        mate_pos: int | None = None,
        tlen: int = 0,
        tags: dict | None = None,
    ) -> pysam.AlignedSegment:
        """Build one record; ``seq`` is always given in STORED orientation."""
        header = pysam.AlignmentHeader.from_dict(self.header_dict())
        record = pysam.AlignedSegment(header)
        record.query_name = qname
        record.flag = flag
        if pos >= 0 and cigar != "*":
            record.reference_id = list(self.ref_seqs).index(contig)
            record.reference_start = pos
        else:
            record.reference_id = -1
            record.reference_start = -1
        if mate_contig is not None and flag & 0x1:
            record.next_reference_id = list(self.ref_seqs).index(mate_contig)
            record.next_reference_start = mate_pos if mate_pos is not None else -1
        else:
            record.next_reference_id = -1
            record.next_reference_start = -1
        record.template_length = tlen
        if cigar != "*":
            record.cigarstring = cigar
        record.query_sequence = seq
        record.query_qualities = pysam.qualitystring_to_array(chr(33 + qual) * len(seq))
        record.mapping_quality = mapq
        for tag, value in (tags or {}).items():
            record.set_tag(tag, value)
        return record

    def pair(
        self,
        qname: str,
        r1_pos: int,
        r2_pos: int,
        read_len: int = 20,
        contig: str = "chr1",
        r1_mutations: dict[int, str] | None = None,
        r2_mutations: dict[int, str] | None = None,
        r1_qual: int = 40,
        r2_qual: int = 40,
        mapq: int = 60,
        r1_cigar: str | None = None,
        r2_cigar: str | None = None,
        r1_lead_clip: int = 0,
        r1_trail_clip: int = 0,
        r2_lead_clip: int = 0,
        r2_trail_clip: int = 0,
        r1_fwd: bool = True,
        tags: dict | None = None,
    ) -> list[pysam.AlignedSegment]:
        """Build an FR pair from the reference.

        Mutations are keyed by reference position (0-based). Clips are given
        in REFERENCE orientation for each mate; reverse mates are stored
        reverse-complemented with mirrored CIGARs. Mutated bases default to
        a base different from both mates' agreement where applicable.
        """
        seq = self.ref_seqs[contig]
        # R1 forward.
        r1_seq = self._build_oriented(
            seq, r1_pos, read_len, r1_mutations, r1_lead_clip, r1_trail_clip, True
        )
        r1_cig = r1_cigar or self._clip_cigar(read_len, r1_lead_clip, r1_trail_clip)
        # R2 reverse covering [r2_pos, r2_pos + read_len) in reference
        # orientation. Per the SAM convention, stored SEQ and CIGAR are in
        # reference orientation for BOTH strands, so R2 stores the
        # reference-orientation string directly and only the strand flag
        # differs.
        r2_stored_seq = self._build_oriented(
            seq, r2_pos, read_len, r2_mutations, r2_lead_clip, r2_trail_clip, True
        )
        r2_cig = r2_cigar or self._clip_cigar(read_len, r2_lead_clip, r2_trail_clip)
        if r1_fwd:
            r1_flag = 0x1 | 0x40
            r2_flag = 0x1 | 0x80 | 0x10 | 0x20
        else:
            r1_flag = 0x1 | 0x40 | 0x10 | 0x20
            r2_flag = 0x1 | 0x80
        r1 = self.make_record(
            qname,
            r1_flag | 0x2,
            r1_pos,
            r1_cig,
            r1_seq,
            qual=r1_qual,
            mapq=mapq,
            contig=contig,
            mate_contig=contig,
            mate_pos=r2_pos,
            tlen=0,
            tags=tags,
        )
        r2 = self.make_record(
            qname,
            r2_flag | 0x2,
            r2_pos,
            r2_cig,
            r2_stored_seq,
            qual=r2_qual,
            mapq=mapq,
            contig=contig,
            mate_contig=contig,
            mate_pos=r1_pos,
            tlen=0,
            tags=tags,
        )
        return [r1, r2]

    def _build_oriented(
        self,
        seq: str,
        pos: int,
        read_len: int,
        mutations: dict[int, str] | None,
        lead_clip: int,
        trail_clip: int,
        forward: bool,
    ) -> str:
        """Reference-orientation read string: clip + aligned + clip."""
        bases = list(seq[pos : pos + read_len])
        for rpos, new_base in (mutations or {}).items():
            if not pos <= rpos < pos + read_len:
                raise ValueError(f"mutation at {rpos} outside read [{pos}, {pos + read_len})")
            bases[rpos - pos] = new_base
        clip_l = "".join("N" for _ in range(lead_clip))
        clip_r = "".join("N" for _ in range(trail_clip))
        text = clip_l + "".join(bases) + clip_r
        return text if forward else rc(text)

    @staticmethod
    def _clip_cigar(read_len: int, lead: int, trail: int) -> str:
        parts = []
        if lead:
            parts.append(f"{lead}S")
        parts.append(f"{read_len}M")
        if trail:
            parts.append(f"{trail}S")
        return "".join(parts)

    # -- running the CLI ----------------------------------------------------

    def run(
        self,
        bam: Path | None = None,
        mode: str = "merge",
        extra: list[str] | None = None,
        force: bool = True,
    ) -> tuple[int, dict[str, Path]]:
        """Run pairmerge in-process; returns (exit_code, outputs)."""
        from pairloom.cli import main

        self._run_count += 1
        tag = f"run{self._run_count}"
        out = self.tmp_path / f"{tag}_out"
        out.mkdir(exist_ok=True)
        outputs: dict[str, Path] = {
            "disagreement_parquet": out / "disagreements.parquet",
            "read_stats_parquet": out / "stats.parquet",
        }
        argv = [
            "--input",
            str(bam or self.last_bam),
            "--reference",
            str(self.ref_path),
            "--mode",
            mode,
            "--disagreement-parquet",
            str(outputs["disagreement_parquet"]),
            "--read-stats-parquet",
            str(outputs["read_stats_parquet"]),
        ]
        if mode == "adjust":
            outputs["adjusted_bam"] = out / "adjusted.bam"
            argv += ["--adjusted-bam", str(outputs["adjusted_bam"])]
        else:
            outputs["merged_bam"] = out / "merged.bam"
            outputs["unmerged_bam"] = out / "unmerged.bam"
            argv += [
                "--merged-bam",
                str(outputs["merged_bam"]),
                "--unmerged-bam",
                str(outputs["unmerged_bam"]),
            ]
        if force:
            argv.append("--force")
        argv += extra or []
        code = main(argv)
        self.last_bam = bam or self.last_bam
        return code, outputs

    last_bam: Path | None = None

    def write_bam_and_set(self, records: list[pysam.AlignedSegment], so: str = "queryname") -> Path:
        path = self.write_bam(records, so=so)
        self.last_bam = path
        return path

    # -- reading outputs ----------------------------------------------------

    def read_bam(self, path: Path) -> list[pysam.AlignedSegment]:
        with pysam.AlignmentFile(str(path)) as f:
            return list(f.fetch(until_eof=True))

    def read_parquet(self, path: Path) -> dict:
        import pyarrow.parquet as pq

        return pq.read_table(str(path)).to_pydict()


@pytest.fixture()
def fx(tmp_path: Path) -> Fixture:
    fixture = Fixture(tmp_path=tmp_path)
    yield fixture
