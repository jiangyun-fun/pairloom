"""Parquet statistical output with explicit schemas and bounded memory.

Schemas (column names, order, and types are stable API):

Disagreement table (``--disagreement-parquet``)::

    seq_name string, pos int64 (1-based), ref string, alt string,
    dp int64, adf int64, adr int64, af float64

Read-stats table (``--read-stats-parquet``)::

    type string, count int64

The disagreement aggregator keeps at most ``max_aggregation_records``
positions in memory; beyond that, sorted partial aggregates spill to Parquet
files that are k-way merged lazily at finalize time (at most one Parquet
batch of rows per spill is resident, never a whole spill file). All counts
are integers, so results are deterministic regardless of worker scheduling.
"""

from __future__ import annotations

import heapq
import os
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from operator import itemgetter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from pairloom.classify import CLASSIFICATION_ORDER
from pairloom.config import PARQUET_COMPRESSIONS, PairLoomConfig
from pairloom.exceptions import DataInvariantError, OutputError

DISAGREEMENT_SCHEMA = pa.schema(
    [
        pa.field("seq_name", pa.string()),
        pa.field("pos", pa.int64()),
        pa.field("ref", pa.string()),
        pa.field("alt", pa.string()),
        pa.field("dp", pa.int64()),
        pa.field("adf", pa.int64()),
        pa.field("adr", pa.int64()),
        pa.field("af", pa.float64()),
    ]
)

READ_STATS_SCHEMA = pa.schema(
    [
        pa.field("type", pa.string()),
        pa.field("count", pa.int64()),
    ]
)

# Row layout: one "position row" per (tid, rpos) with alt == "" carrying the
# position-level dp/has_disg/ref base, plus one "alt row" per alternative
# allele with per-strand counts.
_SPILL_SCHEMA = pa.schema(
    [
        pa.field("tid", pa.int32()),
        pa.field("rpos", pa.int64()),
        pa.field("alt", pa.string()),
        pa.field("ref", pa.string()),
        pa.field("dp", pa.int64()),
        pa.field("has_disg", pa.bool_()),
        pa.field("adf", pa.int64()),
        pa.field("adr", pa.int64()),
    ]
)

SpillRow = tuple[int, int, str, str, int, bool, int, int]

#: Sort key of a spill row: (tid, rpos, alt) — position rows (alt == "")
#: precede their alt rows at the same position.
_ROW_KEY = itemgetter(0, 1, 2)

#: Maximum number of spill streams merged at once. More spills than this are
#: first reduced by cascading merges into intermediate spill files, so
#: finalize memory is bounded by configuration (fan-in x batch size), never
#: by dataset size.
_MERGE_FANIN = 16

SUPPORTED_COMPRESSIONS = PARQUET_COMPRESSIONS


@dataclass(slots=True)
class _Entry:
    ref_base: str
    dp: int = 0
    has_disg: bool = False
    alts: dict[str, list[int]] = field(default_factory=dict)


class DisagreementAggregator:
    """Memory-bounded aggregation of overlapping-base observations.

    Observation unit: one eligible read base from a mate whose partner also
    has an eligible aligned base at the same reference position.
    """

    def __init__(self, cfg: PairLoomConfig) -> None:
        self.cfg = cfg
        self._entries: dict[tuple[int, int], _Entry] = {}
        self._spill_paths: list[str] = []
        self._spill_dir: str | None = None
        self._spill_seq = 0
        self.reported_positions = 0

    def set_spill_dir(self, path: str) -> None:
        self._spill_dir = path

    def add_observation(
        self,
        tid: int,
        rpos0: int,
        ref_base: str,
        base: str,
        strand_reverse: bool,
        pair_disagrees: bool,
    ) -> None:
        key = (tid, rpos0)
        entry = self._entries.get(key)
        if entry is None:
            if len(self._entries) >= self.cfg.max_aggregation_records:
                self.spill()
            entry = _Entry(ref_base)
            self._entries[key] = entry
        elif entry.ref_base != ref_base:  # pragma: no cover - guards a real bug
            raise DataInvariantError(
                f"inconsistent reference base at tid {tid} pos {rpos0 + 1}: "
                f"{entry.ref_base!r} vs {ref_base!r}"
            )
        entry.dp += 1
        if base != ref_base:
            counts = entry.alts.setdefault(base, [0, 0])
            counts[1 if strand_reverse else 0] += 1
        if pair_disagrees:
            entry.has_disg = True

    def spill(self) -> None:
        """Flush in-memory entries to a sorted spill file."""
        if not self._entries:
            return
        if self._spill_dir is None:  # pragma: no cover - pipeline always sets it
            raise OutputError("disagreement spill directory not configured")
        self._spill_seq += 1
        path = os.path.join(self._spill_dir, f"disag_spill_{self._spill_seq:06d}.parquet")
        rows = self._entry_rows(sorted(self._entries.items()))
        # Intermediates stay uncompressed; only the final table honors
        # --parquet-compression.
        _write_spill_rows(path, rows, "none", self.cfg.parquet_batch_size)
        self._spill_paths.append(path)
        self._entries = {}

    def _entry_rows(self, items: list[tuple[tuple[int, int], _Entry]]) -> list[SpillRow]:
        rows: list[SpillRow] = []
        for (tid, rpos), entry in items:
            rows.append((tid, rpos, "", entry.ref_base, entry.dp, entry.has_disg, 0, 0))
            for alt in sorted(entry.alts):
                adf, adr = entry.alts[alt]
                rows.append((tid, rpos, alt, "", 0, False, adf, adr))
        return rows

    def _reduce_spills(self) -> list[str]:
        """Cascading merge: reduce spill files to at most _MERGE_FANIN files.

        Each pass merges groups of spill files into sorted intermediate spill
        files (consuming the inputs), until few enough remain for one final
        bounded-fan-in merge. Intermediate files reuse the spill schema and
        directory.
        """
        paths = list(self._spill_paths)
        while len(paths) > _MERGE_FANIN:
            reduced: list[str] = []
            for start in range(0, len(paths), _MERGE_FANIN):
                group = paths[start : start + _MERGE_FANIN]
                if len(group) == 1:
                    reduced.append(group[0])
                    continue
                self._spill_seq += 1
                merged_path = os.path.join(
                    self._spill_dir, f"disag_spill_{self._spill_seq:06d}.parquet"
                )
                readers = [pq.ParquetFile(p) for p in group]
                try:
                    streams = [_spill_stream(pf, self.cfg.parquet_batch_size) for pf in readers]
                    _write_spill_rows(
                        merged_path,
                        heapq.merge(*streams, key=_ROW_KEY),
                        "none",
                        self.cfg.parquet_batch_size,
                    )
                finally:
                    for pf in readers:
                        pf.close()
                for consumed in group:
                    os.unlink(consumed)
                reduced.append(merged_path)
            paths = reduced
        self._spill_paths = paths
        return paths

    def finalize(self, out_path: Path, ref_names: list[str]) -> int:
        """Merge spills and memory, write the final Parquet table.

        Returns the number of reported positions. Positions without mate
        disagreement never appear.
        """
        paths = self._reduce_spills()
        streams: list[Iterable[SpillRow]] = [self._entry_rows(sorted(self._entries.items()))]
        self._entries = {}
        readers: list[pq.ParquetFile] = []
        try:
            for spill in paths:
                pf = pq.ParquetFile(spill)
                readers.append(pf)
                streams.append(_spill_stream(pf, self.cfg.parquet_batch_size))
            merged = heapq.merge(*streams, key=_ROW_KEY)
            reported = _write_final_table(
                out_path,
                merged,
                ref_names,
                self.cfg.parquet_compression,
                self.cfg.parquet_batch_size,
            )
        finally:
            for pf in readers:
                pf.close()
        self.reported_positions = reported
        return reported


def _spill_stream(pf: pq.ParquetFile, batch_size: int) -> Iterator[SpillRow]:
    """Yield one spill's rows lazily, a Parquet batch at a time.

    Materializing whole spill files as lists would make finalize's memory
    grow with the entire dataset, defeating the spilling bound; generators
    keep the k-way merge lazy.
    """
    for batch in pf.iter_batches(batch_size=batch_size):
        cols = [batch.column(i).to_pylist() for i in range(8)]
        yield from zip(*cols, strict=True)


def _rows_to_table(rows: list[SpillRow]) -> pa.Table:
    return pa.Table.from_arrays(
        [
            pa.array([r[0] for r in rows], type=pa.int32()),
            pa.array([r[1] for r in rows], type=pa.int64()),
            pa.array([r[2] for r in rows], type=pa.string()),
            pa.array([r[3] for r in rows], type=pa.string()),
            pa.array([r[4] for r in rows], type=pa.int64()),
            pa.array([r[5] for r in rows], type=pa.bool_()),
            pa.array([r[6] for r in rows], type=pa.int64()),
            pa.array([r[7] for r in rows], type=pa.int64()),
        ],
        schema=_SPILL_SCHEMA,
    )


def _write_spill_rows(
    path: str, rows: Iterable[SpillRow], compression: str, batch_size: int
) -> None:
    """Write sorted spill rows in batches; never materializes all rows."""
    writer = None
    batch: list[SpillRow] = []
    try:
        writer = pq.ParquetWriter(path, _SPILL_SCHEMA, compression=_codec(compression))
        for row in rows:
            batch.append(row)
            if len(batch) >= batch_size:
                writer.write_table(_rows_to_table(batch))
                batch.clear()
        if batch:
            writer.write_table(_rows_to_table(batch))
    except OSError as exc:
        raise OutputError(f"cannot write disagreement spill {path}: {exc}") from exc
    finally:
        if writer is not None:
            writer.close()


def _write_final_table(
    out_path: Path,
    merged_rows,
    ref_names: list[str],
    compression: str,
    batch_size: int,
) -> int:
    """Consume the sorted merged row stream and write the final table."""
    batch_cols: dict[str, list] = {f.name: [] for f in DISAGREEMENT_SCHEMA}
    reported = 0
    writer = None
    try:
        writer = pq.ParquetWriter(
            str(out_path), DISAGREEMENT_SCHEMA, compression=_codec(compression)
        )
        current_pos: tuple[int, int] | None = None
        dp = 0
        has_disg = False
        alts: dict[str, list[int]] = {}
        ref_base = ""

        def flush_position() -> None:
            nonlocal reported
            if current_pos is None:
                return
            if not has_disg or not alts:
                return
            reported += 1
            tid, rpos = current_pos
            seq_name = ref_names[tid]
            for alt in sorted(alts):
                adf, adr = alts[alt]
                batch_cols["seq_name"].append(seq_name)
                batch_cols["pos"].append(rpos + 1)
                batch_cols["ref"].append(ref_base)
                batch_cols["alt"].append(alt)
                batch_cols["dp"].append(dp)
                batch_cols["adf"].append(adf)
                batch_cols["adr"].append(adr)
                batch_cols["af"].append((adf + adr) / dp if dp else 0.0)
            if len(batch_cols["pos"]) >= batch_size:
                writer.write_table(_batch_to_table(batch_cols))

        for tid, rpos, alt, row_ref, row_dp, row_flag, adf, adr in merged_rows:
            pos_key = (tid, rpos)
            if pos_key != current_pos:
                flush_position()
                current_pos = pos_key
                dp = 0
                has_disg = False
                alts = {}
                ref_base = ""
            if alt == "":
                dp += row_dp
                has_disg = has_disg or row_flag
                ref_base = row_ref
            else:
                counts = alts.setdefault(alt, [0, 0])
                counts[0] += adf
                counts[1] += adr
        flush_position()
        if batch_cols["pos"]:
            writer.write_table(_batch_to_table(batch_cols))
    except OSError as exc:
        raise OutputError(f"cannot write disagreement Parquet {out_path}: {exc}") from exc
    finally:
        if writer is not None:
            writer.close()
    return reported


def _batch_to_table(cols: dict[str, list]) -> pa.Table:
    table = pa.Table.from_arrays(
        [pa.array(cols[f.name], type=f.type) for f in DISAGREEMENT_SCHEMA],
        schema=DISAGREEMENT_SCHEMA,
    )
    for name in cols:
        cols[name].clear()
    return table


def _codec(compression: str) -> str | None:
    return None if compression == "none" else compression


def write_read_stats_parquet(
    out_path: Path,
    counts: dict[str, int],
    compression: str,
) -> None:
    """Write the six classification counts as a deterministic table."""
    ordered = [(name, int(counts[name])) for name in CLASSIFICATION_ORDER]
    table = pa.Table.from_arrays(
        [
            pa.array([name for name, _ in ordered], type=pa.string()),
            pa.array([count for _, count in ordered], type=pa.int64()),
        ],
        schema=READ_STATS_SCHEMA,
    )
    try:
        pq.write_table(table, str(out_path), compression=_codec(compression))
    except OSError as exc:
        raise OutputError(f"cannot write read-stats Parquet {out_path}: {exc}") from exc
