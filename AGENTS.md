# AGENTS.md

Guidance for coding agents working in this repository. PairLoom is a Python
package (CLI `pairmerge`) for pair-aware BAM realignment, overlap-disagreement
statistics, and paired-end read merging. Correctness, determinism, bounded
memory, and fail-loud error handling are design constraints, not preferences.

## Environment and commands

The project is a pixi workspace (`pixi.toml` + `pixi.lock`). Run everything
through pixi so the workspace env — not the ambient PATH — supplies binaries:

| Task | Command |
| ---- | ------- |
| Set up env (conda deps + editable install) | `pixi install` |
| Full test suite (151 tests) | `pixi run test` |
| Stop at first failure | `pixi run test-fast` |
| Lint (ruff check) | `pixi run lint` |
| Format check / apply | `pixi run format` / `pixi run format-fix` |
| Build wheel | `pixi run build` |
| Benchmark (1/2/4 threads, writes `benchmarks/results/`) | `pixi run bench` |
| CLI smoke | `pixi run cli` |

Python >= 3.12, < 3.14; `pysam` >= 0.24, `pyarrow` >= 25, `numpy` >= 2.5.

Editable install gotcha: hatchling is configured with `dev-mode-dirs = ["src"]`
in `pyproject.toml`. If imports suddenly stop resolving, check for
`_pairloom.pth` in the pixi env's site-packages before debugging anything else.

## Repository layout

- `src/pairloom/` — the package, one concern per module:
  - `cigarops.py` CIGAR math/validation; `alignment.py` placements and
    reference-orientation read structures; `candidates.py` exact-block
    candidate discovery; `pairing.py` joint R1/R2 selection;
  - `mapq.py` the conservative MAPQ policy table; `overlap.py` overlap
    observations; `merge.py` consensus merging; `classify.py` template
    classification and tagging;
  - `records.py` template grouping and record-policy validation;
    `reference.py` cached FASTA access; `parquetout.py` spilling statistical
    output; `validate.py` preflight path checks; `pipeline.py` streaming
    orchestration, parallelism, atomic publish; `cli.py` argparse entry.
- `tests/` — pytest; `conftest.py` builds a canonical reference with
  engineered repeats (constants `REPEAT1`, `REPEAT2`, `LONG_REPEAT`, `UNIQUE`)
  and an in-process CLI runner `fx.run(...)`. All expected values are derived
  from those reference strings.
- `benchmarks/benchmark.py` — reproducible performance measurement.
- `scripts/build_synthetic_dataset.py` — regenerates `data/synthetic/`.
- `docs/plans/` — plan/report pairs, the durable record of substantive work;
  `docs/misc/` — requirements audit and rolling next/recommendation state.
  Follow that pattern for middle/large changes.

## Critical invariants (each was violated once; do not repeat)

1. **SAM orientation**: BAM SEQ/CIGAR/QUAL are stored in REFERENCE orientation
   for both strands. Never revcomp a reverse-strand record before comparing
   to the reference. (Verified against `samtools calmd`.)
2. **pysam API**: `AlignedSegment` has no `remove_tag()` in pysam 0.24. Remove
   a tag with `set_tag(tag, None)`.
3. **MAPQ policy is conservative and split by role**: eligibility filters
   (`--min-mapq`, `--min-base-quality`) always test the ORIGINAL aligner MAPQ;
   output records carry the adjusted values; merged records carry
   `min(adjusted R1, adjusted R2)`. Never let an ambiguous repeat placement
   earn MAPQ 60 — including via truncated searches (candidate cap or
   `--max-pair-evaluations` exhaustion must flag truncation).
4. **TLEN**: opposite signs within a template; on equal start coordinates
   read1 is positive.
5. **Mate fields are cross-record**: after any adjustment, refresh RNEXT/PNEXT
   (including the unmapped mate of an orphan) and reciprocal `MQ`.
6. **Publish semantics in `run_pipeline`**: temp files then `os.replace`;
   pre-existing outputs are backed up and restored on failure; after the last
   publication the run is COMMITTED and no later failure may delete outputs —
   see the `committed` flag before touching the exception handlers.
7. **Fail loud, never swallow**: no `except: pass`, no silent record skipping.
   An AST scan test enforces this. Malformed input raises typed
   `PairLoomError`s with read/reference context.
8. **Determinism**: output must be byte-identical across `--threads` counts.
   Keep statistics as integer aggregations and tie-breaks coordinate-based.
9. **Record parsing edges**: normalize the record CIGAR before walking it
   (adjacent identical ops share a boundary and would overwrite each other);
   expand SAM `=` SEQ shorthand against the reference before scoring,
   statistics, or merging (clipped/inserted `=` stays literal). Note:
   pysam's `query_sequence` setter drops `query_qualities` — re-set them.

## Testing and regression discipline

- Every bug fix lands with a regression test derived from its reproduction;
   the suite exists to pin behavior, not to chase coverage numbers.
- When behavior changes deliberately, update the affected expectations AND the
  committed synthetic outputs (`scripts/build_synthetic_dataset.py`, then run
  `pairmerge` per `data/synthetic/EXPECTED.md`).
- The README is the source of truth for algorithms, policies, schemas, and
  defaults — update it in the same change as any behavior change.

## Benchmarks and performance

Measured, not claimed. After changes on hot paths, run `pixi run bench` and
keep single-thread throughput within the range documented in the README's
Performance section. Tracked results are `benchmarks/results/RESULTS.md` and
`benchmark.parquet`; regenerated run data under `benchmarks/results/` is
gitignored.

## Git conventions

Conventional Commits (`fix:`, `feat:`, `test:`, `docs:`, `chore:`, ...), no
attribution trailers. One logical change per commit; keep the tree green at
every commit. Do not push unless asked.
