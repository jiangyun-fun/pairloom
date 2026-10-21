# PairLoom

Pair-aware BAM realignment, overlap analysis, and paired-end read merging.

PairLoom (`pairmerge` CLI) processes paired-end BAM alignments against a
reference and can:

1. adjust mate placements jointly, deterministically resolving equivalent
   placements in repetitive reference sequence;
2. recompute MAPQ conservatively after adjustment (never overstating
   uniqueness for repeated reads);
3. detect and count overlapping-base disagreements between mates as Parquet;
4. classify every template and tag it in the BAM;
5. merge overlapping mates into a single consensus record with a correct,
   reconstructed CIGAR and preserved outer soft clips.

Correctness, determinism, bounded memory, and fail-loud error handling are
design constraints, not afterthoughts: every unexpected condition aborts the
run with a nonzero exit code, and outputs are published atomically.

## Installation

```bash
git clone <repo> && cd pairloom
pixi install          # conda deps + editable install of the package
pixi run pairmerge --help
```

Requirements: Python ≥ 3.12, `pysam` ≥ 0.22 (HTSlib), `pyarrow` ≥ 16,
`numpy`. The reference FASTA must carry a `.fai` index
(`samtools faidx ref.fa`). Tests additionally need `samtools` and `pytest`
(both included in the pixi dev environment).

## Quick start

```bash
# Adjust (and tag) alignments only; mates stay separate records:
pairmerge --input input.bam --reference ref.fa --mode adjust \
  --adjusted-bam adjusted.bam \
  --disagreement-parquet disagreements.parquet \
  --read-stats-parquet stats.parquet

# Merge overlapping pairs into consensus reads:
pairmerge --input input.bam --reference ref.fa --mode merge \
  --merged-bam merged.bam --unmerged-bam unmerged.bam \
  --disagreement-parquet disagreements.parquet \
  --read-stats-parquet stats.parquet
```

Every output file has its own explicit path argument; no filename is ever
inferred from a prefix or directory, and outputs never overwrite the inputs.
Existing output files are refused unless `--force` is given.

## Record policy

- **Primary alignments** drive adjustment, classification, and merging.
- **Secondary (0x100) and supplementary (0x800) records** pass through to
  the outputs unchanged and untagged. A template carrying any non-primary
  record is never merged (reason `non_primary_records`): merging replaces
  the template's primary records, so this policy discards nothing and
  duplicates nothing.
- **Unmapped records** pass through and are tagged with their template's
  classification.
- **Templates with only non-primary records** (possible in filtered BAMs)
  have no primary alignment evidence: they pass through unchanged and
  untagged and classify as `neither_aligned`.
- **Orphans** (a mapped primary whose mate is missing or unmapped) receive
  single-read leftmost normalization (spec §3.1). With an unmapped mate
  present, the unmapped record's locus follows the mapped read.
- A template with more than two primary records violates pairing invariants
  and fails loudly (`DataInvariantError`).
- Coordinate-sorted inputs are query-name-sorted to a temporary file for
  streaming mate grouping; query-name-sorted inputs stream directly.
  Outputs are written in template order unless `--sort-output`.

## Algorithms

### Repeat-aware candidate discovery

BAMs carry one primary placement per read even when the read matches
repeated sequence. For every mapped primary, PairLoom verifies the
alignment base by base against the reference (CIGAR `M` is never trusted as
proof of a match) and discovers candidates:

- `--candidate-source bam`: only the placement present in the BAM;
- `--candidate-source reference` (default): additionally slide the read's
  exact-match ungapped block along the reference inside
  `--realign-window` (default 500 bp) per side, capped at
  `--max-candidates-per-read` (default 10).

Only exact-match ungapped blocks are slid, so sliding can never override
alignment evidence: every candidate of such a read is equivalent by
construction. Mismatch- or indel-bearing alignments keep their BAM
placement. Discovery is bounded and never claims global uniqueness; hitting
the candidate cap flags the read *truncated* and triggers a conservative
MAPQ cap.

### Joint R1/R2 optimization

Candidate pairs are evaluated against the documented preference hierarchy:

1. both mates on the same reference (candidates always are; cross-contig
   pairs therefore keep their BAM placements);
2. library orientation (`--library-orientation`, default `fr`:
   leftmost mate forward, rightmost reverse, by reference position) and
   fragment length within `--min-fragment-length`/`--max-fragment-length`
   (default 0–1000);
3. combined alignment score: pairs scoring below
   `best − --alignment-score-tolerance` (default 0) are discarded as
   clearly inferior;
4. largest aligned reference overlap;
5. deterministic tie-breaks: leftmost fragment start, rightmost fragment
   end, then R1/R2 start coordinates.

`--pair-strategy pair-aware` (default) applies the score-tolerance filter
first and then prefers overlap within the surviving set. `--pair-strategy
max-overlap` skips the filter and prefers overlap among all
constraint-satisfying pairs. When no candidate combination satisfies the
constraints, the **original** BAM placements are kept — mates are never
slid independently. The 0x2 (proper pair) flag is set exactly when a valid
candidate combination was chosen.

### MAPQ policy

| Condition                                            | Adjusted MAPQ        |
| ---------------------------------------------------- | -------------------- |
| ≥ 2 equivalent-best placements remain                | `min(orig, 3)`       |
| search truncated (discovery or pair evaluation)      | `min(orig, 10)`      |
| placement unchanged                                  | `orig`               |
| changed, verified exact, original MAPQ ≥ 30          | `60`                 |
| changed, otherwise                                   | `min(orig, 25)`      |

MAPQ 60 is therefore never assigned merely because an alignment is an
exact match: the placement must also be unambiguous in the searched
window, neither candidate discovery nor pair evaluation
(`--max-pair-evaluations`) may have been truncated, and the original
aligner must have considered the read high quality.
`--force-mapq60-exact` (compatibility mode) assigns 60 to every exact
match and prints an explicit warning about overstating mapping
confidence.

Note: the `--min-mapq` eligibility filters (for overlap observations and
merging) apply to the **original** aligner MAPQ, not to the conservative
values PairLoom writes, so capping our own output MAPQ cannot disqualify a
pair whose input evidence was strong.

Adjusted records keep their strand, flags (0x2 recomputed), mate
coordinates, and TLEN consistently updated; `NM`/`MD` are recomputed for
moved (always exact) placements as `NM:i:0` / `MD:<span>`; stale scoring
tags that cannot be reproduced faithfully (`AS`, `XS`, `SA`, `MC`, `MQ`,
`XA`) are removed rather than guessed at.

### Overlap analysis and disagreement statistics

Two bases overlap when they align to the same reference coordinate
(CIGAR-aware; soft clips, insertions, deletions, and reference skips never
produce observations). Only A/C/G/T observations passing
`--min-base-quality` (default 20) with `--min-mapq` (default 20) are
eligible; strand is taken from the actual alignment orientation. BAM SEQ
is stored in reference orientation for both strands, so comparisons are
strand-trivial and `adf`/`adr` count forward/reverse observations directly.

Only positions where at least one eligible pair contains **different**
overlapping nucleotides are reported. At a reported position, counts use
all eligible overlapping observations, including agreeing ones from other
pairs. `dp` counts read-base observations (not molecules); a merged
consensus base is never an additional observation; statistics are computed
from the two original adjusted mates, before any consensus.

### Consensus merging

Merging walks the reference span (never string concatenation):

- outer soft clips are preserved with qualities: the leftmost read's
  leading clip, and the trailing clip of the mate reaching the fragment's
  rightmost reference position (a contained mate's trailing clip is
  internal and dropped), both in reference orientation; when mates share
  an aligned start or end, the mate supporting the longer clip owns it;
- internal soft clips are dropped (the other mate's aligned bases cover
  those positions);
- overlapping aligned bases appear exactly once;
- deletions observed by one or both mates survive; one mate deleting where
  the other aligns a base is `indel_conflict` (unmergeable);
- insertions merge when both mates insert the identical sequence
  (per-base quality maximum); a single-mate insertion inside the other
  mate's span is `insertion_conflict`;
- uncovered middle positions are never synthesized (`fragment_gap`);
- hard-clipped mates are not merged (`hard_clipped_mate`: their full
  sequence is not recoverable from the BAM);
- mates must pass `--min-mapq`, share one read group (a conflict is
  `read_group_conflict`), satisfy the library orientation model, and
  overlap on at least `--min-overlap` aligned reference bases (default 10).

Quality model (sequencing errors assumed independent):

- agreeing bases: `Q = min(q1 + q2, --max-consensus-quality)`;
- disagreeing bases resolve toward the higher-quality base when the gap
  reaches `--min-quality-difference` (default 5), keeping `Q = min(q_hi,
  --max-consensus-quality)`;
- unresolved disagreements emit `N` at Q2 under `--conflict-policy
  quality` (default) or make the pair unmergeable under `strict`.

Merged records are synthetic single fragments: flag 0, the template's
QNAME, `POS` = leftmost aligned reference coordinate, `MAPQ = min of the
two adjusted MAPQs` (the `--min-mapq` eligibility gate still uses the
original aligner MAPQ), no mate coordinates, TLEN 0, `NM`/`MD` recomputed,
`RG` preserved, the status tag set to `mergeable`. Unmerged outputs never
duplicate a successfully merged template.

## Parquet schemas

`--disagreement-parquet` (ZSTD by default, `--parquet-compression`):

| Column     | Type    | Meaning                                    |
| ---------- | ------- | ------------------------------------------ |
| `seq_name` | string  | reference sequence name                    |
| `pos`      | int64   | 1-based genomic position                   |
| `ref`      | string  | reference nucleotide                       |
| `alt`      | string  | alternative nucleotide                     |
| `dp`       | int64   | eligible overlapping read-base depth       |
| `adf`      | int64   | forward-strand observations supporting ALT |
| `adr`      | int64   | reverse-strand observations supporting ALT |
| `af`       | float64 | `(adf + adr) / dp`                         |

Rows are ordered by reference-header order, position, then alt allele;
positions without mate disagreement never appear; multiple alternative
alleles at one position occupy separate rows. Tables are written in row
batches (`--parquet-batch-size`, default 65536) and stay correct for zero
rows.

`--read-stats-parquet`:

| Column  | Type    | Meaning            |
| ------- | ------- | ------------------ |
| `type`  | string  | classification     |
| `count` | int64   | template count     |

Categories (always all six rows, in this order, mutually exclusive,
`total` = sum of the rest): `total`, `mergeable`, `left_aligned_only`,
`right_aligned_only`, `both_aligned`, `neither_aligned`. Left/right mean
R1/R2 (read number), never physical reference order. `left_aligned_only`
also covers non-paired (0x1-unset) primary records.

## Classification tagging

`--status-tag` (default `ZC`) writes `ZC:Z:<classification>` on every
primary record of a template (both mates share the value; merged records
carry `mergeable`). A pre-existing tag value that is not a PairLoom
classification string is treated as unrelated data: PairLoom refuses to
overwrite it and fails loudly with a suggestion to pick another tag.

## Output integrity and failure policy

- **Fail loudly**: no `except: pass`, no logged-and-ignored errors, no
  silent record skipping, no substituted defaults for inconsistent data.
  Malformed CIGARs, templates with >2 primary records, missing reference
  sequences, corrupt BAMs, and internal invariant violations raise typed
  exceptions (`PairLoomError` hierarchy) with read/reference context.
- Exit codes: 0 success, 2 usage/validation errors, 1 processing
  failures. `--debug` prints full tracebacks; the default output is a
  concise actionable message.
- **Atomic publish**: all BAM/Parquet outputs are written to temp files on
  the destination filesystem and `os.replace`d into place only after every
  writer succeeded. On any failure the temps and any already-published
  files are removed (best effort) — incomplete outputs are never presented
  as complete. Under `--force`, pre-existing outputs are moved aside before
  publication and restored if the run fails, so a failed run never destroys
  the previous results. Index files (`.bai`, only with `--index-output`, which
  implies `--sort-output`) are published with the same discipline and are
  validated like every other output: an existing index requires `--force`,
  and an index path colliding with another requested output is rejected
  before processing starts.
- **Worker failures propagate**: multiprocessing workers' exceptions abort
  the whole run; apparently-successful-but-incomplete output is impossible.
- Expected biological outcomes (no overlap, unmapped mates, conflicting
  indels, quality-gate failures, orientation incompatibilities) are
  **not** errors: they produce classifications and documented merge-rejection
  reasons (summarized on stderr), and the records are retained.

## Performance

Measured, not claimed — run `pixi run bench` (or
`python benchmarks/benchmark.py --pairs 20000 --threads 1,2,4`) to
reproduce. Results below: 20,000 pairs, 2×100 bp reads, 2 Mbp reference
with a 1 kb period-4 repeat island every 50 kb, 20% of pairs drawn from
repeat islands, 2% both-unmapped templates, 0.2%/base SNVs, mode `merge`
(see `benchmarks/results/RESULTS.md` for the freshly generated table):

- 1 thread: ≈ 2.6–4.2k pairs/s, ≈ 120–220 MiB peak RSS;
- 4 threads: ≈ 5.8–10.1k pairs/s (2.2–2.4× speedup), ≈ 260 MiB peak RSS;
- results are byte-identical across thread counts (tested in the suite).

Memory stays bounded: templates stream in query-name order, chunked work
units bound parallel IPC, disagreement aggregation holds at most
`--max-aggregation-records` positions in memory before spilling sorted
Parquet partials, and finalize merges those partials lazily with a bounded
fan-in (cascading merges when spills exceed it) — each open stream holds at
most one Parquet batch of rows, so finalize memory is a function of
configuration, never of dataset size. Overlap analysis fetches only the
intersection of the mates' aligned spans, so distant mates never allocate
the span between them. Reference contigs are cached under a
base budget (contigs larger than the budget stream through HTSlib). The
whole BAM is never loaded.

## Selected defaults for spec ambiguities (§15)

1. Indel-bearing alignments are not slid (only exact ungapped blocks are);
   indel re-anchoring is a documented limitation, not a hidden heuristic.
2. Adjust mode evaluates mergeability for classification/tagging but emits
   no consensus records.
3. Cross-contig and constraint-violating pairs keep their original
   placements and have 0x2 cleared.
4. NM/MD are recomputed only where PairLoom changed the record; unchanged
   records keep their input tags. Merged/adjusted exact placements are
   orientation-independent (`NM:i:0`, `MD:<span>`).
5. `--max-fragment-length` defaults to 1000; `--min-overlap` to 10;
   `--min-base-quality`/`--min-mapq` to 20; window 500; candidates 10;
   batch 65536; aggregation spill threshold 10^6 positions.
6. Non-paired primary records classify as `left_aligned_only`.

## Known limitations

- Only exact-match ungapped blocks slide; equivalent-placement
  normalization for gapped alignments is future work.
- Uniqueness claims are bounded by the search window; MAPQ 60 additionally
  requires original-MAPQ corroboration.
- RNA-seq reference skips (N) inside a mate make the pair unmergeable
  (`reference_skip_unsupported`).
- `H`-clipped mates cannot merge (full sequence unavailable in BAM).
- Disagreement statistics cover SNV-style base observations; insertions,
  deletions, and skips are excluded by design.

## Development

```bash
pixi run test          # full pytest suite (137 tests, all 35 spec scenarios)
pixi run lint          # ruff check
pixi run format        # ruff format --check
pixi run bench         # reproducible benchmark
```

The test suite builds synthetic references with engineered repeats and
asserts exact expected outputs; generated BAMs are validated with
`samtools quickcheck`; Parquet outputs are read back with PyArrow and
compared field by field.

## License

MIT.
