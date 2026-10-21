# Synthetic example dataset — expected outputs

Regenerate with `pixi run python scripts/build_synthetic_dataset.py`, then:

```bash
pairmerge --input data/synthetic/input.bam --reference data/synthetic/ref.fa \
  --mode merge --merged-bam /tmp/m.bam --unmerged-bam /tmp/u.bam \
  --disagreement-parquet /tmp/d.parquet --read-stats-parquet /tmp/s.parquet --force
```

The `expected_*` files in this directory are the committed reference
outputs produced by the current version and validated by the test suite.

## Input (5 templates)

| template          | layout                                                      |
| ----------------- | ----------------------------------------------------------- |
| pair_mergeable    | R1 [500,520) fwd with 1 SNV at 516; R2 [510,530) rev        |
| pair_repeat       | both mates 16bp inside the [400,464) tandem repeat          |
| pair_far          | R1 [600,620), R2 [650,670): no overlap                      |
| pair_half         | R1 mapped [620,640); R2 unmapped                            |
| pair_none         | both mates unmapped                                         |

## Expected results

- Classification counts: `total=5, mergeable=2, left_aligned_only=1,
  right_aligned_only=0, both_aligned=1, neither_aligned=1`.
- `pair_mergeable` merges to a 30 bp consensus at chr1:501; the SNV at
  chr1:516 is an equal-quality conflict, so the consensus carries `N` at
  Q2 there; agreeing overlap bases combine to Q60; `NM:i:1`,
  `MD:Z:15A14`.
- `pair_repeat` is jointly adjusted to (400, 400) — the leftmost
  full-overlap candidate pair in the period-4 repeat — and merges to
  `16M`, `NM:i:0`, `MD:Z:16`; the consensus carries the conservative
  MAPQ 3 (min of the two adjusted mate MAPQs — eligibility gates on the
  original aligner MAPQ 60, the output never does).
- `pair_far` stays unmerged with reason `insufficient_overlap`;
  `pair_half` keeps its two records with `ZC:Z:left_aligned_only`;
  `pair_none` records carry `ZC:Z:neither_aligned`.
- Disagreement table: exactly one row — `chr1`, pos 516, ref base `A`,
  alt `C`, `dp=2`, `adf=1`, `adr=0`, `af=0.5` — computed from the two
  original mates, never from the consensus.
