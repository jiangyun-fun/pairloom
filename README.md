# PairLoom

Pair-aware BAM realignment, overlap analysis, and paired-end read merging.

> Work in progress — this README is expanded as implementation completes.

## Quick start

```bash
pixi run pairmerge --input input.bam --reference ref.fa --mode adjust \
  --adjusted-bam adjusted.bam \
  --disagreement-parquet disagreements.parquet \
  --read-stats-parquet stats.parquet

pixi run pairmerge --input input.bam --reference ref.fa --mode merge \
  --merged-bam merged.bam --unmerged-bam unmerged.bam \
  --disagreement-parquet disagreements.parquet \
  --read-stats-parquet stats.parquet
```
