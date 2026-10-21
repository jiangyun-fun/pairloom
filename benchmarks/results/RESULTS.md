# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 7.6 | 2,633 | 224 |
| 2 | 4.8 | 4,170 | 258 |
| 4 | 3.8 | 5,223 | 259 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.22 |
| candidate_discovery | 2.57 |
| pair_optimization | 1.41 |
| overlap_analysis | 2.43 |
| consensus_merge | 4.92 |
| parquet_aggregation | 0.88 |
| record_processing | 13.79 |
