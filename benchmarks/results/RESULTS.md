# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 9.3 | 2,154 | 225 |
| 2 | 5.6 | 3,580 | 259 |
| 4 | 4.1 | 4,900 | 260 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.23 |
| candidate_discovery | 2.71 |
| pair_optimization | 1.56 |
| overlap_analysis | 2.39 |
| consensus_merge | 8.73 |
| parquet_aggregation | 0.93 |
| record_processing | 18.22 |
