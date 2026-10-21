# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 7.6 | 2,643 | 224 |
| 2 | 4.7 | 4,280 | 258 |
| 4 | 3.7 | 5,396 | 259 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.22 |
| candidate_discovery | 2.63 |
| pair_optimization | 1.43 |
| overlap_analysis | 2.48 |
| consensus_merge | 5.07 |
| parquet_aggregation | 0.91 |
| record_processing | 14.09 |
