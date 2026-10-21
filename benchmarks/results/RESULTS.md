# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 10.0 | 1,993 | 225 |
| 2 | 7.5 | 2,679 | 259 |
| 4 | 4.6 | 4,341 | 261 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.27 |
| candidate_discovery | 3.15 |
| pair_optimization | 1.81 |
| overlap_analysis | 2.90 |
| consensus_merge | 6.19 |
| parquet_aggregation | 1.20 |
| record_processing | 17.47 |
