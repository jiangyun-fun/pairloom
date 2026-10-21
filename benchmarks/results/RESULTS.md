# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 7.8 | 2,562 | 225 |
| 2 | 5.0 | 3,985 | 260 |
| 4 | 3.6 | 5,534 | 261 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.24 |
| candidate_discovery | 2.70 |
| pair_optimization | 1.46 |
| overlap_analysis | 2.58 |
| consensus_merge | 5.20 |
| parquet_aggregation | 0.95 |
| record_processing | 14.56 |
