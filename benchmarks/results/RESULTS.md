# PairLoom benchmark results

- pairs: 20,000, reference: 2,000,000 bp, chunk size: 256
- input: 20% repeat-island pairs, 2% both-unmapped templates, 0.2%/base SNVs

## Scaling

| threads | seconds | pairs/s | peak RSS (MiB) |
| --- | --- | --- | --- |
| 1 | 7.7 | 2,590 | 225 |
| 2 | 5.0 | 4,005 | 259 |
| 4 | 3.5 | 5,663 | 260 |

## Stage profile (single thread, cumulative seconds)

| stage | seconds |
| --- | --- |
| bam_read_and_grouping | 0.23 |
| candidate_discovery | 2.63 |
| pair_optimization | 1.43 |
| overlap_analysis | 2.48 |
| consensus_merge | 5.12 |
| parquet_aggregation | 0.92 |
| record_processing | 14.28 |
