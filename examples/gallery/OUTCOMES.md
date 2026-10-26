# Classic gallery — one template per PairLoom outcome

`classic.fa` (5 contigs, deterministic seed 20261026) and `classic.bam`
(48 records, 24 templates) exercise every success path, every
classification, and every merge-failure reason that PairLoom can surface.
Sequences are reference-oriented in storage for both strands; strand
lives in the flags (the PairLoom/bwa convention).

Rebuild:

```bash
pixi run -- python examples/build_classic_bam.py --out-dir examples/gallery
pixi run -- samtools faidx examples/gallery/classic.fa
```

Reference design: `ctgA`/`ctgB` are random 2 kb contigs; `ctgDup1` and
`ctgDup2` are identical (rescue ambiguity); `ctgRep` holds 11 exact
copies of a 40 bp unit (rescue truncation at the default cap of 10).

## Default run (`--mode merge`, no flags)

mergeable=5, both_aligned=15 (15 distinct reasons × 1), left/right=1/1,
neither=2, merged_non_primary=1.

| Template | Outcome | Evidence on the record |
| --- | --- | --- |
| ok_overlap | MERGED, 15 bp overlap | `ZC:Z:mergeable`, MAPQ 60 |
| ok_adjacent | MERGED, zero-overlap adjacent join | `ZC:Z:mergeable`, MAPQ 60 |
| ok_equal_qual_N | MERGED, equal-quality conflict | consensus carries `N` at the conflict (`MD:Z:15C9`) |
| ok_rescued | MERGED via joint rescue | `ZR:i:1`, MAPQ 0, mates were MAPQ 0 |
| ok_nonprimary | MERGED despite supplementary | `ZN:i:1`, MAPQ 0, supplementary dropped |
| fail_cross_contig | different_contig | primaries on ctgA / ctgB |
| fail_chimeric | chimeric_read | 20M10S primaries, clips = other contig prefix |
| fail_gap | fragment_gap | mates 160 bp apart, span ≤ 1000 |
| fail_frag_len | fragment_length_out_of_bounds | span 1520 > 1000 |
| fail_orientation | orientation_incompatible | same-strand mates |
| fail_indel | indel_conflict | 10M2D10M vs exact 20M |
| fail_insertion | insertion_conflict | GG vs CC insertions at one locus |
| fail_hardclip | hard_clipped_mate | 5H20M primary |
| fail_noqual | qualities_missing | QUAL absent |
| fail_refskip | reference_skip_unsupported | 10M5N10M |
| fail_rg | read_group_conflict | RG rgA vs rgB |
| fail_rescue_unsupported | rescue_unsupported_alignment | MAPQ 19 + gapped CIGAR |
| fail_rescue_ambiguous | rescue_ambiguous | blocks on ctgDup1 AND ctgDup2 |
| fail_rescue_no_pair | rescue_no_compatible_pair | unique blocks, cross-contig, MAPQ 0 |
| fail_rescue_truncated | rescue_search_truncated | block occurs 11× on ctgRep (cap 10) |
| cls_left_only / cls_right_only | one mate unmapped | mapped mate `ZC`-tagged, normalized |
| cls_neither | both unmapped | pass-through, `ZC:Z:neither_aligned` |
| cls_only_nonprimary | only a secondary record | pass-through, untagged, neither_aligned |

## Flag runs that surface the policy-gated reasons

| Flag | Extra reasons surfaced | Notes |
| --- | --- | --- |
| `--no-rescue` | `mapq_below_min` ×3 | the three MAPQ-0 pairs skip rescue; cross-contig MAPQ-0 pairs keep `different_contig` |
| `--no-merge-non-primary` | `non_primary_records` ×1 | ok_nonprimary refused instead of merged; `merged_non_primary=0` |
| `--conflict-policy strict` | `strict_conflict` ×1 | ok_equal_qual_N refused |
| `--min-overlap 30` | `insufficient_overlap` ×9 | every gallery overlap is < 30 bp, so the overlap gate fires before the gap/length/indel gates |

## Not exercised

`empty_consensus` is a defensive branch (`merge.py:219-220`: a merged
CIGAR with no M/I/D operations — a fragment consisting only of clips).
No legitimate alignment produces it; the gallery does not attempt a
degenerate construction.
