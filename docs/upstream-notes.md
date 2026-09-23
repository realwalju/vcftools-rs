# Notes on VCFtools 0.1.17 behaviour

vcftools-rs reproduces VCFtools' output exactly, **including the behaviours
below**, so that it can be swapped in without changing results. They are
listed here so users know about them. Items marked *confirmed* were
reproduced against the Bioconda build of VCFtools 0.1.17 with
`scripts/probe_upstream_bugs.sh`.

## Confirmed

### `--site-pi` gives wrong values above ~23,000 diploid samples

`output_per_site_nucleotide_diversity` computes the number of pairwise
comparisons as `int pairs = total_alleles * (total_alleles - 1)`, and the
mismatch count in an `int`. Both overflow 32-bit integers once a site has
more than 46,340 called chromosomes. Example (one site, half `0|0`, half
`1|1`, expected π ≈ 0.5):

| Samples | VCFtools `PI` | Correct |
|---:|---:|---:|
| 20,000 | 0.500013 | 0.500013 |
| 24,000 | -0.578599 | 0.50001 |
| 40,000 | -0.520186 | 0.500006 |

Related unsigned 32-bit products (`--window-pi` comparison counts,
`--TajimaD`'s `n*n` and `ui*ui`) wrap at larger sizes (over ~32,768
diploid samples); not yet separately demonstrated.

### `--non-ref-af-any` / `--max-non-ref-af-any` have no effect on their own

In `filter_sites_by_frequency_and_call_rate`, the "all alternate alleles
failed" test checks `min_non_ref_af`/`max_non_ref_af` instead of the
`_any` variants, so the `_any` thresholds are only applied when a plain
`--non-ref-af`/`--max-non-ref-af` is also given. On the chr22 test
fixture, `--non-ref-af-any 0.9` keeps all 1,500 sites. The allele-count
counterpart (`--non-ref-ac-any`) uses the right variables and works.

## Other faithfully reproduced details

- `--weir-fst-pop` also acts as `--keep` for the listed samples.
- Weir & Cockerham Fst sums per-population sample sizes with
  `std::accumulate(..., 0)` (an `int` seed), truncating after every
  addition. Sample sizes are normally whole numbers, so this rarely matters.
- `--missing-site` treats a phased call with a missing second allele
  (`0|.`) as a haploid call.
- `--from-bp`/`--to-bp` require exactly one `--chr`; repeated `--chr` of
  the same name counts once.
- Calls removed by `--minDP`/`--maxDP`/`--minGQ` count as missing
  chromosomes for `--max-missing`, and calls without DP are removed by
  `--minDP`.
