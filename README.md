# vcftools-rs

A faster reimplementation of the population-genetics statistics in
[VCFtools](https://vcftools.github.io/) (Danecek et al., *Bioinformatics* 2011),
producing **byte-identical output** to VCFtools 0.1.17.

VCFtools is still the de facto tool for Weir & Cockerham Fst, nucleotide
diversity (π), Tajima's D and heterozygosity-based inbreeding coefficients
from VCF files: its paper is cited over 2,000 times a year (2,658 in 2024,
per OpenAlex), and these statistics are not available in bcftools. But
VCFtools is single-threaded and slow on large datasets. vcftools-rs gives
the same results, byte for byte, 6–10× faster on one core and ~45× faster
on 28 cores (1000 Genomes chr22), so existing pipelines and published
analyses can be reproduced unchanged, only faster.

> **Unofficial.** This project is not affiliated with or endorsed by the
> VCFtools authors. All statistical methods are theirs; please cite the
> [VCFtools paper](https://doi.org/10.1093/bioinformatics/btr330) (see
> `CITATION.cff`). Report problems with this reimplementation here, not to
> the VCFtools issue tracker.

## Install

Prebuilt Linux binaries are attached to each
[GitHub release](https://github.com/realwalju/vcftools-rs/releases). To build
from source (Rust 1.80+ and a C compiler):

```bash
cargo install --git https://github.com/realwalju/vcftools-rs vcftools-rs
```

Usage is the same as VCFtools for the supported options, e.g.
`vcftools-rs --gzvcf in.vcf.gz --weir-fst-pop popA.txt --weir-fst-pop popB.txt --out result`.
`vcftools-rs --help` lists everything supported.

## Scope

Statistics (VCFtools option names and output files):

| Option | Output |
|---|---|
| `--freq`, `--counts` | `.frq`, `.frq.count` |
| `--het` | `.het` |
| `--hardy` | `.hwe` |
| `--missing-site`, `--missing-indv` | `.lmiss`, `.imiss` |
| `--site-pi`, `--window-pi N [--window-pi-step N]` | `.sites.pi`, `.windowed.pi` |
| `--TajimaD N` | `.Tajima.D` |
| `--weir-fst-pop F --weir-fst-pop F [--fst-window-size N --fst-window-step N]` | `.weir.fst`, `.windowed.weir.fst` |

Sample filters: `--keep`, `--remove`, `--indv`, `--remove-indv`.

Site filters: `--chr`, `--not-chr`, `--from-bp`, `--to-bp`, `--positions`,
`--exclude-positions`, `--snp`, `--snps`, `--exclude`, `--remove-indels`,
`--keep-only-indels`, `--min-alleles`, `--max-alleles`, `--minQ`,
`--remove-filtered-all`, `--remove-filtered`, `--keep-filtered`, `--phased`,
`--maf`, `--max-maf`, `--non-ref-af[-any]`, `--max-non-ref-af[-any]`,
`--max-missing`, `--mac`, `--max-mac`, `--non-ref-ac[-any]`,
`--max-non-ref-ac[-any]`, `--max-missing-count`, `--hwe`, `--min-meanDP`,
`--max-meanDP`.

Genotype filters (calls are marked filtered, not removed): `--minDP`,
`--maxDP`, `--minGQ`.

Input: plain VCF, BGZF (`bgzip`) or ordinary gzip. `--threads N` sets the
number of worker threads (default: all cores).

Differences from VCFtools: the run log goes to standard error only (no
`.log` file), and only one statistic may be requested per run (as in
VCFtools). Byte-identity is validated on Linux (glibc); number formatting
comes from the platform's C library.

Not implemented: BCF input, `--bed`/`--exclude-bed`, `--thin`, `--mask`,
INFO-flag site filters, genotype FILTER-flag filters
(`--remove-filtered-geno[-all]`), `--derived`, and VCFtools' other output
types.

## Relationship to other tools

- **bcftools** is VCFtools' modern successor for VCF manipulation. Its
  `+fill-tags` plugin computes allele frequencies/counts and Hardy-Weinberg
  statistics (as INFO fields, in a different format), but it has no Fst,
  nucleotide diversity, Tajima's D or inbreeding-coefficient calculations.
  For `--freq`, `--counts`, `--hardy` and missingness, vcftools-rs mainly
  offers VCFtools-format compatibility; for the population-genetics
  statistics it is a faster route to the tool people already use.
- **plink2** computes frequencies, HWE, heterozygosity, missingness and Fst
  quickly, with its own formats and numerics (see the benchmark note below).
- **pixy** (Korunes & Samuk 2021) estimates π and dxy correctly in the
  presence of missing data from all-sites VCFs. VCFtools' windowed π does
  not distinguish missing from invariant sites and can be biased when data
  are missing; vcftools-rs reproduces VCFtools' estimate exactly, bias
  included. For new analyses of π with substantial missing data, consider
  pixy; use vcftools-rs to reproduce or speed up VCFtools-based work.

## How it is faster

1. **Parallel decompression.** BGZF files are independent ~64 KB deflate
   blocks; batches of blocks are decompressed concurrently (libdeflate).
2. **Parallel, order-preserving processing.** Decompressed text is cut into
   chunks of whole lines, processed on all cores, and results are consumed
   in file order. Statistics that sum floating-point values across sites
   (E(HOM) in `--het`, window sums) are accumulated strictly in file order
   so rounding matches VCFtools bit-for-bit.
3. **Cheap parsing.** Rows whose genotypes are all 3-byte `a|b`/`a/b`
   columns are validated in one vectorisable pass and decoded without
   per-column branching; allele counts use a branch-free histogram.
4. **Buffered output** instead of flushing after every line.

## Faithfulness

VCFtools' arithmetic is transcribed operation-for-operation, including its
quirks: `std::accumulate` with an `int` seed truncating Fst sample sizes,
unsigned wrap-around in π, `std::min` NaN behaviour in the MAF filter, the
implicit `--keep` performed by `--weir-fst-pop`, chromosome ordering of
windowed output, and C `printf` number formatting (`-nan` included).

This includes two confirmed VCFtools issues, reproduced deliberately so
results stay identical: `--site-pi` overflows above ~23,000 diploid
samples, and `--non-ref-af-any` has no effect on its own. See
[docs/upstream-notes.md](docs/upstream-notes.md).

## Verification

- `scripts/compare.sh [--threads N] <names>` runs the benchmark commands on
  1000 Genomes chr22 and checks each output with `cmp` against VCFtools'.
- `scripts/fuzz.py CASES SEED` is a differential fuzz test: it generates
  random messy VCFs (missing/half-missing/haploid genotypes, phased/unphased
  mixes, multi-allelic sites, indels, GT not first in FORMAT, DP/GQ values
  including missing and GQ > 99, FILTER/QUAL/ID variety, repeated
  chromosomes; plain, BGZF and gzip) with random sample, site and genotype
  filters, and requires identical output (or both tools failing).
  `scripts/fuzz_many.sh CASES SEED...` runs several seeds.

## Benchmark

1000 Genomes Phase 3 chr22 (1,103,547 variants × 2,504 samples), Intel
i7-14700K (20 cores / 28 threads), WSL2 Ubuntu 24.04. Wall-clock seconds,
median of 3 runs (monotonic clock). Every output is byte-identical to
VCFtools 0.1.17.

| Command | VCFtools | vcftools-rs, 1 thread | vcftools-rs, 28 threads |
|---|---:|---:|---:|
| `--freq` | 85.8 | 10.7 (8.0×) | 1.87 (45.9×) |
| `--counts` | 85.4 | 10.6 (8.1×) | 1.80 (47.4×) |
| `--het` | 96.1 | 14.2 (6.8×) | 2.11 (45.5×) |
| `--hardy` | 93.5 | 15.2 (6.2×) | 2.13 (43.9×) |
| `--missing-site` | 80.7 | 8.5 (9.5×) | 1.77 (45.6×) |
| `--missing-indv` | 78.5 | 7.7 (10.1×) | 1.86 (42.2×) |
| `--site-pi` | 89.3 | 11.9 (7.5×) | 1.90 (47.0×) |
| `--window-pi 10000` | 86.5 | 11.6 (7.4×) | 1.93 (44.8×) |
| `--TajimaD 10000` | 85.9 | 11.7 (7.4×) | 1.89 (45.5×) |
| `--weir-fst-pop` (CEU vs YRI) | 83.1 | 9.3 (9.0×) | 1.77 (46.9×) |
| `--weir-fst-pop` windowed | 83.1 | 9.2 (9.1×) | 1.78 (46.7×) |

For context, plink2 2.0.0a6.9 computing the overlapping statistics
(`--freq`, `--hardy`, `--het`, `--missing`) from the same VCF takes about
10.5 s on 1 thread and 5.7 s on 28 threads (dominated by VCF import); its
output format and numerics differ from VCFtools. See
`results/plink2_context.tsv`.

Raw tables: `results/baseline.tsv`, `results/summary_1thread.txt`,
`results/summary_28threads.txt`. Reproduce: `scripts/fetch_data.sh`,
`scripts/baseline.sh`, `scripts/final.sh`, `scripts/plink2_context.sh`.

### Second dataset: 1001 Genomes *Arabidopsis thaliana*

First 264,200 variants of Chr1 × 1,135 accessions from the 1001 Genomes
release v3.1 VCF (every call is `GT:GQ:DP`, with real missing data), Fst
between the `germany` and `western_europe` ADMIXTURE groups. All outputs
are byte-identical to VCFtools, with and without `--minDP 3 --minGQ 20
--max-missing 0.8`. Median of 3 runs, seconds:

| Command | VCFtools | vcftools-rs, 1 thread | vcftools-rs, 28 threads |
|---|---:|---:|---:|
| `--freq` | 13.5 | 4.9 (2.7×) | 0.74 (18×) |
| `--het` | 13.9 | 6.4 (2.2×) | 0.86 (16×) |
| `--hardy` | 14.1 | 5.8 (2.4×) | 0.81 (17×) |
| `--missing-site` | 12.8 | 4.4 (2.9×) | 0.70 (18×) |
| `--site-pi` | 13.4 | 5.0 (2.7×) | 0.72 (19×) |
| `--TajimaD 10000` | 12.8 | 4.8 (2.7×) | 0.71 (18×) |
| `--weir-fst-pop` | 13.4 | 5.0 (2.7×) | 0.72 (19×) |
| `--freq --minDP 3 --minGQ 20 --max-missing 0.8` | 34.5 | 18.3 (1.9×) | 1.69 (20×) |

Single-thread gains are smaller than on 1000 Genomes because these
genotype columns are not bare `a|b` calls, so the vectorised fast path does
not apply; the genotype filters also re-scan each sample column for DP and
GQ separately (a known optimisation target). Reproduce with
`scripts/fetch_athal.sh` and `scripts/athal_bench.sh`
(`results/athal_benchmark.tsv`).

## Credit

All statistical methods and their exact numerical behaviour are those of
VCFtools by Adam Auton, Petr Danecek and contributors (LGPL-3.0). This
reimplementation was written with AI assistance (Claude) and is validated
by output comparison against VCFtools.
