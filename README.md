# vcftools-rs

A faster reimplementation of the population-genetics statistics in
[VCFtools](https://vcftools.github.io/) (Danecek et al., *Bioinformatics* 2011),
producing **byte-identical output** to VCFtools 0.1.17.

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
`--max-non-ref-ac[-any]`, `--max-missing-count`, `--hwe`.

Input: plain VCF, BGZF (`bgzip`) or ordinary gzip. `--threads N` sets the
number of worker threads (default: all cores).

Not implemented: BCF input, `--bed`/`--exclude-bed`, `--thin`, `--mask`,
INFO-flag and mean-depth site filters, genotype-level filters (`--minDP`,
`--minGQ`, ...), `--derived`, and VCFtools' other output types.

## Build and run

```bash
cargo build --release --manifest-path vcftools-rs/Cargo.toml
vcftools-rs/target/release/vcftools-rs --gzvcf in.vcf.gz --freq --out result
```

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

## Verification

- `scripts/compare.sh [--threads N] <names>` runs the benchmark commands on
  1000 Genomes chr22 and checks each output with `cmp` against VCFtools'.
- `scripts/fuzz.py CASES SEED` is a differential fuzz test: it generates
  random messy VCFs (missing/half-missing/haploid genotypes, phased/unphased
  mixes, multi-allelic sites, indels, GT not first in FORMAT, FILTER/QUAL/ID
  variety, repeated chromosomes; plain, BGZF and gzip) with random sample
  and site filters, and requires identical output (or both tools failing).

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

## Credit

All statistical methods and their exact numerical behaviour are those of
VCFtools by Adam Auton, Petr Danecek and contributors (LGPL-3.0). This
reimplementation was written with AI assistance (Claude) and is validated
by output comparison against VCFtools.
