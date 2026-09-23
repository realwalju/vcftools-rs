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
i7-14700K, WSL2 Ubuntu 24.04. See `results/baseline.tsv` (VCFtools, median
of 3 runs) and `results/summary_1thread.txt`, `results/summary_28threads.txt`.

Reproduce: `scripts/fetch_data.sh`, `scripts/baseline.sh`, `scripts/final.sh`.

## Credit

All statistical methods and their exact numerical behaviour are those of
VCFtools by Adam Auton, Petr Danecek and contributors (LGPL-3.0). This
reimplementation was written with AI assistance (Claude) and is validated
by output comparison against VCFtools.
