# Changelog

## 0.1.0 (unreleased)

First release. Output validated byte-for-byte against VCFtools 0.1.17.

- Statistics: `--freq`, `--counts`, `--het`, `--hardy`, `--missing-site`,
  `--missing-indv`, `--site-pi`, `--window-pi`, `--TajimaD`,
  `--weir-fst-pop` (per-site and windowed).
- Sample filters: `--keep`, `--remove`, `--indv`, `--remove-indv`.
- Site filters: region/position/ID lists, allele type and count, QUAL,
  FILTER flags, phasing, mean depth, frequency, allele count, missingness
  and HWE filters.
- Genotype filters: `--minDP`, `--maxDP`, `--minGQ`.
- Input: plain VCF, BGZF and gzip. Parallel decompression and processing
  (`--threads`).
