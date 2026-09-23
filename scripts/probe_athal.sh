#!/usr/bin/env bash
# Probe the 1001 Genomes VCF: download speed and FORMAT fields (first 2 MB).
set -uo pipefail
URL=https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
T=$(mktemp -d)
start=$(date +%s.%N)
curl -fsS -r 0-2097151 "$URL" -o "$T/head.gz"
end=$(date +%s.%N)
echo "2 MB in $(echo "$end - $start" | bc) s"
# Truncated stream: decompress what we can.
gzip -dc "$T/head.gz" 2>/dev/null > "$T/head.vcf"
grep -m3 '^##FORMAT' "$T/head.vcf"
grep -v '^#' "$T/head.vcf" | head -3 | cut -f1-10 | cut -c1-200
echo "samples: $(grep -m1 '^#CHROM' "$T/head.vcf" | awk -F'\t' '{print NF-9}')"
echo "complete data rows in 2 MB: $(grep -vc '^#' "$T/head.vcf")"
~/micromamba/envs/vcfbench/bin/htsfile "$T/head.gz" 2>&1 | head -1
rm -rf "$T"
