#!/usr/bin/env bash
# Second benchmark dataset: 1001 Genomes Arabidopsis thaliana (1,135
# accessions, FORMAT GT:GQ:DP). The full VCF is 18 GB and has no tabix
# index, so we take the first 400 MB of the BGZF file (start of Chr1),
# keep complete rows, and re-compress. Population lists come from the
# project's ADMIXTURE group table.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
D="$ROOT/data/athal"
mkdir -p "$D"
cd "$D"
URL=https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
[ -s slice.raw.gz ] || curl -fsS -r 0-419430399 "$URL" -o slice.raw.gz

# The last BGZF block is cut off: decompress what is readable, then drop
# the final (partial) line.
{ gzip -dc slice.raw.gz 2>/dev/null || true; } | sed '$d' | bgzip -@8 -c > athal_chr1_slice.vcf.gz
echo "rows: $(bgzip -dc athal_chr1_slice.vcf.gz | grep -vc '^#')"
echo "last: $(bgzip -dc athal_chr1_slice.vcf.gz | tail -1 | cut -f1-2)"

# ADMIXTURE groups (public Google Sheet linked from 1001genomes/admixture-map).
SHEET=https://docs.google.com/spreadsheets/d/1_jp6KKfUC0z1WteS9HUIp4f3ssFZlCV3KXPIvrTW0Fk/export?format=csv
[ -s groups.csv ] || curl -fsSL "$SHEET" -o groups.csv
head -3 groups.csv
