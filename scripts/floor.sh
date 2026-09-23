#!/usr/bin/env bash
# Reference points: raw decompression cost (the floor for any VCF tool) and
# plink2 on the same VCF input, for context against the VCFtools baseline.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
VCF="$ROOT/data/chr22.vcf.gz"
T=$(mktemp -d)

t() { /usr/bin/time -f "%e" -o "$T/t" "$@" > /dev/null 2>&1; printf "%-38s %6ss\n" "$LABEL" "$(cat "$T/t")"; }

LABEL="decompress, 1 thread (bgzip -d)"   t bgzip -dc "$VCF"
LABEL="decompress, 16 threads (bgzip -@16)" t bgzip -dc -@16 "$VCF"
LABEL="plink2 --freq, 1 thread"          t plink2 --vcf "$VCF" --freq --threads 1 --out "$T/p1"
LABEL="plink2 --freq, 16 threads"        t plink2 --vcf "$VCF" --freq --threads 16 --out "$T/p16"
rm -rf "$T"
