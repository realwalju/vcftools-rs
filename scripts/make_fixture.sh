#!/usr/bin/env bash
# Build the small real-data test fixture from 1000 Genomes chr22:
# header + first 1,500 variants (all 2,504 samples), and CEU/YRI lists.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
mkdir -p "$ROOT/tests/data"
cd "$ROOT/tests/data"
# awk stops early, so the decompressor gets SIGPIPE; that is expected.
{ bgzip -dc "$ROOT/data/chr22.vcf.gz" | awk '/^#/ {print; next} {print; if (++n == 1500) exit}' || true; } \
    | bgzip -c > chr22_1500.vcf.gz
cp "$ROOT/data/CEU.txt" "$ROOT/data/YRI.txt" .
ls -la
