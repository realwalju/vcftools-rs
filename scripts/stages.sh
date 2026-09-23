#!/usr/bin/env bash
# Time each pipeline stage in isolation (single thread) to find hot spots.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/vcftools-rs/target/release/vcftools-rs"
VCF="$ROOT/data/chr22.vcf.gz"
for stage in decompress lines genotypes; do
    /usr/bin/time -f "%e" -o /tmp/stage.t "$BIN" --gzvcf "$VCF" --threads "${1:-1}" --x-stage-$stage --out /tmp/stage 2>/dev/null
    printf "%-12s %6ss\n" "$stage" "$(cat /tmp/stage.t)"
done
