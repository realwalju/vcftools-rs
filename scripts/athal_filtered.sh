#!/usr/bin/env bash
# Re-check and re-time the genotype-filtered Arabidopsis case.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
D="$ROOT/data/athal"
VCF="$D/athal_chr1_slice.vcf.gz"
OURS="$ROOT/vcftools-rs/target/release/vcftools-rs"
F="--minDP 3 --minGQ 20 --max-missing 0.8"
bash "$ROOT/scripts/check_real.sh" "$VCF" "$D/germany.txt" "$D/western_europe.txt" $F || exit 1
T=$(mktemp -d)
median() { sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}'; }
for th in 1 28; do
    t=$(for r in 1 2 3; do python3 "$ROOT/scripts/mtime.py" "$T/t" "$OURS" --gzvcf "$VCF" --freq $F --threads $th --out "$T/o" > /dev/null 2>&1; cat "$T/t"; done | median)
    echo "ours --freq $F, $th thread(s): ${t}s  (VCFtools: 34.53s)"
done
rm -rf "$T"
