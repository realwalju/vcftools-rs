#!/usr/bin/env bash
# Context benchmark: plink2 computing overlapping statistics from the same
# VCF (its output format and numerics differ from VCFtools, so this is a
# speed reference only). Median of 3 runs, monotonic clock.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
VCF="$ROOT/data/chr22.vcf.gz"
T=$(mktemp -d)
OUT="$ROOT/results/plink2_context.tsv"
printf "statistic\tthreads\tplink2_s\n" > "$OUT"

median() { sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}'; }

for stat in freq hardy het missing; do
    for th in 1 28; do
        for r in 1 2 3; do
            python3 "$ROOT/scripts/mtime.py" "$T/t" \
                plink2 --vcf "$VCF" --$stat --threads $th --out "$T/p" > /dev/null 2>&1
            cat "$T/t"
        done | median | xargs printf "%s\t%s\t%s\n" "$stat" "$th" | tee -a "$OUT"
    done
done
rm -rf "$T"
