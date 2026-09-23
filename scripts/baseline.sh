#!/usr/bin/env bash
# Time VCFtools on each in-scope command (monotonic clock, RUNS repeats,
# median reported). Writes reference outputs to results/vcftools/<name>.*
# and a timing table to results/baseline.tsv.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUNS=${RUNS:-3}
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench

VCF="$ROOT/data/chr22.vcf.gz"
D="$ROOT/data"
OUT="$ROOT/results/vcftools"
mkdir -p "$OUT"
TSV="$ROOT/results/baseline.tsv"
printf "command\twall_s\truns_s\n" > "$TSV.tmp"

declare -A CMDS=(
    [freq]="--freq"
    [counts]="--counts"
    [het]="--het"
    [hardy]="--hardy"
    [missing_site]="--missing-site"
    [missing_indv]="--missing-indv"
    [site_pi]="--site-pi"
    [window_pi]="--window-pi 10000"
    [tajimad]="--TajimaD 10000"
    [fst_site]="--weir-fst-pop $D/CEU.txt --weir-fst-pop $D/YRI.txt"
    [fst_window]="--weir-fst-pop $D/CEU.txt --weir-fst-pop $D/YRI.txt --fst-window-size 10000 --fst-window-step 5000"
)

for name in freq counts het hardy missing_site missing_indv site_pi window_pi tajimad fst_site fst_window; do
    times=()
    for r in $(seq "$RUNS"); do
        python3 "$ROOT/scripts/mtime.py" "$OUT/$name.time" \
            vcftools --gzvcf "$VCF" ${CMDS[$name]} --out "$OUT/$name" > /dev/null 2>&1
        times+=("$(cat "$OUT/$name.time")")
    done
    median=$(printf "%s\n" "${times[@]}" | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    printf "%s\t%s\t%s\n" "$name" "$median" "$(IFS=,; echo "${times[*]}")" | tee -a "$TSV.tmp"
done
mv "$TSV.tmp" "$TSV"
column -t "$TSV"
