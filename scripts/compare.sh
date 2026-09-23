#!/usr/bin/env bash
# Build vcftools-rs, run it on the benchmark commands given as arguments
# (names from baseline.sh), and check each output is byte-identical to the
# VCFtools reference. Usage: compare.sh [--threads N] freq counts ...
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source ~/.cargo/env
THREADS=()
if [ "${1:-}" = "--threads" ]; then THREADS=(--threads "$2"); shift 2; fi

(cd "$ROOT/vcftools-rs" && cargo build --release -q) || exit 1
BIN="$ROOT/vcftools-rs/target/release/vcftools-rs"
VCF="$ROOT/data/chr22.vcf.gz"
D="$ROOT/data"
REF="$ROOT/results/vcftools"
OUT="$ROOT/results/ours"
mkdir -p "$OUT"

declare -A CMDS=(
    [freq]="--freq"                          [freq_ext]=frq
    [counts]="--counts"                      [counts_ext]=frq.count
    [het]="--het"                            [het_ext]=het
    [hardy]="--hardy"                        [hardy_ext]=hwe
    [missing_site]="--missing-site"          [missing_site_ext]=lmiss
    [missing_indv]="--missing-indv"          [missing_indv_ext]=imiss
    [site_pi]="--site-pi"                    [site_pi_ext]=sites.pi
    [window_pi]="--window-pi 10000"          [window_pi_ext]=windowed.pi
    [tajimad]="--TajimaD 10000"              [tajimad_ext]=Tajima.D
    [fst_site]="--weir-fst-pop $D/CEU.txt --weir-fst-pop $D/YRI.txt"  [fst_site_ext]=weir.fst
    [fst_window]="--weir-fst-pop $D/CEU.txt --weir-fst-pop $D/YRI.txt --fst-window-size 10000 --fst-window-step 5000"  [fst_window_ext]=windowed.weir.fst
)

printf "%-14s %10s %10s %9s  %s\n" command vcftools_s ours_s speedup output
for name in "$@"; do
    ext=${CMDS[${name}_ext]}
    rm -f "$OUT/$name.$ext"   # never compare against a stale file
    python3 "$ROOT/scripts/mtime.py" "$OUT/$name.time" \
        "$BIN" --gzvcf "$VCF" ${CMDS[$name]} "${THREADS[@]}" --out "$OUT/$name" 2> "$OUT/$name.log"
    rc=$?
    ours=$(cat "$OUT/$name.time")
    base=$(awk -v n="$name" '$1==n {print $2}' "$ROOT/results/baseline.tsv")
    if [ $rc -ne 0 ]; then verdict="FAILED (exit $rc)";
    elif cmp -s "$REF/$name.$ext" "$OUT/$name.$ext"; then verdict="IDENTICAL"; else
        verdict="DIFFERS: $(cmp "$REF/$name.$ext" "$OUT/$name.$ext" 2>&1 | head -1)"; fi
    printf "%-14s %10s %10s %8.1fx  %s\n" "$name" "$base" "$ours" "$(echo "$base / $ours" | bc -l)" "$verdict"
done
