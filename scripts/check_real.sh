#!/usr/bin/env bash
# Run every supported statistic with VCFtools and vcftools-rs on a real VCF
# and require byte-identical output. Population files for Fst are optional.
# Usage: check_real.sh VCF [POP_A POP_B] [EXTRA_ARGS...]
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VCF=$1; shift
POPS=()
if [ $# -ge 2 ] && [ -f "$1" ] && [ -f "$2" ]; then POPS=("$1" "$2"); shift 2; fi
EXTRA=("$@")
VCFTOOLS=${VCFTOOLS:-$(command -v vcftools || echo ~/micromamba/envs/vcfbench/bin/vcftools)}
OURS="$ROOT/vcftools-rs/target/release/vcftools-rs"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
flag=--vcf; [[ $VCF == *.gz ]] && flag=--gzvcf

cmds=(
    "--freq:frq" "--counts:frq.count" "--het:het" "--hardy:hwe"
    "--missing-site:lmiss" "--missing-indv:imiss" "--site-pi:sites.pi"
    "--window-pi 10000:windowed.pi" "--window-pi 10000 --window-pi-step 2500:windowed.pi"
    "--TajimaD 10000:Tajima.D"
)
if [ ${#POPS[@]} -eq 2 ]; then
    cmds+=("--weir-fst-pop ${POPS[0]} --weir-fst-pop ${POPS[1]}:weir.fst"
           "--weir-fst-pop ${POPS[0]} --weir-fst-pop ${POPS[1]} --fst-window-size 10000 --fst-window-step 5000:windowed.weir.fst")
fi

fail=0
for c in "${cmds[@]}"; do
    args=${c%:*}; ext=${c##*:}
    "$VCFTOOLS" $flag "$VCF" $args "${EXTRA[@]}" --out "$T/ref" > /dev/null 2>&1; rc1=$?
    "$OURS" $flag "$VCF" $args "${EXTRA[@]}" --out "$T/ours" 2> /dev/null; rc2=$?
    if [ $rc1 -ne 0 ] && [ $rc2 -ne 0 ]; then verdict="both rejected input"
    elif [ $rc1 -ne $rc2 ]; then verdict="EXIT CODE DIFFERS ($rc1 vs $rc2)"; fail=1
    elif cmp -s "$T/ref.$ext" "$T/ours.$ext"; then verdict="identical ($(($(wc -l < "$T/ref.$ext") - 1)) rows)"
    else verdict="DIFFERS"; fail=1; fi
    printf "%-75s %s\n" "$args ${EXTRA[*]}" "$verdict"
    rm -f "$T"/ref.* "$T"/ours.*
done
exit $fail
