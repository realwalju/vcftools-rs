#!/usr/bin/env bash
# Second dataset: 1001 Genomes A. thaliana Chr1 slice (see fetch_athal.sh).
# 1) byte-identity of every statistic, with and without DP/GQ filters;
# 2) benchmark (median of 3): VCFtools vs vcftools-rs at 1 and 28 threads.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
D="$ROOT/data/athal"
VCF="$D/athal_chr1_slice.vcf.gz"
OURS="$ROOT/vcftools-rs/target/release/vcftools-rs"

# Population lists from the ADMIXTURE group column (sample IDs are the
# accession ids in column 1).
python3 - "$D" <<'EOF'
import csv, collections, sys
d = sys.argv[1]
rows = list(csv.DictReader(open(f"{d}/groups.csv")))
counts = collections.Counter(r["group"] for r in rows)
print("groups:", dict(counts))
for g in ("germany", "western_europe"):
    with open(f"{d}/{g}.txt", "w") as f:
        f.write("".join(r["id"] + "\n" for r in rows if r["group"] == g))
EOF
PA="$D/germany.txt"; PB="$D/western_europe.txt"

echo "== correctness"
bash "$ROOT/scripts/check_real.sh" "$VCF" "$PA" "$PB" || exit 1
bash "$ROOT/scripts/check_real.sh" "$VCF" "$PA" "$PB" --minDP 3 --minGQ 20 --max-missing 0.8 || exit 1

echo "== benchmark"
T=$(mktemp -d)
OUT="$ROOT/results/athal_benchmark.tsv"
printf "command\tvcftools_s\tours_1t_s\tours_28t_s\n" > "$OUT"
median() { sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}'; }
time3() { for r in 1 2 3; do python3 "$ROOT/scripts/mtime.py" "$T/t" "$@" > /dev/null 2>&1; cat "$T/t"; done | median; }
cmds=("--freq" "--het" "--hardy" "--missing-site" "--site-pi" "--TajimaD 10000"
      "--weir-fst-pop $PA --weir-fst-pop $PB"
      "--freq --minDP 3 --minGQ 20 --max-missing 0.8")
for c in "${cmds[@]}"; do
    v=$(time3 vcftools --gzvcf "$VCF" $c --out "$T/v")
    o1=$(time3 "$OURS" --gzvcf "$VCF" $c --threads 1 --out "$T/o")
    o28=$(time3 "$OURS" --gzvcf "$VCF" $c --out "$T/o")
    label=${c//$D\//}
    printf "%s\t%s\t%s\t%s\n" "$label" "$v" "$o1" "$o28" | tee -a "$OUT"
done
rm -rf "$T"
