#!/usr/bin/env bash
# Probe two suspected VCFtools 0.1.17 issues found while transcribing it.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VT=~/micromamba/envs/vcfbench/bin/vcftools
T=$(mktemp -d)
FIX="$ROOT/tests/data/chr22_1500.vcf.gz"

echo "== 1. --non-ref-af-any appears to have no effect on its own"
for args in "" "--non-ref-af-any 0.5" "--non-ref-af-any 0.9" "--non-ref-ac-any 500"; do
    $VT --gzvcf "$FIX" --freq $args --out "$T/o" > /dev/null 2>&1
    printf "  %-24s sites kept: %s\n" "${args:-(no filter)}" "$(($(wc -l < "$T/o.frq") - 1))"
done

echo "== 2. --site-pi with many samples (int overflow?)"
# One biallelic site, n diploid samples, half 0|0 and half 1|1:
# 2n chromosomes, n copies of each allele, so pi = n*n*2 / (2n*(2n-1)).
python3 - "$T" <<'EOF'
import sys
d = sys.argv[1]
for n in (1000, 20000, 24000, 40000):
    with open(f"{d}/n{n}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                "\t".join(f"S{i}" for i in range(n)) + "\n")
        f.write("1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t" +
                "\t".join("0|0" if i < n // 2 else "1|1" for i in range(n)) + "\n")
EOF
for n in 1000 20000 24000 40000; do
    $VT --vcf "$T/n$n.vcf" --site-pi --out "$T/p$n" > /dev/null 2>&1
    got=$(tail -1 "$T/p$n.sites.pi" | cut -f3)
    expect=$(python3 -c "n=$n; c=2*n; print('%g' % ((n*n*2)/(c*(c-1))))")
    printf "  %6d samples: VCFtools PI=%-12s correct=%s\n" "$n" "$got" "$expect"
done
rm -rf "$T"
