#!/usr/bin/env bash
# Sanity-check VCFtools reference outputs: log messages, row counts, samples.
cd "$(dirname "$0")/../results/vcftools"
for f in *.log; do
    echo "== $f"
    grep -iE "error|warning|After filtering|Run Time" "$f"
done
echo
wc -l *.frq *.frq.count *.het *.hwe *.lmiss *.imiss *.sites.pi *.windowed.pi *.Tajima.D *.weir.fst *.windowed.weir.fst
echo
head -3 freq.frq het.het hardy.hwe site_pi.sites.pi tajimad.Tajima.D fst_site.weir.fst fst_window.windowed.weir.fst
echo
# Measure peak RSS again with a direct, verbose time call on one command
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
file "$(command -v vcftools)"
/usr/bin/time -v vcftools --gzvcf ../../data/chr22.vcf.gz --chr 22 --to-bp 2000000 --freq --out /tmp/rsscheck 2>&1 \
    | grep -E "Maximum resident|Elapsed"
