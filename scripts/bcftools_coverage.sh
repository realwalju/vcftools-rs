#!/usr/bin/env bash
# What can bcftools compute natively? Plugin list + fill-tags tag list.
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
bcftools --version | head -1
echo "== plugins"
bcftools plugin -l 2>&1 | tr '\n' ' '; echo
echo "== fill-tags: available tags"
bcftools +fill-tags -h 2>&1 | sed -n '/Available tags/,/^$/p'
echo "== search plugin/command help for popgen terms"
for p in $(bcftools plugin -l 2>/dev/null); do
    bcftools +$p -h 2>&1 | grep -qiE 'fst|tajima|nucleotide diversity|\bpi\b|inbreeding|heterozygosity' && echo "plugin mentions popgen term: $p"
done
bcftools stats -h 2>&1 | grep -iE 'fst|tajima|diversity|inbreed' || echo "bcftools stats: no Fst/Tajima/diversity/inbreeding options"
