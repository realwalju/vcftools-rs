#!/usr/bin/env bash
# Download 1000 Genomes Phase 3 benchmark data (public, EBI mirror).
set -euo pipefail
cd "$(dirname "$0")/../data"

B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
F=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

for f in "$F" "$F.tbi" integrated_call_samples_v3.20130502.ALL.panel; do
    [ -s "$f" ] || wget -q --show-progress=off "$B/$f"
done
ln -sf "$F" chr22.vcf.gz
ln -sf "$F.tbi" chr22.vcf.gz.tbi

# Population lists for Fst (tab-separated panel: sample pop super_pop gender)
for pop in CEU YRI CHB; do
    awk -v p="$pop" '$2==p {print $1}' integrated_call_samples_v3.20130502.ALL.panel > "$pop.txt"
done
ls -la
