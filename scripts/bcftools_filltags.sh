#!/usr/bin/env bash
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/bin/micromamba shell hook -s bash)"
micromamba activate vcfbench
bcftools +fill-tags -h 2>&1 | head -45
echo "== smpl-stats"
bcftools +smpl-stats -h 2>&1 | head -20
