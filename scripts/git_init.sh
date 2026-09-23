#!/usr/bin/env bash
# Create the local repository and make the first commit.
set -euo pipefail
cd "$(dirname "$0")/.."
rm -f scripts/refactor_gts.sh
git init -q -b main
git add -A
git status --short
git commit -q -F - <<'EOF'
vcftools-rs: byte-identical, parallel reimplementation of VCFtools stats

Implements --freq, --counts, --het, --hardy, --missing-site,
--missing-indv, --site-pi, --window-pi, --TajimaD and --weir-fst-pop
(per-site and windowed), sample filters and the common site filters,
with output byte-identical to VCFtools 0.1.17. Includes benchmark,
comparison and differential fuzz-testing scripts.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline
