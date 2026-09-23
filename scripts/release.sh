#!/usr/bin/env bash
# Tag and push a release. Usage: release.sh VERSION (e.g. 0.1.0)
# The Release workflow then builds binaries into a draft GitHub release.
set -euo pipefail
cd "$(dirname "$0")/.."
v=$1
grep -q "^version = \"$v\"" vcftools-rs/Cargo.toml || { echo "Cargo.toml version is not $v"; exit 1; }
grep -q "^## $v " CHANGELOG.md || { echo "CHANGELOG has no $v entry"; exit 1; }
if [ -n "$(git status --porcelain)" ]; then
    git add -A
    printf 'Release %s\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n' "$v" | git commit -q -F -
fi
git tag -a "v$v" -m "vcftools-rs $v"
git push -q origin main "v$v"
git log --oneline -1
git tag -l "v$v"
