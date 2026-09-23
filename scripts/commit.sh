#!/usr/bin/env bash
# Commit all tracked changes with the message given on stdin.
set -euo pipefail
cd "$(dirname "$0")/.."
git rm -q --cached --ignore-unmatch results/baseline.tsv.tmp
git add -A
git status --short
git commit -q -F -
git log --oneline | head -5
