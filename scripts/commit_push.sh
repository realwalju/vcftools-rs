#!/usr/bin/env bash
# Commit all changes with the message on stdin and push to origin.
set -euo pipefail
cd "$(dirname "$0")/.."
git add -A
git status --short
git commit -q -F -
git push -q origin main
git log --format='%h %ae %s' -1
git status -sb | head -1
