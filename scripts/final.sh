#!/usr/bin/env bash
# Full benchmark: all statistics at 1 thread and all threads; saves tables.
cd "$(dirname "$0")/.."
ALL="freq counts het hardy missing_site missing_indv site_pi window_pi tajimad fst_site fst_window"
bash scripts/compare.sh --threads 1 $ALL > results/summary_1thread.txt
bash scripts/compare.sh $ALL > results/summary_28threads.txt
echo "== 1 thread";   cat results/summary_1thread.txt
echo "== 28 threads"; cat results/summary_28threads.txt
