#!/usr/bin/env bash
# Run the differential fuzz test over several seeds.
# Usage: fuzz_many.sh CASES SEED...
cases=$1; shift
status=0
for s in "$@"; do
    python3 "$(dirname "$0")/fuzz.py" "$cases" "$s" | tail -5 || status=1
done
exit $status
