#!/usr/bin/env bash
# Build vcftools-rs (release) and run its unit tests.
set -o pipefail
source ~/.cargo/env
cd "$(dirname "$0")/../vcftools-rs"
cargo build --release 2>&1 | grep -E "^(warning|error)|-->|Finished" | head -60
cargo test --release 2>&1 | grep -E "test result|FAILED|panicked|error" | head -20
