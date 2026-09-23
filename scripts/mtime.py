#!/usr/bin/env python3
"""Run a command and print its elapsed time from a monotonic clock.
(/usr/bin/time uses the wall clock, which WSL2 can step during a run.)
Usage: mtime.py OUTFILE cmd args...  -- stdout/stderr of cmd pass through."""
import subprocess
import sys
import time

start = time.monotonic()
rc = subprocess.run(sys.argv[2:]).returncode
with open(sys.argv[1], "w") as f:
    f.write(f"{time.monotonic() - start:.2f}\n")
sys.exit(rc)
