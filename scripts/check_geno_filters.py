#!/usr/bin/env python3
"""Sanity check that genotype/depth filters are exercised: on one random
VCF, show how many calls each filter removes (VCFtools vs vcftools-rs)."""
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fuzz  # noqa: E402

with tempfile.TemporaryDirectory() as d:
    vcf = os.path.join(d, "t.vcf")
    rng = random.Random(42)
    while True:  # want a file whose rows carry DP and GQ
        fuzz.make_vcf(rng, vcf, 30, 2000)
        if sum(":GQ" in l or "GQ:" in l for l in open(vcf)) > 300:
            break
    for args in (["--minDP", "10"], ["--maxDP", "8"], ["--minGQ", "20"], ["--min-meanDP", "20"]):
        totals = []
        for tool in (fuzz.VCFTOOLS, fuzz.OURS):
            prefix = os.path.join(d, "o")
            subprocess.run([tool, "--vcf", vcf, "--missing-indv", "--out", prefix] + args, capture_output=True)
            rows = [l.split("\t") for l in open(prefix + ".imiss").read().splitlines()[1:]]
            totals.append((sum(int(r[1]) for r in rows), sum(int(r[2]) for r in rows)))
        print(f"{' '.join(args):18} VCFtools N_DATA/FILTERED={totals[0]}  ours={totals[1]}")
