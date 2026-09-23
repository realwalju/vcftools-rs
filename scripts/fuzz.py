#!/usr/bin/env python3
"""Differential fuzz test: generate random messy VCFs and require
vcftools-rs output to be byte-identical to VCFtools for every statistic.

Covers: missing and half-missing genotypes, phased/unphased mixes, haploid
calls, multi-allelic sites, "." ALT, lower-case alleles, GT not first in
FORMAT, samples lacking GT, sites-with-no-data, repeated chromosomes,
--keep/--remove, and window step options.

Usage: fuzz.py [N_CASES] [SEED]
"""
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OURS = os.path.join(ROOT, "vcftools-rs/target/release/vcftools-rs")
ENV_BIN = os.path.expanduser("~/micromamba/envs/vcfbench/bin")
VCFTOOLS = os.path.join(ENV_BIN, "vcftools")
BGZIP = os.path.join(ENV_BIN, "bgzip")

BASES = "ACGT"


def genotype(rng, n_alleles, style):
    """One GT string, drawn from a mix of well-formed shapes."""
    r = rng.random()
    allele = lambda: str(rng.randrange(n_alleles))
    if style == "haploid_site" and r < 0.5:
        return rng.choice([allele(), "."])
    if r < 0.08:
        return rng.choice(["./.", ".|."])
    if r < 0.12:
        return rng.choice(["0/.", "./1", ".|0", "1|."])
    if r < 0.14 and style != "clean":
        return rng.choice([allele(), "."])  # stray haploid call
    sep = "|" if (style == "phased" or rng.random() < 0.5) else "/"
    return allele() + sep + allele()


def make_vcf(rng, path, n_samples, n_sites):
    samples = [f"S{i:03d}" for i in range(n_samples)]
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    chroms = rng.choice([["1"], ["1", "2"], ["chrA", "chrB", "chrA"]])
    # Rows without GT make VCFtools reject genotype-based filters, so only
    # some files contain them.
    fmt_weights = [60, 20, 15, 5] if rng.random() < 0.2 else [60, 20, 15, 0]
    sites = []  # (chrom, pos, id) for building filter lists
    for chrom in chroms:
        pos = rng.randrange(1, 500)
        for _ in range(n_sites // len(chroms)):
            pos += rng.choice([1, 2, 5, 17, 60, 250, 900, 4000])
            ref = rng.choice(BASES)
            n_alt = rng.choices([0, 1, 2, 3], weights=[3, 80, 12, 5])[0]
            alts = rng.sample([b for b in BASES if b != ref], n_alt)
            if rng.random() < 0.1:  # indel
                if rng.random() < 0.5:
                    ref += "".join(rng.choices(BASES, k=rng.randrange(1, 4)))
                elif alts:
                    alts[0] += "".join(rng.choices(BASES, k=rng.randrange(1, 4)))
            alt = ",".join(alts) if alts else "."
            vid = f"rs{rng.randrange(10**6)}" if rng.random() < 0.8 else "."
            qual = rng.choice([".", "0", "12.5", "29.99", "30", "55", "1e2"])
            flt = rng.choice(["PASS", "PASS", ".", "LowQual", "q10", "q10;s50", "LowQual;PASS", ""])
            sites.append((chrom, pos, vid))
            if rng.random() < 0.05:
                ref, alt = ref.lower(), alt.lower()
            n_alleles = 1 + len(alts)
            style = rng.choices(["mixed", "phased", "clean", "haploid_site"], weights=[50, 25, 15, 10])[0]
            fmt = rng.choices(["GT", "GT:DP", "DP:GT", "DP"], weights=fmt_weights)[0]
            cols = []
            for _ in samples:
                gt = genotype(rng, n_alleles, style)
                dp = str(rng.randrange(0, 60))
                if fmt == "GT":
                    cols.append(gt)
                elif fmt == "GT:DP":
                    cols.append(gt if rng.random() < 0.05 else f"{gt}:{dp}")
                elif fmt == "DP:GT":
                    cols.append(dp if rng.random() < 0.05 else f"{dp}:{gt}")
                else:
                    cols.append(dp)
            lines.append(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{flt}\t.\t{fmt}\t" + "\t".join(cols))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return samples, sites


def site_filters(rng, d, sites):
    """0-3 random site filters (VCFtools' commonly used ones)."""
    chroms = sorted({s[0] for s in sites})
    pos = sorted(s[1] for s in sites)

    def listfile(name, rows):
        p = os.path.join(d, name)
        with open(p, "w") as f:
            f.write("\n".join(rows) + "\n")
        return p

    pool = [
        lambda: ["--chr", rng.choice(chroms)],
        lambda: ["--chr", rng.choice(chroms), "--from-bp", str(rng.choice(pos)), "--to-bp", str(rng.choice(pos) + 5000)],
        lambda: ["--not-chr", rng.choice(chroms)],
        lambda: ["--maf", rng.choice(["0.01", "0.05", "0.2"])],
        lambda: ["--max-maf", rng.choice(["0.1", "0.3"])],
        lambda: ["--max-missing", rng.choice(["0.5", "0.9", "1"])],
        lambda: ["--mac", rng.choice(["1", "2", "5"])],
        lambda: ["--max-mac", rng.choice(["3", "20"])],
        lambda: ["--max-missing-count", rng.choice(["0", "3", "10"])],
        lambda: ["--min-alleles", "2", "--max-alleles", "2"],
        lambda: ["--max-alleles", "3"],
        lambda: [rng.choice(["--remove-indels", "--keep-only-indels"])],
        lambda: ["--minQ", rng.choice(["0", "20", "30"])],
        lambda: ["--remove-filtered-all"],
        lambda: ["--remove-filtered", rng.choice(["LowQual", "s50"])],
        lambda: ["--keep-filtered", rng.choice(["q10", "PASS"])],
        lambda: ["--hwe", rng.choice(["0.001", "0.05"])],
        lambda: ["--phased"],
        lambda: ["--non-ref-af", "0.1"],
        lambda: ["--max-non-ref-af", "0.5"],
        lambda: ["--non-ref-af-any", "0.1"],
        lambda: ["--non-ref-ac", "2"],
        lambda: ["--max-non-ref-ac-any", "10"],
        lambda: ["--snp", rng.choice(sites)[2]],
        lambda: ["--snps", listfile("snps.txt", [s[2] for s in rng.sample(sites, len(sites) // 2)])],
        lambda: ["--exclude", listfile("excl.txt", [s[2] for s in rng.sample(sites, len(sites) // 4)])],
        lambda: ["--positions", listfile("pos.txt", ["#CHROM\tPOS"] + [f"{c}\t{p}" for c, p, _ in rng.sample(sites, len(sites) // 2)])],
        lambda: ["--exclude-positions", listfile("xpos.txt", [f"{c}\t{p}" for c, p, _ in rng.sample(sites, len(sites) // 4)])],
    ]
    args = []
    for f in rng.sample(pool, rng.choice([0, 0, 1, 1, 2, 3])):
        args += f()
    return args


def commands(rng, d, samples, sites):
    def pop(name, k):
        p = os.path.join(d, name)
        with open(p, "w") as f:
            f.write("\n".join(rng.sample(samples, k)) + "\n")
        return p

    k = max(2, len(samples) // 3)
    popA, popB = pop("popA.txt", k), pop("popB.txt", k)
    keep, remove = pop("keep.txt", len(samples) // 2), pop("remove.txt", len(samples) // 5)
    win = str(rng.choice([100, 1000, 5000]))
    step = str(rng.choice([50, 100, 500, 1000]))
    filt = rng.choice([[], ["--keep", keep], ["--remove", remove]]) + site_filters(rng, d, sites)
    return [
        (["--freq"] + filt, "frq"),
        (["--counts"] + filt, "frq.count"),
        (["--het"] + filt, "het"),
        (["--hardy"] + filt, "hwe"),
        (["--missing-site"] + filt, "lmiss"),
        (["--missing-indv"] + filt, "imiss"),
        (["--site-pi"] + filt, "sites.pi"),
        (["--window-pi", win] + filt, "windowed.pi"),
        (["--window-pi", win, "--window-pi-step", step] + filt, "windowed.pi"),
        (["--TajimaD", win] + filt, "Tajima.D"),
        (["--weir-fst-pop", popA, "--weir-fst-pop", popB] + site_filters(rng, d, sites), "weir.fst"),
        (["--weir-fst-pop", popA, "--weir-fst-pop", popB, "--fst-window-size", win, "--fst-window-step", step]
         + site_filters(rng, d, sites), "windowed.weir.fst"),
    ]


def run(tool, vcf, args, prefix):
    flag = "--gzvcf" if vcf.endswith(".gz") else "--vcf"
    cmd = [tool, flag, vcf, "--out", prefix] + args
    if tool == OURS:
        cmd += ["--threads", "3"]
    return subprocess.run(cmd, capture_output=True, text=True).returncode


def main():
    n_cases = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    failures = 0
    checks = 0
    both_failed = 0
    rows_compared = 0
    for case in range(n_cases):
        rng = random.Random(seed * 100003 + case)
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, "t.vcf")
            samples, sites = make_vcf(rng, vcf, rng.randrange(4, 40), rng.randrange(20, 3000))
            if case % 3 == 1:
                subprocess.run([BGZIP, "-f", vcf], check=True)
                vcf += ".gz"
            elif case % 3 == 2:
                subprocess.run(["gzip", "-f", vcf], check=True)
                vcf += ".gz"
            for args, ext in commands(rng, d, samples, sites):
                rc_ref = run(VCFTOOLS, vcf, args, os.path.join(d, "ref"))
                rc_ours = run(OURS, vcf, args, os.path.join(d, "ours"))
                ref_f, ours_f = os.path.join(d, f"ref.{ext}"), os.path.join(d, f"ours.{ext}")
                checks += 1
                ref_data = open(ref_f, "rb").read() if os.path.exists(ref_f) else None
                ours_data = open(ours_f, "rb").read() if os.path.exists(ours_f) else None
                if rc_ref != 0:
                    both_failed += rc_ours != 0
                elif ref_data is not None:
                    rows_compared += ref_data.count(b"\n") - 1
                # When both tools reject the input, partial output of the
                # aborted run is not meaningful; agreement on failure suffices.
                both_error = rc_ref != 0 and rc_ours != 0
                if (rc_ref == 0) != (rc_ours == 0) or (not both_error and ref_data != ours_data):
                    failures += 1
                    keep_dir = os.path.join(ROOT, "results", "fuzz_failures", f"case{case}_{ext}")
                    os.makedirs(keep_dir, exist_ok=True)
                    subprocess.run(["cp", "-r", d + "/.", keep_dir])
                    print(f"FAIL case {case} {' '.join(args)} (rc {rc_ref} vs {rc_ours}) -> {keep_dir}")
    print(f"{checks - failures}/{checks} checks identical over {n_cases} random VCFs (seed {seed}); "
          f"{rows_compared} output rows compared; {both_failed} checks where both tools exited with an error")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
