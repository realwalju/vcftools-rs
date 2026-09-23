//! --weir-fst-pop (output_weir_and_cockerham_fst and the windowed variant
//! with --fst-window-size). Weir & Cockerham (1984), multi-allelic form,
//! transcribed operation-for-operation from VCFtools.

use super::{for_each_site, window_range, ChromBins, Runs, Scratch};
use crate::fmt;
use crate::record::Gt;
use crate::Ctx;

pub const SITE_HEADER: &[u8] = b"CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\n";
pub const WINDOW_HEADER: &[u8] = b"CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST\n";

/// Membership of each sample in each population (kept samples only).
pub struct Pops {
    pub members: Vec<Vec<bool>>,
}

/// Reads population files the way VCFtools does: the first
/// whitespace-separated token of each line names a sample.
pub fn read_pop_samples(path: &str) -> Vec<String> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|_| crate::fatal(&format!("Could not open Individual file: {path}")));
    text.lines().filter_map(|l| l.split_whitespace().next().map(str::to_string)).collect()
}

impl Pops {
    pub fn new(ctx: &Ctx, files: &[String]) -> Self {
        let members = files
            .iter()
            .map(|f| {
                let names: std::collections::HashSet<String> = read_pop_samples(f).into_iter().collect();
                (0..ctx.n_indv).map(|i| ctx.include[i] && names.contains(&ctx.names[i])).collect()
            })
            .collect();
        Pops { members }
    }
}

/// Per-site result: (sum_a, sum_all, fst).
fn site_fst(pops: &Pops, gts: &[Gt], n_alleles: usize) -> (f64, f64, f64) {
    let n_pops = pops.members.len();
    let mut n = vec![0.0f64; n_pops];
    let mut p = vec![vec![0.0f64; n_alleles]; n_pops];
    let mut pbar = vec![0.0f64; n_alleles];
    let mut hbar = vec![0.0f64; n_alleles];
    let mut ssqr = vec![0.0f64; n_alleles];
    let mut sum_nsqr = 0.0f64;
    let mut n_hom = vec![0u32; n_alleles];
    let mut n_het = vec![0u32; n_alleles];

    for i in 0..n_pops {
        // entry::get_multiple_genotype_counts
        n_hom.iter_mut().for_each(|x| *x = 0);
        n_het.iter_mut().for_each(|x| *x = 0);
        for (g, &member) in gts.iter().zip(&pops.members[i]) {
            if !member || g.excluded {
                continue;
            }
            for uj in 0..n_alleles as i32 {
                if g.a == uj && g.b == uj {
                    n_hom[uj as usize] += 1;
                } else if (g.a == uj || g.b == uj) && g.a != -1 && g.b != -1 {
                    n_het[uj as usize] += 1;
                }
            }
        }
        for j in 0..n_alleles {
            n[i] += n_hom[j] as f64 + 0.5 * n_het[j] as f64;
            p[i][j] = n_het[j].wrapping_add(n_hom[j].wrapping_mul(2)) as f64;
            pbar[j] += p[i][j];
            hbar[j] += n_het[j] as f64;
        }
        for j in 0..n_alleles {
            p[i][j] /= 2.0 * n[i];
        }
        sum_nsqr += n[i] * n[i];
    }
    // std::accumulate with an int initial value truncates after every add.
    let n_sum = n.iter().fold(0i32, |acc, &x| (acc as f64 + x) as i32) as f64;
    let nbar = n_sum / n_pops as f64;

    for j in 0..n_alleles {
        pbar[j] /= n_sum * 2.0;
        hbar[j] /= n_sum;
    }
    for j in 0..n_alleles {
        for i in 0..n_pops {
            ssqr[j] += n[i] * (p[i][j] - pbar[j]) * (p[i][j] - pbar[j]);
        }
        ssqr[j] /= (n_pops as f64 - 1.0) * nbar;
    }
    let nc = (n_sum - (sum_nsqr / n_sum)) / (n_pops as f64 - 1.0);
    let r = n_pops as f64;
    let (mut sum_a, mut sum_all) = (0.0f64, 0.0f64);
    for j in 0..n_alleles {
        let a = (ssqr[j] - (pbar[j] * (1.0 - pbar[j]) - (((r - 1.0) * ssqr[j]) / r) - (hbar[j] / 4.0)) / (nbar - 1.0))
            * nbar
            / nc;
        let b = (pbar[j] * (1.0 - pbar[j])
            - (ssqr[j] * (r - 1.0) / r)
            - hbar[j] * (((2.0 * nbar) - 1.0) / (4.0 * nbar)))
            * nbar
            / (nbar - 1.0);
        let c = hbar[j] / 2.0;
        if !a.is_nan() && !b.is_nan() && !c.is_nan() {
            sum_a += a;
            sum_all += a + b + c;
        }
    }
    (sum_a, sum_all, sum_a / sum_all)
}

/// Running genome-wide sums for the log (mean and weighted Fst).
#[derive(Default)]
pub struct Totals {
    sum1: f64,
    sum2: f64,
    sum3: f64,
    count: f64,
}

impl Totals {
    fn add(&mut self, (sum_a, sum_all, fst): (f64, f64, f64)) {
        self.sum1 += sum_a;
        self.sum2 += sum_all;
        self.sum3 += fst;
        self.count += 1.0;
    }

    pub fn log(&self) -> String {
        format!(
            "Weir and Cockerham mean Fst estimate: {:.5}\nWeir and Cockerham weighted Fst estimate: {:.5}\n",
            self.sum3 / self.count,
            self.sum1 / self.sum2
        )
    }
}

/// Per-site pass. Returns output text plus the non-NaN site values (in
/// order) so totals and windows can be summed sequentially.
#[derive(Default)]
pub struct Chunk {
    pub text: Vec<u8>,
    pub sites: Runs<(i32, (f64, f64, f64))>,
}

pub fn process(ctx: &Ctx, pops: &Pops, per_site_output: bool, text: &[u8], c: &mut Chunk) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, _| {
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        let v = site_fst(pops, gts, alleles.len());
        if per_site_output {
            let out = &mut c.text;
            out.extend_from_slice(site.chrom);
            out.push(b'\t');
            fmt::int(out, site.pos);
            out.push(b'\t');
            fmt::g(out, v.2);
            out.push(b'\n');
        }
        if !v.2.is_nan() {
            c.sites.push(site.chrom, (site.pos, v));
        }
    });
}

pub struct Windows {
    size: i32,
    step: i32,
    pub totals: Totals,
    bins: ChromBins<[f64; 4]>,
}

impl Windows {
    pub fn new(size: i32, step: i32) -> Self {
        let step = if step <= 0 || step > size { size } else { step };
        Windows { size, step, totals: Totals::default(), bins: ChromBins::default() }
    }

    /// Adds sites to the totals and (if windowed) to the bins, in order.
    pub fn add(&mut self, sites: Runs<(i32, (f64, f64, f64))>, windowed: bool) {
        for (chrom, recs) in sites.runs {
            for (pos, v) in recs {
                self.totals.add(v);
                if !windowed {
                    continue;
                }
                self.bins.visit(&chrom);
                let (first, last) = window_range(pos, self.size, self.step);
                let bins = self.bins.get(&chrom);
                for idx in first..last {
                    if idx as usize >= bins.len() {
                        bins.resize(idx as usize + 1, [0.0; 4]);
                    }
                    let b = &mut bins[idx as usize];
                    b[0] += v.0;
                    b[1] += v.1;
                    b[2] += v.2;
                    b[3] += 1.0;
                }
            }
        }
    }

    pub fn write(&mut self, out: &mut Vec<u8>) {
        let chrs = self.bins.chrs.clone();
        let (size, step) = (self.size, self.step);
        for chrom in chrs {
            for (s, b) in self.bins.get(&chrom).iter().enumerate() {
                if b[1] != 0.0 && !b[0].is_nan() && !b[1].is_nan() && b[3] > 0.0 {
                    let start = (s as u32).wrapping_mul(step as u32);
                    out.extend_from_slice(&chrom);
                    out.push(b'\t');
                    fmt::int(out, start.wrapping_add(1));
                    out.push(b'\t');
                    fmt::int(out, start.wrapping_add(size as u32));
                    out.push(b'\t');
                    fmt::g(out, b[3]);
                    out.push(b'\t');
                    fmt::g(out, b[0] / b[1]);
                    out.push(b'\t');
                    fmt::g(out, b[2] / b[3]);
                    out.push(b'\n');
                }
            }
        }
    }
}
