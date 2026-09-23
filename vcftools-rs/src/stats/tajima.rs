//! --TajimaD (variant_file::output_Tajima_D), after Carlson et al. 2005.

use super::{for_each_site, ChromBins, Runs, Scratch};
use crate::fmt;
use crate::Ctx;

pub const HEADER: &[u8] = b"CHROM\tBIN_START\tN_SNPS\tTajimaD\n";

pub struct SiteRec {
    idx: u32,
    p: f64,
}

pub fn sites(ctx: &Ctx, window: i32, text: &[u8], recs: &mut Runs<SiteRec>) {
    let c = 1.0 / window as f64;
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, counts| {
        if alleles.len() != 2 {
            return;
        }
        let idx = (site.pos as f64 * c) as u32;
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        let n = ctx.allele_counts(gts, 2, counts);
        recs.push(site.chrom, SiteRec { idx, p: counts[0] as f64 / n as f64 });
    });
}

pub struct Windows {
    window: i32,
    bins: ChromBins<(i32, f64)>,
}

impl Windows {
    pub fn new(window: i32) -> Self {
        Windows { window, bins: ChromBins::default() }
    }

    pub fn add(&mut self, recs: Runs<SiteRec>) {
        for (chrom, sites) in recs.runs {
            for r in sites {
                let bins = self.bins.get(&chrom);
                if r.idx as usize >= bins.len() {
                    bins.resize(r.idx as usize + 1, (0, 0.0));
                }
                self.bins.visit(&chrom);
                if r.p > 0.0 && r.p < 1.0 {
                    let b = &mut self.bins.get(&chrom)[r.idx as usize];
                    b.0 += 1;
                    b.1 += r.p * (1.0 - r.p);
                }
            }
        }
    }

    pub fn write(&mut self, ctx: &Ctx, out: &mut Vec<u8>) {
        let n = (ctx.n_kept as u32).wrapping_mul(2);
        if n < 2 {
            crate::fatal("Require at least two chromosomes!");
        }
        let (mut a1, mut a2) = (0.0f64, 0.0f64);
        for ui in 1..n {
            a1 += 1.0 / ui as f64;
            a2 += 1.0 / ui.wrapping_mul(ui) as f64;
        }
        let b1 = n.wrapping_add(1) as f64 / 3.0 / (n - 1) as f64;
        let b2 = 2.0 * n.wrapping_mul(n).wrapping_add(n).wrapping_add(3) as f64 / 9.0 / n as f64 / (n - 1) as f64;
        let c1 = b1 - (1.0 / a1);
        let c2 = b2 - ((n + 2) as f64 / (a1 * n as f64)) + (a2 / a1 / a1);
        let e1 = c1 / a1;
        let e2 = c2 / ((a1 * a1) + a2);

        let chrs = self.bins.chrs.clone();
        for chrom in chrs {
            let mut output = false;
            for (s, &(big_s, sum)) in self.bins.get(&chrom).iter().enumerate() {
                let mut d = f64::NAN;
                if big_s > 0 {
                    let pi = 2.0 * sum * n as f64 / (n - 1) as f64;
                    let tw = big_s as f64 / a1;
                    let var = (e1 * big_s as f64) + e2 * big_s as f64 * (big_s - 1) as f64;
                    d = (pi - tw) / var.sqrt();
                    output = true;
                }
                if output {
                    out.extend_from_slice(&chrom);
                    out.push(b'\t');
                    fmt::int(out, (s as u32).wrapping_mul(self.window as u32));
                    out.push(b'\t');
                    fmt::int(out, big_s);
                    out.push(b'\t');
                    fmt::g(out, d);
                    out.push(b'\n');
                }
            }
        }
    }
}
