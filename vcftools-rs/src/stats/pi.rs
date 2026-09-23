//! --site-pi (output_per_site_nucleotide_diversity) and --window-pi
//! (output_windowed_nucleotide_diversity, without --mask). Integer
//! arithmetic follows the C++ types, including unsigned wrap-around.

use super::{for_each_site, window_range, ChromBins, Runs, Scratch};
use crate::fmt;
use crate::Ctx;

pub const SITE_HEADER: &[u8] = b"CHROM\tPOS\tPI\n";
pub const WINDOW_HEADER: &[u8] = b"CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n";

pub fn site(ctx: &Ctx, text: &[u8], out: &mut Vec<u8>) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, counts| {
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        ctx.allele_counts(gts, alleles.len(), counts);
        let total = counts.iter().fold(0i32, |a, &c| a.wrapping_add(c)) as u32;
        let mut mismatches: i32 = 0;
        for &c in counts.iter() {
            let other = total.wrapping_sub(c as u32) as i32;
            mismatches = mismatches.wrapping_add(c.wrapping_mul(other));
        }
        let pairs = total.wrapping_mul(total.wrapping_sub(1)) as i32;
        let pi = mismatches as f64 / pairs as f64;

        out.extend_from_slice(site.chrom);
        out.push(b'\t');
        fmt::int(out, site.pos);
        out.push(b'\t');
        fmt::g(out, pi);
        out.push(b'\n');
    });
}

pub struct SiteRec {
    pos: i32,
    mismatches: u32,
    comparisons: u64,
    polymorphic: bool,
}

pub fn window_sites(ctx: &Ctx, text: &[u8], recs: &mut Runs<SiteRec>) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, counts| {
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        let n = ctx.allele_counts(gts, alleles.len(), counts);
        let mut mismatches = 0u32;
        for &c in counts.iter() {
            mismatches = mismatches.wrapping_add((c as u32).wrapping_mul(n.wrapping_sub(c as u32)));
        }
        if mismatches == 0 {
            return; // Site is actually fixed.
        }
        recs.push(
            site.chrom,
            SiteRec {
                pos: site.pos,
                mismatches,
                comparisons: n.wrapping_mul(n.wrapping_sub(1)) as u64,
                polymorphic: counts[0] < n as i32,
            },
        );
    });
}

/// [N_variant_sites, N_variant_site_pairs, N_mismatches, N_polymorphic_sites]
type Bin = [u64; 4];

pub struct Windows {
    size: i32,
    step: i32,
    bins: ChromBins<Bin>,
}

impl Windows {
    pub fn new(size: i32, step: i32) -> Self {
        let step = if step <= 0 || step > size { size } else { step };
        Windows { size, step, bins: ChromBins::default() }
    }

    pub fn add(&mut self, recs: Runs<SiteRec>) {
        for (chrom, sites) in recs.runs {
            for r in sites {
                let (first, last) = window_range(r.pos, self.size, self.step);
                if self.bins.visit(&chrom) {
                    self.bins.get(&chrom).resize(1, [0; 4]);
                }
                let bins = self.bins.get(&chrom);
                if last >= bins.len() as i32 {
                    bins.resize((last + 1) as usize, [0; 4]);
                }
                for idx in first..last {
                    let b = &mut bins[idx as usize];
                    b[0] += 1;
                    b[1] = b[1].wrapping_add(r.comparisons);
                    b[2] = b[2].wrapping_add(r.mismatches as u64);
                    if r.polymorphic {
                        b[3] += 1;
                    }
                }
            }
        }
    }

    pub fn write(&mut self, ctx: &Ctx, out: &mut Vec<u8>) {
        let n_kept_chr = (2 * ctx.n_kept) as i32;
        let n_comparisons = n_kept_chr.wrapping_mul(n_kept_chr - 1) as i64 as u64;
        let chrs = self.bins.chrs.clone();
        for chrom in chrs {
            let (size, step) = (self.size, self.step);
            for (s, b) in self.bins.get(&chrom).iter().enumerate() {
                if b[3] == 0 && b[2] == 0 {
                    continue;
                }
                let n_mono = (size as i64 as u64).wrapping_sub(b[0]);
                let n_pairs = b[1].wrapping_add(n_mono.wrapping_mul(n_comparisons));
                let pi = b[2] as f64 / n_pairs as f64;
                let start = (s as u32).wrapping_mul(step as u32);
                out.extend_from_slice(&chrom);
                out.push(b'\t');
                fmt::int(out, start.wrapping_add(1));
                out.push(b'\t');
                fmt::int(out, start.wrapping_add(size as u32));
                out.push(b'\t');
                fmt::int(out, b[3]);
                out.push(b'\t');
                fmt::int(out, n_mono);
                out.push(b'\t');
                fmt::g(out, pi);
                out.push(b'\n');
            }
        }
    }
}
