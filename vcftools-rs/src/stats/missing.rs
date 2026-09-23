//! --missing-site (output_site_missingness) and --missing-indv
//! (output_indv_missingness). Calls removed by genotype filters are
//! reported as N_GENOTYPE(S)_FILTERED and excluded from N_DATA.

use super::{for_each_site, Scratch};
use crate::fmt;
use crate::Ctx;

pub const SITE_HEADER: &[u8] = b"CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS\n";
pub const INDV_HEADER: &[u8] = b"INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS\n";

pub fn site(ctx: &Ctx, text: &[u8], out: &mut Vec<u8>) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, _, gts, _| {
        let gts = gts.get(site, ctx.n_indv);
        let (mut miss, mut tot, mut filtered) = (0u32, 0u32, 0u32);
        for (g, &inc) in gts.iter().zip(&ctx.include) {
            if !inc {
                continue;
            }
            if g.excluded {
                filtered += 1;
                continue;
            }
            if g.a == -1 {
                miss += 1;
            }
            if g.b == -1 {
                miss += 1;
            }
            tot += 2;
            // Phased missing second allele = haploid genome.
            if (g.b == -1 && g.phase == b'|') || g.b == -2 {
                tot = tot.wrapping_sub(1);
                miss = miss.wrapping_sub(1);
            }
        }
        out.extend_from_slice(site.chrom);
        out.push(b'\t');
        fmt::int(out, site.pos);
        out.push(b'\t');
        fmt::int(out, tot);
        out.push(b'\t');
        fmt::int(out, filtered);
        out.push(b'\t');
        fmt::int(out, miss);
        out.push(b'\t');
        fmt::g(out, miss as f64 / tot as f64);
        out.push(b'\n');
    });
}

/// Per-individual totals for one chunk.
#[derive(Default)]
pub struct IndvCounts {
    pub tot: Vec<u32>,
    pub miss: Vec<u32>,
    pub filtered: Vec<u32>,
}

impl IndvCounts {
    pub fn add(&mut self, other: IndvCounts) {
        if self.tot.is_empty() {
            *self = other;
            return;
        }
        let pairs = [(&mut self.tot, other.tot), (&mut self.miss, other.miss), (&mut self.filtered, other.filtered)];
        for (acc, add) in pairs {
            for (a, b) in acc.iter_mut().zip(add) {
                *a += b;
            }
        }
    }
}

pub fn indv(ctx: &Ctx, text: &[u8], acc: &mut IndvCounts) {
    if acc.tot.is_empty() {
        acc.tot = vec![0; ctx.n_indv];
        acc.miss = vec![0; ctx.n_indv];
        acc.filtered = vec![0; ctx.n_indv];
    }
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, _, gts, _| {
        let gts = gts.get(site, ctx.n_indv);
        for (i, g) in gts.iter().enumerate() {
            if g.excluded {
                acc.filtered[i] += 1;
                continue;
            }
            if g.a == -1 {
                acc.miss[i] += 1;
            }
            acc.tot[i] += 1;
        }
    });
}

pub fn write_indv(ctx: &Ctx, acc: &IndvCounts, out: &mut Vec<u8>) {
    let get = |v: &Vec<u32>, i: usize| v.get(i).copied().unwrap_or(0);
    for i in 0..ctx.n_indv {
        if !ctx.include[i] {
            continue;
        }
        let (tot, miss) = (get(&acc.tot, i), get(&acc.miss, i));
        out.extend_from_slice(ctx.names[i].as_bytes());
        out.push(b'\t');
        fmt::int(out, tot);
        out.push(b'\t');
        fmt::int(out, get(&acc.filtered, i));
        out.push(b'\t');
        fmt::int(out, miss);
        out.push(b'\t');
        fmt::g(out, miss as f64 / tot as f64);
        out.push(b'\n');
    }
}
