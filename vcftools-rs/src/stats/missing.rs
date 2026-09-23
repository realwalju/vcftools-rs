//! --missing-site (output_site_missingness) and --missing-indv
//! (output_indv_missingness). Genotype filters are not supported, so
//! N_GENOTYPE_FILTERED is always 0.

use super::{for_each_site, Scratch};
use crate::fmt;
use crate::Ctx;

pub const SITE_HEADER: &[u8] = b"CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS\n";
pub const INDV_HEADER: &[u8] = b"INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS\n";

pub fn site(ctx: &Ctx, text: &[u8], out: &mut Vec<u8>) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, _, gts, _| {
        let gts = gts.get(site, ctx.n_indv);
        let (mut miss, mut tot) = (0u32, 0u32);
        for (g, &inc) in gts.iter().zip(&ctx.include) {
            if !inc {
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
        out.extend_from_slice(b"\t0\t");
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
}

impl IndvCounts {
    pub fn add(&mut self, other: IndvCounts) {
        if self.tot.is_empty() {
            *self = other;
            return;
        }
        for (a, b) in self.tot.iter_mut().zip(other.tot) {
            *a += b;
        }
        for (a, b) in self.miss.iter_mut().zip(other.miss) {
            *a += b;
        }
    }
}

pub fn indv(ctx: &Ctx, text: &[u8], acc: &mut IndvCounts) {
    if acc.tot.is_empty() {
        acc.tot = vec![0; ctx.n_indv];
        acc.miss = vec![0; ctx.n_indv];
    }
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, _, gts, _| {
        let gts = gts.get(site, ctx.n_indv);
        for (i, g) in gts.iter().enumerate() {
            if g.a == -1 {
                acc.miss[i] += 1;
            }
            acc.tot[i] += 1;
        }
    });
}

pub fn write_indv(ctx: &Ctx, acc: &IndvCounts, out: &mut Vec<u8>) {
    for i in 0..ctx.n_indv {
        if !ctx.include[i] {
            continue;
        }
        let (tot, miss) = (acc.tot.get(i).copied().unwrap_or(0), acc.miss.get(i).copied().unwrap_or(0));
        out.extend_from_slice(ctx.names[i].as_bytes());
        out.push(b'\t');
        fmt::int(out, tot);
        out.extend_from_slice(b"\t0\t");
        fmt::int(out, miss);
        out.push(b'\t');
        fmt::g(out, miss as f64 / tot as f64);
        out.push(b'\n');
    }
}
