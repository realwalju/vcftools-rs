//! --het (variant_file::output_het): method-of-moments inbreeding F.
//!
//! E(HOM) is a floating-point sum over sites. To reproduce VCFtools'
//! rounding exactly, each individual's sum is accumulated in site order;
//! chunks only compute the per-site terms, and terms are applied in file
//! order (parallelised across individuals, which keeps each sum's order).

use rayon::prelude::*;

use super::{for_each_site, Scratch};
use crate::fmt;
use crate::Ctx;

pub const HEADER: &[u8] = b"INDV\tO(HOM)\tE(HOM)\tN_SITES\tF\n";

/// Expected-homozygosity term for one site, and the kept individuals that
/// are missing at that site (and so do not receive it).
pub struct Term {
    value: f64,
    missing: Vec<u32>,
}

#[derive(Default)]
pub struct Chunk {
    obs_hom: Vec<u32>,
    n_sites: Vec<u32>,
    terms: Vec<Term>,
}

pub fn process(ctx: &Ctx, text: &[u8], c: &mut Chunk) {
    if c.obs_hom.is_empty() {
        c.obs_hom = vec![0; ctx.n_indv];
        c.n_sites = vec![0; ctx.n_indv];
    }
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, counts| {
        if alleles.len() != 2 {
            return;
        }
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        let n = ctx.allele_counts(gts, 2, counts);
        let freq = if n > 0 { counts[1] as f64 / n as f64 } else { -1.0 };
        if freq <= f64::EPSILON || 1.0 - freq <= f64::EPSILON {
            return;
        }
        let n = n as f64;
        let value = 1.0 - (2.0 * freq * (1.0 - freq) * (n / (n - 1.0)));
        let mut missing = Vec::new();
        for (i, (g, &inc)) in gts.iter().zip(&ctx.include).enumerate() {
            if !inc {
                continue;
            }
            if !g.excluded && g.a > -1 && g.b > -1 {
                c.n_sites[i] += 1;
                if g.a == g.b {
                    c.obs_hom[i] += 1;
                }
            } else {
                missing.push(i as u32);
            }
        }
        c.terms.push(Term { value, missing });
    });
}

#[derive(Default)]
pub struct Acc {
    obs_hom: Vec<u32>,
    n_sites: Vec<u32>,
    expected: Vec<f64>,
    /// Terms not yet applied, in file order. Applying in large batches
    /// keeps the parallel step coarse-grained.
    pending: Vec<Term>,
}

const BLOCK: usize = 64;
const BATCH: usize = 1 << 16;

impl Acc {
    pub fn add(&mut self, mut c: Chunk) {
        if c.obs_hom.is_empty() {
            return;
        }
        if self.obs_hom.is_empty() {
            self.obs_hom = vec![0; c.obs_hom.len()];
            self.n_sites = vec![0; c.obs_hom.len()];
            self.expected = vec![0.0; c.obs_hom.len()];
        }
        for (a, b) in self.obs_hom.iter_mut().zip(&c.obs_hom) {
            *a += b;
        }
        for (a, b) in self.n_sites.iter_mut().zip(&c.n_sites) {
            *a += b;
        }
        self.pending.append(&mut c.terms);
        if self.pending.len() >= BATCH {
            self.apply();
        }
    }

    /// Adds pending terms to each individual's E(HOM) in file order,
    /// parallelised across individuals.
    fn apply(&mut self) {
        let terms = std::mem::take(&mut self.pending);
        let terms = &terms;
        self.expected.par_chunks_mut(BLOCK).enumerate().for_each(|(blk, e)| {
            let base = (blk * BLOCK) as u32;
            let end = base + e.len() as u32;
            for t in terms {
                if t.missing.is_empty() {
                    e.iter_mut().for_each(|x| *x += t.value);
                } else {
                    for (k, x) in e.iter_mut().enumerate() {
                        let i = base + k as u32;
                        if i < end && t.missing.binary_search(&i).is_err() {
                            *x += t.value;
                        }
                    }
                }
            }
        });
    }

    pub fn write(&mut self, ctx: &Ctx, out: &mut Vec<u8>) {
        self.apply();
        for i in 0..ctx.n_indv {
            if !ctx.include[i] || self.n_sites.get(i).copied().unwrap_or(0) == 0 {
                continue;
            }
            let (o, e, n) = (self.obs_hom[i], self.expected[i], self.n_sites[i]);
            let f = (o as f64 - e) / (n as f64 - e);
            out.extend_from_slice(ctx.names[i].as_bytes());
            out.push(b'\t');
            fmt::int(out, o);
            out.push(b'\t');
            fmt::fixed(out, 1, e);
            out.push(b'\t');
            fmt::int(out, n);
            out.push(b'\t');
            fmt::fixed(out, 5, f);
            out.push(b'\n');
        }
    }
}
