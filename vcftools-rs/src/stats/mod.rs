//! Statistics. Each module mirrors one `variant_file::output_*` function
//! in VCFtools' variant_file_output.cpp. Per-site work runs in parallel;
//! anything that sums floating-point values across sites is accumulated
//! sequentially in file order so rounding matches VCFtools exactly.

pub mod freq;
pub mod fst;
pub mod het;
pub mod hwe;
pub mod missing;
pub mod pi;
pub mod tajima;

use crate::pipeline::for_each_line;
use crate::record::{Alleles, Gt, Site};
use crate::Ctx;

/// Genotypes of the current site, decoded at most once (by a filter or by
/// the statistic, whichever asks first).
#[derive(Default)]
pub struct Genotypes {
    buf: Vec<Gt>,
    ready: bool,
}

impl Genotypes {
    pub fn get(&mut self, site: &Site, n_indv: usize) -> &[Gt] {
        if !self.ready {
            site.genotypes(n_indv, &mut self.buf);
            self.ready = true;
        }
        &self.buf
    }
}

/// Per-chunk reusable buffers.
#[derive(Default)]
pub struct Scratch {
    pub alleles: Alleles,
    pub gts: Genotypes,
    pub counts: Vec<i32>,
}

/// Runs `f` for each data line in `text` that passes the site filters.
/// Genotypes are decoded lazily via `Genotypes::get`.
pub fn for_each_site(
    ctx: &Ctx,
    text: &[u8],
    s: &mut Scratch,
    mut f: impl FnMut(&Site, &Alleles, &mut Genotypes, &mut Vec<i32>),
) {
    for_each_line(text, |line| {
        let site = Site::parse(line, &mut s.alleles);
        s.gts.ready = false;
        if ctx.filter.passes(ctx, &site, &s.alleles, &mut s.gts, &mut s.counts) {
            f(&site, &s.alleles, &mut s.gts, &mut s.counts);
        }
    });
}

impl Ctx {
    /// entry::get_N_chr: total ploidy of kept individuals.
    pub fn n_chr(&self, gts: &[Gt]) -> u32 {
        gts.iter().zip(&self.include).filter(|(_, &inc)| inc).map(|(g, _)| g.ploidy as u32).sum()
    }

    /// entry::get_genotype_counts for a biallelic site: (hom1, het, hom2).
    pub fn genotype_counts(&self, gts: &[Gt]) -> (u32, u32, u32) {
        let (mut b11, mut b12, mut b22) = (0u32, 0u32, 0u32);
        for (g, &inc) in gts.iter().zip(&self.include) {
            if inc && g.a > -1 && g.b > -1 {
                if g.a != g.b {
                    b12 += 1;
                } else if g.a == 0 {
                    b11 += 1;
                } else if g.a == 1 {
                    b22 += 1;
                } else {
                    crate::fatal("Unknown allele in genotype");
                }
            }
        }
        (b11, b12, b22)
    }

    /// entry::is_diploid: every kept individual has a diploid call.
    pub fn is_diploid(&self, gts: &[Gt]) -> bool {
        gts.iter().zip(&self.include).all(|(g, &inc)| !inc || g.ploidy == 2)
    }

    /// entry::get_allele_counts over kept individuals; returns the number
    /// of non-missing chromosomes. Out-of-range allele indices (undefined
    /// behaviour upstream) are counted as present but not tallied.
    pub fn allele_counts(&self, gts: &[Gt], n_alleles: usize, counts: &mut Vec<i32>) -> u32 {
        // Branch-free histogram. Slot 0: missing (< 0); slots 1..=n: allele
        // a-1; slot n+1: out-of-range index (non-missing but not tallied).
        // Two tallies break the store-to-load dependency between a and b.
        let slot = |a: i32| ((a + 1).max(0) as usize).min(n_alleles + 1);
        let mut ta = [0u32; 8];
        let mut tb = [0u32; 8];
        let (ta, tb): (&mut [u32], &mut [u32]) = if n_alleles + 2 <= 8 {
            (&mut ta[..n_alleles + 2], &mut tb[..n_alleles + 2])
        } else {
            counts.clear();
            counts.resize(2 * (n_alleles + 2), 0);
            return self.allele_counts_slow(gts, n_alleles, counts);
        };
        if self.n_kept == self.n_indv {
            for g in gts {
                ta[slot(g.a)] += 1;
                tb[slot(g.b)] += 1;
            }
        } else {
            for (g, &inc) in gts.iter().zip(&self.include) {
                if inc {
                    ta[slot(g.a)] += 1;
                    tb[slot(g.b)] += 1;
                }
            }
        }
        counts.clear();
        counts.extend((1..=n_alleles).map(|i| (ta[i] + tb[i]) as i32));
        ta[1..].iter().chain(&tb[1..]).sum()
    }

    fn allele_counts_slow(&self, gts: &[Gt], n_alleles: usize, counts: &mut Vec<i32>) -> u32 {
        counts.clear();
        counts.resize(n_alleles, 0);
        let mut n = 0u32;
        for (g, &inc) in gts.iter().zip(&self.include) {
            if !inc {
                continue;
            }
            for a in [g.a, g.b] {
                if a > -1 {
                    if let Some(c) = counts.get_mut(a as usize) {
                        *c += 1;
                    }
                    n += 1;
                }
            }
        }
        n
    }
}

/// Per-site records grouped into runs of the same chromosome, so chunk
/// results carry each chromosome name once.
pub struct Runs<T> {
    pub runs: Vec<(Vec<u8>, Vec<T>)>,
}

impl<T> Default for Runs<T> {
    fn default() -> Self {
        Runs { runs: Vec::new() }
    }
}

impl<T> Runs<T> {
    pub fn push(&mut self, chrom: &[u8], t: T) {
        match self.runs.last_mut() {
            Some((c, v)) if c.as_slice() == chrom => v.push(t),
            _ => self.runs.push((chrom.to_vec(), vec![t])),
        }
    }
}

/// The window indices [first, last) that a site at `pos` falls into, as
/// computed by VCFtools' windowed pi and Fst.
pub fn window_range(pos: i32, size: i32, step: i32) -> (i32, i32) {
    let first = ((pos.wrapping_sub(size)) as f64 / step as f64).ceil() as i32;
    let last = (pos as f64 / step as f64).ceil() as i32;
    (first.max(0), last)
}

/// std::map<string, vector<Bin>> plus the order chromosomes were first
/// seen (VCFtools prints in that order, repeating a chromosome if the
/// input revisits it).
pub struct ChromBins<B> {
    pub bins: std::collections::HashMap<Vec<u8>, Vec<B>>,
    pub chrs: Vec<Vec<u8>>,
}

impl<B> Default for ChromBins<B> {
    fn default() -> Self {
        ChromBins { bins: Default::default(), chrs: Vec::new() }
    }
}

impl<B: Clone + Default> ChromBins<B> {
    /// Records `chrom` in the output order if it differs from the last one.
    /// Returns true when a new run started.
    pub fn visit(&mut self, chrom: &[u8]) -> bool {
        if self.chrs.last().map(Vec::as_slice) != Some(chrom) {
            self.chrs.push(chrom.to_vec());
            true
        } else {
            false
        }
    }

    pub fn get(&mut self, chrom: &[u8]) -> &mut Vec<B> {
        if !self.bins.contains_key(chrom) {
            self.bins.insert(chrom.to_vec(), Vec::new());
        }
        self.bins.get_mut(chrom).unwrap()
    }
}
