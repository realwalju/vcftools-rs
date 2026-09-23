//! --hardy (variant_file::output_hwe) with the exact SNP-HWE test of
//! Wigginton, Cutler & Abecasis (2005), transcribed from entry::SNPHWE.

use super::{for_each_site, Scratch};
use crate::fmt;
use crate::Ctx;

pub const HEADER: &[u8] =
    b"CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS\n";

pub fn process(ctx: &Ctx, text: &[u8], out: &mut Vec<u8>) {
    let mut s = Scratch::default();
    let mut probs = Vec::new();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, counts| {
        if alleles.len() != 2 {
            return;
        }
        let gts = gts.get(site, ctx.n_indv);
        if !ctx.is_diploid(gts) {
            return;
        }
        let n = ctx.allele_counts(gts, 2, counts);
        let freq = counts[0] as f64 / n as f64;

        let (b11, b12, b22) = ctx.genotype_counts(gts);
        let tot = b11.wrapping_add(b12).wrapping_add(b22) as f64;
        let exp_11 = freq * freq * tot;
        let exp_12 = 2.0 * freq * (1.0 - freq) * tot;
        let exp_22 = (1.0 - freq) * (1.0 - freq) * tot;
        let d = |b: u32, e: f64| (b as f64 - e) * (b as f64 - e) / e;
        let chisq = d(b11, exp_11) + d(b12, exp_12) + d(b22, exp_22);
        let (p_hwe, p_lo, p_hi) = snphwe(b12 as i32, b11 as i32, b22 as i32, &mut probs);

        out.extend_from_slice(site.chrom);
        out.push(b'\t');
        fmt::int(out, site.pos);
        out.push(b'\t');
        fmt::int(out, b11);
        out.push(b'/');
        fmt::int(out, b12);
        out.push(b'/');
        fmt::int(out, b22);
        out.push(b'\t');
        fmt::fixed(out, 2, exp_11);
        out.push(b'/');
        fmt::fixed(out, 2, exp_12);
        out.push(b'/');
        fmt::fixed(out, 2, exp_22);
        for v in [chisq, p_hwe, p_lo, p_hi] {
            out.push(b'\t');
            fmt::sci(out, 6, v);
        }
        out.push(b'\n');
    });
}

/// Returns (p_hwe, p_lo, p_hi). `het_probs` is scratch space.
pub fn snphwe(obs_hets: i32, obs_hom1: i32, obs_hom2: i32, het_probs: &mut Vec<f64>) -> (f64, f64, f64) {
    if obs_hom1 + obs_hom2 + obs_hets == 0 {
        return (1.0, 1.0, 1.0);
    }
    let obs_homc = if obs_hom1 < obs_hom2 { obs_hom2 } else { obs_hom1 };
    let obs_homr = if obs_hom1 < obs_hom2 { obs_hom1 } else { obs_hom2 };
    let rare_copies = 2 * obs_homr + obs_hets;
    let genotypes = obs_hets + obs_homc + obs_homr;

    het_probs.clear();
    het_probs.resize(rare_copies as usize + 1, 0.0);

    let mut mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
    if (rare_copies & 1) ^ (mid & 1) != 0 {
        mid += 1;
    }

    let mut curr_homr = (rare_copies - mid) / 2;
    let mut curr_homc = genotypes - mid - curr_homr;
    het_probs[mid as usize] = 1.0;
    let mut sum = het_probs[mid as usize];
    let mut curr_hets = mid;
    while curr_hets > 1 {
        let h = curr_hets as usize;
        het_probs[h - 2] = het_probs[h] * curr_hets as f64 * (curr_hets as f64 - 1.0)
            / (4.0 * (curr_homr as f64 + 1.0) * (curr_homc as f64 + 1.0));
        sum += het_probs[h - 2];
        curr_homr += 1;
        curr_homc += 1;
        curr_hets -= 2;
    }

    curr_homr = (rare_copies - mid) / 2;
    curr_homc = genotypes - mid - curr_homr;
    curr_hets = mid;
    while curr_hets <= rare_copies - 2 {
        let h = curr_hets as usize;
        het_probs[h + 2] = het_probs[h] * 4.0 * curr_homr as f64 * curr_homc as f64
            / ((curr_hets as f64 + 2.0) * (curr_hets as f64 + 1.0));
        sum += het_probs[h + 2];
        curr_homr -= 1;
        curr_homc -= 1;
        curr_hets += 2;
    }

    for p in het_probs.iter_mut() {
        *p /= sum;
    }

    let obs = obs_hets as usize;
    let mut p_hi = het_probs[obs];
    for p in &het_probs[obs + 1..] {
        p_hi += p;
    }
    let mut p_lo = het_probs[obs];
    for p in het_probs[..obs].iter().rev() {
        p_lo += p;
    }
    let mut p_hwe = 0.0;
    for &p in het_probs.iter() {
        if p > het_probs[obs] {
            continue;
        }
        p_hwe += p;
    }
    if p_hwe > 1.0 {
        p_hwe = 1.0;
    }
    (p_hwe, p_lo, p_hi)
}
