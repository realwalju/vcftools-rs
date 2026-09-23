//! --freq / --counts (variant_file::output_frequency), without --derived.

use super::{for_each_site, Scratch};
use crate::fmt;
use crate::Ctx;

pub fn header(counts: bool) -> &'static [u8] {
    if counts {
        b"CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n"
    } else {
        b"CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n"
    }
}

pub fn process(ctx: &Ctx, counts: bool, text: &[u8], out: &mut Vec<u8>) {
    let mut s = Scratch::default();
    for_each_site(ctx, text, &mut s, |site, alleles, gts, allele_counts| {
        let gts = gts.get(site, ctx.n_indv);
        let n_alleles = alleles.len();
        let n_chr = ctx.allele_counts(gts, n_alleles, allele_counts);

        out.extend_from_slice(site.chrom);
        out.push(b'\t');
        fmt::int(out, site.pos);
        out.push(b'\t');
        fmt::int(out, n_alleles as u32);
        out.push(b'\t');
        fmt::int(out, n_chr);
        for (i, &c) in allele_counts.iter().enumerate() {
            out.push(b'\t');
            out.extend_from_slice(alleles.get(i));
            out.push(b':');
            if counts {
                fmt::int(out, c);
            } else {
                fmt::g(out, c as f64 / n_chr as f64);
            }
        }
        out.push(b'\n');
    });
}
