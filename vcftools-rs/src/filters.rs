//! Site filters (entry::apply_filters in VCFtools' entry_filters.cpp),
//! applied in the same order and with the same comparison semantics.
//! Genotype filters (--minDP, --maxDP, --minGQ) mark individual calls as
//! excluded rather than removing the site.
//! Not supported: --bed/--exclude-bed, --thin, --mask, INFO-flag filters,
//! genotype FILTER-flag filters (--remove-filtered-geno*).

use std::collections::{HashMap, HashSet};
use std::io::Read;

use crate::record::{str2double, str2int, Alleles, Site};
use crate::stats::Genotypes;
use crate::{fatal, Ctx};

const INT_MAX: f64 = i32::MAX as f64;

/// Raw filter options as given on the command line (VCFtools defaults).
pub struct FilterArgs {
    pub keep_only_indels: bool,
    pub remove_indels: bool,
    pub snps: Vec<String>,
    pub snps_file: Option<String>,
    pub exclude_file: Option<String>,
    pub filter_keep: Vec<String>,
    pub filter_remove: Vec<String>,
    pub remove_all_filtered: bool,
    pub chrs: Vec<String>,
    pub not_chrs: Vec<String>,
    pub from_bp: i32,
    pub to_bp: i32,
    pub positions: Option<String>,
    pub exclude_positions: Option<String>,
    pub min_alleles: i32,
    pub max_alleles: i32,
    pub min_quality: f64,
    pub min_mean_depth: f64,
    pub max_mean_depth: f64,
    pub phased: bool,
    pub min_gq: f64,
    pub min_dp: i32,
    pub max_dp: i32,
    pub min_maf: f64,
    pub max_maf: f64,
    pub min_nraf: f64,
    pub max_nraf: f64,
    pub min_nraf_any: f64,
    pub max_nraf_any: f64,
    pub min_call_rate: f64,
    pub min_mac: f64,
    pub max_mac: f64,
    pub min_nrac: f64,
    pub max_nrac: f64,
    pub min_nrac_any: f64,
    pub max_nrac_any: f64,
    pub max_missing_count: f64,
    pub min_hwe: f64,
}

impl Default for FilterArgs {
    fn default() -> Self {
        FilterArgs {
            keep_only_indels: false,
            remove_indels: false,
            snps: Vec::new(),
            snps_file: None,
            exclude_file: None,
            filter_keep: Vec::new(),
            filter_remove: Vec::new(),
            remove_all_filtered: false,
            chrs: Vec::new(),
            not_chrs: Vec::new(),
            from_bp: -1,
            to_bp: i32::MAX,
            positions: None,
            exclude_positions: None,
            min_alleles: -1,
            max_alleles: i32::MAX,
            min_quality: -1.0,
            min_mean_depth: -1.0,
            max_mean_depth: f64::MAX,
            phased: false,
            min_gq: -1.0,
            min_dp: -1,
            max_dp: i32::MAX,
            min_maf: -1.0,
            max_maf: f64::MAX,
            min_nraf: -1.0,
            max_nraf: f64::MAX,
            min_nraf_any: -1.0,
            max_nraf_any: f64::MAX,
            min_call_rate: 0.0,
            min_mac: -1.0,
            max_mac: INT_MAX,
            min_nrac: -1.0,
            max_nrac: INT_MAX,
            min_nrac_any: -1.0,
            max_nrac_any: INT_MAX,
            max_missing_count: INT_MAX,
            min_hwe: -1.0,
        }
    }
}

impl FilterArgs {
    /// Consumes `flag` (and its value via `value`) if it is a filter
    /// option. Numeric values are parsed with atoi/atof like VCFtools.
    pub fn parse(&mut self, flag: &str, value: &mut dyn FnMut() -> String) -> bool {
        let i = |s: String| crate::record::atoi(s.as_bytes());
        let f = |s: String| crate::record::atof(s.as_bytes());
        match flag {
            "--keep-only-indels" => self.keep_only_indels = true,
            "--remove-indels" => self.remove_indels = true,
            "--snp" => self.snps.push(value()),
            "--snps" => self.snps_file = Some(value()),
            "--exclude" => self.exclude_file = Some(value()),
            "--keep-filtered" => self.filter_keep.push(value()),
            "--remove-filtered" => self.filter_remove.push(value()),
            "--remove-filtered-all" => self.remove_all_filtered = true,
            // std::set upstream: repeated names count once.
            "--chr" => push_unique(&mut self.chrs, value()),
            "--not-chr" => push_unique(&mut self.not_chrs, value()),
            "--from-bp" => self.from_bp = i(value()),
            "--to-bp" => self.to_bp = i(value()),
            "--positions" => self.positions = Some(value()),
            "--exclude-positions" => self.exclude_positions = Some(value()),
            "--min-alleles" => self.min_alleles = i(value()),
            "--max-alleles" => self.max_alleles = i(value()),
            "--minQ" => self.min_quality = f(value()),
            "--min-meanDP" => self.min_mean_depth = f(value()),
            "--max-meanDP" => self.max_mean_depth = f(value()),
            "--phased" => self.phased = true,
            "--minGQ" => self.min_gq = f(value()),
            "--minDP" => self.min_dp = i(value()),
            "--maxDP" => self.max_dp = i(value()),
            "--maf" => self.min_maf = f(value()),
            "--max-maf" => self.max_maf = f(value()),
            "--non-ref-af" => self.min_nraf = f(value()),
            "--max-non-ref-af" => self.max_nraf = f(value()),
            "--non-ref-af-any" => self.min_nraf_any = f(value()),
            "--max-non-ref-af-any" => self.max_nraf_any = f(value()),
            "--max-missing" => self.min_call_rate = f(value()),
            "--mac" => self.min_mac = i(value()) as f64,
            "--max-mac" => self.max_mac = i(value()) as f64,
            "--non-ref-ac" => self.min_nrac = i(value()) as f64,
            "--max-non-ref-ac" => self.max_nrac = i(value()) as f64,
            "--non-ref-ac-any" => self.min_nrac_any = i(value()) as f64,
            "--max-non-ref-ac-any" => self.max_nrac_any = i(value()) as f64,
            "--max-missing-count" => self.max_missing_count = i(value()) as f64,
            "--hwe" => {
                self.max_alleles = 2;
                self.min_hwe = f(value());
            }
            _ => return false,
        }
        true
    }
}

impl FilterArgs {
    /// parameters::check_parameters (filter-related checks, same exit codes).
    pub fn validate(&self) {
        let err = |msg: &str, code: i32| -> ! {
            eprintln!("\n\nError: {msg}\n");
            std::process::exit(code);
        };
        let d = FilterArgs::default();
        if !self.chrs.is_empty() && !self.not_chrs.is_empty() {
            err("Cannot specify chromosomes to keep and to exclude", 1);
        }
        if self.to_bp < self.from_bp {
            err("End position must be greater than Start position.", 1);
        }
        if (self.to_bp != i32::MAX || self.from_bp != -1) && self.chrs.len() != 1 {
            err("Require a single chromosome when specifying a range.", 2);
        }
        if self.max_maf < self.min_maf {
            err("Maximum MAF must be not be less than Minimum MAF.", 4);
        }
        if self.max_mac < self.min_mac {
            err("Maximum MAC must be not be less than Minimum MAC.", 4);
        }
        if self.min_maf != d.min_maf && !(0.0..=1.0).contains(&self.min_maf) {
            err("MAF must be between 0 and 1.", 4);
        }
        if self.max_maf != d.max_maf && !(0.0..=1.0).contains(&self.max_maf) {
            err("Maximum MAF must be between 0 and 1.", 4);
        }
        if (self.min_nraf != d.min_nraf && !(0.0..=1.0).contains(&self.min_nraf))
            || (self.min_nraf_any != d.min_nraf_any && !(0.0..=1.0).contains(&self.min_nraf_any))
        {
            err("Non-Ref Allele Frequency must be between 0 and 1.", 4);
        }
        if self.max_nraf < self.min_nraf || self.max_nraf_any < self.min_nraf_any {
            err("Maximum Non-Ref Allele Frequency must not be less that Minimum Non-Ref AF.", 4);
        }
        if self.max_nrac < self.min_nrac || self.max_nrac_any < self.min_nrac_any {
            err("Maximum Non-Ref Allele Count must not be less that Minimum Non-Ref AC.", 4);
        }
        if self.min_call_rate > 1.0 {
            err("Minimum Call rate cannot be greater than 1.", 5);
        }
        if self.max_alleles < self.min_alleles {
            err("Max Number of Alleles must be greater than Min Number of Alleles.", 6);
        }
        if self.max_mean_depth < self.min_mean_depth {
            err("Max Mean Depth must be greater the Min Mean Depth.", 7);
        }
        if self.max_dp < self.min_dp {
            err("Max Genotype Depth must be greater than Min Genotype Depth.", 9);
        }
    }
}

type PositionSet = HashMap<Vec<u8>, HashSet<i32>>;

pub struct SiteFilter {
    a: FilterArgs,
    snps_keep: Option<HashSet<Vec<u8>>>,
    snps_exclude: Option<HashSet<Vec<u8>>>,
    filter_keep: HashSet<Vec<u8>>,
    filter_remove: HashSet<Vec<u8>>,
    region: Option<(Vec<u8>, i32, i32)>,
    positions_keep: Option<PositionSet>,
    positions_exclude: Option<PositionSet>,
    chr_keep: HashSet<Vec<u8>>,
    chr_exclude: HashSet<Vec<u8>>,
    freq_active: bool,
    mac_active: bool,
}

impl SiteFilter {
    pub fn new(a: FilterArgs) -> Self {
        if a.keep_only_indels && a.remove_indels {
            fatal("Can't both keep and remove all indels!");
        }
        let bytes_set = |v: &[String]| v.iter().map(|s| s.as_bytes().to_vec()).collect::<HashSet<_>>();
        let snps_keep = if a.snps.is_empty() && a.snps_file.is_none() {
            None
        } else {
            let mut s = bytes_set(&a.snps);
            if let Some(f) = &a.snps_file {
                s.extend(first_tokens(&read_text(f, "SNPs to Keep")));
            }
            Some(s)
        };
        let snps_exclude = a.exclude_file.as_ref().map(|f| first_tokens(&read_text(f, "SNPs to Exclude")).collect());
        // --from-bp/--to-bp only apply together with exactly one --chr.
        let region = (a.chrs.len() == 1 && !(a.from_bp == -1 && a.to_bp == i32::MAX))
            .then(|| (a.chrs[0].as_bytes().to_vec(), a.from_bp, a.to_bp));
        let freq_active = !(a.min_maf <= 0.0
            && a.max_maf >= 1.0
            && a.min_call_rate <= 0.0
            && a.min_nraf <= 0.0
            && a.max_nraf >= 1.0
            && a.min_nraf_any <= 0.0
            && a.max_nraf_any >= 1.0);
        let mac_active = !(a.min_mac <= 0.0
            && a.max_mac == INT_MAX
            && a.min_nrac <= 0.0
            && a.max_nrac == INT_MAX
            && a.min_nrac_any <= 0.0
            && a.max_nrac_any == INT_MAX
            && a.max_missing_count == INT_MAX);
        SiteFilter {
            snps_keep,
            snps_exclude,
            filter_keep: bytes_set(&a.filter_keep),
            filter_remove: bytes_set(&a.filter_remove),
            region,
            positions_keep: a.positions.as_ref().map(|f| read_positions(f)),
            positions_exclude: a.exclude_positions.as_ref().map(|f| read_positions(f)),
            chr_keep: bytes_set(&a.chrs),
            chr_exclude: bytes_set(&a.not_chrs),
            freq_active,
            mac_active,
            a,
        }
    }

    /// entry::apply_filters. `g` caches the decoded genotypes for reuse by
    /// the statistic; `counts` is scratch space.
    pub fn passes(&self, ctx: &Ctx, site: &Site, alleles: &Alleles, g: &mut Genotypes, counts: &mut Vec<i32>) -> bool {
        let a = &self.a;
        let na = alleles.len();

        if a.keep_only_indels || a.remove_indels {
            let ref_len = alleles.get(0).len();
            let is_indel = ref_len != 1 || (1..na).any(|i| alleles.get(i).len() != ref_len);
            if (a.keep_only_indels && !is_indel) || (a.remove_indels && is_indel) {
                return false;
            }
        }
        if self.snps_exclude.as_ref().is_some_and(|s| s.contains(site.id)) {
            return false;
        }
        if self.snps_keep.as_ref().is_some_and(|s| !s.contains(site.id)) {
            return false;
        }
        if a.remove_all_filtered || !self.filter_remove.is_empty() || !self.filter_keep.is_empty() {
            let mut flags: Vec<&[u8]> =
                if site.filter == b"." { Vec::new() } else { site.filter.split(|&c| c == b';').collect() };
            flags.sort();
            if !self.filter_keep.is_empty() && !flags.iter().any(|f| self.filter_keep.contains(*f)) {
                return false;
            }
            if flags.first() != Some(&&b"PASS"[..]) {
                if a.remove_all_filtered && !flags.is_empty() {
                    return false;
                }
                if !a.remove_all_filtered && flags.iter().any(|f| self.filter_remove.contains(*f)) {
                    return false;
                }
            }
        }
        if let Some((chr, start, end)) = &self.region {
            if site.chrom != chr.as_slice() || site.pos < *start || site.pos > *end {
                return false;
            }
        }
        if let Some(keep) = &self.positions_keep {
            if !keep.get(site.chrom).is_some_and(|s| s.contains(&site.pos)) {
                return false;
            }
        }
        if let Some(ex) = &self.positions_exclude {
            if ex.get(site.chrom).is_some_and(|s| s.contains(&site.pos)) {
                return false;
            }
        }
        if !self.chr_keep.is_empty() {
            if !self.chr_keep.contains(site.chrom) {
                return false;
            }
        } else if self.chr_exclude.contains(site.chrom) {
            return false;
        }
        if (na as i32) < a.min_alleles || (na as i32) > a.max_alleles {
            return false;
        }
        if a.min_quality >= 0.0 && str2double(site.qual) < a.min_quality {
            return false;
        }
        if a.min_mean_depth > 0.0 || a.max_mean_depth != f64::MAX {
            // Kept samples without a depth still count in the denominator.
            let (mut sum, mut n) = (0.0f64, 0u32);
            site.for_each_subfield(ctx.n_indv, site.dp_idx, |i, dp| {
                if ctx.include[i] {
                    let depth = dp.map_or(-1, str2int);
                    if depth >= 0 {
                        sum += depth as f64;
                    }
                    n += 1;
                }
            });
            let mean = sum / n as f64;
            if mean < a.min_mean_depth || mean > a.max_mean_depth {
                return false;
            }
        }
        if a.phased {
            let gts = g.get(site, ctx.n_indv);
            if gts.iter().zip(&ctx.include).any(|(g, &inc)| inc && g.phase != b'|') {
                return false;
            }
        }
        // Genotype filters (GQ, then DP upstream; both only mark exclusions,
        // so one pass over the sample columns reads both fields).
        let gq_idx = if a.min_gq > 0.0 { site.gq_idx } else { -1 };
        let dp_idx = if a.min_dp > 0 || a.max_dp != i32::MAX { site.dp_idx } else { -1 };
        if gq_idx != -1 || dp_idx != -1 {
            g.get(site, ctx.n_indv);
            let gts = g.get_mut();
            site.for_each_subfield2(ctx.n_indv, gq_idx, dp_idx, |i, gq, dp| {
                if gq_idx != -1 {
                    // set_indv_GQUALITY: missing is -1, values above 99 are capped.
                    let mut q = gq.map_or(-1.0, str2double);
                    if q != -1.0 && q > 99.0 {
                        q = 99.0;
                    }
                    if q < a.min_gq {
                        gts[i].excluded = true;
                    }
                }
                if dp_idx != -1 {
                    let depth = dp.map_or(-1, str2int);
                    if depth < a.min_dp || depth > a.max_dp {
                        gts[i].excluded = true;
                    }
                }
            });
        }
        if self.freq_active {
            if site.gt_idx == -1 {
                fatal("Require Genotypes in variant file to filter by frequency and/or call rate");
            }
            let gts = g.get(site, ctx.n_indv);
            let n = ctx.allele_counts(gts, na, counts);
            let mut pass = true;
            let mut maf = f64::MAX;
            let mut n_failed = 0u32;
            for (ui, &c) in counts.iter().enumerate() {
                let freq = c as f64 / n as f64;
                maf = cmin(maf, cmin(freq, 1.0 - freq));
                if ui > 0 && (freq < a.min_nraf || freq > a.max_nraf) {
                    pass = false;
                }
                if ui > 0 && (freq < a.min_nraf_any || freq > a.max_nraf_any) {
                    n_failed += 1;
                }
            }
            // Upstream tests min/max_non_ref_af (not the _any variants) here.
            if (a.min_nraf > 0.0 || a.max_nraf < 1.0) && n_failed == (na as u32).wrapping_sub(1) {
                pass = false;
            }
            if maf < a.min_maf || maf > a.max_maf {
                pass = false;
            }
            if (n as f64 / ctx.n_chr(gts) as f64) < a.min_call_rate {
                pass = false;
            }
            if !pass {
                return false;
            }
        }
        if self.mac_active {
            let gts = g.get(site, ctx.n_indv);
            let n = ctx.allele_counts(gts, na, counts);
            let n_chr = ctx.n_chr(gts);
            let mut pass = !(na <= 1 && a.min_mac > 0.0);
            let mut mac = i32::MAX;
            let mut n_failed = 0u32;
            for (ui, &c) in counts.iter().enumerate() {
                mac = mac.min(c);
                if ui > 0 && ((c as f64) < a.min_nrac || (c as f64) > a.max_nrac) {
                    pass = false;
                }
                if ui > 0 && ((c as f64) < a.min_nrac_any || (c as f64) > a.max_nrac_any) {
                    n_failed += 1;
                }
            }
            if (a.min_nrac_any > 0.0 || a.max_nrac_any < INT_MAX) && n_failed == (na as u32).wrapping_sub(1) {
                pass = false;
            }
            if (mac as f64) < a.min_mac || (mac as f64) > a.max_mac {
                pass = false;
            }
            if n_chr.wrapping_sub(n) as f64 > a.max_missing_count {
                pass = false;
            }
            if !pass {
                return false;
            }
        }
        if a.min_hwe > 0.0 {
            if na > 2 {
                fatal("Tried to return the genotype counts of a non-biallelic SNP");
            }
            let gts = g.get(site, ctx.n_indv);
            let (b11, b12, b22) = ctx.genotype_counts(gts);
            let (p_hwe, _, _) = crate::stats::hwe::snphwe(b12 as i32, b11 as i32, b22 as i32, &mut Vec::new());
            if p_hwe < a.min_hwe {
                return false;
            }
        }
        true
    }
}

fn push_unique(v: &mut Vec<String>, s: String) {
    if !v.contains(&s) {
        v.push(s);
    }
}

/// std::min(a, b): returns a unless b < a (so NaN in `b` is ignored).
fn cmin(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

fn read_text(path: &str, what: &str) -> String {
    let mut f = std::fs::File::open(path).unwrap_or_else(|_| fatal(&format!("Could not open {what} file {path}")));
    let mut raw = Vec::new();
    f.read_to_end(&mut raw).unwrap_or_else(|e| fatal(&e.to_string()));
    if raw.starts_with(&[31, 139]) {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .unwrap_or_else(|e| fatal(&e.to_string()));
        raw = out;
    }
    String::from_utf8_lossy(&raw).into_owned()
}

/// First whitespace-separated token of each non-blank line (`in >> tmp;
/// in.ignore(max, '\n')`).
fn first_tokens(text: &str) -> impl Iterator<Item = Vec<u8>> + '_ {
    text.lines().filter_map(|l| l.split_whitespace().next().map(|t| t.as_bytes().to_vec()))
}

/// --positions / --exclude-positions: "CHROM POS" per line, '#' lines skipped.
fn read_positions(path: &str) -> PositionSet {
    let mut out: PositionSet = HashMap::new();
    for line in read_text(path, "Positions").lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut it = line.split_whitespace();
        if let (Some(chr), Some(pos)) = (it.next(), it.next()) {
            out.entry(chr.as_bytes().to_vec()).or_default().insert(crate::record::atoi(pos.as_bytes()));
        }
    }
    out
}
