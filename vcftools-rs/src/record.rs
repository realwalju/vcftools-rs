//! Parsing of one VCF data line, following VCFtools' `vcf_entry` rules
//! exactly (including its quirks) so downstream statistics match.

use memchr::memchr;

/// A genotype as VCFtools represents it: two allele indices (-1 = missing)
/// and a ploidy of 1 or 2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gt {
    pub a: i32,
    pub b: i32,
    pub ploidy: u8,
    /// b'|' or b'/'; haploid calls count as phased.
    pub phase: u8,
    /// Removed by a genotype filter (VCFtools' include_genotype == false).
    pub excluded: bool,
}

impl Gt {
    const fn new(a: i32, b: i32, ploidy: u8, phase: u8) -> Self {
        Gt { a, b, ploidy, phase, excluded: false }
    }
}

const MISSING_DIPLOID: Gt = Gt::new(-1, -1, 2, b'/');

/// REF followed by ALT alleles, upper-cased, "." ALTs dropped. Reused
/// across lines to avoid per-line allocation.
#[derive(Default)]
pub struct Alleles {
    bytes: Vec<u8>,
    ends: Vec<usize>,
}

impl Alleles {
    pub fn len(&self) -> usize {
        self.ends.len()
    }

    pub fn get(&self, i: usize) -> &[u8] {
        let start = if i == 0 { 0 } else { self.ends[i - 1] };
        &self.bytes[start..self.ends[i]]
    }

    fn push(&mut self, a: &[u8]) {
        self.bytes.extend(a.iter().map(u8::to_ascii_uppercase));
        self.ends.push(self.bytes.len());
    }
}

/// The fixed columns of a data line plus the raw sample text.
pub struct Site<'a> {
    pub chrom: &'a [u8],
    pub pos: i32,
    pub id: &'a [u8],
    pub qual: &'a [u8],
    pub filter: &'a [u8],
    /// Indices of GT, DP and GQ within FORMAT, or -1.
    pub gt_idx: i32,
    pub dp_idx: i32,
    pub gq_idx: i32,
    pub samples: &'a [u8],
}

impl<'a> Site<'a> {
    /// Parses the fixed columns of `line`; alleles are written to `alleles`.
    pub fn parse(line: &'a [u8], alleles: &mut Alleles) -> Self {
        let mut rest = line;
        let mut next = || {
            let end = memchr(b'\t', rest).unwrap_or(rest.len());
            let f = &rest[..end];
            rest = if end < rest.len() { &rest[end + 1..] } else { b"" };
            f
        };
        let chrom = next();
        let pos = atoi(next());
        let id = next();
        let reference = next();
        let alt = next();
        let qual = next();
        let filter = next();
        let _info = next();
        let format = next();

        alleles.bytes.clear();
        alleles.ends.clear();
        alleles.push(reference);
        // VCFtools splits ALT with getline(',') until eof; "." entries are skipped.
        for a in alt.split(|&c| c == b',') {
            if a != b"." {
                alleles.push(a);
            }
        }

        // FORMAT_to_idx is a map, so a repeated key keeps its last position.
        let (mut gt_idx, mut dp_idx, mut gq_idx) = (-1, -1, -1);
        if !format.is_empty() {
            for (i, key) in format.split(|&c| c == b':').enumerate() {
                match key {
                    b"GT" => gt_idx = i as i32,
                    b"DP" => dp_idx = i as i32,
                    b"GQ" => gq_idx = i as i32,
                    _ => {}
                }
            }
        }
        Site { chrom, pos, id, qual, filter, gt_idx, dp_idx, gq_idx, samples: rest }
    }

    /// Calls `f(i, sub)` with the FORMAT sub-field `idx` of each of the
    /// `n_indv` samples (None when absent), as parse_genotype_entry would
    /// find it.
    pub fn for_each_subfield(&self, n_indv: usize, idx: i32, mut f: impl FnMut(usize, Option<&[u8]>)) {
        self.for_each_subfield2(n_indv, idx, -1, |i, a, _| f(i, a));
    }

    /// Like `for_each_subfield`, but extracts two sub-fields in one pass.
    pub fn for_each_subfield2(
        &self,
        n_indv: usize,
        idx_a: i32,
        idx_b: i32,
        mut f: impl FnMut(usize, Option<&[u8]>, Option<&[u8]>),
    ) {
        let mut rest = self.samples;
        for i in 0..n_indv {
            let end = memchr(b'\t', rest).unwrap_or(rest.len());
            let (a, b) = subfields2(&rest[..end], idx_a, idx_b);
            f(i, a, b);
            rest = if end < rest.len() { &rest[end + 1..] } else { b"" };
        }
    }

    /// Decodes all `n_indv` sample genotypes into `out`.
    pub fn genotypes(&self, n_indv: usize, out: &mut Vec<Gt>) {
        out.clear();
        out.reserve(n_indv);
        if self.gt_idx == 0 && n_indv > 0 && self.samples.len() + 1 == 4 * n_indv && uniform_row(self.samples) {
            // Every column is "a|b" or "a/b": decode without per-column checks.
            let (body, last) = self.samples.split_at(4 * (n_indv - 1));
            out.extend(body.chunks_exact(4).map(|c| decode3(c[0], c[1], c[2])));
            out.push(decode3(last[0], last[1], last[2]));
            return;
        }
        let s = self.samples;
        let mut i = 0;
        for _ in 0..n_indv {
            // Fast path: the whole column is a 3-byte diploid GT ("0|1").
            if self.gt_idx == 0
                && i + 3 <= s.len()
                && (s[i + 1] == b'|' || s[i + 1] == b'/')
                && (i + 3 == s.len() || s[i + 3] == b'\t')
            {
                out.push(decode_gt(&s[i..i + 3]));
                i += 4;
                continue;
            }
            let rest = if i < s.len() { &s[i..] } else { b"" };
            let end = memchr(b'\t', rest).unwrap_or(rest.len());
            out.push(genotype(&rest[..end], self.gt_idx));
            i += end + 1;
        }
    }
}

/// Extracts and decodes the GT sub-field of one sample column
/// (vcf_entry::parse_genotype_entry + set_indv_GENOTYPE_and_PHASE).
#[inline]
pub fn genotype(field: &[u8], gt_idx: i32) -> Gt {
    if gt_idx == 0 && field.len() == 3 && (field[1] == b'|' || field[1] == b'/') {
        return decode_gt(field);
    }
    match subfield(field, gt_idx) {
        Some(s) => decode_gt(s),
        None => MISSING_DIPLOID,
    }
}

/// Sub-fields `a` and `b` (either may be -1) of a sample column, found in
/// a single scan; each matches what `subfield` would return.
#[inline]
fn subfields2(field: &[u8], a: i32, b: i32) -> (Option<&[u8]>, Option<&[u8]>) {
    let (mut ra, mut rb) = (None, None);
    let last = a.max(b);
    if last < 0 {
        return (ra, rb);
    }
    let (mut i, mut start) = (0, 0);
    loop {
        let end = memchr(b':', &field[start..]).map(|o| start + o);
        let sub = &field[start..end.unwrap_or(field.len())];
        if i == a {
            ra = Some(sub);
        }
        if i == b {
            rb = Some(sub);
        }
        match end {
            Some(e) if i < last => {
                i += 1;
                start = e + 1;
            }
            _ => return (ra, rb),
        }
    }
}

/// The `idx`-th ':'-separated sub-field of a sample column.
#[inline]
fn subfield(field: &[u8], idx: i32) -> Option<&[u8]> {
    if idx < 0 {
        return None;
    }
    let (mut i, mut start) = (0, 0);
    while let Some(off) = memchr(b':', &field[start..]) {
        if i == idx {
            return Some(&field[start..start + off]);
        }
        i += 1;
        start += off + 1;
    }
    (i == idx).then(|| &field[start..])
}

/// True if `s` (length 4n-1) is n columns of 3-byte diploid GTs separated
/// by tabs. Non-short-circuiting so the compiler can vectorize it.
#[inline]
fn uniform_row(s: &[u8]) -> bool {
    let (body, last) = s.split_at(s.len() - 3);
    let sep_ok = |c: u8| (c == b'|') | (c == b'/');
    body.chunks_exact(4).fold(true, |ok, c| ok & sep_ok(c[1]) & (c[3] == b'\t')) & sep_ok(last[1])
}

#[inline(always)]
fn decode3(a: u8, sep: u8, b: u8) -> Gt {
    let dec = |c: u8| if c == b'.' { -1 } else { c as i32 - b'0' as i32 };
    Gt::new(dec(a), dec(b), 2, sep)
}

#[inline]
fn decode_gt(s: &[u8]) -> Gt {
    if s.len() == 3 && (s[1] == b'/' || s[1] == b'|') {
        let a = if s[0] == b'.' { -1 } else { s[0] as i32 - b'0' as i32 };
        let b = if s[2] == b'.' { -1 } else { s[2] as i32 - b'0' as i32 };
        return Gt::new(a, b, 2, s[1]);
    }
    let sep = |c: &u8| *c == b'/' || *c == b'|';
    match s.iter().position(sep) {
        Some(p) => {
            if s.iter().rposition(sep) != Some(p) {
                crate::fatal("Polyploidy found, and not supported by vcftools");
            }
            Gt::new(str2int(&s[..p]), str2int(&s[p + 1..]), 2, s[p])
        }
        None => Gt::new(str2int(s), -1, 1, b'|'),
    }
}

/// header::str2int: "" or "." is missing (-1), otherwise C `atoi`.
pub fn str2int(s: &[u8]) -> i32 {
    if s.is_empty() || s == b"." {
        -1
    } else {
        atoi(s)
    }
}

/// C `atof`. Plain decimals (`12`, `3.5`, `.5`) are parsed in Rust, whose
/// float parsing is correctly rounded like glibc's strtod, so results are
/// bit-identical; anything else (sign, exponent, inf/nan, whitespace,
/// trailing text) goes to the C library.
pub fn atof(s: &[u8]) -> f64 {
    if is_plain_decimal(s) {
        // SAFETY: plain decimals are ASCII.
        if let Ok(v) = unsafe { std::str::from_utf8_unchecked(s) }.parse::<f64>() {
            return v;
        }
    }
    let mut stack = [0u8; 64];
    let mut heap;
    let buf: &mut [u8] = if s.len() < stack.len() {
        &mut stack[..s.len() + 1]
    } else {
        heap = vec![0u8; s.len() + 1];
        &mut heap
    };
    buf[..s.len()].copy_from_slice(s);
    buf[s.len()] = 0;
    unsafe { libc::atof(buf.as_ptr() as *const libc::c_char) }
}

/// `[0-9]+` or `[0-9]*\.[0-9]+` or `[0-9]+\.`
fn is_plain_decimal(s: &[u8]) -> bool {
    let mut digits = 0;
    let mut dots = 0;
    for &c in s {
        match c {
            b'0'..=b'9' => digits += 1,
            b'.' => dots += 1,
            _ => return false,
        }
    }
    digits > 0 && dots <= 1
}

/// header::str2double: "" or "." is missing (-1), otherwise C `atof`.
pub fn str2double(s: &[u8]) -> f64 {
    if s.is_empty() || s == b"." {
        -1.0
    } else {
        atof(s)
    }
}

/// C `atoi`: optional leading whitespace and sign, then digits; 0 if none.
pub fn atoi(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && matches!(s[i], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
        i += 1;
    }
    let neg = i < s.len() && s[i] == b'-';
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        i += 1;
    }
    let mut v: i64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        v = v.wrapping_mul(10).wrapping_add((s[i] - b'0') as i64);
        i += 1;
    }
    (if neg { -v } else { v }) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotypes() {
        assert_eq!(genotype(b"0|1", 0), Gt::new(0, 1, 2, b'|'));
        assert_eq!(genotype(b"./.", 0), Gt::new(-1, -1, 2, b'/'));
        assert_eq!(genotype(b"1", 0), Gt::new(1, -1, 1, b'|'));
        assert_eq!(genotype(b".", 0), Gt::new(-1, -1, 1, b'|'));
        assert_eq!(genotype(b"10/2:35", 0), Gt::new(10, 2, 2, b'/'));
        assert_eq!(genotype(b"35:0/1", 1), Gt::new(0, 1, 2, b'/'));
        assert_eq!(genotype(b"35", 1), MISSING_DIPLOID);
        assert_eq!(genotype(b"0|1", -1), MISSING_DIPLOID);
    }

    #[test]
    fn subfields_match_single_lookup() {
        let cols: [&[u8]; 6] = [b"0/1:35:12", b"0/1:35", b"0/1", b"", b"::", b"a:b:c:d"];
        for f in cols {
            for a in -1..5 {
                for b in -1..5 {
                    assert_eq!(subfields2(f, a, b), (subfield(f, a), subfield(f, b)), "{f:?} {a} {b}");
                }
            }
        }
    }

    #[test]
    fn atof_matches_libc() {
        let cases: [&[u8]; 16] = [
            b"0", b"12", b"35.5", b".5", b"7.", b"99.99", b"0.1", b"123456789.123456789", b"1e2", b"-3", b"+4",
            b" 5", b"nan", b"inf", b"12abc", b"",
        ];
        for s in cases {
            let mut buf = s.to_vec();
            buf.push(0);
            let want = unsafe { libc::atof(buf.as_ptr() as *const libc::c_char) };
            let got = atof(s);
            assert!(got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()), "{s:?}: {got} vs {want}");
        }
        // Random plain decimals, including long ones that stress rounding.
        let mut x: u64 = 0x9E3779B97F4A7C15;
        let mut next = |m: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % m
        };
        for _ in 0..200_000 {
            let mut s = Vec::new();
            for _ in 0..next(20) {
                s.push(b'0' + next(10) as u8);
            }
            if next(2) == 0 {
                s.push(b'.');
                for _ in 0..next(25) {
                    s.push(b'0' + next(10) as u8);
                }
            }
            let mut buf = s.clone();
            buf.push(0);
            let want = unsafe { libc::atof(buf.as_ptr() as *const libc::c_char) };
            assert_eq!(atof(&s).to_bits(), want.to_bits(), "{:?}", String::from_utf8_lossy(&s));
        }
    }
}
