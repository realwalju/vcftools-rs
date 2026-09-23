//! Number formatting that matches C++ iostream output byte-for-byte.
//! The default `ostream << double` is printf("%g") with precision 6;
//! `setprecision`/`fixed` map to "%.*g" / "%.*f". Delegating to the C
//! library's snprintf guarantees identical rounding and nan/inf spelling.

use std::os::raw::c_char;

#[inline]
fn printf_f64(out: &mut Vec<u8>, fmt: &[u8], prec: i32, x: f64) {
    let mut buf = [0u8; 64];
    let n = unsafe {
        libc::snprintf(buf.as_mut_ptr() as *mut c_char, buf.len(), fmt.as_ptr() as *const c_char, prec, x)
    };
    if n >= 0 && (n as usize) < buf.len() {
        out.extend_from_slice(&buf[..n as usize]);
    } else {
        // Very large fixed-point values; fall back to an allocation.
        let mut big = vec![0u8; n.max(0) as usize + 1];
        let m = unsafe {
            libc::snprintf(big.as_mut_ptr() as *mut c_char, big.len(), fmt.as_ptr() as *const c_char, prec, x)
        };
        out.extend_from_slice(&big[..m.max(0) as usize]);
    }
}

/// `out << x` with default stream flags.
#[inline]
pub fn g(out: &mut Vec<u8>, x: f64) {
    printf_f64(out, b"%.*g\0", 6, x);
}

/// `out << scientific << setprecision(p) << x`.
#[inline]
pub fn sci(out: &mut Vec<u8>, p: i32, x: f64) {
    printf_f64(out, b"%.*e\0", p, x);
}

/// `out << fixed << setprecision(p) << x`.
#[inline]
pub fn fixed(out: &mut Vec<u8>, p: i32, x: f64) {
    printf_f64(out, b"%.*f\0", p, x);
}

#[inline]
pub fn int(out: &mut Vec<u8>, x: impl itoa::Integer) {
    out.extend_from_slice(itoa::Buffer::new().format(x).as_bytes());
}

#[cfg(test)]
mod tests {
    #[test]
    fn matches_iostream() {
        let mut v = Vec::new();
        super::g(&mut v, 0.000199681469648562);
        assert_eq!(v, b"0.000199681");
        v.clear();
        let zero = std::hint::black_box(0.0f64);
        super::g(&mut v, zero / zero);
        assert_eq!(v, b"-nan");
    }
}
