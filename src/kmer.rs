//! 2-bit packing and the rolling canonical scan the adapter search runs on.
//!
//! Adapter k-mers are matched in canonical form because adscrub runs before any
//! strand normalisation: a read carrying the adapter on either strand must be
//! caught.

/// Two-bit code of a base, or `None` for anything that is not ACGT.
#[inline]
pub fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
pub fn revcomp_packed(mut code: u64, k: usize) -> u64 {
    let mut out = 0u64;
    for _ in 0..k {
        out = (out << 2) | (3 - (code & 3));
        code >>= 2;
    }
    out
}

#[inline]
pub fn canonical_packed(code: u64, k: usize) -> u64 {
    let reverse = revcomp_packed(code, k);
    if reverse < code {
        reverse
    } else {
        code
    }
}

/// Canonical packing of one k-length window, or `None` if it holds a non-ACGT base.
pub fn canonical_packed_of(seq: &[u8], k: usize) -> Option<u64> {
    if seq.len() < k {
        return None;
    }
    let mut code = 0u64;
    for &base in &seq[..k] {
        code = (code << 2) | base_code(base)?;
    }
    Some(canonical_packed(code, k))
}

/// Visit every k-window as `(start, Some(canonical))`, or `(start, None)` where the
/// window holds a non-ACGT base. Rolling, so no per-position allocation.
pub fn for_each_canonical<F: FnMut(usize, Option<u64>)>(seq: &[u8], k: usize, mut visit: F) {
    if seq.len() < k || k == 0 || k > 32 {
        return;
    }
    let mask = if k * 2 >= 64 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut code = 0u64;
    let mut filled = 0usize;
    for (index, &base) in seq.iter().enumerate() {
        match base_code(base) {
            Some(value) => {
                code = ((code << 2) | value) & mask;
                filled += 1;
                if filled >= k {
                    visit(index + 1 - k, Some(canonical_packed(code, k)));
                }
            }
            None => {
                filled = 0;
                code = 0;
                if index + 1 >= k {
                    visit(index + 1 - k, None);
                }
            }
        }
    }
}
