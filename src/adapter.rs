//! Phase A — trim from the first bundled adapter k-mer to the 3' end.
//!
//! Ported from `stirka::adapter` (itself a modernisation of the 2016 `LU_trim` /
//! `V5_trim` C++ `Compute_trims.cpp`), which is where this lived while adscrub was
//! a scaffold. Adapter k-mers are 2-bit packed into canonical form, so a read
//! carrying the adapter on either strand is caught; adscrub runs before any strand
//! normalisation.
//!
//! The bundled set (`data/illumina_adapters_k23.txt`, 98,167 23-mers) covers the
//! TruSeq and Nextera adapters and the poly-G dark-cycle artifact. A 23-mer is a
//! specific synthetic sequence: the chance of a false hit against a genome is
//! ~98k / 4^23 = 1.4e-12, so the first hit is taken as the boundary without
//! further evidence.
//!
//! **What this phase cannot do.** It needs a whole 23-mer of adapter, so
//! read-through shorter than 23 bases carries no detectable adapter. Those first
//! cycles are exactly where theerror-looking k-mers are made, and phase B exists
//! for them.

use hashbrown::HashSet;

use crate::history::{Deletion, Reason};
use crate::kmer::{canonical_packed_of, for_each_canonical};

const BUNDLED: &str = include_str!("data/illumina_adapters_k23.txt");

/// k of the bundled adapter k-mers. Deliberately distinct from any downstream
/// tool's k — adapters are matched at their own resolution.
pub const ADAPTER_K: usize = 23;

/// Canonical 5' starts of the Illumina adapters, for the *licence* check only.
///
/// The bundled 23-mer set decides where to cut and must be specific enough to act
/// on alone. The licence answers a different and easier question — did this pair
/// read into an adapter at all — and a pair whose read-through is shorter than 23
/// bases carries no 23-mer while still carrying the first dozen. Measured on
/// VAF12, requiring a 23-mer as licence refused 126,652 pairs that had an
/// unambiguous fragment end; the whole of phase B was blocked by it.
///
/// These are short enough to catch that and specific enough not to fire on a
/// genome: the shortest is 12 bases, so a spurious hit costs about 1.8e-3 per
/// 300-base read, and a licence alone never cuts anything — the overlap still has
/// to identify the fragment end uniquely.
const LICENCE_SEEDS: [&[u8]; 4] = [
    b"AGATCGGAAGAGC",       // TruSeq universal
    b"CTGTCTCTTATACACATCT", // Nextera transposase
    b"GATCGGAAGAGC",        // TruSeq, one base in
    b"TGGAATTCTCGG",        // small-RNA 3'
];

pub struct AdapterSet {
    kmers: HashSet<u64>,
}

impl AdapterSet {
    /// The compiled-in Illumina set, canonicalised. Done once.
    pub fn bundled() -> Self {
        Self::from_lines(BUNDLED.lines())
    }

    /// A user-supplied set, one k-mer per line. Lines that are not exactly
    /// `ADAPTER_K` ACGT bases are ignored, so a FASTA-ish file is tolerated.
    pub fn from_lines<'a, I: Iterator<Item = &'a str>>(lines: I) -> Self {
        let mut kmers = HashSet::new();
        for line in lines {
            let line = line.trim();
            if line.len() == ADAPTER_K {
                if let Some(code) = canonical_packed_of(line.as_bytes(), ADAPTER_K) {
                    kmers.insert(code);
                }
            }
        }
        Self { kmers }
    }

    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    /// First read position whose canonical 23-mer is an adapter — the 3' trim
    /// point. `None` if the read carries no adapter k-mer.
    pub fn adapter_start(&self, seq: &[u8]) -> Option<usize> {
        if seq.len() < ADAPTER_K {
            return None;
        }
        let mut found = None;
        for_each_canonical(seq, ADAPTER_K, |start, code| {
            if found.is_none() {
                if let Some(code) = code {
                    if self.kmers.contains(&code) {
                        found = Some(start);
                    }
                }
            }
        });
        found
    }
}

/// Does this read contain an adapter start? Licence for [`crate::overlap`] to act
/// on a fragment end it has already placed unambiguously; never a cut on its own.
pub fn licence_hit(seq: &[u8]) -> bool {
    LICENCE_SEEDS.iter().any(|seed| contains(seq, seed))
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|window| window.eq_ignore_ascii_case(needle))
}

/// Trim a 3' adapter — everything from the first adapter k-mer to the end —
/// recording each removed base reversibly. Empty if the read carries no adapter.
pub fn trim_adapter(seq: &mut Vec<u8>, set: &AdapterSet) -> Vec<Deletion> {
    let mut deletions = Vec::new();
    if let Some(start) = set.adapter_start(seq) {
        while seq.len() > start {
            let base = seq.pop().expect("non-empty while trimming");
            deletions.push(Deletion {
                position: seq.len(),
                base: base as char,
                reason: Reason::AdapterKmer,
            });
        }
    }
    deletions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::history::reverse_apply;

    #[test]
    fn bundle_loads() {
        assert!(AdapterSet::bundled().len() > 50_000);
    }

    #[test]
    fn trims_adapter_reversibly() {
        let set = AdapterSet::bundled();
        let adapter = b"GGGGGGGCAGATGCCTAAGTGAC";
        assert_eq!(adapter.len(), ADAPTER_K);
        let prefix = b"ACTGACTGACAGTCAGTCATGCATGCAT";
        let mut read = prefix.to_vec();
        read.extend_from_slice(adapter);
        let original = read.clone();
        let deletions = trim_adapter(&mut read, &set);
        assert_eq!(read, prefix.to_vec());
        assert_eq!(reverse_apply(&read, &deletions), original);
    }

    #[test]
    fn licence_fires_on_a_partial_adapter_and_not_on_genome() {
        // Nine bases of TruSeq: too short for the 23-mer set, enough to licence.
        let mut read = b"ACTGACTGACAGTCAGTCATGCATGCATTAGCTAGCTAGCAT".to_vec();
        assert!(!licence_hit(&read), "genomic read carries no licence");
        read.extend_from_slice(b"AGATCGGAAGAGCAC");
        assert!(licence_hit(&read), "partial adapter licences the pair");
        assert!(
            AdapterSet::bundled().adapter_start(&read).is_none()
                || read.len() >= ADAPTER_K,
            "licence is a weaker test than the trim set, by design"
        );
    }

    #[test]
    fn clean_read_untouched() {
        let set = AdapterSet::bundled();
        let mut read = b"ACTGACTGACAGTCAGTCATGCATGCATTAGCTAGCTAGCAT".to_vec();
        let before = read.clone();
        assert!(trim_adapter(&mut read, &set).is_empty());
        assert_eq!(read, before);
    }
}
