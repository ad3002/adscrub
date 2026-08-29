//! Phase B — place the fragment end from the mate overlap, for read-through too
//! short to carry a detectable adapter k-mer.
//!
//! When the insert is shorter than the read, read 1 and the reverse complement of
//! read 2 align over the insert at a positive shift `i`; the leading `i` bases of
//! RC(read 2) and the matching tail of read 1 are read-through. Ported from
//! `stirka::overlap` (2016 `Compute_trims.cpp`, "Overlap too short"), with two
//! changes that are the whole point of the port.
//!
//! # Why this phase cannot stand on its own
//! Sequence similarity between mates is not evidence that the fragment ended. In a
//! tandem array of period `p` the two mates match at every shift `i`, `i+p`,
//! `i+2p`, … and a pair lying *entirely inside* an array — no read-through at all —
//! matches at a small shift by construction. Trimming on that cuts real array
//! sequence, consistently, from every read over the array: a collapsed tandem
//! manufactured by the trimmer, in the sequence class an assembler cares about
//! most. The ported code took the smallest matching shift, which is conservative
//! when a read-through is real and wrong when it is not.
//!
//! Two guards, and both are refusals rather than tie-breaks:
//!
//! 1. **The shift must be unique.** Every shift is scanned, not just the first;
//!    more than one match means the alignment does not identify a fragment end,
//!    which is exactly the tandem signature. Refuse.
//! 2. **An adapter licenses the cut.** By default the pair must carry a bundled
//!    adapter k-mer somewhere — the evidence that read-through happened at all.
//!    The overlap then says *where* it started, which is earlier than the first
//!    detectable adapter k-mer, and that is what phase B adds over phase A.
//!
//! Disabling the licence widens the catch to read-through with no detectable
//! adapter, and re-opens the tandem risk to whatever the uniqueness guard does not
//! catch. It is a flag, not a default.

use rayon::prelude::*;

use crate::adapter::AdapterSet;
use crate::history::{Deletion, Reason, ReadHistory};

/// Shorter overlaps are too chance-prone to identify a fragment end.
pub const MIN_OVERLAP: usize = 30;
/// Allowed mismatch fraction over the overlap — sequencing error in the insert.
pub const MAX_MISMATCH_RATE: f64 = 0.1;

#[derive(Default, Debug, Clone, Copy)]
pub struct OverlapStats {
    pub pairs_examined: usize,
    pub pairs_trimmed: usize,
    pub bases_trimmed: usize,
    /// Pairs where more than one shift matched — the tandem signature.
    pub refused_ambiguous: usize,
    /// Pairs with a unique shift but no adapter anywhere in the pair.
    pub refused_unlicensed: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapConfig {
    pub min_overlap: usize,
    pub max_mismatch_rate: f64,
    /// Require a bundled adapter k-mer in the pair before trimming. Turning this
    /// off trades the tandem guarantee for reach; see the module comment.
    pub require_adapter_licence: bool,
}

impl Default for OverlapConfig {
    fn default() -> Self {
        Self {
            min_overlap: MIN_OVERLAP,
            max_mismatch_rate: MAX_MISMATCH_RATE,
            require_adapter_licence: true,
        }
    }
}

#[inline]
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn matches_at(first: &[u8], second_rc: &[u8], shift: usize, overlap: usize, budget: usize) -> bool {
    let mut mismatches = 0usize;
    for offset in 0..overlap {
        if first[offset] != second_rc[shift + offset] {
            mismatches += 1;
            if mismatches > budget {
                return false;
            }
        }
    }
    true
}

/// Outcome of looking for a fragment end in one pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shift {
    /// Exactly one shift explains the pair: `(shift, overlap)`.
    Unique(usize, usize),
    /// No shift explains it — the ordinary long-insert case.
    None,
    /// Several shifts explain it equally well, so the fragment end is not
    /// identifiable. A tandem array does this.
    Ambiguous,
}

/// Every shift `i > 0` at which `first[0..ov]` matches `second_rc[i..i+ov]` within
/// the mismatch budget, with `ov = min(first.len, second_rc.len - i) >= min_overlap`.
/// The overlap only shrinks as `i` grows, so the scan stops at the floor.
pub fn find_fragment_end(first: &[u8], second_rc: &[u8], config: &OverlapConfig) -> Shift {
    let mut found: Option<(usize, usize)> = None;
    let mut shift = 1usize;
    while shift < second_rc.len() {
        let overlap = (second_rc.len() - shift).min(first.len());
        if overlap < config.min_overlap {
            break;
        }
        let budget = (overlap as f64 * config.max_mismatch_rate) as usize;
        if matches_at(first, second_rc, shift, overlap, budget) {
            if found.is_some() {
                return Shift::Ambiguous;
            }
            found = Some((shift, overlap));
        }
        shift += 1;
    }
    match found {
        Some((shift, overlap)) => Shift::Unique(shift, overlap),
        None => Shift::None,
    }
}

fn pop_tail(seq: &mut Vec<u8>, count: usize, history: &mut ReadHistory) -> usize {
    for _ in 0..count {
        let Some(base) = seq.pop() else { break };
        history.deletions.push(Deletion {
            position: seq.len(),
            base: base as char,
            reason: Reason::PairReadThrough,
        });
    }
    count
}

/// Trim one pair. Returns bases removed, 0 if the pair was left alone.
pub fn trim_pair(
    first: &mut Vec<u8>,
    first_history: &mut ReadHistory,
    second: &mut Vec<u8>,
    second_history: &mut ReadHistory,
    adapters: &AdapterSet,
    config: &OverlapConfig,
    stats: &mut OverlapStats,
) -> usize {
    stats.pairs_examined += 1;
    if first.len() < config.min_overlap || second.len() < config.min_overlap {
        return 0;
    }
    let second_rc = revcomp(second);
    let (shift, overlap) = match find_fragment_end(first, &second_rc, config) {
        Shift::Unique(shift, overlap) => (shift, overlap),
        Shift::Ambiguous => {
            stats.refused_ambiguous += 1;
            return 0;
        }
        Shift::None => return 0,
    };
    // The licence asks only whether this pair read into an adapter at all, so it
    // uses the short seeds rather than the 23-mer trim set: read-through shorter
    // than 23 bases carries no 23-mer, and that is exactly the case phase B is for.
    if config.require_adapter_licence
        && !crate::adapter::licence_hit(first)
        && !crate::adapter::licence_hit(second)
    {
        stats.refused_unlicensed += 1;
        return 0;
    }
    let mut bases = pop_tail(first, first.len() - overlap, first_history);
    bases += pop_tail(second, shift, second_history);
    stats.pairs_trimmed += 1;
    stats.bases_trimmed += bases;
    bases
}

/// Run phase B over consecutive mates `(2i, 2i+1)`. A trailing unpaired read is
/// left untouched.
pub fn process(
    sequences: &mut [Vec<u8>],
    histories: &mut [ReadHistory],
    adapters: &AdapterSet,
    config: &OverlapConfig,
) -> OverlapStats {
    let pairs = sequences.len() / 2;
    let stats: Vec<OverlapStats> = sequences
        .par_chunks_mut(2)
        .zip(histories.par_chunks_mut(2))
        .take(pairs)
        .map(|(seq_pair, history_pair)| {
            let mut local = OverlapStats::default();
            let (left_seq, right_seq) = seq_pair.split_at_mut(1);
            let (left_history, right_history) = history_pair.split_at_mut(1);
            trim_pair(
                &mut left_seq[0],
                &mut left_history[0],
                &mut right_seq[0],
                &mut right_history[0],
                adapters,
                config,
                &mut local,
            );
            local
        })
        .collect();
    stats.into_iter().fold(OverlapStats::default(), |mut total, one| {
        total.pairs_examined += one.pairs_examined;
        total.pairs_trimmed += one.pairs_trimmed;
        total.bases_trimmed += one.bases_trimmed;
        total.refused_ambiguous += one.refused_ambiguous;
        total.refused_unlicensed += one.refused_unlicensed;
        total
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::history::reverse_apply;

    fn adapters() -> AdapterSet {
        AdapterSet::bundled()
    }

    /// A short insert with real adapter on both mates: the case phase B is for.
    #[test]
    fn trims_read_through_reversibly() {
        let insert = b"ACGTTGCAACCGGTTAACCGGTTACGTACGTTGCAACCGGTTAACCGGTTACGTA";
        // A real adapter start, so the pair is licensed the way a real one is.
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCC";
        let mut first = insert.to_vec();
        first.extend_from_slice(adapter);
        let mut second = revcomp(insert);
        second.extend_from_slice(adapter);
        let (before_first, before_second) = (first.clone(), second.clone());
        let mut first_history = ReadHistory::new("r1");
        let mut second_history = ReadHistory::new("r2");
        let mut stats = OverlapStats::default();
        let bases = trim_pair(
            &mut first,
            &mut first_history,
            &mut second,
            &mut second_history,
            &adapters(),
            &OverlapConfig::default(),
            &mut stats,
        );
        assert!(bases > 0, "read-through pair trimmed");
        assert_eq!(stats.pairs_trimmed, 1);
        assert_eq!(reverse_apply(&first, &first_history.deletions), before_first);
        assert_eq!(reverse_apply(&second, &second_history.deletions), before_second);
    }

    /// The guard this port exists for: a pair lying entirely inside a tandem array,
    /// with no read-through at all, must not lose a single base.
    #[test]
    fn tandem_array_is_refused() {
        let unit = b"ACGTTGCAATCG"; // period 12
        let mut array = Vec::new();
        while array.len() < 200 {
            array.extend_from_slice(unit);
        }
        let mut first = array[0..100].to_vec();
        let mut second = revcomp(&array[40..140]);
        let (before_first, before_second) = (first.clone(), second.clone());
        let mut first_history = ReadHistory::new("t1");
        let mut second_history = ReadHistory::new("t2");
        let mut stats = OverlapStats::default();
        let bases = trim_pair(
            &mut first,
            &mut first_history,
            &mut second,
            &mut second_history,
            &adapters(),
            &OverlapConfig::default(),
            &mut stats,
        );
        assert_eq!(bases, 0, "a tandem pair must not be trimmed");
        assert_eq!(first, before_first);
        assert_eq!(second, before_second);
        assert!(
            stats.refused_ambiguous > 0 || stats.refused_unlicensed > 0,
            "and the refusal must be counted, not silent"
        );
    }

    /// A long insert: the mates do not overlap and nothing happens.
    #[test]
    fn long_insert_untouched() {
        let mut first = b"ACGTACGTTTGCAAGGCCTTAACCGGATCGATCGGGCCTTAAGGCCATATCGCG".to_vec();
        let mut second = b"TTTTGCGCGCATATCCCGGGAAATTTCCCGGGTTTAAACCCGGGTTTAAACCCG".to_vec();
        let (before_first, before_second) = (first.clone(), second.clone());
        let mut first_history = ReadHistory::new("l1");
        let mut second_history = ReadHistory::new("l2");
        let mut stats = OverlapStats::default();
        trim_pair(
            &mut first,
            &mut first_history,
            &mut second,
            &mut second_history,
            &adapters(),
            &OverlapConfig::default(),
            &mut stats,
        );
        assert_eq!(first, before_first);
        assert_eq!(second, before_second);
    }

    /// Unique shift, but no adapter anywhere: refused by default, and counted.
    #[test]
    fn unlicensed_overlap_is_refused_by_default() {
        let insert = b"ACGTTGCAACCGGTTAACCGGTTACGTACGTTGCAACCGGTTAACCGGTTACGTA";
        let mut first = insert.to_vec();
        first.extend_from_slice(b"TTTTTTTT");
        let mut second = revcomp(insert);
        second.extend_from_slice(b"AAAAAAAA");
        assert!(!crate::adapter::licence_hit(&first) && !crate::adapter::licence_hit(&second));
        let mut first_history = ReadHistory::new("u1");
        let mut second_history = ReadHistory::new("u2");
        let mut stats = OverlapStats::default();
        let bases = trim_pair(
            &mut first,
            &mut first_history,
            &mut second,
            &mut second_history,
            &adapters(),
            &OverlapConfig::default(),
            &mut stats,
        );
        assert_eq!(bases, 0);
        assert_eq!(stats.refused_unlicensed, 1);
    }
}
