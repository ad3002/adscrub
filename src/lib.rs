//! adscrub — adapter and read-through scrubbing for short-read pairs.
//!
//! The first stage of the genemayers pipeline: physical contamination removal,
//! before anything algorithmic. What the sequencer wrote that is not DNA has to
//! leave here, because downstream every base is treated as evidence — adapter
//! k-mers become "real" high-frequency k-mers, and read-through bases become
//! minority variants that no k-mer-spectrum argument can recognise as artefacts.
//!
//! Two phases, in order:
//!
//! * [`adapter`] — a bundled Illumina adapter k-mer matched anywhere in the read
//!   trims from there to the 3' end. Unambiguous, per-read, and blind to
//!   read-through shorter than the 23-base k-mer.
//! * [`overlap`] — the mate overlap places the fragment end for exactly that
//!   blind spot. It refuses whenever the placement is not identifiable, which is
//!   what keeps it from shortening tandem arrays.
//!
//! Every cut is recorded reversibly ([`history`]); a trim that cannot be undone
//! is a deletion of evidence at the one stage that runs before anything could
//! re-derive it.

pub mod adapter;
pub mod history;
pub mod kmer;
pub mod overlap;

pub use adapter::{trim_adapter, AdapterSet, ADAPTER_K};
pub use history::{reverse_apply, Deletion, ReadHistory, Reason};
pub use overlap::{OverlapConfig, OverlapStats, Shift};

/// What one scrubbing run did.
#[derive(Debug, Default, Clone, Copy)]
pub struct ScrubStats {
    pub reads: usize,
    pub reads_adapter_trimmed: usize,
    pub bases_adapter_trimmed: usize,
    pub overlap: OverlapStats,
}

/// Scrub a set of reads in place. `paired` treats consecutive entries as mates.
pub fn scrub(
    sequences: &mut [Vec<u8>],
    histories: &mut [ReadHistory],
    adapters: &AdapterSet,
    config: &OverlapConfig,
    paired: bool,
) -> ScrubStats {
    let mut stats = ScrubStats { reads: sequences.len(), ..Default::default() };
    for (sequence, history) in sequences.iter_mut().zip(histories.iter_mut()) {
        let deletions = trim_adapter(sequence, adapters);
        if !deletions.is_empty() {
            stats.reads_adapter_trimmed += 1;
            stats.bases_adapter_trimmed += deletions.len();
            history.deletions.extend(deletions);
        }
    }
    if paired {
        stats.overlap = overlap::process(sequences, histories, adapters, config);
    }
    stats
}
