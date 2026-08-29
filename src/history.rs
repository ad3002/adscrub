//! Every cut is reversible. A trim that cannot be undone is a deletion of
//! evidence, and this stage runs before anything that could re-derive it.

use serde::{Deserialize, Serialize};

/// One removed base, with the position it occupied at the time of removal.
/// Trims pop from the 3' end, so replaying in reverse order restores the read.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Deletion {
    pub position: usize,
    pub base: char,
    /// Which rule removed it, so a history file says why and not only what.
    pub reason: Reason,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Reason {
    /// A bundled adapter k-mer matched at this position or later.
    AdapterKmer,
    /// The mate overlap placed the fragment end before this position.
    PairReadThrough,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadHistory {
    pub read_id: String,
    pub deletions: Vec<Deletion>,
}

impl ReadHistory {
    pub fn new(read_id: impl Into<String>) -> Self {
        Self { read_id: read_id.into(), deletions: Vec::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.deletions.is_empty()
    }
}

/// Restore a trimmed read from its history.
pub fn reverse_apply(trimmed: &[u8], deletions: &[Deletion]) -> Vec<u8> {
    let mut out = trimmed.to_vec();
    for deletion in deletions.iter().rev() {
        out.insert(deletion.position.min(out.len()), deletion.base as u8);
    }
    out
}
