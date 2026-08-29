//! adscrub CLI — scrub a tilde-pair read file and write the cleaned stream plus a
//! reversible history.
//!
//! Input is the pipeline's `mate1~mate2` seq-per-line format, gzip optional, as
//! sequenced: each mate written in the direction the sequencer read it. A file
//! whose mates have been reoriented cannot be scrubbed by the overlap phase,
//! because the overlap is what places the fragment end.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use adscrub::{scrub, AdapterSet, OverlapConfig, ReadHistory};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

fn open(path: &str) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::with_capacity(1 << 22, MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 22, file)))
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let mut reads = String::new();
    let mut out = String::new();
    let mut history_path = String::new();
    let mut adapters_path = String::new();
    let mut config = OverlapConfig::default();
    let mut index = 1;
    while index + 1 < args.len() {
        match args[index].as_str() {
            "--reads" => reads = args[index + 1].clone(),
            "--out" => out = args[index + 1].clone(),
            "--history" => history_path = args[index + 1].clone(),
            "--adapters" => adapters_path = args[index + 1].clone(),
            "--min-overlap" => config.min_overlap = args[index + 1].parse().expect("number"),
            "--max-mismatch-rate" => {
                config.max_mismatch_rate = args[index + 1].parse().expect("number")
            }
            // Widens the catch to read-through carrying no detectable adapter, and
            // gives up the guarantee that a tandem array is never shortened.
            "--no-adapter-licence" => {
                config.require_adapter_licence = args[index + 1].parse::<bool>().map(|v| !v).unwrap_or(false)
            }
            other => panic!("unknown option {other}"),
        }
        index += 2;
    }
    assert!(!reads.is_empty() && !out.is_empty(), "--reads and --out are required");

    let adapters = if adapters_path.is_empty() {
        AdapterSet::bundled()
    } else {
        let text = std::fs::read_to_string(&adapters_path)?;
        AdapterSet::from_lines(text.lines())
    };
    eprintln!("adapter set: {} canonical {}-mers", adapters.len(), adscrub::ADAPTER_K);

    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut histories: Vec<ReadHistory> = Vec::new();
    let mut reader = open(&reads)?;
    let mut line = String::new();
    let mut pair_index = 0usize;
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end();
        let mut mates = trimmed.split('~');
        let first = mates.next().unwrap_or("");
        let second = mates.next().unwrap_or("");
        sequences.push(first.as_bytes().to_vec());
        histories.push(ReadHistory::new(format!("{pair_index}/1")));
        sequences.push(second.as_bytes().to_vec());
        histories.push(ReadHistory::new(format!("{pair_index}/2")));
        pair_index += 1;
        line.clear();
    }
    eprintln!("read {pair_index} pairs");

    let stats = scrub(&mut sequences, &mut histories, &adapters, &config, true);

    let mut writer = BufWriter::new(GzEncoder::new(File::create(&out)?, Compression::new(1)));
    for pair in sequences.chunks(2) {
        writeln!(
            writer,
            "{}~{}",
            String::from_utf8_lossy(&pair[0]),
            String::from_utf8_lossy(pair.get(1).map(|v| v.as_slice()).unwrap_or(b""))
        )?;
    }
    writer.flush()?;

    if !history_path.is_empty() {
        let mut writer = BufWriter::new(File::create(&history_path)?);
        for history in histories.iter().filter(|h| !h.is_empty()) {
            writeln!(writer, "{}", serde_json::to_string(history).expect("serialise"))?;
        }
        writer.flush()?;
    }

    eprintln!(
        "adapter phase: {} reads trimmed, {} bases\n\
         overlap phase: {} pairs examined, {} trimmed, {} bases; \
         refused {} ambiguous (tandem-like), {} unlicensed",
        stats.reads_adapter_trimmed,
        stats.bases_adapter_trimmed,
        stats.overlap.pairs_examined,
        stats.overlap.pairs_trimmed,
        stats.overlap.bases_trimmed,
        stats.overlap.refused_ambiguous,
        stats.overlap.refused_unlicensed
    );
    Ok(())
}
