# adscrub — engine kickoff instructions

You (the AI engine) are implementing **adscrub v0.1**, a standalone tool for PacBio adapter detection and trimming.

**Status:** spec is complete and reviewed; **implementation is currently DEFERRED** per project owner direction (v0.1 ARIADNA prototype proceeds without adapter scrubbing initially — adapter contamination ~0.01-0.5% is non-blocking for first proof-of-concept). Build this when assembly pipeline maturity demands it (v0.5+) or when user explicitly requests.

When you DO start, follow this guide.

## Read these first (in order)

1. **`../CLAUDE.md`** — umbrella guide for all tools in this folder. Gene Myers style, dependency policy, workflow, what NOT to do. **READ THIS FIRST.**
2. **`../../wiki/specifications/adscrub/README.md`** — directory navigation for adscrub specs.
3. **`../../wiki/specifications/adscrub/versions/v0.1.md`** — version spec (CLI, I/O, output formats, build artifacts).
4. **`../../wiki/specifications/adscrub/components/spec_adapter_search.md`** — the central algorithm.
5. **`../../wiki/specifications/glossary.md`** — shared vocabulary.
6. **`../../wiki/specifications/spec_format.md`** + **`../../wiki/specifications/spec_process.md`** — meta-specs about how specs work + TEST-FIRST + mutation-test discipline.

## What you are building

A Rust **library + binary** (per user direction "это библиотека + быстрый поиск в стиле майерса в много потоков"):

- **Library:** `adscrub::scan_read`, `adscrub::classify`, `adscrub::run` — for in-process invocation by other tools.
- **Binary:** `adscrub` CLI — standalone use.

What it does:
1. Takes long reads (FASTA/FASTQ, gzip-optional) + adapter database (default embedded PacBio SMRTbell, user-overridable).
2. Multi-threaded k-mer-based search of adapter sequences against each read.
3. Classifies per read: `Clean | TerminalTrim | ChimericDefer`.
4. Emits cleaned reads (working stream) + chimeric reads (deferred bucket) + per-read JSONL history.

## Pipeline position

```
raw HiFi → adscrub → stirka → syncphase → ariadna
```

**adscrub is the FIRST stage** — physical contamination removal before any algorithmic processing.

## Implementation kickoff

### Step 1 — read all the specs above

Don't start coding before you've read every doc in the "Read these first" list.

### Step 2 — set up Cargo project (lib + bin)

```toml
# Cargo.toml
[package]
name = "adscrub"
version = "0.1.0"
edition = "2021"
authors = ["..."]
license = "MIT OR Apache-2.0"

[lib]
name = "adscrub"
path = "src/lib.rs"

[[bin]]
name = "adscrub"
path = "src/main.rs"

[dependencies]
clap = { version = "4", features = ["derive"] }
noodles-fastq = "..."
flate2 = "..."
serde = { version = "1", features = ["derive"] }
serde_json = "..."
hashbrown = "..."
rayon = "..."
log = "..."
env_logger = "..."

[dev-dependencies]
proptest = "..."
criterion = "..."
```

Use latest stable Rust edition. Pin specific dep versions in Cargo.lock at first commit.

### Step 3 — set up Makefile

```makefile
.PHONY: build test bench install clean ecoli-roundtrip synthetic-tests mutants ci

build:
	cargo build --release

test:
	cargo test

bench:
	cargo bench

install:
	cargo install --path .

clean:
	cargo clean

ecoli-roundtrip:
	# Run on E.coli fixture; verify acceptance criteria from spec
	# Fixture path: /mnt/data/claude/2026-02-03_ariadna/data/samples/ecoli_e2348_69_30x_hifi.fastq.gz
	# md5: 99e3702405d27db8aa0ca21765de49fd
	./target/release/adscrub -o /tmp/adscrub_ecoli /mnt/data/claude/2026-02-03_ariadna/data/samples/ecoli_e2348_69_30x_hifi.fastq.gz
	# Verify: clean_passthrough_fraction >= 0.99; chimeric_deferred_fraction <= 0.01

synthetic-tests:
	cargo test --test integration_synthetic_terminal
	cargo test --test integration_synthetic_chimera
	cargo test --test integration_passthrough
	cargo test --test integration_reversibility

mutants:
	cargo install cargo-mutants 2>/dev/null || true
	cargo mutants

ci: build test
```

### Step 4 — set up source layout

Per spec `versions/v0.1.md` §"Build artifacts":

```
src/
├── lib.rs                           # public library API
├── main.rs                          # CLI entry, depends on lib
├── adapter_db.rs                    # adapter database loader (file + embedded defaults)
├── adapter_index.rs                 # k-mer index over adapter sequences
├── search.rs                        # multi-threaded per-read scan
├── classify.rs                      # rule-based no-adapter/terminal-trim/chimera classification
├── history.rs                       # JSONL output + reversibility for terminal trims
└── output.rs                        # write all output files (clean.fa.gz, chimeric.fa.gz, ...)

data/
└── adapters_pacbio_smrtbell.fa     # embedded default adapter DB (committed; ~5 sequences)

tests/
├── integration_ecoli.rs
├── integration_synthetic_terminal.rs
├── integration_synthetic_chimera.rs
├── integration_passthrough.rs
├── integration_reversibility.rs
├── property_determinism.rs
├── gen_synthetic_adapters.py        # test fixture generator
└── fixtures/
    ├── synthetic_terminal_adapters.fa
    ├── synthetic_chimeras.fa
    ├── clean_synthetic_reads.fa
    └── custom_adapters.fa

benches/
└── ecoli_benchmark.rs
```

### Step 5 — TEST-FIRST

Per `../CLAUDE.md` §"Step 6 — build the test suite first":

1. Generate test fixtures (or use placeholders that fail). Engine creates real fixtures via `tests/gen_synthetic_adapters.py`.
2. Write test functions for every Style C concrete test case in the algorithm spec (initially `#[ignore]` or panicking).
3. Write property tests:
   - Determinism (P1)
   - Conservation (P2)
   - Reversibility for terminal trims (P3)
4. THEN implement modules.

Test names match spec test_id field (e.g., `aas_synthetic_left_terminal`, `aas_synthetic_internal_chimera`, `adscrub_ecoli_basic`, etc.).

### Step 6 — implement in dependency order

Suggested order:
1. **`adapter_db.rs`** — load adapter sequences from FASTA file OR embedded default. Tests: round-trip load.
2. **`adapter_index.rs`** — build hashbrown::HashMap<u64, AdapterMetadata> from adapter sequences (forward + revcomp). Tests: index size, lookup correctness.
3. **`search.rs`** — per-read scan: rolling Wang hash, lookup in index, collect hits. Tests: synthetic input with known adapter at known position; verify hits returned.
4. **`classify.rs`** — cluster_hits + classification rules (Clean / TerminalTrim / ChimericDefer). Tests: synthetic terminal vs chimera fixtures.
5. **`history.rs`** — JSONL serialization + reverse_apply for reversibility. Tests: terminal trim → serialize → reconstruct → match original.
6. **`output.rs`** — orchestrate emission of all output files (clean.fa.gz, chimeric.fa.gz, history.jsonl.gz, summary.tsv).
7. **`lib.rs`** — public API: `scan_read`, `classify`, `run`. Tests: end-to-end small synthetic.
8. **`main.rs`** — clap CLI parsing → invoke lib::run.

### Step 7 — mutation testing

Run `cargo mutants` and verify each spec'd mutation (M1-M8) is caught by at least one test. Per spec coverage table in `components/spec_adapter_search.md` §"Mutation testing expectations".

### Step 8 — acceptance run

```bash
make build
make test
make ecoli-roundtrip
make mutants
```

If all pass and acceptance criteria met, v0.1 is `accepted`-ready. Update spec status in wiki from `review` to `accepted`.

## Open questions to resolve

Per spec `versions/v0.1.md` and `components/spec_adapter_search.md` §"Open questions":

| Q | Spec | Resolution approach |
|---|---|---|
| -k default | versions/v0.1 Q1 | Engine: 16 (per spec default); tune empirically |
| --min-adapter-hits default | versions/v0.1 Q2 | Engine: 5 (per spec default); tune |
| Library API surface | versions/v0.1 Q3 | Engine: per ENGINE_CHOICE in spec — `scan_read`, `classify`, `run` |
| Adapter DB defaults | versions/v0.1 Q5 | Engine: 3-5 PacBio SMRTbell variants from canonical PacBio docs |
| Reverse-complement strand | components Q1 | Engine: include both forward + revcomp in index (recommended) |
| Mutation tolerance | components Q2 | Engine: strict exact match (per spec default) |

Pre-resolved (per user direction in spec):
- ✅ Library + binary structure
- ✅ Multi-threaded via rayon
- ✅ Defer chimeric reads whole (no split)
- ✅ JSONL history with reversibility for terminal trims
- ✅ Embed PacBio SMRTbell adapters by default

## Reference implementations to consult

- **HiFiAdapterFilt** at https://github.com/sheinasim/HiFiAdapterFilt — reference for SMRTbell adapter sequences + workflow logic
- **PacBio lima** at https://github.com/PacificBiosciences/barcoding — vendor official adapter trimmer (heavy stack; reference only)
- **MBG** at `../../external/MBG/src/FastHasher.cpp` — Wang invertible hash (~30 LOC vendored verbatim)
- **rukki** at `../../external/rukki/src/` — NHGRI Rust style precedent (clap, hashbrown, minimal deps)

When uncertain about algorithm choice:

- **`../../wiki/projects/DASCRUBBER.md`** — Gene's predecessor read-scrubbing tool (architectural ancestor)
- **`../../wiki/projects/DEXTRACTOR.md`** — Gene's PacBio raw-data extractor
- **`../../wiki/algorithms/heng_li_optimizations.md`** §3 — Wang invertible hash details

## Coding style reminders

Per `../CLAUDE.md`:
- Single-purpose binary + library.
- Minimal dependencies.
- Result + ? error handling; no silent swallow.
- Determinism (fixed seed; sorted reductions if parallel).
- TEST-FIRST.
- No comments explaining WHAT.

## When to ask the user (project owner @ad3002)

1. First check spec, wiki, reference repos.
2. If still unclear: leave a `// QUESTION:` TODO + use a documented default.
3. Surface in status report.

## Summary

Standalone Rust library + binary for PacBio adapter detection. Multi-threaded k-mer search. Defers chimeric reads (whole) instead of discarding. Reversibility for terminal trims via JSONL history. First stage of ARIADNA pipeline.

**Status: implementation deferred until needed.** Spec is engine-consumable when the time comes.

When in doubt, read what Gene Myers wrote.
