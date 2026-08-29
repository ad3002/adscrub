# adscrub

Adapter and read-through scrubbing for short-read pairs — the first stage of the
genemayers pipeline.

```
raw reads → adscrub → stirka|sverka → syncphase → aindexrs → clew|ariadna
```

Physical contamination removal, before anything algorithmic. What the sequencer
wrote that is not DNA has to leave here, because downstream every base is treated
as evidence: adapter k-mers become "real" high-frequency k-mers, and read-through
bases become minority variants that no k-mer-spectrum argument can recognise as
artefacts.

## The two phases, and why there are two

**Phase A — the adapter k-mer.** A bundled Illumina adapter 23-mer matched anywhere
in the read trims from there to the 3' end. Adapter k-mers are canonical, so a read
carrying the adapter on either strand is caught; adscrub runs before any strand
normalisation. A 23-mer is a specific synthetic sequence — the chance of a false hit
against a genome is ~98k / 4²³ ≈ 1.4e-12 — so the first hit is acted on without
further evidence.

Its blind spot is structural: read-through shorter than 23 bases carries no 23-mer,
and those first cycles are exactly where the artefact k-mers are made.

**Phase B — the mate overlap.** When the insert is shorter than the read, read 1 and
the reverse complement of read 2 align over the insert at a positive shift; the shift
*is* the fragment end, at single-base resolution, with no adapter sequence needed.
That covers phase A's blind spot.

## The hazard phase B has to be built around

Sequence similarity between mates is not evidence that the fragment ended.

In a tandem array of period `p` the two mates match at every shift `i`, `i+p`,
`i+2p`, … A pair lying **entirely inside an array**, with no read-through at all,
matches at a small shift by construction. Trimming on that cuts real array sequence
— consistently, from every read over the array. The result is a collapsed tandem
manufactured by the trimmer, in the sequence class an assembler cares about most.

The 2016 ancestor of this code took the smallest matching shift. That is conservative
when a read-through is real and simply wrong when it is not.

Two guards, and both are refusals rather than tie-breaks:

1. **The shift must be unique.** Every shift is scanned, not just the first. More
   than one match means the alignment does not identify a fragment end — the tandem
   signature. Refuse, and count the refusal.
2. **An adapter licenses the cut.** The pair must carry an adapter start somewhere:
   the evidence that read-through happened at all. The overlap then says *where* it
   started, earlier than the first detectable 23-mer, and that is what phase B adds.

The licence deliberately uses **short seeds** (`AGATCGGAAGAGC` and friends, 12–19
bases) rather than the 23-mer trim set. They answer an easier question — did this
pair read into an adapter — and a pair whose read-through is under 23 bases carries
no 23-mer while still carrying the first dozen. Measured: requiring a 23-mer as
licence refused 126,652 VAF12 pairs that had an unambiguous fragment end, blocking
the whole of phase B. A licence never cuts anything on its own; the overlap still has
to identify the end uniquely.

`--no-adapter-licence true` widens the catch to read-through with no detectable
adapter and gives up the tandem guarantee to whatever the uniqueness guard misses.
It is a flag, not a default.

## Measured

| | K-12, 95 b mates, long insert | VAF12, MiSeq 2×301, short insert |
|---|--:|--:|
| pairs whose mates overlap | 0.23 % | 94.7 % |
| **adapter phase** | 0 reads, **0 bases** | 10,030 reads, 860,179 bases |
| **overlap phase** | 0 pairs, **0 bases** | 19,210 pairs, 1,213,668 bases |
| refused, ambiguous (tandem-like) | 0 | 20 |
| refused, unlicensed | 40 of 2,000,000 | 107,442 |

**Inert on a long-insert library**, which is the acceptance criterion that matters:
a trimmer that fires there destroys real 3' sequence.

Downstream, on VAF12, scrubbing removes 13.5 % of the strand-locked minority k-mer
class and 32.6 % of the paired-allele class, and moves both toward where real
sequence sits in the read. It does not remove the class — read-through accounts for
97 % of it only among the extreme pairs whose fragment is shorter than a single
read, which are 0.75 % of that library.

## Reversibility

Every cut is recorded as a `Deletion` with its position, base and reason, one JSONL
record per read. `reverse_apply` restores the read. A trim that cannot be undone is
a deletion of evidence at the one stage that runs before anything could re-derive
it.

## Use

```bash
adscrub --reads pairs.seq.gz --out scrubbed.seq.gz --history scrub.jsonl
```

Input is the pipeline's `mate1~mate2` seq-per-line format, gzip optional, **as
sequenced** — each mate in the direction the sequencer read it. A file whose mates
have been reoriented cannot be scrubbed by phase B, because the overlap is what
places the fragment end.

| flag | default | |
|---|---|---|
| `--adapters` | bundled | one k-mer per line, k = 23 |
| `--min-overlap` | 30 | shorter overlaps cannot identify a fragment end |
| `--max-mismatch-rate` | 0.1 | sequencing error inside the insert |
| `--no-adapter-licence` | false | see the hazard section before setting it |

As a library: `adscrub::scrub`, or `AdapterSet` / `trim_adapter` / `overlap::trim_pair`
for finer control.

## Provenance

Phase A and the overlap geometry are ported from `stirka::adapter` and
`stirka::overlap`, which are themselves a modernisation of the 2016 `LU_trim` /
`V5_trim` C++ (`Compute_trims.cpp`). They lived in stirka while adscrub was a
scaffold. The uniqueness guard, the adapter licence and the short licence seeds are
new here; the tandem hazard they exist for was raised by @ad3002.

The read-through measurement behind phase B is in the research state as
`RES-refree-cross-library-error-classes` (research_compiler,
`research/2026_04_25_genemayers/`).

## Status

v0.1. Illumina paired short reads. The PacBio SMRTbell case in
`wiki/specifications/adscrub/` — long reads, chimera deferral, `ChimericDefer` — is
specified and **not** implemented here.
