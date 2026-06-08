# Overview

## Motivation

ssign's proximity-based substrate filtering (±N genes around each SS
component) is conservative by design: it keeps the substrate set small
and high-confidence at the cost of missing substrates that are encoded
distantly from their secretion machinery. We don't know how big that
miss is. No published bacterial tool quantifies recall delta between
proximity-based and whole-genome substrate prediction.

If the gap is small (<5% of true substrates miss the proximity window),
ssign as-is is defensible and this project doesn't have legs.

If the gap is large (>20%), a meaningful fraction of real biology is
hidden, and a sequence-aware classifier that uses ssign's tool outputs
as features could recover those substrates without the per-protein
network calls or false-positive blow-up that comes from running
DeepSecE / PLM-Effector genome-wide.

## Decision tree

```
Step 1: Quantify the gap (see 01_gap_validation.md)
    ├── Gap <5%      → stop. ssign as-is. Document the result.
    ├── Gap 5-20%    → marginal. Decide based on which substrates are missed
    │                   (T6SS effectors? niche pathogen-specific?)
    └── Gap >20%     → proceed to Step 2.

Step 2: Proof-of-concept multimodal classifier
    ├── Train on existing DeepSecE/PLM-E data (no new ssign runs)
    ├── Sequence-only baseline vs multimodal (with tool-output features)
    ├── If multimodal lift <2% AUROC over PLM-E head-to-head → stop
    └── If multimodal lift ≥2% AUROC → proceed to Step 3

Step 3: Full dataset (~5000 reference bacterial genomes)
    ├── Run ssign-lite (Bakta + MacSyFinder + DLP + DSE + SignalP + PLM-E)
    ├── ~1700 GPU-h, ~1 week wallclock at CX3 current throughput
    └── Train final model(s) + benchmark
```

## How it works (plain English)

The model's job is: given a whole bacterial genome, and given a
specific secretion system that MacSyFinder has detected in that
genome, predict which proteins in the genome are substrates of THAT
system.

This is different from "is this protein an effector somewhere." A
genome can have several secretion systems (e.g., P. aeruginosa has
three distinct T6SS clusters), and different proteins are substrates
of different systems. The model has to match protein to system, not
just label proteins.

Step by step at inference time:

1. User gives ssign a genome.
2. ssign runs the existing pipeline (Bakta, MacSyFinder, DeepLocPro,
   SignalP, DeepSecE, PLM-Effector) on the whole genome. This produces
   per-protein features: localization probabilities, signal-peptide
   class, per-type effector scores.
3. For each MacSyFinder-detected secretion system in the genome:
   - For each protein in the genome:
     - Build a feature vector that combines protein features +
       system features + pair features (e.g., distance from this
       protein to this system on the chromosome).
     - The model scores this (protein, system) pair: how likely is
       this protein a substrate of this specific system?
   - Sort proteins by score, return the top-ranked as predicted
     substrates of that system.
4. The user gets a per-system list of predicted substrates.

This replaces (or augments) the current proximity-window filter. The
proximity filter is essentially saying "if a protein is within ±3
genes of a system component, score it 1, else 0." The model learns a
smarter version of that score using all the available features.

Why this should help: proximity is one signal among many. Some real
substrates encode distantly from their machinery (especially T4
effectors in Legionella, ~330 distributed across the genome). A model
that knows "DLP says extracellular AND DSE says T4 AND it has a Sec
signal peptide AND it sits in a known SS operon" can pick up those
distant substrates without false-positive blow-up.

## Architecture (technical summary)

Frozen PLM backbone (ESM-2 650M and ESM-C 300M, both tested) producing
per-protein embedding, with two pooling strategies (whole-sequence
mean and N-terminal mean) both tested. Concatenated with tabular
features per (protein, system) pair:

- Protein features: DLP probabilities, SignalP probabilities, DSE
  per-type scores, PLM-Effector per-type scores
- System features: SS type, component count, chromosomal span
- Pair features: distance to system, in-operon, strand match

Multi-task head over T1 / T2 / T3 / T4 / T6, conditioned on the input
system's type. See `03_architecture.md` for full detail.

## Open questions (pending agent reports)

1. **PLM backbone**: ESM-2-650M (MIT, mature) vs ESM-3 (Cambrian
   non-commercial, multimodal) vs ESM-C (Cambrian, ESM-2 successor).
   Decision criterion: ablation lift over ESM-2 baseline; license
   compatibility with shipping inside ssign.
2. **Negative-set construction**: PU loss vs weighted CE vs focal;
   genome-context negatives vs DeepSecE curated negatives; hard-negative
   mining criterion.
3. **Data sourcing pipeline**: MMseqs2 vs CD-HIT dedup threshold;
   GraphPart vs leave-one-genus-out for held-out splits; re-BLAST
   strategy for anonymized PLM-E / DeepSecE FASTAs.

## Out of scope (for now)

- T5SS / T7SS substrate prediction. No curated training data exists for
  T5SS effectors as a class (architecture HMMs only). T7SS is Gram+,
  outside ssign's scope.
- ESMFold structure prediction as an extra feature. Adds a heavy
  compute step; revisit if ESM-3 with structure shows a clear win.
- Active learning / hard-negative mining beyond a static recipe.

## What this is NOT

- Not a replacement for ssign's proximity filtering. If it ships, it's
  an ssign module that runs alongside (or instead of) proximity, with
  the user choosing the mode.
- Not a publication commitment. Each decision-tree branch can exit
  early. We only commit to the next step after the previous one
  produces evidence.
