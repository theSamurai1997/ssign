# Architecture decision

## Backbone PLM: ESM-2 650M (with ESM-C 300M as a head-to-head comparison)

Train every model variant twice: once with **ESM-2 650M** as the
backbone, once with **ESM-C 300M**. Pick the winner per task on the
held-out test set.

ESM-2 is the safer default (MIT license, two published bacterial-
effector classifiers already use it, ~3-4k cited downstream papers).
ESM-C 300M is a 2024 alternative that may or may not lift performance.
We don't know which wins on our task because no published benchmark
covers this. Training the second backbone costs ~3 GPU-h per model
variant: negligible.

ESM-3 is rejected (Cambrian Non-Commercial license blocks commercial
downstream use; structure modality is useless without crystal
structures; no independent benchmark shows it beating ESM-2 sequence-
only).

## Pooling strategy: test both

Also train every variant with two pooling strategies:

- **Whole-sequence mean-pool**: average ESM embedding across all amino
  acids.
- **N-terminal mean-pool**: average ESM embedding across only the
  first ~50 amino acids.

For T1/T3/T4 substrates the "secrete me" signal lives at the N-terminus
and N-terminal pool captures it more cleanly. For T2/T6 the signal is
more distributed. TXSelect tested both; we test both on our data too.

## What the model actually predicts (the important change)

**The model does NOT predict "is this protein an effector?".** It
predicts "is this protein a substrate of THIS specific secretion-system
instance in this genome?". This matters because:

- A genome can have multiple instances of the same SS type
  (P. aeruginosa PAO1 has 3 distinct T6SS clusters: H1, H2, H3, each
  with its own effectors).
- A protein can be a substrate of one T6SS but not another in the same
  genome.
- The user wants to know which substrate goes with which detected
  system, not just "is this an effector somewhere."

At inference, the loop is:

```
for genome in input_genomes:
    detected_systems = MacSyFinder(genome)           # e.g., [Dot/Icm, H1-T6SS, H2-T6SS]
    proteins = extract_proteins(genome)
    for system in detected_systems:
        for protein in proteins:
            features = (
                protein_features(protein),           # ESM + DLP + DSE + SignalP
                system_features(system),             # type, component count, location
                pair_features(protein, system),      # distance, co-operon, etc.
            )
            score = model.predict(features, system_type=system.type)
        report_top_K_by_score(system, proteins)
```

Output: per detected system, a ranked list of predicted substrates.

## Classifier architecture (revised for system-instance prediction)

```
Inputs per (protein, system) pair:
  Protein features:
    - ESM-2 (or ESM-C) embedding, pooled (whole or N-terminal)  [1280-d or 960-d]
    - DeepLocPro probabilities (5 classes)
    - SignalP probabilities (6 classes)
    - DeepSecE per-type probabilities (5)
    - PLM-Effector per-type score (5)

  System features:
    - System type (T1 / T2 / T3 / T4 / T6) — one-hot                [5]
    - Component count (number of MacSyFinder hits in this system)   [1]
    - System chromosomal span (start, end, length)                   [3]

  Pair features:
    - Distance from protein to nearest system component (in genes)  [1]
    - In same operon as system?                                      [1]
    - Strand match with system?                                      [1]

Head:
  - Tabular features (system + pair + scalar protein features) -> dense -> 64-d
  - Concatenate with pooled PLM embedding -> [PLM_dim + 64]
  - MLP: -> 256 -> 128 -> 5 (multi-task heads, one per SS type)
  - Sigmoid (multi-label per pair)

At training time:
  - For each (protein, system) pair in a labeled genome:
      target = 1 if protein is a known substrate of THIS system, else 0
  - One genome with one T4SS produces N_proteins training pairs
  - One genome with three T6SS produces 3 × N_proteins training pairs
```

The multi-task head conditions on the system type so the same model
serves all SS types. Only the head corresponding to the input system's
type contributes loss for that pair.

## Multi-task vs per-system models

Per Teo's request, benchmark both:

| Model variant | What |
|---|---|
| seq-only multi-task | PLM only, all 5 SS types, multi-task head |
| seq-only per-type | PLM only, one model per SS type (×5) |
| multimodal multi-task | PLM + tabular + pair features, multi-task |
| multimodal per-type | PLM + tabular + pair features, per-type (×5) |
| baseline (PLM-Effector) | published tool, no retraining |

12 variants (1+5+1+5 = 12, plus the baseline as 0-cost reference).
Each variant trained with 2 backbones (ESM-2, ESM-C) × 2 pooling
strategies (whole, N-terminal) = **48 trained models**.

Training one model: ~1 GPU-h on RTX 6000 (frozen backbone, small MLP
head). 48 models: ~48 GPU-h, ~2 days at 12 concurrent on CX3. Cheap
relative to the dataset build.

## Embedding cache

- ESM-2 650M, 1280-d, fp16, mean-pooled: ~2.56 KB / protein
- ESM-C 300M, 960-d, fp16, mean-pooled: ~1.92 KB / protein
- For 5000 genomes × 5000 proteins = 25M proteins:
  - ESM-2: ~64 GB
  - ESM-C: ~48 GB
- For both backbones × both pooling strategies: ~225 GB total.
- HPC scratch is fine. Don't ship cached embeddings; compute on demand
  per user genome at inference (~30-90 s GPU).

## Open ablations (after baselines work)

- Continue-pretrain ESM-2 on validation-genome proteins (1-2 epochs
  MLM, lr=1e-5) before freezing. Cheap, sometimes lifts AUPRC.
- With and without genomic-context pair features (isolated lift
  estimate).
- CLS-token pooling as a third strategy.
