## Why

The secretion-classifier model (scoped in the `reidmat/secretion-classifier` repo, docs 00-05) needs a labelled training set of (protein, secretion-system-instance) examples with per-protein tool features. The effector-recovery benchmark already produced the hard part, instance-assigned, citation-verified effectors with their tool signals and gene-distance to machinery, so the training set should be assembled from those tables rather than rebuilt. This change builds the dataset; model training is a later change. The label-side work is unblocked now; the feature side waits on the Phase 2 panel ssign runs (pilot queued on CX3).

## What Changes

- Audit the 347 `predicted`-evidence corpus rows to the same standard as the validated gold set (DOI resolves + UniProt/locus cross-check); keep all that survive, tagged with an `evidence_tier` so training weights validated over predicted (no hard cutoff).
- Assign predicted effectors to a specific system instance: auto-assign when the genome has exactly one system of that SS type; for multi-instance genomes, a literature-audit step reads the sourcing DOI to resolve the instance; unresolved rows remain instance-unknown type-level positives.
- Source T5SS effectors (absent from the corpus) via agents, by subtype, under the machinery-answer-key anti-hallucination contract (verbatim quote + resolvable DOI + real locus_tag): T5bSS as normal (protein, instance) TpsA substrates; T5aSS/T5cSS as self-secreted autotransporters (their own positive, `self_secreted=true`), usable as a MacSyFinder-detection confidence check.
- Assemble the labelled feature matrix: join per-protein tool signals (DeepLocPro / DeepSecE / SignalP / PLM-Effector + ESM embedding) from ssign run output, pair-features (gene-distance to the assigned system), and system-features, into positive (protein, instance) rows plus a PU-learning candidate/unlabeled set (EffectorP-style biological hard-negative filter, nnPU prior).
- Out of scope (recorded, not built here): the model training loop (ESM-2/ESM-C backbones, pooling, nnPU); and the T5aSS-neighbor DLP/DSE observation, which stays an exploratory side-study in the benchmark (task 6b), not a training label.

## Capabilities

### New Capabilities
- `predicted-effector-audit`: verify the predicted corpus rows (citation + identity) and assign an evidence tier so they can be added as lower-weight positives.
- `predicted-instance-assignment`: resolve each predicted effector to a specific system instance (unique-genome auto-assign, else literature audit, else instance-unknown).
- `t5ss-effector-sourcing`: agent-sourced T5SS effector set by subtype with the anti-hallucination contract and a `self_secreted` flag.
- `training-feature-matrix`: assemble labelled (protein, instance) feature rows + the PU candidate/unlabeled set from ssign run output and the benchmark gene-order/machinery tables.

### Modified Capabilities
<!-- none: this is a new analysis project; it reads the benchmark outputs but changes no existing spec -->

## Impact

- New analysis + data under `validation_sweeps/benchmark/` (or a sibling dataset dir); reuses `scripts/bench_index.py` and `scripts/bench_runout.py`.
- Reads the effector-recovery benchmark outputs (gold set, instances, machinery answer key, gene-order index, Phase 2 `results_raw`/`actual_per_effector`). Depends on the Phase 2 panel runs for the feature side.
- Consumes external literature via agents (PubMed/Crossref) for the predicted-instance audit and T5SS sourcing; no ssign pipeline code changes.
- Produces the training dataset the future `secretion-classifier` model change will consume; model code may live in the `reidmat/secretion-classifier` repo.
