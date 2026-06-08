# Negative-set construction and PU learning

Captured from the background research agent on 2026-06-08. Concrete
recipe for handling positive-unlabeled training data with genome-context
negatives.

## Recipe (pseudocode)

```python
# ---- Pre-train (optional but high-value) ----
# Continue-pretrain ESM-2 on the unlabeled proteins from the ~10
# validation genomes (masked language modeling, 1-2 epochs, lr=1e-5).
# Frozen for downstream.

# ---- Per-fold training set (leave-one-genome-out) ----
positives_P    = known_effectors                          # ~1-2k
candidates_U   = [p for p in all_other_proteins           # ~30k -> ~5k
                  if (DeepLocPro.extracellular_p(p) > 0.3
                      or DeepSecE.score(p)         > 0.3
                      or SignalP.has_sp(p)
                      or p.has_TM_helix)]
# Homology filter: drop U items >= 40% identity to any P
candidates_U   = mmseqs_filter(candidates_U, positives_P,
                               min_id=0.40, action="drop_from_U")

# ---- Estimate class prior pi on U ----
pi_hat = tice(positives_P, candidates_U)   # TIcE (fast) or KM2 (more accurate, smaller U)

# ---- Loss: nnPU + focal weighting on negative-risk term ----
loss = nnPU(pi=pi_hat,
            base_loss=focal(gamma=2.0, alpha=pi_hat),
            beta=0.0, gamma_nn=1.0)

# ---- Batching: stratified + online hard-negative mining ----
# Per minibatch: B/2 positives, B/2 from U.
# Of the U half: 25% uniform, 75% top-loss from a 4x overdraw (OHEM).
# Warmup: first 2 epochs uniform U sampling only (curriculum).

# ---- Decision threshold ----
# Tune on held-out genome by maximizing F1 vs reliable-negatives
# (reliable-neg = U items with score in bottom decile after epoch 1, "spy" method).
```

## Why these choices

**nnPU loss** (Kiryo et al. 2017): the deep-learning standard. The
unbiased Elkan-Noto estimator overfits with neural nets; nnPU clips the
negative-risk term at zero. PU-GO (Oxford Bioinformatics 2024) confirms
nnPU works directly on protein-function prediction.

**TIcE for class-prior estimation**: decision-tree induction, fast,
second only to KM2 in accuracy but stable. KM2 (Ramaswamy 2016) is more
accurate but doesn't scale beyond ~10k samples. For ~5k filtered U,
TIcE is the sweet spot. DEDPUL (Ivanov 2020) outperforms both if we can
afford a density estimator. PULSNAR (PeerJ 2024) handles SCAR-violation
if we suspect biased positive sampling.

**Focal loss with `gamma=2.0, alpha=pi_hat`**: canonical default for
imbalanced classification. Inverse-frequency CE alone is weaker than
focal on biology data. SMOTE in PLM-embedding space mostly fails past
~256 dims (nearest-neighbor interpolation breaks); skip it, let nnPU +
focal handle imbalance.

**Biological hard-negative filter**: the EffectorP 3.0 recipe (restrict
negatives to "secreted non-effectors"). Drop U items with predicted
cytoplasmic localization, no signal peptide, no TM helix, no DSE
signal. Cuts ~30k raw candidates to ~5k informative ones. Cite
EffectorP 3.0 as precedent for SS effectors.

**OHEM (online hard-example mining)**: standard in detection;
demonstrated for PLM embeddings in arXiv:2405.17902 (2024). With a
warmup of uniform sampling, then escalating to top-loss negatives,
training stays stable while focusing capacity on hard cases.

## Alternative framings (rejected for now)

- **MIL (multi-instance learning)**: wrong fit. MIL assumes we don't
  know which protein is positive; we do. Useful only with unannotated
  genomes where bag label = "has Dot/Icm system."
- **One-class classification**: trains on positives only. Underperforms
  PU when U is informative (PMC10326160 review). Skip.
- **Self-supervised pre-training on unlabeled genome proteins**: yes,
  included in the recipe above as an optional first step.
- **Active learning**: pilot only after the first model converges.
  PULSNAR + active labeling pairs naturally.

## Reference implementations

- nnPU in PyTorch: https://github.com/kiryor/nnPUlearning
- nnPU + variants: https://github.com/cimeister/pu-learning
- sklearn wrappers: https://github.com/pulearn/pulearn
- VPU (variational PU): https://github.com/HC-Feynman/vpu

## Key citations

- Kiryo et al. 2017, "Positive-Unlabeled Learning with Non-Negative
  Risk Estimator." arXiv:1703.00593
- Lin et al. 2017, "Focal Loss for Dense Object Detection."
- EffectorP 3.0, bioRxiv 2021, doi:10.1101/2021.07.28.454080
- Effectidor, Oxford Bioinformatics 2022, doi:10.1093/bioinformatics/btac162
- PU-GO, Oxford Bioinformatics 2024, doi:10.1093/bioinformatics/btae404
- ProtTucker (contrastive PLM), PMC9188115
