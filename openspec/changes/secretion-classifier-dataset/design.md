## Context

The `reidmat/secretion-classifier` repo (docs 00-05) designs a multimodal classifier that predicts, per detected system instance, which proteins are its substrates, using a frozen PLM embedding plus tabular features from MacSyFinder/DeepLocPro/DeepSecE/SignalP and PU learning (nnPU) over a biologically-filtered unlabeled set. The effector-recovery benchmark (`validation_sweeps/benchmark/`) already produced the instance-assigned, citation-verified positive effectors with their per-tool signals and gene-distance to machinery. This change assembles the training dataset from those tables. It does not train the model.

The benchmark provides the reusable substrate: `data/phase1/effector_gold_set_phase1.tsv` (582 validated, instance-tagged), `data/machinery/instances.tsv` + `machinery_answer_key.tsv` (per-instance machinery loci), `data/phase1/{gene_order_index,ceiling_per_effector}.tsv` (gene-distance), and Phase 2 `actual_per_effector.<tag>.tsv` + `results_raw` (per-protein tool signals). Shared libs `scripts/bench_index.py` (drift-tolerant matching) and `scripts/bench_runout.py` (run-output reader, sequence/locus bridge) are reused.

## Goals / Non-Goals

**Goals:**
- Add the predicted corpus rows as verified, evidence-tiered, instance-resolved-where-possible positives.
- Source T5SS (absent from the corpus) by subtype with the same anti-hallucination contract used for the machinery answer key.
- Emit one labelled feature row per (protein, instance) plus a PU candidate/unlabeled set ready for nnPU training.
- Keep the label side fully runnable now; cleanly gate the feature side on the Phase 2 panel runs.

**Non-Goals:**
- Model training (backbones, pooling, nnPU loss) — a later change.
- The T5aSS-neighbor DLP/DSE statistical observation — stays exploratory in benchmark task 6b, never a training label.
- Re-deriving machinery via MacSyFinder for labels (circular; literature-only, as in the benchmark).

## Decisions

- **Reuse benchmark tables as the positive backbone, not a fresh build.** The instance assignment + citation verification are the expensive parts and are already done for the 582 validated rows. Alternative (rebuild from PLM-E/DeepSecE/SecReT4/6 FASTAs per data-audit doc 02) is deferred; those anonymized sets lack instance and genome context, which the benchmark already resolved.
- **Predicted rows: verify, keep, weight, don't cut.** Per Teo: train on validated + predicted with heavier weight on validated, rather than a hard high-confidence cutoff. `evidence_tier` is a column consumed by the training loss later; the audit only removes rows with broken citations / mismatched identifiers (the same bar as the gold set).
- **Instance assignment: auto then literature, never guess.** Single-same-type-instance genomes auto-assign (Phase 1 precedent). Multi-instance genomes get a literature-audit agent that reads the sourcing DOI; resolved rows carry a verbatim quote. Unresolved rows are kept as instance-unknown type-level positives (protein features only, pair-features null) rather than dropped or guessed.
- **T5SS by subtype with distinct label conventions.** T5bSS TpsA = normal (protein, instance) substrate; T5aSS/T5cSS = self-secreted autotransporter as its own positive (`self_secreted=true`), used to score MacSyFinder's T5a/c detection confidence. Agents seed on PF03797 (AT β-domain), PF03895 (YadA/trimeric), and TpsA/TpsB(POTRA) families; output the gold-set schema so rows drop straight in.
- **Feature side gated on runs.** Label-side artifacts (audited positives, instances, T5SS set) are produced now. The feature matrix joins `results_raw` signals via `bench_runout` once the Phase 2 panel runs land; genomes without output are reported as pending-run, never emitted as complete.
- **PU unlabeled set via the EffectorP-style filter (doc 04).** Candidates = non-effector proteins with an extracellular/secretion signal (DLP/DSE/SignalP/TM), minus those highly similar to a positive; labelled unlabeled, not negative.

## Risks / Trade-offs

- **Predicted-instance literature audit is low-yield** → many multi-instance papers won't name the specific system; mitigate by keeping unresolved rows as type-level positives instead of losing them.
- **T5SS sourcing has no curated DB** → rely on agents + Pfam seeds + literature; mitigate hallucination with the verbatim-quote + resolvable-DOI + real-locus_tag contract and an independent verification pass (mirror benchmark task 3.6).
- **T5aSS/T5cSS self-secreted positives are weaker labels** → they teach "is this autotransporter really secreted," not substrate selection; keep them flagged so the model/eval can down-weight or hold them out.
- **Feature side blocked on CX3 runs** → the pilot is queued; mitigate by finishing all label-side tasks first so only the join remains when runs return.

## Open Questions

- Exact training weight ratio validated:predicted — deferred to the model-training change (this change only records `evidence_tier`).
- Whether instance-unknown type-level positives are used in the instance model or only a separate type-level head — deferred to model design.
- T5SS subtype coverage depth (T5dSS/T5eSS are rare) — sourcing agents report what exists; depth decided after the first pass.
