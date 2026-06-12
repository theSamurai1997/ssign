# Secretion-classifier training dataset

Labelled positive set for the multimodal secretion-system substrate classifier
(`reidmat/secretion-classifier`, docs 00-05). This directory holds the **label side** of
the dataset: verified, instance-assigned, evidence-tiered positive effectors. The feature
side (per-protein tool signals + ESM embedding + PU unlabeled set) is added in group 4 and
is gated on the Phase-2 ssign panel runs (see *Run dependency* below).

**Deliverable:** `positives_all.tsv` — 925 positive (protein, system-instance) examples (930 from the
build chain, minus 5 dropped by the group-6 citation-audit overlay; see "Citation-audit overlay" below)
across all five secreted-effector secretion-system types, built from the effector-recovery
benchmark tables.

```
            validated   predicted   total
  T1SS            28         177       205
  T2SS            83           7        90
  T3SS           237         116       353
  T4SS           111          17       128
  T5SS            23           0        23
  T6SS           123           8       131
  ----------------------------------------
  all            605         325       930
```

## How it is built

The benchmark already did the expensive parts (instance assignment + citation verification)
for the 582 validated effectors. This dataset reuses those tables as the backbone and adds
two things: the audited predicted corpus, and sourced T5SS. Build order, each script reading
the previous outputs (run with the repo `.venv`, from `scripts/`):

| Script | Produces | Role |
|---|---|---|
| `30_extract_predicted.py` | `30_predicted_*.tsv` | pull the 347 predicted corpus rows, apply the gold-set biology/status filters |
| `31_verify_predicted.py` | `predicted_audited.tsv` | hold each to the gold-set bar (DOI resolves + UniProt/locus cross-check); 325 survive |
| `32_tiered_positives.py` | `positives_tiered.tsv` | union validated (582) + predicted (325) with `evidence_tier` |
| `33_assign_predicted_instances.py` | `predicted_instanced.tsv` | auto-assign single-instance genomes; flag multi-instance as ambiguous |
| `34_fold_recitation.py` | `doi_recite.jsonl` | fold the broken-citation re-sourcing (verbatim-quote + Crossref contract) |
| `35_fold_instances.py` | `positives_instanced.tsv` | apply the literature-audit instance resolutions |
| `36_resolve_verify_t5ss.py` | `t5ss_effectors.tsv` | resolve + verify agent-sourced T5SS to RefSeq loci |
| `37_fold_t5ss.py` | **`positives_all.tsv`** | final union: validated + predicted-instanced + T5SS, one canonical instance column (930 rows) |
| `40_pair_features.py` | `pair_features.tsv` | gene-distance + reachability + system-features per positive (task 4.2; no runs needed) |
| `41_citation_consistency.py` | `citation_consistency_found.tsv` | group-6 pass-1: deterministic CrossRef check of the ssign-found effectors' sourcing DOIs |
| `42_pass2_verify.py` | `pass2_results.tsv` | group-6 pass-2: merge the agent re-audit + re-verify every returned DOI |
| `43_apply_citation_corrections.py` | **`positives_all.tsv`** (overlay) | apply the audit: fix DOIs, relabel BopA/BopE T3SS, drop 5 unsupported/duplicate → **925 rows** |

Shared helpers: `scripts/bench_io.py` (tsv read/write, DOI normalize), `scripts/bench_index.py`
(drift-tolerant genome/locus matching). Re-running `37_fold_t5ss.py` rebuilds the deliverable
from the three branch tables; the full chain reproduces from the benchmark Phase-1 gold set
plus the source corpus.

**Citation-audit overlay (group 6).** `43_apply_citation_corrections.py` is a post-build overlay:
it edits `positives_all.tsv` in place (930 → **925**) after `37` produces it, so the canonical
training table is the corrected 925-row version. **Run order matters:** a rebuild via `37` reverts
to the uncorrected 930, so `43` must be re-applied after any `37` re-run. The pristine pre-overlay
table is snapshotted to `positives_all.pre_citation_audit.tsv`; removed rows + reasons in
`positives_removed_citation.tsv`; every field change in `citation_corrections_log.tsv`. The overlay
changed 2 SS-type labels (BopA/BopE T6SS→T3SS) and dropped 5 rows, so `40_pair_features.py` should be
re-run before the feature matrix is built.

## Label conventions (load-bearing)

These columns are what the training loss and evaluation read. Get them right.

**`evidence_tier`** — `validated` | `predicted`. Validated = the benchmark Phase-1 gold set
(582) plus the 23 verified T5SS rows. Predicted = corpus rows that passed the same audit bar
(DOI resolves + identifiers cross-check) but come from prediction-derived evidence. Train on
both; weight validated heavier. The exact weight ratio is deferred to the model-training
change — this dataset only records the tier, it does not cut predicted rows.

**`instance_source`** — how the row's system instance was assigned:

| value | meaning | `type_level` |
|---|---|---|
| `gold` | validated effector, instance from Phase-1 curation | `no` if an instance was assigned, else `yes` |
| `auto` | predicted, genome had exactly one system of that SS type | `no` |
| `literature` | predicted, multi-instance genome resolved by a literature-audit agent (verbatim quote in `instance_quote`) | `no` |
| `none` | no specific instance assigned (unresolved predicted, or a type-level T5b substrate) | `yes` |
| `self` | T5a/c/d/e autotransporter — its own system, no separate instance | `no` |

Counts: gold 582, none 214, auto 95, literature 27, self 12.

**`type_level`** — `yes` when no specific instance is assigned (669 `no`, 261 `yes`). Type-level
positives carry protein features only; their pair-features (gene-distance to machinery) are
null because there is no assigned instance to measure against. Whether the instance model uses
them or only a separate type-level head is a model-design decision (deferred).

**`sys_instance_id`** — the single canonical benchmark instance id (e.g. `T1SS_R06`, `T3SS_20`),
unified across tiers and strictly canonical-or-blank: validated rows take it from
`ceiling_per_effector.instance_id` (the benchmark's authoritative assignment, not the gold set's
messy literature label), predicted rows from the group-2 instance assignment. Blank for
type-level rows. Group 4's feature join keys on this column; 657 rows carry an instance.

**`self_secreted`** / **`subtype`** — T5SS only. `self_secreted=true` for T5a/c/d/e (the
autotransporter passenger secretes itself, the row is a confidence label on MacSyFinder's
autotransporter call, *not* a substrate-selection example); `false` for T5b, whose TpsA is a
genuine substrate of its separate TpsB pore. `subtype` ∈ {T5aSS(7), T5bSS(11), T5cSS(3),
T5dSS(1), T5eSS(1)}. Keep the self-secreted rows flagged so the model/eval can down-weight or
hold them out. See `t5ss_sourcing_brief.md` for the full subtype rationale.

**`citation_status`** — `RESOLVED` (323) | `UNRESOLVED` (25) for predicted + T5SS rows; blank
for validated gold rows (the gold set's verification lives in `verification_status`). The 25
UNRESOLVED are predicted rows whose primary citation did not resolve even after the re-citation
pass; they are kept (the protein/locus identity is verified) but flagged.

## Schema (46 columns)

- **Identity** (1-6): `gene`, `uniprot`, `locus_tag`, `organism`, `refseq_genome`, `ss_type`.
- **Instance** (7, 30-34): `sys_instance_id`, `instance_source`, `type_level`,
  `instance_candidates`, `instance_quote`, `instance_source_doi`.
- **Evidence & citation** (8-14, 29): `evidence_level`, `evidence_tier`, `primary_ref`,
  `verification_status`, `verification_notes`, `citation_status`, `audit_tier`, `uniprot_status`.
- **Benchmark provenance** (15-28): `family`, `length`, `proteome`, `source`, `gate`,
  `ceiling_source`, `placement_tier`, `species_match`, `placement_start/stop/strand`,
  `placement_effector_locus`, `testable`, `testable_reason`. Carried from the benchmark gold-set
  build; populated for validated/predicted, blank for T5SS.
- **T5SS sourcing** (35-46): `subtype`, `self_secreted`, `contig`, `start`, `stop`, `strand`,
  `quote`, `uniprot_note`, `locus_match`, `doi_resolves`, `gene_in_quote`, `verified`.
  T5SS rows only; blank for the other tiers.

## PU (positive-unlabeled) set semantics

`positives_all.tsv` is the **P** set. The classifier trains with nnPU over P plus an
**unlabeled** set U (doc 04), built in group 4: non-effector proteins that pass an
EffectorP-style secretion filter (extracellular/secretion signal from DLP/DSE/SignalP/TM),
minus proteins highly similar to a positive. U is labelled *unlabeled, not negative* — some
unknown fraction are real effectors not yet in the literature, which is exactly what nnPU
handles. Building U needs the per-protein tool signals, so it is gated on the runs below.

## Pair-features (group 4, task 4.2 — done now, no runs needed)

`40_pair_features.py` → `pair_features.tsv`, one row per positive (keyed on `protein_id`). The
gene-distance to the assigned instance's nearest machinery locus plus the proximity-window flags
are a function of the labels + the benchmark's own gene-order index and machinery answer key, so
this is the one feature task that does **not** wait on the ssign runs. Columns: `nearest_dist`,
`nearest_tier`, `nearest_locus`, `reachable_n3/n5/n7`, `n_machinery`, `pair_source`.

`pair_source` records provenance: `ceiling` (456, validated, distance precomputed by the
benchmark) | `computed` (91, predicted, distance computed here from the gene-order index) |
`none` (261, type-level, no instance) | `self` (12, T5a/c/d/e autotransporter) | `unreachable`
(110, instance assigned but effector/machinery not co-located on a resolved replicon — reported,
not guessed). 547/930 carry a distance; only 95 sit within ±3 genes of machinery, the proximity
rule's recall ceiling. `protein_id` doubles as the stable key the ESM embedding step (4.3) caches on.

## Run dependency (rest of feature side, group 4)

The remaining feature columns — DLP/DSE/SignalP/PLM-E scores (4.1), the ESM embedding itself
(4.3), and the PU unlabeled set (4.4) — are **not** in `positives_all.tsv` and need the per-protein
tool signals from the Phase-2 ssign panel runs. They are joined in `training_dataset.tsv` (4.5) via
`bench_runout`. Those runs are on Imperial CX3. Until they land, `positives_all.tsv` +
`pair_features.tsv` are the run-independent backbone; positives whose genome has no run output yet
are reported as pending-run, never emitted as feature-complete. The T5a-neighbor DLP/DSE
observation stays an exploratory side-study in the benchmark (task 6b), it is never a training label.

## Files

- `positives_all.tsv` — **the deliverable** (930 positives, label-complete).
- `pair_features.tsv` — gene-distance pair-features + system-features per positive (task 4.2).
- `positives_tiered.tsv`, `positives_instanced.tsv`, `t5ss_effectors.tsv` — the three branch
  tables `37` unions; re-runnable inputs.
- `predicted_audited.tsv` + `predicted_audit_provenance.tsv` — per-row fate of the predicted audit.
- `t5ss_unplaceable.tsv` — 6 sourced T5SS examples with no resolvable locus_tag (excluded, not fabricated).
- `t5ss_sourcing_brief.md` — subtype label conventions + Pfam seeds + anti-hallucination contract.
- `predicted_instances_ambiguous.tsv`, `instance_audit_results.json`, `doi_recite.jsonl` — instance/citation audit provenance.
