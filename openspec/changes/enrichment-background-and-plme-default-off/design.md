## Context

The enrichment test (`enrichment_testing.py`) compares the fraction of tool-positive proteins in a secretion system's ±3 neighborhood against a genome background `p_bg`, estimated from a random null sample of non-neighborhood proteins (`sample_null_proteins.py`, default 200). Validation on PAO1 showed the 200-sample background runs low (sampling noise on a ~1.5% base rate) and that PLM-Effector, recently wired in as a third predictor, is unfit for this test: it calls 25% of the proteome positive and shows no real per-system enrichment. The predictions themselves already run whole-genome on the benchmark fleet, so an exact background is available for free there. This change corrects the background estimate and demotes PLM-Effector to opt-in.

## Goals / Non-Goals

**Goals:**
- More accurate enrichment significance via a larger (or exact) background.
- PLM-Effector off by default; no PLM-Effector in the enrichment test.
- Consistent ≥ 0.8 positivity convention across DLP, DSE, and (when used) PLME.
- No change to the proximity / cross-validation substrate-finding logic beyond the PLME gate.

**Non-Goals:**
- Removing PLM-Effector from the install (vendored code, weights, extra). Deferred.
- Re-tuning DLP/DSE thresholds or the proximity window.
- Changing the binomial test or BH FDR method.

## Decisions

- **n_null default 200 → 1000.** Empirically (null-size sweep) 1000 converges to the exact whole-genome background while 200 over-calls. Alternative considered: a Laplace/pseudocount smoothing of the 200-sample rate — rejected as a band-aid that still discards available signal. 1000 is a one-line default change with a known compute cost (~+6 min/genome only in neighborhood mode; free in whole-genome mode).
- **Exact background when predictors ran whole-genome.** When the whole-genome predictions exist (fleet, any `--*-whole-genome` run), compute `p_bg` from ALL non-neighborhood, non-component proteins rather than sampling. Zero extra compute, removes sampling noise entirely. The null-sample path stays for default neighborhood-mode runs (where only sampled proteins have predictions).
- **PLM-Effector off by default.** `skip_plm_effector` resolves to True unless explicitly enabled. Keeps it available for opt-in and for the secretion-classifier feature matrix. Alternative (full removal) rejected for now: the classifier work consumes PLME as a feature.
- **Drop PLME from the enrichment test.** `score_scope` emits DLP/DSE only; the runner stops passing `--plme`. Removes a 25%-background predictor whose inclusion is misleading and saves the (already-spent) compute of feeding it.
- **PLME 0.8 max_prob gate where it remains a binary call.** Every per-type native threshold is ≤ 0.8, so a max_prob ≥ 0.8 gate is strictly stricter than the native `passes_threshold` and a clean, uniform rule. Applies in `enrichment_testing.is_plme_positive` (for the opt-in case) and `cross_validate_predictions._plm_effector_flag`.

## Risks / Trade-offs

- [Changing the PLME cross-validation gate shifts substrate calls in existing opt-in runs] → Document in the change; the gate is stricter (fewer, higher-confidence PLME calls), which is the intended direction. Covered by updated unit tests.
- [n_null=1000 adds compute in default neighborhood-mode runs] → ~+6 min/genome worst case (per the effort fit); acceptable, and free in whole-genome mode where most fleet runs operate.
- [Users with scripts that rely on PLME running by default see different output] → It is a deliberate default-behavior change, called out as BREAKING in the proposal and in docs/CHANGELOG.

## Migration Plan

- Pure default + logic change; no data migration. Existing outputs are not rewritten.
- Rollback: revert the commit; defaults return to 200 / PLME-on.
- Fleet: CX3 checkout must `git pull` before launch to pick this up (also fixes the pre-existing PLME-enrichment-wiring gap).

## Open Questions

- Full removal of PLM-Effector from the install — revisit after the secretion-classifier-dataset work decides whether PLME features are worth keeping.
- Whether to expose `n_null` as a user-tunable CLI/GUI option (consistent with the user-controllable-thresholds direction) or leave it a constant for now.
