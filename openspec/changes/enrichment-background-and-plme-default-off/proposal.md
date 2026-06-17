## Why

Validation of the enrichment statistical test on the PAO1 smoke run (job 3013556; analysis in `validation_sweeps/benchmark/analysis/enrichment_validation/`) surfaced two problems. First, the 200-protein null sample underestimates the genome background positive rate (DLP 0.5% sampled vs ~1.3-1.65% true; DSE 1.0% vs ~1.4-1.7%), which inflates significance; a null-size sweep shows 200 over-calls while 1000 matches the all-proteins exact value. Second, PLM-Effector calls ~25% of the PAO1 proteome as secreted effectors (18% even gated at max_prob ≥ 0.8) and contributes no reliable enrichment signal (2/18 systems significant, the real T3SS depleted). Genome-scale over-prediction is a documented, expected failure mode for sequence-only effector predictors trained on balanced sets (McDermott 2010, PMID 20974833; T3SEpp, PMC7406222; the PLM-Effector paper itself, Zheng 2026 bbag143, reports no genome-scale validation). PAO1 truly has ~4 T3SS effectors, so 25% is biologically implausible.

## What Changes

- Bump the default null-sample size from 200 to 1000 for the enrichment background estimate.
- When predictors ran whole-genome, estimate the background from ALL non-neighborhood proteins (exact, zero extra compute) instead of a random subsample.
- **BREAKING (default behavior):** PLM-Effector runs OFF by default for all ssign runs. It remains installable and opt-in (extended tier); it is NOT removed from the install.
- Drop PLM-Effector from the enrichment test entirely (DLP/DSE only), removing both the misleading signal and the extra compute.
- Gate PLM-Effector positivity at max_prob ≥ 0.8 wherever it remains a binary call (enrichment + cross-validation), for consistency with the DLP/DSE 0.8 convention.

## Capabilities

### New Capabilities
- `enrichment-stats`: the per-system/per-broad-type binomial enrichment test — which predictors it tests (DLP/DSE only), how the genome background is estimated (null-sample size, exact whole-genome path), and the BH FDR call.
- `plme-prediction`: PLM-Effector's role in the pipeline — off by default, opt-in only, max_prob ≥ 0.8 positivity gate, excluded from the enrichment test, retained as a proximity/cross-validation screen and a classifier-training feature.

### Modified Capabilities
<!-- none: fresh OpenSpec repo, no existing specs -->

## Impact

- `src/ssign_app/core/runner.py`: `n_null_proteins` default 200→1000; `skip_plm_effector` default resolves to skipped; enrichment-step wiring (drop `--plme`); whole-genome exact-background path.
- `src/ssign_app/scripts/enrichment_testing.py`: remove PLME from `score_scope`/CLI; `is_plme_positive` gains the 0.8 gate; background source.
- `src/ssign_app/scripts/sample_null_proteins.py`: default sample size; all-non-neighborhood mode.
- `src/ssign_app/scripts/cross_validate_predictions.py`: `_plm_effector_flag` 0.8 gate.
- `src/ssign_app/scripts/ssign_lib/constants.py`: any shared default/threshold constants.
- CLI (`cli.py`) + Streamlit GUI defaults; docs (README, configure, output_files); unit tests (`test_enrichment_testing.py`, runner/CLI default tests).
- Out of scope (deferred open question): fully removing PLM-Effector from the install (vendored `scripts/plm_effector/`, weights download, `pyproject` extra) — pending the secretion-classifier-dataset decision on whether PLM-E features earn their place.
