## Why

ssign runs span minutes (base tier, small genome) to hours (extended tier with whole-genome predictors), and the user gets no estimate of how long a run will take until it finishes. We now have 196 cleaned per-tool calibration points across four machines, enough to fit a per-tool effort model and project a run's remaining time. The hard part is that the user's machine speed is unknown at run start, so the estimate must infer it online from the first tools that complete.

## What Changes

- Add an **offline effort-model fit**: one machine-agnostic `effort_T(size)` per tool, fit from the cleaned calibration data, with the whole-genome predictor regime tagged separately from the neighborhood regime. Leave-one-out validation reports per-tool error so we know which fits to trust.
- Add an **online ETA estimator**: at t=0 it knows the tool set and `n_proteins`, emits a wide prior ETA, then infers this machine's `cpu_rate` / `gpu_rate` / `io_factor` as each tool finishes and re-projects the remaining tools with a confidence interval.
- Add a **replay harness**: feed a finished run's per-tool wallclocks through the estimator step by step to measure how fast the ETA converges to the truth.
- The estimator is **tool-set-agnostic**: it composes per-tool models, so any tier or cherry-picked subset (even a never-before-run combination) is handled by summing whichever tools are active.
- Out of scope (follow-up change): wiring the estimator into `runner.py` as a live user-facing in-run ETA. This change ships the model + validation only, so the fits are proven before they surface to users.

## Capabilities

### New Capabilities
- `runtime-effort-model`: machine-agnostic per-tool effort fitting from cleaned calibration data, regime-aware (whole-genome vs neighborhood vs substrate), with leave-one-out error reporting and serializable fitted coefficients.
- `runtime-eta-online`: online machine-rate inference + remaining-time projection from a partially-complete run, plus the offline replay harness that validates convergence.

### Modified Capabilities

(none — no existing spec-level behavior changes; the pipeline's tool execution is untouched)

## Impact

- **New code (ssign package)**: a `runtime/` module holding the effort model, the online estimator, and bundled fitted coefficients (JSON). No changes to `runner.py` in this change.
- **Calibration data (local memory, not git)**: consumes `clean.py` output from `calibration/runs.jsonl`; the fit step reads the live cleaned set. The fitting script lives with the data; only the resulting coefficients ship in the package.
- **Dependencies**: numpy/scipy only (already core deps). No new packages.
- **No runtime/behavior change** for existing users until the follow-up wiring change.
