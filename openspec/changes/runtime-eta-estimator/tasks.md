## 1. Offline effort fit (local, with the calibration data)

- [x] 1.1 Add the per-tool limiting-factor table (tool → cpu/gpu/io) as a dict, mirrored from the README, in the shipped `effort_model.py` (single source of truth for both fit and online).
- [x] 1.2 Write `calibration/fit.py`: import `clean_rows` from `clean.py`, group clean points by (tool, regime), fit linear `effort = a*size + b` per group (numpy lstsq); for groups with <3 points or a single size, fall back to mean wallclock and tag low-confidence.
- [x] 1.3 Add leave-one-out CV in `fit.py`: per (tool, regime) report median absolute error (s) and median percentage error; flag groups whose median % error exceeds a stated threshold.
- [x] 1.4 Serialize the fit to `src/ssign_app/runtime/coefficients.json` (per tool/regime: a, b, n_points, loo_pct_err, confidence flag) + a header recording calibration row-count and date.
- [x] 1.5 Print the LOO report and show Teo which fits are trustworthy vs thin (expect whole-genome predictor regime flagged low-confidence).

## 2. Effort model (shipped package)

- [x] 2.1 Create `src/ssign_app/runtime/__init__.py` and `effort_model.py`: load `coefficients.json`, expose `effort(tool, size, regime)` returning reference-machine seconds; raise/clear-fallback for unknown tools.
- [x] 2.2 Add regime resolution: given a tool + the run's flags (whole-genome vs neighborhood), pick the right (tool, regime) coefficients and the right size field.
- [x] 2.3 Unit-test effort_model: known coefficients round-trip; regime selection picks the right curve; low-confidence tools return effort + a wide flag.

## 3. Online estimator (shipped package)

- [x] 3.1 Implement `estimator.py` state: known tool set, `n_proteins`, per-limiting-factor rate estimates (start unset), list of (tool, regime, size) still pending.
- [x] 3.2 Implement the t=0 prior: sum per-tool efforts (max over the parallel predictor block, sum elsewhere) with a wide CI from the historical machine-rate spread; sanity-check against `_pipeline_total` rows at similar size/tier.
- [x] 3.3 Implement `on_tool_complete(tool, size, wallclock)`: infer `rate_obs = effort/wallclock` for that tool's limiting factor, update the factor's running mean (down-weight high-LOO tools), re-project pending tools, return ETA + CI.
- [x] 3.4 Handle the parallel predictor block correctly (contribution = max of projected member times, not sum).
- [x] 3.5 Unit-test the estimator: prior is wide; first CPU completion narrows CPU-bound ETAs; first GPU completion narrows GPU-bound; a synthetic NFS-slow EggNOG raises io_factor; a cherry-picked never-seen subset still returns an ETA.

## 4. Replay / convergence harness (local)

- [x] 4.1 Write `calibration/replay.py`: take a finished run's per-tool wallclocks, feed them to the estimator in completion order, record ETA + CI + error after each step (using only data available up to that step).
- [x] 4.2 Run replay on the existing clean runs + the PAO1 smoke run; produce a numbered convergence figure (ETA-vs-truth as completions accrue) per the publication-plots defaults.
- [x] 4.3 Report convergence summary: after how many completions the ETA lands within X% of the observed total, per tier.

## 5. Wire-up readiness + docs

- [x] 5.1 Document the runtime module + how to refit (re-run `clean.py` → `fit.py` → commit new coefficients.json) in the calibration README and a short module docstring.
- [x] 5.2 Note the follow-up: where `runner.py` would call the estimator (t=0 + per-step-callback hooks) without implementing it, so the next change has a clear seam.
- [x] 5.3 Run the simplify skill over the new module + scripts; fix findings.
