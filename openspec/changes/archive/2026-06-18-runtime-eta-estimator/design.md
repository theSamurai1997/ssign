## Context

The calibration system (`calibration/runs.jsonl` + `clean.py`, in local memory, not git) now yields 196 clean per-tool points across CX3-A40/RTX6000/L40S/A100. The README already specifies the model shape: `wallclock_T = effort_T(size) / rate_lim(T, machine)`, with `effort_T` machine-agnostic (fit offline, CX3-A40 = reference rate 1.0) and `rate_lim` latent at run start, inferred online from completions. This change implements that model and validates it; it does **not** wire it into `runner.py` yet.

Key facts the design must respect:
- Three size axes, one per tool class: `n_proteins` (whole-genome tools + predictors with a whole-genome flag), `n_neighborhood` (default predictors), `n_substrates` (annotation tools). `clean.py` already assigns the correct axis per row.
- Predictors have two regimes with very different slopes (whole-genome PLM-E ~40 min vs neighborhood ~minutes); they must be fit separately.
- Most real users run the low-n regime (neighborhood/substrate sizes), so accuracy there matters more than at whole-genome scale.
- The whole-genome predictor regime is currently thin (2-3 points); the 67-genome benchmark fleet will fill it. The model must degrade gracefully when a regime is data-poor.

## Goals / Non-Goals

**Goals:**
- A per-tool offline effort fit, regime-aware, with leave-one-out error so we know which fits to trust.
- An online estimator that emits a prior at t=0 and tightens as cpu/gpu/io rates are inferred, for any tool subset.
- A replay harness that quantifies ETA convergence against finished runs.
- Bundle fitted coefficients as a JSON artifact the package ships; keep the raw data + fitter local.

**Non-Goals:**
- Wiring into `runner.py` / live progress UI (follow-up change).
- Non-linear effort models (start linear; revisit only if LOO error demands it).
- Predicting DTU-remote (network-bound) tool times; those are queue-dominated and tracked separately.
- Per-protein-length effort (size = protein count, not residue count) for v1.

## Decisions

**1. Linear effort `effort_T = a*size + b`, per (tool, regime).** Rationale: the README calls for it, the data is roughly linear in protein count, and it's interpretable. Alternatives: log-linear or piecewise (rejected for v1 — not enough points per tool to justify; revisit if LOO error is high). The intercept `b` absorbs fixed overhead (model load, process spawn), which dominates at low n and is exactly the regime we care most about.

**2. Three latent machine rates (`cpu_rate`, `gpu_rate`, `io_factor`), assigned by a per-tool limiting-factor table.** The table already exists in the README. Each tool maps to one limiting factor; its observed time updates only that rate. Alternative: a single global speed scalar (rejected — a fast-GPU/slow-NFS machine like CX3 on NFS would be mis-modeled; EggNOG's 0-8557s spread proves IO must be separable).

**3. Reference machine = CX3-A40 (rate 1.0).** Effort is fit in A40-seconds. A new machine with `gpu_rate=2.0` runs GPU tools twice as fast. Since `clean.py` does not yet normalize across machines, v1 fits effort on the pooled clean points and treats cross-machine scatter as noise folded into the CI; the online step is what corrects for the actual machine. A later refinement can jointly fit effort + per-machine rates (matrix factorization), but that needs more machines.

**4. Rate update = recency-light running mean.** Each completion contributes `rate_obs = effort_T(size)/wallclock`; the machine rate is the mean of observations for that limiting factor, optionally down-weighting tools whose LOO error is high (their `rate_obs` is noisier). Alternative: Kalman filter (rejected — overkill for ~5-10 completions per run).

**5. Module layout.** Shipped package: `src/ssign_app/runtime/effort_model.py` (load coefficients, `effort(tool, size, regime)`), `src/ssign_app/runtime/estimator.py` (online state machine), `src/ssign_app/runtime/coefficients.json` (bundled fit). Local-only (with the data): `calibration/fit.py` (reads `clean.py`, writes coefficients.json + LOO report), `calibration/replay.py` (convergence harness). This keeps raw data and the fitter out of the package while shipping only the small JSON.

**6. Limiting-factor table is the contract between fit and online.** It lives in one place (a dict in `effort_model.py`, mirrored from the README) so the online estimator and the fit agree on which rate each tool keys to.

## Risks / Trade-offs

- [Thin whole-genome predictor regime (2-3 pts)] → Mark those fits low-confidence; widen their CI; refit when the benchmark fleet lands ~67 more whole-genome predictor points. The estimator must not present a confident whole-genome ETA yet.
- [Cross-machine scatter folded into effort noise in v1] → Acceptable because the online step corrects per-machine; document that effort CIs are inflated until joint fitting lands. Revisit when ≥4 machines have ≥5 shared tools.
- [EggNOG IO bimodality (NFS vs `--dbmem`)] → Handle via `io_factor`, not a separate effort curve; the online step detects the slow branch from EggNOG's own completion. Risk: EggNOG is often late in the run, so IO is inferred late — note this limitation.
- [Stale coefficients] → coefficients.json records the calibration row-count + date; the fit step warns if the live clean set has grown materially since the bundled fit.
- [Mis-sized rows silently re-entering] → The fit only ever consumes `clean.py`'s clean set, never raw `runs.jsonl`; this is enforced by importing `clean_rows`, not re-reading the file.

## Migration Plan

Additive only — new `runtime/` module + local fit/replay scripts. No existing code paths change, so no rollback concern. The follow-up change introduces the `runner.py` call site behind a flag.

## Open Questions

- CI method: bootstrap over the clean points, or analytic from the linear fit's residual std? Lean bootstrap (small n, non-normal). Decide at implementation.
- Should `_pipeline_total` rows feed the t=0 prior directly, or should the prior be the per-tool effort sum with a tier-level correction? Lean per-tool sum (composes for cherry-picked subsets); use `_pipeline_total` only to sanity-check the sum.
