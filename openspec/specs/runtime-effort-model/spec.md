# runtime-effort-model Specification

## Purpose
TBD - created by archiving change runtime-eta-estimator. Update Purpose after archive.
## Requirements
### Requirement: Per-tool effort model fit from cleaned calibration data

The system SHALL fit one machine-agnostic effort function `effort_T(size)` per tool from the cleaned calibration fit-set (the output of `clean.py`). The model SHALL start as linear, `effort_T = a_T * size + b_T`, and SHALL be expressed in reference-machine units (CX3-A40 rate = 1.0). The fit SHALL only consume rows that `clean.py` marks as clean fit points; dropped, mis-sized, aggregate, and buggy rows MUST be excluded.

#### Scenario: Fit a tool with sufficient clean points
- **WHEN** a tool has at least 3 clean fit points spanning more than one size
- **THEN** the system produces `a_T`, `b_T`, and a goodness-of-fit metric for that tool

#### Scenario: Tool with too few points
- **WHEN** a tool has fewer than 3 clean fit points
- **THEN** the system marks the tool's effort model as low-confidence and falls back to the mean wallclock of the available points rather than extrapolating a slope

### Requirement: Regime separation for predictor tools

For the predictor tools (DeepLocPro, SignalP, DeepSecE, PLM-Effector), the system SHALL fit the whole-genome regime (`n_proteins`) separately from the neighborhood regime (`n_neighborhood`), because they are different effort curves. The system MUST select the regime matching how the tool will actually run (whole-genome flag vs default neighborhood).

#### Scenario: Neighborhood-regime prediction for a default run
- **WHEN** a predictor will run on the +/-3 neighborhood (no whole-genome flag)
- **THEN** the effort model used is the neighborhood-regime fit keyed on `n_neighborhood`

#### Scenario: Whole-genome-regime prediction
- **WHEN** a predictor will run with its whole-genome flag set
- **THEN** the effort model used is the whole-genome-regime fit keyed on `n_proteins`

### Requirement: Leave-one-out validation report

The system SHALL report leave-one-out cross-validation error per tool, so the trustworthiness of each fit is visible before the coefficients ship. The report SHALL express error in both absolute seconds and percentage of observed wallclock.

#### Scenario: Validation surfaces an unreliable fit
- **WHEN** a tool's leave-one-out median percentage error exceeds a stated threshold
- **THEN** the report flags that tool so its ETA contribution is treated as wide / low-confidence

### Requirement: Serializable fitted coefficients

The system SHALL serialize the fitted coefficients (per tool, per regime) to a JSON artifact that the online estimator loads at runtime. The artifact SHALL record the calibration data version (row count and date) it was fit from, so a stale fit is detectable.

#### Scenario: Coefficients round-trip
- **WHEN** the fit is serialized and reloaded by the estimator
- **THEN** the reloaded coefficients reproduce the same `effort_T(size)` values

