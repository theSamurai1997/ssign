## 1. PLM-Effector off by default

- [x] 1.1 Make `skip_plm_effector` resolve to True by default (`constants.py`: moved `plm_effector` out of `_BASE_ENABLED` into a new `_OPT_IN` set kept in `_ALL_TOOLS`, so the tier resolver emits skip=True at every tier)
- [x] 1.2 CLI opt-in via the existing `--no-skip-plm-effector` (BooleanOptionalAction); help text updated. Streamlit default off (the GUI checkbox now enables PLME when checked)
- [x] 1.3 PLM-Effector stays installable: no change to `pyproject` extras, vendored `scripts/plm_effector/`, or weights handling
- [x] 1.4 Unit tests: `test_plm_effector_off_by_default_every_tier`, `test_plm_effector_opt_in_respected`, updated `test_skip_flags_align_with_install_tier`

## 2. Drop PLM-Effector from the enrichment test

- [x] 2.1 `enrichment_testing.py`: removed PLME from `score_scope` (DLP/DSE only) and dropped the `--plme`/`p_plme`/`plme` plumbing from `main`
- [x] 2.2 `runner.py` `_step_enrichment`: stopped passing `--plme`
- [x] 2.3 `test_enrichment_testing.py`: `test_plme_never_emitted`; score_scope tests updated to the DLP/DSE-only signature

## 3. PLM-Effector 0.8 positivity gate (where still a binary call)

- [x] 3.1 N/A as written: `is_plme_positive` was REMOVED from `enrichment_testing.py` (PLME is fully dropped from the test, so a gated-but-unused function would be dead code). The gate lives only where PLME remains a binary call — see 3.2.
- [x] 3.2 `cross_validate_predictions._plm_effector_flag(plm_row, conf_threshold)`: require `max_stacking >= conf` (0.8) in addition to `passes_threshold`; caller threads `conf_threshold`
- [x] 3.3 Unit tests: `test_plm_effector_below_confidence_gate_is_negative` (0.6 → negative), `test_plm_effector_at_confidence_gate_is_positive` (0.8 → positive); `_plm_e_row` helper carries `max_stacking`

## 4. Enrichment background sizing

- [x] 4.1 Default null-sample size 200 → 1000 (`runner.py` `n_null_proteins`, `sample_null_proteins.py` `--n` default)
- [x] 4.2 Exact-background path: `runner._step_sample_null_proteins` passes `--n -1` (= all) when DLP and DSE both ran whole-genome; `sample_null` treats `n <= 0` as "entire complement"
- [x] 4.3 Small-pool case handled (`n >= len(pool)` returns the whole pool; existing `test_small_pool_returns_all`)
- [x] 4.4 Unit tests: `test_n_null_proteins_default_is_1000`, `test_nonpositive_n_returns_all_candidates`

## 5. Docs + consistency sweep

- [x] 5.1 Updated README, CLAUDE.md key-params, design_decisions.md (§3.1 revision), CHANGELOG (default-behaviour change). `docs/reference/cli.md` already documented PLME off-by-default
- [x] 5.2 Grepped docs/CLI/GUI for stale PLME-on / "equal predictor" / n_null=200 references and fixed the load-bearing ones
- [x] 5.3 Full unit suite green: 1371 passed

## 6. Validation

- [x] 6.1 Verified the updated enrichment on real PAO1 data (phase2_work): DLP/DSE-only output, exact whole-genome background (5282 non-neighborhood proteins, p_DLP=0.0134, p_DSE=0.0144). Note: the `02_*` analysis driver still passes `--plme` and would need that removed to re-run as-is; production behaviour is covered by unit tests + this smoke
- [x] 6.2 NOTES.md already records that CX3 must `git pull` before the fleet launch
