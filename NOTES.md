# Deferred work

Tracks items skipped during tasks. One bullet per item: what, why, trigger to revisit.

## CX3 environment

- **signalp6 + deeplocpro PATH on fresh nodes** (task #22). Now auto-discovered via `_find_in_conda_envs` in `core/runner.py`. Binaries execute directly without `conda activate`, which works for pip-installed Python entry points because the shebang pins the env's Python. **Open risk:** if either tool's torch build needs system CUDA libs (rather than its own bundled CUDA wheels), `LD_LIBRARY_PATH` won't include `<env>/lib/` and the subprocess will crash with `libcudnn.so.X not found`. Mitigation if it ever fires: prepend `<env>/lib` to `LD_LIBRARY_PATH` in `run_signalp.py` / `run_deeplocpro.py` when the binary path lives under a conda env. Trigger to revisit: a CUDA-related ImportError from either wrapper.

## PLM-Effector performance

- **Ensemble model checkpoint cache** (task #16). `run_ensemble` now called 85× instead of 5×; re-reads ~150 MB of small files per call. Page cache makes real cost ~150 MB total, so lower priority than the raw number suggests. Revisit: if PLM-E step wallclock is still the long pole after `0445d94` cross-type caching is validated.
- **FP16/BF16 + `--batch-size 32`**. Could drop PLM-E ~74m → ~12-15m. Needs validation that FP16 doesn't shift predictions. A40 has 48 GB VRAM, room to grow batch size. Either alone gives ~2×; combined ~5×. Revisit: after cross-type caching speedup is confirmed and we have a baseline to compare against.

## Torch.load safety

- **Migrate `run_deepsece.py` to `weights_only=True` + `add_safe_globals`.** DeepSecE's checkpoint is a state_dict (already loaded via `model.load_state_dict(...)`), not a whole-module pickle like PLM-Effector's, so it's a safe candidate for the stricter loader PyTorch 2.6 introduced. PLM-E itself can't migrate (whole-module saves; would need an upstream refactor). Revisit: any time a new ssign dep bumps torch and `weights_only=False` triggers a deprecation warning.

## EggNOG

- **`--dbmem` blocked on GPU node** (task #2). EggNOG cache may still hang on SQLite mmap even with cache pre-warmed; `--dbmem` would fix that but needs 50 GB RAM. GPU nodes are capped at 32 GB. Workaround: run EggNOG step on a 64 GB CPU-only node, rejoin outputs. Revisit: once cross-type caching run completes and we attempt a full tier-2 again.
