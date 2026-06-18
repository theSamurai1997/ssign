## Why

A fleet genome (BX470251, Photorhabdus laumondii TT01) failed because DeepLocPro crashed (GPU OOM, deterministic) on **plu2670 / kolossin synthetase — 16,367 aa, the largest known prokaryotic protein** (UniProt Q7N3P5). DeepLocPro is a core step, so a single mega-protein cascaded the whole genome to failure (8/24 steps). The wrapper has no length guard. Any toxin / secondary-metabolite-rich genome (Photorhabdus, Xenorhabdus, giant adhesins, NRPS/PKS megasynthases) carries such proteins, so this will hit real users and undercuts the publication / zero-maintenance goals.

## What Changes

- Partition the DeepLocPro input by sequence length before running: sequences over a max length are set aside, the rest go to DeepLocPro as normal.
- Over-length proteins are emitted in the output as explicitly *not predicted* (sentinel localization, zero probabilities, a "skipped: length > N aa" note) and a warning lists them, so a single mega-protein can no longer crash the run, and the skip is visible rather than silent.
- Add a configurable max-length constant (default 5000 aa).
- Apply to both the local and remote (DTU) DeepLocPro paths via a shared partition step.

## Capabilities

### New Capabilities
- `deeplocpro-length-guard`: how the DeepLocPro wrapper handles sequences too long for the model — the length threshold, that over-length proteins are skipped (not crashed on), and how they appear in the output.

### Modified Capabilities
<!-- none: fresh OpenSpec repo, no existing predictor spec -->

## Impact

- `src/ssign_app/scripts/run_deeplocpro.py`: partition input before dispatch (local + remote); append skipped proteins to the parsed output.
- `src/ssign_app/scripts/ssign_lib/constants.py`: add `DEEPLOCPRO_MAX_AA` (default 5000).
- Possibly `ssign_lib/fasta_io.py`: a small reusable `partition_by_length` helper (or keep it local to the wrapper).
- Tests: `tests/unit/test_run_deeplocpro.py` — partition logic, skipped-row emission, all-short passthrough, an over-length protein is excluded from the DeepLocPro call and present-as-skipped in output.
- Downstream: a skipped protein has no extracellular call → cross_validate treats it as non-secreted, which is correct for a cytoplasmic NRPS megasynthase. No schema change.
- Out of scope (note for later): the same defensive guard for DeepSecE / SignalP / PLM-Effector. They survived BX470251 (DeepSecE handled the full proteome; SignalP ran neighborhood-only), so only DeepLocPro is fixed here.
