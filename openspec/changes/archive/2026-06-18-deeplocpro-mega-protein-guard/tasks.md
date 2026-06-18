## 1. Length-partition helper + constant

- [x] 1.1 Add `DEEPLOCPRO_MAX_AA = 5000` to `ssign_lib/constants.py` with a comment (kolossin/mega-protein rationale) and an env override (e.g. `SSIGN_DEEPLOCPRO_MAX_AA`)
- [x] 1.2 Add a pure `partition_by_length(records, max_aa) -> (kept: dict, skipped: list[(id, length)])` helper (in `run_deeplocpro.py` or `ssign_lib/fasta_io.py`), reusing `read_fasta_records`

## 2. Wire the guard into the DeepLocPro wrapper

- [x] 2.1 In `main()` (before dispatch to local/remote): read input, partition, write the kept records to a temp FASTA, pass that to `run_local_deeplocpro` / `run_remote_deeplocpro`; carry the skipped list forward
- [x] 2.2 After `parse_deeplocpro_output`, append one sentinel row per skipped protein: `predicted_localization="Not predicted (too long)"`, all probs `0.0`, `product` noting `skipped: length > {MAX} aa`
- [x] 2.3 Log a warning listing skipped protein IDs + lengths
- [x] 2.4 Confirm the 0-sequence and all-kept-empty edge cases still behave (no crash when everything is skipped or nothing is)

## 3. Tests

- [x] 3.1 `partition_by_length`: mixed input splits correctly at the threshold; boundary (== max is kept); all-short passthrough; all-long → empty kept
- [x] 3.2 Sentinel-row emission: a skipped protein yields exactly one not-predicted row with zero probs and the schema's columns
- [x] 3.3 Integration (no real DeepLocPro binary): with a stubbed/short-circuited model call, an over-length protein is excluded from what the model sees and present-as-skipped in the final output; a short protein passes through
- [x] 3.4 Downstream safety: a cross_validate test that a protein absent/skipped from DLP output is treated as non-extracellular without error

## 4. Validate

- [x] 4.1 Run the unit suite (`pytest tests/unit/ -v`) green
- [x] 4.2 Rerun BX470251 on CX3 — DeepLocPro completed, genome reached 24/24 (25 secreted proteins), plu2670 (16367 aa) + plu3123 (5457 aa) emitted as "Not predicted (too long)" → fleet at 67/67
- [x] 4.3 Note in NOTES.md that the mega-protein guard is shipped; clear the deferred item
