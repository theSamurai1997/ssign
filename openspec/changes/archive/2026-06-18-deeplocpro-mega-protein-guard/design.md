## Context

`run_deeplocpro.run_local_deeplocpro` invokes `deeplocpro -f input.faa -o out -g negative -d cuda` over the whole input with no length handling (the only batching is a 500-sequence-COUNT cap on the remote DTU path). DeepLocPro embeds each sequence with an ESM-style transformer whose attention memory scales ~O(L²); a 16,367-aa protein OOMs a 24 GB GPU and exits non-zero, and since DeepLocPro is a core step the whole genome fails. `ssign_lib/fasta_io.py` already provides `read_fasta_records` and `write_fasta`. The output schema (locus_tag, predicted_localization, extracellular_prob, periplasmic_prob, outer_membrane_prob, cytoplasmic_prob, cytoplasmic_membrane_prob, product) is produced by `parse_deeplocpro_output`.

## Goals / Non-Goals

**Goals:**
- A single over-length protein can never crash DeepLocPro (and thus the genome).
- Over-length proteins are visibly skipped (logged + a sentinel output row), not silently dropped.
- Works for both local and remote DeepLocPro paths with one shared step.

**Non-Goals:**
- Predicting localization for mega-proteins (truncating a 16k NRPS gives a meaningless call; skipping is correct since these are cytoplasmic machinery, not secretion substrates).
- Guarding DeepSecE / SignalP / PLM-Effector (they survived BX470251; defensive follow-up only).
- Changing the output schema or downstream consumers.

## Decisions

- **Skip, don't truncate.** A skipped protein simply has no DeepLocPro row's positive call; cross_validate's `dlp.get(L, {})` → `{}` → `is_dlp_positive` False, i.e. treated as non-extracellular. That is the correct answer for a cytoplasmic NRPS megasynthase, and truncation would fabricate a localization. Alternative (truncate to N-terminal) rejected: semantically wrong and still risks OOM at large caps.
- **Emit a sentinel row** for each skipped protein: `predicted_localization` = "Not predicted (too long)", all probs 0.0, `product` carrying "skipped: length > {MAX} aa". Keeps the skip explicit in the output and the figures/report, rather than the protein silently vanishing.
- **Partition once, before dispatch.** A pure helper `partition_by_length(records, max_aa) -> (kept, skipped_ids)` runs in `main()` (or a shared pre-step); the kept records are written to a temp FASTA passed to local/remote; skipped IDs are appended after parsing. This avoids duplicating the guard in both paths and keeps the dispatch functions unchanged.
- **Threshold = 5000 aa, configurable** via `DEEPLOCPRO_MAX_AA` in constants.py (and an env/CLI override). 5000 keeps essentially all real substrates (most are < 1500 aa; even large autotransporters rarely exceed it) while cutting the genuine megas. Validated by the BX470251 rerun; lower if a ≤5000 aa protein still OOMs.

## Risks / Trade-offs

- [A genuine large secretion substrate (>5000 aa, e.g. a giant adhesin) gets skipped] → It still appears in the output (flagged skipped) and other predictors (DSE/SignalP) still run on it; only its DLP localization is absent. Threshold is tunable. Acceptable: such proteins are rare and the alternative is a crashed run.
- [Threshold too high, a ≤5000 aa protein still OOMs on a smaller GPU] → tune `DEEPLOCPRO_MAX_AA` down; the mechanism is unaffected.

## Migration Plan

- Pure additive behavior; no schema or data migration. Existing runs unaffected (no genome in the prior fleet except BX470251 has a >5000 aa protein that mattered).
- Rerun BX470251 with the fix to confirm 67/67.

## Open Questions

- Extend the same partition guard to DeepSecE / SignalP / PLM-Effector defensively? Deferred; revisit if any of them is observed to OOM on a mega-protein.
- Expose `DEEPLOCPRO_MAX_AA` as a user-facing CLI flag now, or keep it a constant/env until someone needs it? Leaning constant + env for this change.
