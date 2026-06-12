# Effector-recovery benchmark — Phase 2 results (full panel, 2026-06-12)

What fraction of literature-curated secretion-system effectors does as-shipped ssign actually
emit, and how does that compare to the ceiling that ssign's proximity rule could in principle
reach? Run on all 67 panel genomes, both with T3SS excluded (the default exclusion) and included.

Reproduce: `scripts/24_actual_call.py` → `25_actual_vs_ceiling.py` → `26_found_systems.py` →
`27_phase2_figures.py`, each `--run-tag panel_genbank_default` and `panel_genbank_t3ss`.
Figures in `figures/01..04`. Per-effector tables in `actual_per_effector.*.tsv`,
`actual_vs_ceiling.*.tsv`, `found_systems.*.tsv`.

## Dataset under test

| | T3SS excluded | T3SS included |
|---|---|---|
| Secretion systems (curated, with ≥1 testable effector) | 52 | 83 |
| Secreted-protein effectors (testable) | 272 | 499 |

"Testable" = the effector's genome was staged and its ORF is in ssign's input. 83 of the 582
curated effectors are non-testable (no genome accession / not in the panel) and are excluded from
every count below.

## Headline: a large ceiling-to-actual gap (Figure 01, 02)

Recall as a fraction of testable effectors (n=499):

| SS | testable | findable @±3 | @±5 | @±7 | found (excl) | found (incl) |
|----|---|---|---|---|---|---|
| T1SS | 25 | — | — | 20 | 15 | 15 |
| T2SS | 77 | 1 | 2 | 3 | 4 | 4 |
| T3SS | 227 | 48 | 66 | 75 | 3 | 15 |
| T4SS | 98 | 6 | 8 | 10 | 0 | 0 |
| T6SS | 72 | 16 | 17 | 19 | 17 | 17 |
| **ALL** | **499** | **91** | **113** | **127** | **39** | **51** |

- ssign runs proximity at **±3** (its default). The honest apples-to-apples is found vs findable@±3:
  **39 of 91** (excluded) / **51 of 91** (included). Even at the generous ±7 ceiling, proximity
  could only ever reach **127/499 (25%)**; ssign emits **8–10%**.
- **T4SS is the starkest gap:** ceiling 10, found **0**. ssign emits nothing for T4SS effectors.
- **T3SS** default=3 is by design (T3SS systems excluded by default → effectors auto-miss); the
  included pass lifts it to 15, still a fraction of the 75 ceiling.
- **T1SS / T6SS** are the near-ceiling types (15/20, 17/19): where proximity already does most of
  the reachable work, ssign captures it.

The gap between the ceiling (what proximity could reach) and the actual (what ssign emits) is the
core motivation for the secretion-classifier: a learned model has headroom the proximity rule is
leaving on the table.

## Why each found effector was emitted (Figure 03)

Mapping each found effector to the specific ssign-detected system it sits within ±3 of (the
systems are recorded under unique `sys_id`s in each run's results.csv):

| | found | distinct detected systems | own-type nearby (legit) | only different-type nearby (accidental) |
|---|---|---|---|---|
| T3SS excluded | 39 | 32 | 33 | 6 |
| T3SS included | 51 | 39 | 45 | 6 |

The large majority of found effectors sit next to a system of **their own type** (right protein,
right reason). Only **6** are accidental: emitted purely because an unrelated system happened to be
within ±3 genes (celA/plaA next to a T4 pilus; BopA/BopE next to a T5c autotransporter; VirA/
ChlaDub1 next to a T5a autotransporter). Those 6 overlap the audit's misassigned/cross-type rows.

## Answer-key citation integrity (Figure 04)

Of the 19 effectors ssign emitted but whose curated machinery sits >7 genes away
(`discordant_audit.md`, one literature agent per sourcing DOI), only **5/19 are fully sound**:

- **6** have a wrong or non-existent sourcing DOI (resolve to unrelated papers — soil ecology,
  mouse-tooth development, lung regeneration — or 404), though the effector biology is often correct.
- **4** are unidentifiable/unsupported rows (EFF00142, EFF00150 placeholder IDs; TseA_T6SS1 never
  named in its cited paper; ChlaDub1 has no T3SS evidence).
- **2** are genuine SS-type misassignments (BopA, BopE are Bsa **T3SS** effectors, mislabelled T6SS
  via a "Bop" name collision).
- **2** are a duplicate row (Tle4 = TplE_alias_Tle4 = PA1510).

This does **not** affect the recall numbers (recall is scored on machinery gene-position, not the
effector's citation), but it is a real training-set quality problem, and the same defect class
almost certainly extends beyond these 19 rows. Repair list tracked in `../../../NOTES.md`.

## Caveats

- Recall is measured against literature-curated machinery loci; where MacSyFinder detects a real
  same-type system at a different instance than the answer key anchored, the effector reads as
  "unreachable@7" yet is correctly emitted (13 of the 19 discordant rows are this benign case).
- GenBank input with `--use-input-annotations` (locus_tags preserved); the pilot confirmed FASTA
  input gives identical recall, so only GenBank was run at panel scale.
