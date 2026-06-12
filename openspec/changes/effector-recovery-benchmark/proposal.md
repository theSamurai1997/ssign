## Why

ssign keeps its substrate set small by only predicting on proteins within +/-N genes of a secretion-system component. We don't actually know how much real biology that proximity rule misses. We need a defensible, fully-cited benchmark that says, per secretion-system type: of the proteins a system is experimentally known to secrete, how many could the proximity rule ever reach, how many are structurally out of reach, and how many ssign actually recovers when run.

A previous attempt produced figures (analysis scripts 05-11) that conflated two different questions and reported artifacts: it judged ssign against an incomplete reference list and used ssign's own detector to define the answer key, making the test partly circular. That work is being discarded and rebuilt correctly.

## What Changes

- **Build a verified effector gold set** from the existing literature-curated corpus (`secretion-classifier/data/verified/`, 949 rows, T1-T6SS). Scope to experimentally `validated` rows, repair broken citations and UniProt IDs on the recoverable `PARTIAL` rows, drop the ~50 biology errors (T6SS apparatus-as-effector, immunity proteins, eukaryotic-host artifacts) and `FAIL` rows. Keep only rows that are both `validated` and `VERIFIED`.
- **Reconcile against externally-curated effector databases** to control the corpus's own selection bias: SecReT4 (T4SS), SecReT6 (T6SS), SecretEPDB and EffectiveDB (T3/T4/T6), and BastionHub (the only resource covering T1SS/T2SS depth). Pull each, map identifiers to our genomes, add net-new experimentally-supported effectors, and verify every addition to the same evidence + citation bar as the corpus. (BastionHub and SecretEPDB are currently unreachable, a Monash host outage; fall back to the BastionHub NAR-2021 supplementary tables / Wayback snapshot until access returns.)
- **T5SS is excluded from the main benchmark.** Its passenger stays attached to the cell surface (not released) and no curated database tracks it as a secreted effector, so the ceiling/actual numbers cover T1/T2/T3/T4/T6 only. The separate observation that ssign flags predicted secreted proteins near T5aSS is logged as an exploratory side study, not part of the headline result.
- **Build an ssign-independent machinery answer key**: for each of the ~98 distinct system instances (across ~65 genomes), the complete set of apparatus gene locus_tags and genome coordinates. Sourced from primary literature only (one agent per instance, strict anti-hallucination output schema), then RefSeq used solely to resolve confirmed gene names to coordinates, then an independent verification pass. **MacSyFinder must not be used here**, since it is ssign's core detector and would make the benchmark circular.
- **Compute the proximity ceiling**: per verified effector, gene-distance to the nearest component of its own system instance; classify reachable vs structurally-impossible at N = 3, 5, 7. Report per SS type and per genome. No ssign run required.
- **Measure actual recall**: run the chosen genome panel through ssign, bridge Bakta locus_tags back to RefSeq by coordinate overlap, count what ssign emits as secreted, and compare against the ceiling.
- **Document everything** as a new, reproducible analysis section under `validation_sweeps/`, with methods, citations, strengths, and limitations.
- **BREAKING (analysis only, no product code):** remove the previous benchmark scripts (`analysis/05_*`–`analysis/11_*`) and figures (`09`, `10`, `11`). The proximity/fragmentation/copies sweep figures (`01`-`04`, `07`) are unaffected and stay.

## Capabilities

### New Capabilities
- `effector-gold-set`: the curated, experimentally-validated, citation-verified effector reference dataset, its repair/verification procedure, and its reconciliation with externally-curated effector databases.
- `machinery-answer-key`: the ssign-independent, literature-sourced apparatus gene location dataset (per system instance) and the agent + verification procedure that builds it.
- `proximity-ceiling-analysis`: the per-effector reachability classification at N = 3, 5, 7 and its per-type / per-genome reporting.
- `actual-recall-benchmark`: the ssign run, Bakta-to-RefSeq coordinate bridge, and ceiling-vs-actual comparison.

### Modified Capabilities
<!-- None. openspec/specs/ is empty; this change introduces analysis capabilities only and touches no existing ssign product spec. -->

## Impact

- **New files** under `validation_sweeps/`: a rebuilt `analysis/` benchmark pipeline, a `benchmark_data/` set of verified effector + machinery answer-key tables, figures, and a docs section.
- **Removed**: `validation_sweeps/analysis/05_*`–`11_*` scripts and figures `09`/`10`/`11` (gitignored research one-offs; the committed sweep figures stay).
- **Reads from** the `secretion-classifier` repo's `data/verified/` corpus as the effector starting point (copied in, not modified there); `predicted`-evidence rows are left untouched for later use.
- **External dependencies**: UniProt REST and NCBI RefSeq (citation/ID/coordinate resolution), PubMed/web access for the literature-curation and verification agents, and the externally-curated effector databases (SecReT4, SecReT6, SecretEPDB, EffectiveDB, BastionHub). No change to ssign's runtime code or dependencies. Note: BastionHub/SecretEPDB access is currently down; the BastionHub NAR-2021 supplement is the fallback for T1SS/T2SS.
- **Compute**: Phase 2 requires ssign runs on the genome panel (CX3 GPU jobs); Phases 0-1 are local + network only.
