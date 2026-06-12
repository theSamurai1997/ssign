# Effector-recovery benchmark

How well can ssign's proximity-based substrate prediction recover experimentally-verified secretion-system effectors, and how well does it actually do it? Per SS type (T1, T2, T3, T4, T6) we report three numbers:

1. **Ceiling** — fraction of verified effectors within +/-N genes (N=3,5,7) of a component of their own system instance. The most the proximity rule could ever recover.
2. **Impossible** — the complement. Effectors structurally out of reach of the proximity rule.
3. **Actual** — what ssign emits as secreted when the genomes are run, bridged back to RefSeq coordinates.

This is the rebuilt benchmark. It replaces the earlier `../analysis/05_*`-`11_*` attempt, which conflated the ceiling and actual questions and used MacSyFinder (ssign's own detector) to place the apparatus, making the test circular. T5SS is excluded from the headline numbers (its product stays cell-surface-attached; no curated database tracks it).

Plan of record: OpenSpec change `effector-recovery-benchmark` (`openspec/changes/effector-recovery-benchmark/`).

## Layout

```
data/
  source_corpus/    read-only copy of the secretion-classifier verified corpus + audit (see PROVENANCE.txt)
  refseq_cache/      RefSeq GenBank per genome (coordinate source)
  ...                gold set, machinery answer key, and result tables land here as phases complete
scripts/             analysis pipeline (numbered)
figures/             numbered output figures
```

## Status

Phase 0 and Phase 1 complete; Phase 2 (actual recall) pending. Tasks tracked in the OpenSpec
change's `tasks.md`.

- **Phase 1 ceiling result + method:** [`docs/phase1_ceiling.md`](docs/phase1_ceiling.md).
  Headline: of 499 testable verified effectors, the ±3-gene rule could reach 18% (T1SS 80%,
  T3SS 21%, T6SS 22%, T4SS 6%, T2SS 1%). Tables in `data/phase1/`, figures in `figures/`.
