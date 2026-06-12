# Phase 2 actual-recall: CX3 run instructions

Runs ssign on the 67-genome benchmark panel to measure what it actually emits as secreted,
for comparison against the Phase 1 ceiling. You submit; I analyse the results.

## Config (fixed in run_benchmark_batch.pbs)

- `--tier extended`, secretion predictors ON (DeepLocPro + SignalP + DeepSecE + PLM-Effector),
  run **local** (GPU) because the compute nodes are offline and the remote DTU fallback can't
  reach the internet.
- annotation OFF (`--skip-eggnog/interproscan/plmblast/protparam`) — runs after the secreted
  call, can't change recall, saves hours.
- DB env auto-set to the flat CX3 layout (`$EPHEMERAL/ssign-databases/{bakta,plm_effector_weights}`),
  local predictors from `~/.conda/envs/{signalp6,deeplocpro}`, DeepSecE checkpoint from
  `~/.ssign/models/`. The job prints a WARN (not fatal) if any is missing.

## Two input modes (we compare them on the pilot)

- `INPUT_MODE=genbank` (default): `inputs_gb/<unit>.gbff` + `--use-input-annotations`. ssign uses
  the RefSeq CDS, so its locus_tags ARE the gold-set tags and recall is in the same gene order
  as the ceiling. No Bakta, no bridge.
- `INPUT_MODE=fasta`: `inputs/<unit>.fasta`. ssign re-annotates with Bakta; results bridged back
  to RefSeq by coordinate overlap (task 6.3).

## One-time copy laptop -> CX3

```bash
rsync -av validation_sweeps/benchmark/inputs/    ttr25@login.cx3.hpc.imperial.ac.uk:bench/inputs/
rsync -av validation_sweeps/benchmark/inputs_gb/ ttr25@login.cx3.hpc.imperial.ac.uk:bench/inputs_gb/
rsync -av validation_sweeps/benchmark/scripts/cx3/ ttr25@login.cx3.hpc.imperial.ac.uk:bench/
```
After this CX3 has `~/bench/{inputs,inputs_gb,batches}/` and `~/bench/run_benchmark_batch.pbs`.

## Submit (gpu72 ~2-job cap -> pairs)

Pilot first — see `batches/SUBMIT.sh`. Three jobs: GenBank+default, FASTA+default,
GenBank+T3SS-included. Send back `~/runs/benchmark_phase2/pilot_*/`. We compare:
GenBank vs Bakta (does re-annotation change recall?) and default vs T3SS-included (real
effectors recovered vs DeepSecE flagellar false positives). Then we lock both and run the
full panel (the commented lines in SUBMIT.sh).

## What to hand back

Per genome: `~/runs/benchmark_phase2/<RUN_TAG>/<unit>/results/<unit>_results.csv` (emitted
secreted) + `<unit>_results_raw.csv` (every protein, locus_tag+contig+start+end+strand) +
`ssign.run.log`. Simplest: `rsync` the whole `benchmark_phase2/` tree back.

## Notes

- DBs live in `ephemeral` (purges after ~30 days idle); fine for a run this week. They are all
  present as of 2026-06-11 (bakta 4 GB, plm_effector_weights 26 GB, taxdump intact).
- If jobs pend or land on the wrong GPU, pin it: add `:gpu_type=RTX6000` to the qsub
  `-l select=...` line (gpu72 sometimes needs it).
