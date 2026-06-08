# Validation sweeps

Behavioral testing runs that probe how ssign's substrate predictions
respond to inputs and parameter changes. Each subdirectory holds one
experimental axis: the raw outputs collected per run, plus the analysis
that aggregates across runs.

## Experiments (2026-06-07 / 2026-06-08 batch)

| Axis | What varies | Values tested | N runs | Genomes |
|---|---|---|---|---|
| `proximity_window/` | `proximity_window` setting (genes flanking each SS component scanned for substrates) | 1, 3, 5, 7, 10, 15 | 30 | E. coli K-12, H. pylori, L. pneumophila, B. pertussis, P. aeruginosa PAO1 |
| `fragmentation/` | Number of artificial contigs the input genome is chopped into via `scripts/chop_genome.py` | 1, 3, 10, 30, 70, 200, 1000 | 35 | same 5 genomes |
| `copies/` | Number of identical copies of the genome submitted in one batched-mode job | 1, 2, 3, 4, 6, 8, 12, 16 | 40 | same 5 genomes |
| `benchmarking/` | Default proximity filter vs. all 4 whole-genome flags enabled, scored against literature-curated effector lists | 2 conditions | 12 | C. burnetii RSA493, L. pneumophila, S. Typhimurium LT2, P. aeruginosa PAO1, V. cholerae N16961, Y. pestis CO92 |

**Total: 117 runs.**

## What each experiment tests

### proximity_window
Default is `proximity_window=3`. Does ssign's recall and precision
shift smoothly across +/-N gene windows, or is there an optimum that
trades a small precision loss for a meaningful recall gain? The
upper end of the sweep (N=15) approximates "permissive proximity";
N=1 is "almost no proximity filter at all."

### fragmentation
Real-world genomes from short-read assemblers come in many contigs,
not as single chromosomes. ssign's proximity rule is "same contig
only," so heavy fragmentation should erode substrate calls (substrates
get separated from their secretion-system components by contig
breaks). This sweep quantifies that erosion across 7 fragmentation
levels using the controlled chop utility in `scripts/chop_genome.py`.
N=1 is the intact chromosome (baseline); N=1000 is heavy
fragmentation (~5 kb median contig for K-12).

### copies
ssign's batched mode pools per-protein predictions (DeepSecE,
PLM-Effector, SignalP, DeepLocPro, etc.) across the input genome set.
Submitting N identical copies of the same genome should produce N
identical per-genome substrate lists. Any drift across copies reveals
non-determinism in the pooling logic, a cross-genome leakage bug, or
an unintended scaling effect.

### benchmarking
Six well-characterized genomes have literature-curated effector lists
that act as ground truth (see `validation_ground_truth/`). For each,
ssign runs twice: once with the default per-component proximity
filter, once with all four whole-genome flags enabled
(`--dlp-whole-genome --dse-whole-genome --signalp-whole-genome
--plme-whole-genome`). The two outputs are scored against the
ground-truth substrate lists. The recall delta quantifies how much
real biology the proximity filter is missing.

## How runs were submitted

All runs went to Imperial CX3 v1_gpu72, RTX6000 GPUs, 32 cores /
64 GB RAM. Standard placement recipe (see project memory
`cx3_gpu72_placement.md`).

Run-directory naming on CX3:

- Single-axis sweeps:
  `<genome>_<axis><value>_<gpu>_<datetime>_<jobid>/`
  e.g. `pseudomonas_pao1_pw5_RTX6000_20260608_011334_2944821/`.
- Copies sweep: `batched_RTX6000_<datetime>_<jobid>/`. Inside each
  batched dir the per-copy outputs live at `<genome>_cal1/`,
  `cal2/`, ..., `cal{N}/`. The "cal" suffix carries no semantic
  meaning (it just numbers the duplicated input copies).
- Benchmarking: `<genome>_<gpu>_<datetime>_<jobid>/` paired (lower
  jobid = proximity filter; higher jobid = whole-genome flags
  enabled). PBS jobids 2946811-2946822.

## What lands in each subdir

Each per-axis subdir contains:
- `<run_dir>/<sample>_substrates_filtered.tsv` (the headline result)
- `<run_dir>/<sample>_results.csv` (full result table)
- `<run_dir>/<sample>_summary.txt` (plain-text per-genome summary)
- `<run_dir>/ssign.run.log` (run trace for sanity checks)
- For batched runs: `<run_dir>/combined_results.csv` and `_pool/`
  pooled-annotation outputs.

The multi-GB intermediates (Bakta, EggNOG, HHpred, runtime_data,
.ssign) are not collected here, only the small result files.

## Analysis

`analysis/` holds the aggregation scripts and figures derived from the
raw outputs above. The first script reads each axis subdir, joins
substrate calls across runs of the same genome, and emits per-axis
sensitivity plots. Figure index lands in `analysis/figure_index.md`
once generated.

## Raw archive

The initial result archive captured on 2026-06-08 is
`mg_sweep_2026-06-08.tar.gz` (downloaded from
`ttr25@login.cx3.hpc.imperial.ac.uk:~/runs/`). The tar holds only the
small files; multi-GB intermediates stay on CX3.

## Reproducing the runs

The submission scripts that fired this batch were one-offs (not
committed). The reproducible submission path going forward:

1. Edit the genome list and parameter range at the top of a new
   `scripts/cx3/submit_validation_sweep.sh` (to be added when we next
   run the sweep).
2. The script loops over (genome, value) pairs, exporting parameter
   overrides via `SSIGN_EXTRA_ARGS`, and submits one PBS job per pair
   using the existing `scripts/cx3/run_k12_validation.pbs` template.
3. Each job writes to its own RUN_DIR per the naming convention above.

Until that script lands, the design recorded in this README is the
authoritative reference for what the 2026-06-08 batch contained.
