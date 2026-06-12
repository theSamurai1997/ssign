#!/bin/bash
# Phase 2 submit plan. Defaults (INPUT_DIR, INPUT_DIR_GB, SSIGN_DB, SSIGN_VENV) come
# from run_benchmark_batch.pbs; override with -v if your CX3 layout differs.
# gpu72 enforces ~2 concurrent jobs/user -> submit in pairs, not all at once.
set -eu
B="$HOME/bench/batches"

# 1) PILOT (4 genomes) — DONE 2026-06-11 (harness validated; GenBank==FASTA recall;
#    T3SS pass recovers real T3SS effectors). Kept for reference; commented out.
# qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=genbank,RUN_TAG=pilot_genbank_default $HOME/bench/run_benchmark_batch.pbs
# qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=fasta,RUN_TAG=pilot_fasta_default $HOME/bench/run_benchmark_batch.pbs
# qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=pilot_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs

# 2) FULL PANEL (67 genomes, 6 batches). Pilot conclusions (2026-06-11):
#    - GenBank vs FASTA: identical emitted set -> run GenBank only (FASTA bridge verified but
#      adds no recall signal; skip it to save GPU).
#    - T3SS: default mode EXCLUDES T3SS, so all 227 testable T3SS effectors auto-miss. The
#      T3SS-included pass is REQUIRED to measure T3SS recall; the default pass is the
#      conservative headline (and shows the DeepSecE-flagellar-FP cost of including T3SS).
#    gpu72 fits ~12 GPUs/user, so each 6-batch wave can run concurrently (~2-3h/batch, 8h cap).
#
#    Wave A — GenBank, T3SS excluded (default, conservative headline):
qsub -v BATCH_FILE=$B/batch_01.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_02.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_03.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_04.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_05.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_06.txt,INPUT_MODE=genbank,RUN_TAG=panel_genbank_default $HOME/bench/run_benchmark_batch.pbs

#    Wave B — GenBank, T3SS included (needed for T3SS recall). Submit after Wave A lands,
#    or together if the queue allows (12 jobs).
qsub -v BATCH_FILE=$B/batch_01.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_02.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_03.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_04.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_05.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
qsub -v BATCH_FILE=$B/batch_06.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=panel_genbank_t3ss $HOME/bench/run_benchmark_batch.pbs
