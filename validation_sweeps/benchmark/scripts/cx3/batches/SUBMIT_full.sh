#!/bin/bash
# Submit plan for the 2026-06-16 "full" panel run:
#   T3SS ON + whole-genome DLP/DSE/PLM-E (SignalP neighborhood only) + enrichment stats + tier-2 annotation
#   (InterProScan, EggNOG, pLM-BLAST, ProtParam). BLASTp + HH-suite are tier-full, left OFF.
#
# gpu72: explicit gpu_type=RTX6000 is mandatory (else "Placement set too small", queues forever).
# 12 GPUs/user cap -> all 6 panel batches can run concurrently. 72h walltime covers the heavy
# whole-genome PLM-E pass + annotation; the smoke run calibrates the real per-genome wallclock.
set -eu
B="$HOME/bench/batches"
PBS="$HOME/bench/run_benchmark_batch.pbs"
RES="select=1:ncpus=16:mem=64gb:ngpus=1:gpu_type=RTX6000"
VARS="INCLUDE_T3SS=1,WHOLE_GENOME=1,ENRICH=1,ANNOT=1,RUN_TAG=full_t3ss"

# ---- 0. PRE-FLIGHT: confirm every local install + DB path exists before burning GPU hours ----
#   bash $B/../preflight_full.sh     # (writes WARN lines for anything missing)

# ---- 1. SMOKE: one annotation-rich genome (PAO1) first. Inspect the log end-to-end before step 2.
qsub -l "$RES" -l walltime=12:00:00 -v "BATCH_FILE=$B/smoke.txt,$VARS" "$PBS"

# ---- 2. FULL PANEL: 6 batches (67 genomes), concurrent. Run ONLY after the smoke log is clean.
# for b in 01 02 03 04 05 06; do
#     qsub -l "$RES" -l walltime=72:00:00 -v "BATCH_FILE=$B/batch_$b.txt,$VARS" "$PBS"
# done
