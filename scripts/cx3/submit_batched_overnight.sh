#!/bin/bash
# Submit ONE PBS job that runs ssign in batched multi-genome mode
# (pooled predictions + annotations across N genomes via MultiGenomeRunner).
#
# Usage:
#   # Default — 2-genome smoke test (K-12 + PAO1) on RTX6000:
#   bash scripts/cx3/submit_batched_overnight.sh
#
#   # All 4 tutorial genomes (K-12 + PAO1 + Vc + Salm):
#   bash scripts/cx3/submit_batched_overnight.sh --tutorial-all
#
#   # Specific genomes:
#   bash scripts/cx3/submit_batched_overnight.sh /path/g1.gbff /path/g2.gbff
#
#   # Override GPU type:
#   bash scripts/cx3/submit_batched_overnight.sh --gpu L40S
#
# The 2-genome default is sized for fast verification runs — K-12 and
# PAO1 differ enough (4314 vs 5680 proteins, different SS profiles) to
# exercise per-genome / pool / split routing without paying the
# annotation cost of 4 genomes.
#
# RTX6000 is the default because it places reliably on Imperial CX3
# v1_gpu72 at the 64-core / 120-GB spec — L40S and others have failed
# the placement-set check at that size as of 2026-06-04.

set -eu

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PBS_SCRIPT="$SCRIPT_DIR/run_batched_multi.pbs"
test -f "$PBS_SCRIPT" || { echo "FATAL: $PBS_SCRIPT not found"; exit 1; }

GPU="RTX6000"
WALLTIME="8:00:00"
USE_TUTORIAL_ALL="0"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --gpu) GPU="$2"; shift 2 ;;
        --walltime) WALLTIME="$2"; shift 2 ;;
        --tutorial-all) USE_TUTORIAL_ALL="1"; shift ;;
        -h|--help) sed -n '2,25p' "$0" | sed 's/^# \?//'; exit 0 ;;
        --) shift; break ;;
        -*) echo "Unknown flag: $1"; exit 1 ;;
        *) break ;;
    esac
done

if [ "$#" -gt 0 ]; then
    GENOMES=("$@")
elif [ "$USE_TUTORIAL_ALL" = "1" ]; then
    GENOMES=(
        "$HOME/ssign-tutorial/ecoli_k12.gbff"
        "$HOME/ssign-tutorial/pseudomonas_pao1.gbff"
        "$HOME/ssign-tutorial/salmonella_typhimurium_lt2.gbff"
        "$HOME/ssign-tutorial/vibrio_cholerae_n16961.gbff"
    )
else
    GENOMES=(
        "$HOME/ssign-tutorial/ecoli_k12.gbff"
        "$HOME/ssign-tutorial/pseudomonas_pao1.gbff"
    )
fi

for g in "${GENOMES[@]}"; do
    test -f "$g" || { echo "FATAL: genome not found: $g"; exit 1; }
done

# Build colon-separated INPUT_GBFFS in one safe assignment — no terminal
# wrap during the qsub line gets to mangle a long quoted string.
GBFFS=""
for g in "${GENOMES[@]}"; do
    GBFFS="${GBFFS}${GBFFS:+:}${g}"
done

echo "Submitting batched ssign job:"
echo "  GPU: $GPU"
echo "  walltime: $WALLTIME"
echo "  ${#GENOMES[@]} genomes:"
for g in "${GENOMES[@]}"; do echo "    $g"; done
echo

jid=$(qsub \
    -l "select=1:ncpus=64:mem=120gb:ngpus=1:gpu_type=$GPU" \
    -l "walltime=$WALLTIME" \
    -N "ssign_batched_${#GENOMES[@]}genomes" \
    -v "INPUT_GBFFS=${GBFFS},GPU_TYPE=${GPU},SSIGN_EXTRA_ARGS=${SSIGN_EXTRA_ARGS:-}" \
    "$PBS_SCRIPT")
echo "Submitted: $jid"
echo
echo "Watch with:"
echo "  qstat -T -u \$USER"
echo "Output dir (after start):"
echo "  \$HOME/runs/batched_${GPU}_<timestamp>/"
