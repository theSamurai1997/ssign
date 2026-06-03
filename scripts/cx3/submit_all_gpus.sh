#!/bin/bash
# Submit one ssign K-12 validation job per available GPU type on CX3.
#
# Whoever starts first wins — they're racing for the same input, the
# rest can be qdel'd once one is running. Each job writes to its own
# RUN_DIR (k12_<GPU_TYPE>_<timestamp>/) so outputs don't clobber.
#
# Side benefit: we get real per-GPU runtime calibration data for
# future scheduling decisions.
#
# Usage:
#   bash scripts/cx3/submit_all_gpus.sh
#   bash scripts/cx3/submit_all_gpus.sh L40S A40 RTX6000   # explicit list
#
# Without arguments, auto-detects every gpu_type currently visible to
# pbsnodes. With arguments, submits exactly those types.

set -eu

PBS_SCRIPT="$(dirname "$(readlink -f "$0")")/run_k12_validation.pbs"
test -f "$PBS_SCRIPT" || { echo "FATAL: $PBS_SCRIPT not found"; exit 1; }

if [ "$#" -gt 0 ]; then
    GPU_TYPES=("$@")
else
    # Auto-detect: every distinct gpu_type advertised by pbsnodes.
    mapfile -t GPU_TYPES < <(pbsnodes -a 2>/dev/null \
        | awk -F' = ' '/gpu_type/ {print $2}' \
        | sort -u)
fi

if [ "${#GPU_TYPES[@]}" -eq 0 ]; then
    echo "FATAL: no GPU types detected. Pass them explicitly:"
    echo "    bash $0 L40S A40 RTX6000"
    exit 1
fi

echo "Submitting one job per GPU type: ${GPU_TYPES[*]}"
echo

for gpu in "${GPU_TYPES[@]}"; do
    job_name="ssign_k12_${gpu}"
    # -l select=... on the qsub command line overrides the script's
    #   #PBS -l select=... directive.
    # -N renames the job so qstat output is readable.
    # -v passes GPU_TYPE into the job environment — the script reads it
    #   to build a per-GPU RUN_DIR.
    jid=$(qsub \
        -l "select=1:ncpus=16:mem=80gb:ngpus=1:gpu_type=${gpu}" \
        -N "${job_name}" \
        -v "GPU_TYPE=${gpu}" \
        "$PBS_SCRIPT")
    echo "  ${gpu}: ${jid}"
done

echo
echo "Check status with:"
echo "    qstat -T -u \$USER"
