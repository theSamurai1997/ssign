#!/bin/bash
# Submit overnight ssign validation runs.
#
# Default behaviour: one job per genome, NO gpu_type constraint — the PBS
# scheduler picks whichever GPU is free first. Per Imperial RCS guidance
# (https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/queues/gpu-jobs/):
# "It is recommended leaving this option empty unless you have a specific
# need for one or the other." Constraining gpu_type made our 2026-06-03
# overnight queue stall because 11 of 16 jobs sat waiting for specific GPU
# types that never freed up.
#
# Walltime defaults to 16h so pLM-BLAST has room to finish on
# autotransporter-heavy genomes (see tasks #46, #4, #77). Each job writes
# to its own RUN_DIR so outputs don't clobber.
#
# Usage:
#   # One job per .gbff in ~/ssign-tutorial/, scheduler picks any GPU:
#   bash scripts/cx3/submit_overnight.sh
#
#   # Explicit genome list:
#   bash scripts/cx3/submit_overnight.sh ~/genomes/pao1.gbff ~/genomes/typhi.gbff
#
#   # Race across specific GPU types (one job per genome per GPU type) —
#   # use for calibration/benchmarking, not for routine submissions:
#   bash scripts/cx3/submit_overnight.sh --gpu L40S,A100
#
#   # Race across every visible GPU type:
#   bash scripts/cx3/submit_overnight.sh --gpu auto
#
#   # Different walltime:
#   bash scripts/cx3/submit_overnight.sh --walltime 24:00:00
#
# After submission:
#   qstat -T -u $USER

set -eu

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PBS_SCRIPT="$SCRIPT_DIR/run_k12_validation.pbs"
AWK_SCRIPT="$SCRIPT_DIR/_list_usable_gpu_types.awk"
test -f "$PBS_SCRIPT" || { echo "FATAL: $PBS_SCRIPT not found"; exit 1; }
test -f "$AWK_SCRIPT" || { echo "FATAL: $AWK_SCRIPT not found"; exit 1; }

WALLTIME="16:00:00"
GPU_FILTER=""  # empty = no gpu_type spec (scheduler picks any)

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --gpu)
            GPU_FILTER="$2"
            shift 2
            ;;
        --walltime)
            WALLTIME="$2"
            shift 2
            ;;
        -h|--help)
            sed -n '2,32p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown flag: $1"
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

# Genomes: positional args, or auto-discover from ~/ssign-tutorial/.
if [ "$#" -gt 0 ]; then
    GENOMES=("$@")
else
    mapfile -t GENOMES < <(
        ls -1 \
            "$HOME"/ssign-tutorial/*.gbff \
            "$HOME"/ssign-tutorial/*.gbk \
            "$HOME"/ssign-tutorial/*.gb \
            2>/dev/null
    )
fi

if [ "${#GENOMES[@]}" -eq 0 ]; then
    echo "FATAL: no genomes to submit."
    echo "Pass them explicitly:   bash $0 path/to/genome.gbff [...]"
    echo "Or add .gbff/.gbk/.gb files under ~/ssign-tutorial/"
    exit 1
fi

for g in "${GENOMES[@]}"; do
    test -f "$g" || { echo "FATAL: genome not found: $g"; exit 1; }
done

# GPU selection logic:
# - empty GPU_FILTER → no gpu_type spec (scheduler picks any visible GPU)
# - "auto"          → discover every visible gpu_type via pbsnodes
# - comma-separated → use those exact types
if [ -z "$GPU_FILTER" ]; then
    GPU_TYPES=("")  # one job, no gpu_type
elif [ "$GPU_FILTER" = "auto" ]; then
    mapfile -t GPU_TYPES < <(pbsnodes -a 2>/dev/null | awk -f "$AWK_SCRIPT" | sort -u)
    if [ "${#GPU_TYPES[@]}" -eq 0 ]; then
        echo "FATAL: --gpu auto but no usable GPU types visible to pbsnodes."
        echo "Pass them explicitly:   bash $0 --gpu L40S,A100"
        exit 1
    fi
else
    IFS=',' read -ra GPU_TYPES <<<"$GPU_FILTER"
fi

N_GENOMES=${#GENOMES[@]}
N_GPUS=${#GPU_TYPES[@]}
N_JOBS=$((N_GENOMES * N_GPUS))
echo "Submitting: $N_GENOMES genome(s) x $N_GPUS GPU type(s) = $N_JOBS jobs"
echo "Walltime: $WALLTIME"
echo "Genomes:"
for g in "${GENOMES[@]}"; do echo "  $g"; done
if [ -z "$GPU_FILTER" ]; then
    echo "GPU: scheduler-assigned (no gpu_type constraint)"
else
    echo "GPUs: ${GPU_TYPES[*]}"
fi
echo

for genome in "${GENOMES[@]}"; do
    tag=$(basename "$genome" | sed 's/\.[^.]*$//')
    for gpu in "${GPU_TYPES[@]}"; do
        if [ -z "$gpu" ]; then
            select_spec="select=1:ncpus=16:mem=80gb:ngpus=1"
            job_name="ssign_${tag}"
            env_vars="INPUT_GBFF=${genome}"
            label="any GPU"
        else
            select_spec="select=1:ncpus=16:mem=80gb:ngpus=1:gpu_type=${gpu}"
            job_name="ssign_${tag}_${gpu}"
            env_vars="GPU_TYPE=${gpu},INPUT_GBFF=${genome}"
            label="${gpu}"
        fi
        jid=$(qsub \
            -l "$select_spec" \
            -l "walltime=${WALLTIME}" \
            -N "$job_name" \
            -v "$env_vars" \
            "$PBS_SCRIPT")
        echo "  ${tag} on ${label}: ${jid}"
    done
done

echo
echo "Check status:"
echo "    qstat -T -u \$USER"
