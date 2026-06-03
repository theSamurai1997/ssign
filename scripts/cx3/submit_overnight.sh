#!/bin/bash
# Submit overnight ssign validation runs across genomes × GPU types.
#
# Each job is given 16 h walltime so pLM-BLAST has room to finish on
# autotransporter-heavy genomes (it scales linearly with query length,
# and T5SS-rich genomes can need 2-4 h+ for the ECOD70 search step —
# see tasks #46 and #4). Each job writes to its own RUN_DIR so outputs
# don't clobber.
#
# Usage:
#   # All .gbff/.gbk/.gb files in ~/ssign-tutorial/, every visible GPU type:
#   bash scripts/cx3/submit_overnight.sh
#
#   # Explicit genome list:
#   bash scripts/cx3/submit_overnight.sh ~/genomes/pao1.gbff ~/genomes/typhi.gbff
#
#   # Restrict to a subset of GPUs:
#   bash scripts/cx3/submit_overnight.sh --gpu L40S,A100
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
GPU_FILTER=""

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
            sed -n '2,30p' "$0" | sed 's/^# \?//'
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

# GPU types: --gpu filter (comma-separated) or auto-discover via pbsnodes.
if [ -n "$GPU_FILTER" ]; then
    IFS=',' read -ra GPU_TYPES <<<"$GPU_FILTER"
else
    mapfile -t GPU_TYPES < <(pbsnodes -a 2>/dev/null | awk -f "$AWK_SCRIPT" | sort -u)
fi

if [ "${#GPU_TYPES[@]}" -eq 0 ]; then
    echo "FATAL: no GPU types detected."
    echo "Pass them explicitly:   bash $0 --gpu L40S,A100"
    exit 1
fi

N_GENOMES=${#GENOMES[@]}
N_GPUS=${#GPU_TYPES[@]}
echo "Overnight submission: $N_GENOMES genome(s) x $N_GPUS GPU type(s) = $((N_GENOMES * N_GPUS)) jobs"
echo "Walltime: $WALLTIME"
echo "Genomes:"
for g in "${GENOMES[@]}"; do echo "  $g"; done
echo "GPUs: ${GPU_TYPES[*]}"
echo

for genome in "${GENOMES[@]}"; do
    tag=$(basename "$genome" | sed 's/\.[^.]*$//')
    for gpu in "${GPU_TYPES[@]}"; do
        jid=$(qsub \
            -l "select=1:ncpus=16:mem=80gb:ngpus=1:gpu_type=${gpu}" \
            -l "walltime=${WALLTIME}" \
            -N "ssign_${tag}_${gpu}" \
            -v "GPU_TYPE=${gpu},INPUT_GBFF=${genome}" \
            "$PBS_SCRIPT")
        echo "  ${tag} on ${gpu}: ${jid}"
    done
done

echo
echo "Check status:"
echo "    qstat -T -u \$USER"
