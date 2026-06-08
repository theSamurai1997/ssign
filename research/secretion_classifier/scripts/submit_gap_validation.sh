#!/bin/bash
# Gap-validation experiment: run ssign on 6 well-characterized genomes
# with proximity filter (default) and with all four whole-genome flags
# enabled. Submit on CX3 v1_gpu72.
#
# Run from $HOME on CX3 after pulling the latest ssign main:
#   bash submit_gap_validation.sh
#
# Output: 12 jobs queued, run dirs at $HOME/runs/gap_<sample>_<config>_<jobid>/

set -eo pipefail

TUT="$HOME/ssign-tutorial"
GENOME_DIR="$HOME/gap_validation_genomes"
mkdir -p "$GENOME_DIR"

# ---- 1. Stage genomes ----
# Use what we already have where possible; otherwise fetch from NCBI.
declare -A GENOMES=(
    [legionella_pneumophila]="$TUT/legionella_pneumophila.gbff"
    [pseudomonas_pao1]="$TUT/pseudomonas_pao1.gbff"
    [salmonella_lt2]="$GENOME_DIR/salmonella_lt2.gbff"
    [coxiella_rsa493]="$GENOME_DIR/coxiella_rsa493.gbff"
    [yersinia_pestis_co92]="$GENOME_DIR/yersinia_pestis_co92.gbff"
    [vibrio_cholerae_n16961]="$GENOME_DIR/vibrio_cholerae_n16961.gbff"
)

declare -A REFSEQ=(
    [salmonella_lt2]="NC_003197"
    [coxiella_rsa493]="NC_002971"
    [yersinia_pestis_co92]="NC_003143"
    [vibrio_cholerae_n16961]="NC_002505,NC_002506"
)

# Fetch missing genomes via NCBI E-utilities
for sample in "${!REFSEQ[@]}"; do
    out="${GENOMES[$sample]}"
    if [[ ! -f "$out" ]]; then
        echo "Fetching $sample ($(echo ${REFSEQ[$sample]}))"
        ids="${REFSEQ[$sample]}"
        # If multi-replicon, fetch each and concatenate
        > "$out"
        for id in ${ids//,/ }; do
            curl -fsSL \
              "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id}&rettype=gbwithparts&retmode=text" \
              >> "$out"
            sleep 1  # be nice to NCBI
        done
        test -s "$out" || { echo "FATAL: failed to fetch $sample"; exit 1; }
    fi
done

# Verify all files exist
for sample in "${!GENOMES[@]}"; do
    test -f "${GENOMES[$sample]}" || { echo "FATAL: missing ${GENOMES[$sample]}"; exit 1; }
done
echo "All 6 genomes ready."

# ---- 2. Submit 12 jobs (6 genomes × 2 configs) ----
PBS_SCRIPT="$HOME/ssign/scripts/cx3/run_k12_validation.pbs"
test -f "$PBS_SCRIPT" || { echo "FATAL: ssign repo not at \$HOME/ssign"; exit 1; }

# Skip every optional annotation tool: we only need the substrate-prediction
# core (Bakta + MacSyFinder + DLP + DSE + SignalP + PLM-E + proximity filter)
# to measure recall against ground-truth substrate lists. Cuts ~45 min runs
# to ~15 min and frees pLM-BLAST GPU contention.
SKIP_ARGS="--skip-blastp --skip-interproscan --skip-hhsuite --skip-plmblast --skip-eggnog --skip-protparam"

for sample in "${!GENOMES[@]}"; do
    input="${GENOMES[$sample]}"

    # Config A: default proximity filter
    export SSIGN_EXTRA_ARGS="$SKIP_ARGS"
    qsub \
        -v "INPUT_GBFF=$input,SSIGN_EXTRA_ARGS" \
        -N "gap_${sample}_proxy" \
        -l select=1:ncpus=32:mem=64gb:ngpus=1:gpu_type=RTX6000 \
        -l walltime=04:00:00 \
        "$PBS_SCRIPT"

    # Config B: all four whole-genome flags
    export SSIGN_EXTRA_ARGS="$SKIP_ARGS --dlp-whole-genome --dse-whole-genome --sp-whole-genome --plme-whole-genome"
    qsub \
        -v "INPUT_GBFF=$input,SSIGN_EXTRA_ARGS" \
        -N "gap_${sample}_whole" \
        -l select=1:ncpus=32:mem=64gb:ngpus=1:gpu_type=RTX6000 \
        -l walltime=04:00:00 \
        "$PBS_SCRIPT"

    echo "Submitted: $sample (proxy + whole)"
done

echo "All 12 jobs queued. Check with: qstat -u \$USER"
