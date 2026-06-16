#!/bin/bash
# Run ON CX3 (login node) before submitting the full panel. Verifies every install + DB the
# WHOLE_GENOME + ANNOT run needs. Pure checks, no GPU, no submission. Override any path via env
# (same names as run_benchmark_batch.pbs). Exits 0 if all green, 1 if anything is missing.
set -u
: "${SSIGN_DB:=$EPHEMERAL/ssign-databases}"
VENV="${SSIGN_VENV:-$HOME/blastp_t5a/ssign/.venv}"
: "${INPUT_DIR_GB:=$HOME/bench/inputs_gb}"
: "${IPS_DIR:=$SSIGN_DB/interproscan/interproscan-5.77-108.0}"
: "${ECOD70_DB:=$SSIGN_DB/plm_blast/ECOD70}"
: "${EGGNOG_DIR:=$SSIGN_DB/eggnog}"
: "${PLMBLAST_SCRIPT:=$HOME/build/pLM-BLAST/scripts/plmblast.py}"

bad=0
ck() {  # ck <test-expr> <path> <label>
    if eval "$1 \"$2\""; then printf '  OK   %-26s %s\n' "$3" "$2"; else printf '  MISS %-26s %s\n' "$3" "$2"; bad=1; fi
}
echo "=== ssign + predictors ==="
ck "test -f" "$VENV/bin/activate" "ssign venv"
ck "test -x" "$HOME/.conda/envs/signalp6/bin/signalp6" "SignalP6 (local)"
ck "test -x" "$HOME/.conda/envs/deeplocpro/bin/deeplocpro" "DeepLocPro (local)"
ck "test -f" "$HOME/.ssign/models/deepsece_checkpoint.pt" "DeepSecE checkpoint"
ck "test -d" "$SSIGN_DB/plm_effector_weights" "PLM-Effector weights"
echo "=== tier-2 annotation ==="
ck "test -x" "$IPS_DIR/interproscan.sh" "InterProScan"
ck "test -d" "$EGGNOG_DIR" "EggNOG data dir"
ck "test -x" "$HOME/.conda/envs/eggnog/bin/emapper.py" "EggNOG emapper"
ck "test -d" "$ECOD70_DB" "pLM-BLAST ECOD70"
ck "test -f" "$PLMBLAST_SCRIPT" "pLM-BLAST script"
echo "=== inputs ==="
ck "test -d" "$INPUT_DIR_GB" "GenBank inputs dir"
n=$(ls "$INPUT_DIR_GB"/*.gbff 2>/dev/null | wc -l); echo "  ..  $n .gbff inputs found (expect 67)"
echo "=== GPU + torch (run inside the venv) ==="
echo "  (on a GPU node: source $VENV/bin/activate; python -c 'import torch;print(torch.cuda.is_available())')"

[ "$bad" = "0" ] && echo "ALL GREEN — safe to submit." || echo "MISSING ITEMS ABOVE — fix before submitting."
exit $bad
