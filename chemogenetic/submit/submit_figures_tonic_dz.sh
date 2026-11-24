#!/bin/bash
TIME="${ZWBA_TIME:-3:00:00}"   # walltime; override per call: ZWBA_TIME=48:00:00 bash ...
# submit_figures_tonic_dz.sh
# Usage (from repo root):
#     bash chemogenetic/submit/submit_figures_tonic_dz.sh --config config_hcrt_trpv1_csn_120min

CONFIG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."; exit 1
fi

PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CONFIG}.py (run from the repo root; check the config name)"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

LOG_DIR="logs/figures_tonic_dz"
mkdir -p "$LOG_DIR"

echo "Submitting tonic ΔZ figures for config: $CONFIG"

sbatch \
    --job-name="figdz_${CONFIG}" \
    --cpus-per-task=4 \
    --mem=32G \
    --output="${LOG_DIR}/${CONFIG}.log" \
    --time="$TIME" \
    --wrap="$PYTHON chemogenetic/run/run_figures_tonic_dz.py --config $CONFIG" || { echo "ERROR: sbatch failed (message above). Nothing was submitted."; exit 1; }

echo "Submitted."
echo "Log: ${LOG_DIR}/${CONFIG}.log"
