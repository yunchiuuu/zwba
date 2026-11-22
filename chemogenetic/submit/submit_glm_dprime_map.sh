#!/bin/bash
# submit_glm_dprime_map.sh
# ========================
# Usage (from repo root):
#     bash chemogenetic/submit/submit_glm_dprime_map.sh --config config_hcrt_trpv1_csn_120min

SCRIPT="chemogenetic/run/run_glm_dprime_map.py"
LOG_DIR="logs/brainmap"

CONFIG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."
    echo "Usage: bash chemogenetic/submit/submit_glm_dprime_map.sh --config config_hcrt_trpv1_csn_120min"
    exit 1
fi

# PYTHON_BIN comes from the config (same as every other submit script)
PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CONFIG}.py (run from the repo root; check the config name)"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

mkdir -p "$LOG_DIR"

sbatch \
    --job-name="glm_dprime_map_${CONFIG}" \
    --ntasks=1 \
    --cpus-per-task=2 \
    --mem=16G \
    --time=1:00:00 \
    --output="${LOG_DIR}/${CONFIG}.log" \
    --wrap="$PYTHON $SCRIPT --config $CONFIG" || { echo "ERROR: sbatch failed (message above). Nothing was submitted."; exit 1; }

echo "Submitted."
echo "Monitor: squeue -u \$USER"
echo "Log:     ${LOG_DIR}/${CONFIG}.log"
