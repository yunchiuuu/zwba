#!/bin/bash
# submit_glm.sh
# =============
# Submit GLM pipeline for one config.
# Runs: refit → ablation → IAAFT null → responder idx
# Uses empirical HCRT regressor if hcrt_all.csv files exist.
#
# Usage (from repo root):
#     bash chemogenetic/submit/submit_glm.sh --config config_hcrt_trpv1_csn_120min

CONFIG=""
EXTRA=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        *) EXTRA+=("$1"); shift ;;   # passed through, e.g. --stages responders --overwrite null
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."
    echo "Usage: bash chemogenetic/submit/submit_glm.sh --config config_hcrt_trpv1_csn_120min"
    exit 1
fi

PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CONFIG}.py (run from the repo root; check the config name)"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

LOG_DIR="logs/glm"
mkdir -p "$LOG_DIR"

echo "Submitting GLM job for config: $CONFIG"

sbatch \
    --job-name="glm_${CONFIG}" \
    --cpus-per-task=32 \
    --mem=256G \
    --output="${LOG_DIR}/${CONFIG}.log" \
    --wrap="$PYTHON chemogenetic/run/run_glm.py --config $CONFIG ${EXTRA[*]}"

echo "Submitted."
echo "Monitor: squeue -u \$USER"
echo "Log:     ${LOG_DIR}/${CONFIG}.log"
