#!/bin/bash
# submit_survey_by_regions.sh
# ============================
# Submits one sbatch job per config for the regional survey pipeline.
#
# Usage (from repo root):
#     bash chemogenetic/submit/submit_survey_by_regions.sh --config config_<exp>     # one config
#     bash chemogenetic/submit/submit_survey_by_regions.sh                           # the default list below
#     bash chemogenetic/submit/submit_survey_by_regions.sh --config config_<exp> --overwrite --top_n 20
# (any other arguments are passed through to run_survey_by_regions.py)

SCRIPT="chemogenetic/run/run_survey_by_regions.py"
LOG_DIR="logs/chemogenetic"

CONFIG=""
EXTRA=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        *) EXTRA+=("$1"); shift ;;
    esac
done

if [ -n "$CONFIG" ]; then
    CONFIGS=("$CONFIG")
else
    # default list when --config is not given: edit for your project
    CONFIGS=(
        "config_hcrt_trpv1_csn_120min"
        "config_hcrt_trpv1_inj_csn_120min"
        "config_ynt185_120min"
    )
fi

mkdir -p "$LOG_DIR"

for CFG in "${CONFIGS[@]}"; do
    PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CFG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CFG}.py (run from the repo root; check the config name)"; exit 1; }
    [ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

    sbatch \
        --job-name="survey_by_regions_${CFG}" \
        --ntasks=1 \
        --cpus-per-task=25 \
        --mem=256G \
        --output="${LOG_DIR}/survey_by_regions_${CFG}.log" \
        --wrap="$PYTHON $SCRIPT --config $CFG ${EXTRA[*]}"

    echo "  Submitted: $CFG"
done

echo ""
echo "All ${#CONFIGS[@]} jobs submitted."
echo "Monitor: squeue -u \$USER"
echo "Logs in: $LOG_DIR/"
