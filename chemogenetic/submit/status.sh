#!/bin/bash
# status.sh
# =========
# Show which pipeline stages are done for every fish in a config. Read-only: runs on the
# login node in seconds, submits nothing, creates nothing.
#
# Usage (from repo root):
#     bash chemogenetic/submit/status.sh --config config_hcrt_trpv1_csn_120min
#     bash chemogenetic/submit/status.sh --config config_<exp> --fish fish3     # filter by expt_ID text
#     bash chemogenetic/submit/status.sh --config config_<exp> --strict         # exit 1 if anything missing/stale

CONFIG=""
EXTRA=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        *) EXTRA+=("$1"); shift ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."
    echo "Usage: bash chemogenetic/submit/status.sh --config config_hcrt_trpv1_csn_120min"
    exit 1
fi

PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CONFIG}.py (run from the repo root; check the config name)"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

"$PYTHON" chemogenetic/run/run_status.py --config "$CONFIG" "${EXTRA[@]}"
