#!/bin/bash
TIME="${ZWBA_TIME:-2:00:00}"   # walltime; override per call: ZWBA_TIME=48:00:00 bash ...
# submit_temporal_zplanes.sh
# ==========================
# Temporal intensity maps for single z planes (one PDF page per plane), tonic ΔZ and phasic d′.
# Needs the epoch arrays from submit_temporal_intensity_map.sh (current epochs / cell filter).
#
# Usage (from repo root):
#     bash chemogenetic/submit/submit_temporal_zplanes.sh --config config_hcrt_trpv1_csn_120min
#     ... --z-bins 3,6,9,12,15,18 --slab 1        choose the planes (coarse-z bins)
#     ... --z-bins all                            every z bin that holds cells (one page each)
#     ... --n-planes 8                            number of evenly spaced planes (default 6)
#     ... --plot-only                             redraw from the saved volumes (other planes, no recompute)

# ── Parse --config argument ───────────────────────────────────────────────────
CONFIG=""
EXTRA=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        --z-bins|--n-planes|--slab|--n-jobs) EXTRA="$EXTRA $1 $2"; shift 2 ;;
        --plot-only) EXTRA="$EXTRA --plot-only"; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."
    echo "Usage: bash chemogenetic/submit/submit_temporal_zplanes.sh --config config_hcrt_trpv1_csn_120min"
    exit 1
fi

# ── Bootstrap PYTHON_BIN from config ─────────────────────────────────────────
PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'chemogenetic/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not read PYTHON_BIN from chemogenetic/config/${CONFIG}.py (run from the repo root; check the config name)"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

LOG_DIR="logs/temporal_zplanes"
mkdir -p "$LOG_DIR"

echo "Submitting temporal z-plane maps job for config: $CONFIG"

sbatch \
    --job-name="tzp_${CONFIG}" \
    --cpus-per-task=8 \
    --mem=64G \
    --output="${LOG_DIR}/${CONFIG}.log" \
    --time="$TIME" \
    --wrap="$PYTHON chemogenetic/run/run_temporal_zplanes.py --config $CONFIG $EXTRA" || { echo "ERROR: sbatch failed (message above). Nothing was submitted."; exit 1; }

echo "Submitted."
echo "Monitor: squeue -u \$USER"
echo "Log:     ${LOG_DIR}/${CONFIG}.log"
