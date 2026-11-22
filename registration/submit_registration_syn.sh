#!/bin/bash
TIME="${ZWBA_TIME:-6:00:00}"   # walltime; override per call: ZWBA_TIME=48:00:00 bash ...
# submit_registration_syn.sh
# ==========================
# Submit one sbatch job per fish for SyN registration to the shared template.
# Fish list and Python path are read from the config passed with --config.
# The template (generate_template.py) must already exist; its versioned
# filename comes from the same config (MEAN_BRAIN_PATH).
#
# A fish is SKIPPED when it already has a complete registration made with the
# CURRENT template:
#     expt_to_mean_warp.nii.gz + expt_to_mean_affine.mat + expt_to_mean_registered_final.nii.gz
#     + registration_meta.json whose template sha256 equals the sha256 of MEAN_BRAIN_PATH.
# (registration_meta.json is the last file the job writes, so it marks a finished run.)
# A fish is (re)submitted when files are missing, the template changed, or --overwrite is given.
#
# --config is a MODULE NAME under registration/config/ (no folder, no .py),
# the same convention as the chemogenetic submit scripts.
#
# Usage (from repo root):
#     bash registration/submit_registration_syn.sh --config config_registration
# Options:
#     --dry-run          show what would be submitted / skipped, submit nothing
#     --overwrite        re-register every fish in all_fish, even if complete
#     --trust-existing   also skip fish that have the three output files but no
#                        registration_meta.json (made before provenance was recorded;
#                        their template cannot be verified)

SCRIPT="registration/run_registration_syn.py"
CPUS=32
MEM="24G"
LOG_DIR="logs/registration"

# ── Parse arguments ───────────────────────────────────────────────────────────
CONFIG=""
DRY_RUN=0
OVERWRITE=0
TRUST_EXISTING=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --overwrite) OVERWRITE=1; shift ;;
        --trust-existing) TRUST_EXISTING=1; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "ERROR: --config is required."
    echo "Usage: bash registration/submit_registration_syn.sh --config config_registration [--dry-run] [--overwrite] [--trust-existing]"
    exit 1
fi

if [[ "$CONFIG" == *.py || "$CONFIG" == */* ]]; then
    echo "ERROR: --config must be a module name (e.g. config_registration), not a path: $CONFIG"
    exit 1
fi

mkdir -p "$LOG_DIR"

# Bootstrap: use system python3 to read PYTHON_BIN from config
PYTHON=$(python3 -c "
import sys
sys.path.insert(0, 'registration/config')
from ${CONFIG} import PYTHON_BIN
print(PYTHON_BIN)
") || { echo "ERROR: could not import ${CONFIG} from registration/config/"; exit 1; }
[ -x "$PYTHON" ] || { echo "ERROR: PYTHON_BIN does not exist or is not executable: $PYTHON  (set PYTHON_BIN in your config)"; exit 1; }

# Plan: one line per fish, "SUBMIT<TAB>expt_ID<TAB>reason" or "SKIP<TAB>expt_ID<TAB>reason"
PLAN=$(REG_OVERWRITE=$OVERWRITE REG_TRUST=$TRUST_EXISTING $PYTHON - "$CONFIG" <<'PYEOF'
import hashlib, json, os, sys
sys.path.insert(0, "registration/config")
cfg = __import__(sys.argv[1])
overwrite = os.environ.get("REG_OVERWRITE") == "1"
trust = os.environ.get("REG_TRUST") == "1"
template = str(cfg.MEAN_BRAIN_PATH)
if not os.path.exists(template):
    print(f"ERROR\t-\ttemplate not found: {template}")
    sys.exit(0)
h = hashlib.sha256()
with open(template, "rb") as f:
    for block in iter(lambda: f.read(1 << 20), b""):
        h.update(block)
sha = h.hexdigest()
need = ["expt_to_mean_warp.nii.gz", "expt_to_mean_affine.mat", "expt_to_mean_registered_final.nii.gz"]
for proj_id, expt_id in cfg.all_fish:
    d = os.path.join(str(cfg.dir_registration), proj_id, expt_id)
    if overwrite:
        print(f"SUBMIT\t{expt_id}\t--overwrite"); continue
    missing = [n for n in need if not os.path.exists(os.path.join(d, n))]
    if missing:
        print(f"SUBMIT\t{expt_id}\tnot registered yet (missing {', '.join(missing)})"); continue
    mp = os.path.join(d, "registration_meta.json")
    if not os.path.exists(mp):
        if trust:
            print(f"SKIP\t{expt_id}\tfiles present, no registration_meta.json (--trust-existing)")
        else:
            print(f"SUBMIT\t{expt_id}\tfiles present but no registration_meta.json (template unverified; use --trust-existing to skip)")
        continue
    try:
        old = json.load(open(mp))["template"]["sha256"]
    except Exception as e:
        print(f"SUBMIT\t{expt_id}\tregistration_meta.json unreadable ({type(e).__name__})"); continue
    if old != sha:
        print(f"SUBMIT\t{expt_id}\ttemplate changed since this registration")
    else:
        print(f"SKIP\t{expt_id}\tcomplete, same template")
PYEOF
) || { echo "ERROR: could not read the fish list / template from ${CONFIG}"; exit 1; }

if [ -z "$PLAN" ]; then
    echo "ERROR: No fish found in ${CONFIG}. Check that fish are uncommented in all_fish."
    exit 1
fi
if echo "$PLAN" | grep -q "^ERROR"; then
    echo "$PLAN" | grep "^ERROR" | cut -f3-
    exit 1
fi

N_TOTAL=$(echo "$PLAN" | wc -l)
N_SKIP=$(echo "$PLAN" | grep -c "^SKIP")
N=$((N_TOTAL - N_SKIP))
echo "$N_TOTAL fish in ${CONFIG}: $N to submit, $N_SKIP skipped."
echo "$PLAN" | while IFS=$'\t' read -r ACTION EXPT_ID REASON; do
    printf "  %-6s %s   (%s)\n" "$ACTION" "$EXPT_ID" "$REASON"
done

if [ "$DRY_RUN" -eq 1 ]; then
    echo ""
    echo "Dry run: nothing submitted."
    exit 0
fi
if [ "$N" -eq 0 ]; then
    echo "Nothing to do."
    exit 0
fi

echo ""
echo "Submitting $N registration jobs..."
while IFS=$'\t' read -r ACTION EXPT_ID REASON; do
    [ "$ACTION" = "SUBMIT" ] || continue
    JOB_NAME="reg_${EXPT_ID}"
    LOG="${LOG_DIR}/${EXPT_ID}.log"

    sbatch \
        --job-name="$JOB_NAME" \
        --cpus-per-task=$CPUS \
        --mem=$MEM \
        --output="$LOG" \
        --time="$TIME" \
        --wrap="$PYTHON $SCRIPT --config $CONFIG --expt_ID $EXPT_ID" || { echo "ERROR: sbatch failed (message above). Stopping."; exit 1; }

    echo "  Submitted: $EXPT_ID"
done <<< "$PLAN"

echo ""
echo "$N jobs submitted ($N_SKIP skipped)."
echo "Monitor with: squeue -u \$USER"
echo "Logs in:      $LOG_DIR/"
