"""
EXPERIMENT_NAME.py
==================
Config for [DRUG] [DURATION]-minute chemogenetic experiment.

Copy this file, rename it to match the experiment
(e.g. ynt185_120min.py, cnqx_csn_90min.py), and fill in every
field marked TODO. Fields marked with a default value are safe to
leave unless this experiment differs.

Location:
    ~/Zebrafish-whole-brain-analysis/chemogenetic/config/EXPERIMENT_NAME.py
"""

import os

# ============================================================
# BASE PATHS  (shared across all chemogenetic experiments)
# ============================================================
dir_voluseg      = "/resnick/groups/Proberlab/yun/lightsheet/"
dir_registration = "/resnick/groups/Proberlab/yun/lightsheet/analysis_output/registration/"
dir_analysis     = "/resnick/groups/Proberlab/yun/lightsheet/analysis_output/chemogenetic/"

# Python interpreter used by every submit_*.sh (they read PYTHON_BIN from this config).
PYTHON_BIN = os.path.expanduser("~/miniconda3/envs/voluseg/bin/python")   # <- your conda env


# ============================================================
# PROJECT FOLDERS  — must match directory names on disk exactly
# ============================================================
CTRL_PROJ = "TODO: e.g. huc-h2b-g8m_csn_120min"
EXPT_PROJ = "TODO: e.g. huc-h2b-g8m_ynt185_120min"


# ============================================================
# FISH LISTS  — each fish is a (proj_ID, expt_ID) tuple
# ============================================================
ctrl_fish = [
    # (CTRL_PROJ, "YYMMDD_expt_ID_fishN"),
]

expt_fish = [
    # (EXPT_PROJ, "YYMMDD_expt_ID_fishN"),
]

all_fish = ctrl_fish + expt_fish


# ============================================================
# IMAGING SPECS  (defaults match light-sheet setup; change if needed)
# ============================================================
sec_per_volume = 1
volume_per_sec = 1
sampling_rate_hz = volume_per_sec

n_slices = 40
depth    = 250                         # microns
binning  = 1
res_x    = 1.52 * binning
res_y    = 1.52 * binning
res_z    = depth / n_slices
rotation_k = 2                         # 180° rotation for mapzebrain alignment


# ============================================================
# DRUG PERFUSION  (recorded in RUN_META only — the GLM regressor is regressor_files below)
# ============================================================
drug_uM  = None    # TODO: e.g. 10.0 (may stay None if you do not use a drug concentration)
V_ml     = 15.0
Q_ml_min = 4.5

# Phase boundaries in FRAMES (minutes * 60 * volume_per_sec); end of a phase is exclusive.
# TODO: adjust for your experiment. Phases must be contiguous and wash_end must equal the
# number of volumes in the recording (run_decompose checks this for every fish).
baseline_start = 0  * 60 * volume_per_sec
baseline_end   = 45 * 60 * volume_per_sec

drug_start = 45 * 60 * volume_per_sec
drug_end   = 90 * 60 * volume_per_sec

wash_start = 90  * 60 * volume_per_sec
wash_end   = 120 * 60 * volume_per_sec  # = total number of volumes at 1 vol/s


# ============================================================
# ANALYSIS EPOCHS  — one rule for every effect-size summary (ΔZ, d′, figures, responder sign)
# ============================================================
# Discard the first EPOCH_DISCARD_MIN of each phase (fish settling / habituation):
#   baseline summary window = [baseline_start + discard, baseline_end)
#   drug     summary window = [drug_start     + discard, drug_end)
# The seven 15-min temporal windows (epoch_dz / epoch_dprime maps) start at drug onset and
# are compared with the same baseline window. Defined in chemogenetic/epochs.py.
# The GLM fit window follows the same rule: INCLUDED_BASELINE (below) is the baseline minutes before
# drug onset that enter the fit, so the fit starts where the baseline summary window starts
# (45 - EPOCH_DISCARD_MIN = 25 for this design). The fit runs to the end of the recording.
EPOCH_DISCARD_MIN = 20


# ============================================================
# DECOMPOSITION  (F → F_tonic + F_phasic)
# ============================================================
f_tonic_window_size = 600              # seconds; default
f_tonic_percentile  = 20               # default


# ============================================================
# TONIC GLM
# ============================================================
input_tag    = "HCRT"                  # label of the user-defined regressor (names glm/inHCRT_K...). Keep "HCRT".
K_global     = 600                     # kernel length (frames); default
drift_global = 1                       # polynomial drift order; default
lam_global   = 0.5                     # ridge penalty; default
lag_global   = 0                       # causal lag (frames); default

# Minutes of baseline in the GLM fit window. Must reach back at least to the start of the
# baseline summary window (baseline_start + EPOCH_DISCARD_MIN), because the responder sign
# is read from that window. Rule: baseline length - EPOCH_DISCARD_MIN (45 - 20 = 25 here).
INCLUDED_BASELINE     = 25.0   # minutes of baseline before drug onset in the GLM fit (= baseline length - EPOCH_DISCARD_MIN)
CLIP_ABS_DZ           = 50.0           # |ΔZ| clipping threshold
NULL_TAG              = "iaaft"        # "iaaft" or "shift"
RESPONDER_NULL_THRESH = 99             # null percentile for responder cutoff (all consumers read this)
PHASIC_RESPONDER_NULL_THRESH = 99          # IAAFT null percentile for the phasic d' responder validation


# ============================================================
# CELL FILTER  (valid-cell mask; F_tonic / F_phasic are never modified)
# ============================================================
# run_decompose.py writes cell_valid.npy per fish; every later stage sets the excluded cells' ΔZ, d',
# and GLM ΔR² to NaN (so they are never responders and drop out of the maps) and samples the IAAFT
# nulls from valid cells only. Criteria use the BASELINE epoch only (independent of the drug effect):
CELL_MAX_ZERO_FRAC     = 0.10     # exclude: > this fraction of baseline samples with F_tonic <= 0 (clipped dim cells)
CELL_MIN_BASELINE_SD   = 1e-3     # exclude: baseline SD of F_tonic below this (ΔZ divides by it; counts)
CELL_MIN_BASELINE_MEAN = 1.0      # exclude: baseline mean of F_tonic below this (counts)
CELL_MAX_ABS_DZ        = 200.0    # backstop on cells that passed the above: max |tonic ΔZ| over the windows; None = off


# ============================================================
# TONIC RESPONDERS  (stage F1: ΔZ of F_tonic; replaces the GLM ΔR² as the definition. The GLM is optional stage F2)
# ============================================================
# Rank-based, within each fish: top p % of |max ΔZ| (max over the seven temporal windows, sign from that window),
# ranked inside baseline-σ strata. No null: the responder FRACTION is fixed by p and not comparable between groups.
TONIC_DZ_TOP_PERCENT      = 5.0          # main set
TONIC_DZ_TOP_PERCENT_SENS = [2.5, 10.0]  # sensitivity sets, saved next to the main one
TONIC_DZ_N_STRATA         = 5            # baseline-σ strata per fish (1 = no stratification)
TONIC_RESPONDER_SOURCE    = "dz"         # which tonic set presence maps read: "dz" (F1) or "glm" (F2)

# Derived: subfolder name for all GLM outputs (do not edit manually)
param_folder_name = (
    f"in{input_tag}_K{K_global}_drift{drift_global}_lam{lam_global}_lag{lag_global}"
)


# ============================================================
# TONIC GLM REGRESSOR  (optional stage F2 only) — defined by YOU; REQUIRED for run_glm.py (no fallback)
# ============================================================
# The GLM kernel is fit to this trace. run_glm.py stops with an error if it cannot be
# read. One CSV per reference recording; `regressor_column` is the raw trace (one value
# per volume, at least wash_end volumes). With "tonic_zscore" each trace becomes F_tonic
# (centered rolling percentile), is z-scored on [0, drug_start) and averaged across files
# (as in hcrt_regressor.ipynb); "none" uses the column as it is (e.g. a measured trace).
# input_tag (above) is only the label that names the GLM output folder (inHCRT_K600_...).
regressor_files = [
    # TODO: "/resnick/groups/Proberlab/<you>/lightsheet/<proj>/<expt>/hcrt_all.csv",
]
regressor_column         = "Mean"
regressor_preprocess     = "tonic_zscore"   # or "none"
regressor_tonic_window   = 600              # volumes
regressor_tonic_quantile = 0.20


# ============================================================
# GROUP TAGS + PLOT META
# ============================================================
EXPT_TAG = "TODO: e.g. YNT185"
CTRL_TAG = "CTRL"

PLOT_META = {
    EXPT_TAG: {"label": EXPT_TAG, "color": "TODO: e.g. steelblue", "alpha": 0.7},
    CTRL_TAG: {"label": CTRL_TAG, "color": "grey",                  "alpha": 0.7},
}

GROUP_OF_FISH = {fish: CTRL_TAG for fish in ctrl_fish}
GROUP_OF_FISH.update({fish: EXPT_TAG for fish in expt_fish})


def get_group_meta(fish):
    """Return (group_tag, plot_meta_dict) for a fish tuple."""
    tag = GROUP_OF_FISH[fish]
    return tag, PLOT_META[tag]


# ============================================================
# COMPARISON TAG + FIGURE PATH
# ============================================================
# TODO: update if EXPT_TAG / CTRL_TAG changed above
COMPARISON_TAG = f"{EXPT_TAG}_vs_{CTRL_TAG}"   # e.g. "YNT185_vs_CTRL"


def comparison_fig_dir(dir_analysis, comparison_tag=None):
    """
    Return the figure directory for a group-level comparison.

    Parameters
    ----------
    dir_analysis : str
    comparison_tag : str or None
        Defaults to COMPARISON_TAG from this config.

    Returns
    -------
    Path
        dir_analysis / comparisons / comparison_tag / figures
    """
    from pathlib import Path
    tag = comparison_tag if comparison_tag is not None else COMPARISON_TAG
    p = Path(dir_analysis) / "comparisons" / tag / "figures"
    p.mkdir(parents=True, exist_ok=True)
    return p
