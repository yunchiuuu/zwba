"""
config_ynt185_120min.py
=======================
Config for the YNT185 HCRTR-agonist 120-minute experiment.

Fish: pan-neuronal huc-h2b-g8m only (no HCRT:TRPV1).
Drug: YNT185 — direct pharmacological HCRT receptor agonist (perfused).
N = 6 fish.

Comparison note:
    No within-config ctrl group. For EXPT vs CTRL comparisons, import
    ctrl_fish_csn from config_hcrt_trpv1_csn_120min.py as the ctrl group.

Location:
    ~/zwba/chemogenetic/config/config_ynt185_120min.py
"""

import os

# ============================================================
# BASE PATHS
# ============================================================
dir_voluseg      = "/resnick/groups/Proberlab/yun/lightsheet/"
dir_registration = "/resnick/groups/Proberlab/yun/lightsheet/analysis_output/registration/"
dir_analysis     = "/resnick/groups/Proberlab/yun/lightsheet/analysis_output/chemogenetic/"

PYTHON_BIN = os.path.expanduser("~/miniconda3/envs/voluseg/bin/python")   # <- your conda env

# ============================================================
# PROJECT FOLDER
# ============================================================
YNT_PROJ = "huc-h2b-g8m_ynt185_120min"

# ============================================================
# FISH LISTS
# ============================================================
ynt_fish = [
    (YNT_PROJ, "260413_huc-h2b-g8m_ynt_10uM_fish1"),
    (YNT_PROJ, "260413_huc-h2b-g8m_ynt_10uM_fish2"),
    (YNT_PROJ, "260413_huc-h2b-g8m_ynt_10uM_fish3"),
    (YNT_PROJ, "260413_huc-h2b-g8m_ynt_10uM_fish4"),
    (YNT_PROJ, "260414_huc-h2b-g8m_ynt_10uM_fish1"),
    (YNT_PROJ, "260414_huc-h2b-g8m_ynt_10uM_fish2"),
]  # N = 6

# alias expected by run_temporal_intensity_map.py
expt_fish = ynt_fish

# ctrl group — same huc-h2b-g8m fish used across all CSN comparisons
CTRL_PROJ = "huc-h2b-g8m_csn_120min"
ctrl_fish = [
    (CTRL_PROJ, "251021_huc-h2b-g8m_csn_10uM_fish1"),
    (CTRL_PROJ, "251118_huc-h2b-g8m_csn_10uM_fish1"),
    (CTRL_PROJ, "251118_huc-h2b-g8m_csn_10uM_fish2"),
    (CTRL_PROJ, "251118_huc-h2b-g8m_csn_10uM_fish3"),
    (CTRL_PROJ, "251118_huc-h2b-g8m_csn_10uM_fish4"),
    (CTRL_PROJ, "251118_huc-h2b-g8m_csn_10uM_fish5"),
    (CTRL_PROJ, "251126_huc-h2b-g8m_csn_10uM_fish1"),
]  # N = 7

all_fish = ynt_fish + ctrl_fish


# ============================================================
# IMAGING SPECS
# ============================================================
sec_per_volume   = 1
volume_per_sec   = 1
sampling_rate_hz = volume_per_sec

n_slices  = 40
depth     = 250             # µm
binning   = 1
res_x     = 1.52 * binning
res_y     = 1.52 * binning
res_z     = depth / n_slices
rotation_k = 2              # 180° rotation for MapZebrain alignment


# ============================================================
# DRUG PERFUSION  (YNT185)
# ============================================================
drug_uM   = 10.0
V_ml      = 15.0
Q_ml_min  = 4

# Epoch windows in FRAMES (minutes * 60 * volume_per_sec)
baseline_start = 0  * 60 * volume_per_sec
baseline_end   = 45 * 60 * volume_per_sec

drug_start = 45 * 60 * volume_per_sec
drug_end   = 90 * 60 * volume_per_sec

wash_start = 90  * 60 * volume_per_sec
wash_end   = 120 * 60 * volume_per_sec


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
# DECOMPOSITION
# ============================================================
df_f_percentile     = 20
f_tonic_window_size = 600   # seconds
f_tonic_percentile  = 20


# ============================================================
# PERMUTATION TEST + BH-FDR
# ============================================================
p_thresh_permutation   = 0.005
n_resample_permutation = 500
BH_Q = 0.05


# ============================================================
# TONIC GLM
# ============================================================
input_tag    = "HCRT"                  # label of the user-defined regressor (names glm/inHCRT_K...). Keep "HCRT".
K_global     = 600                     # kernel length (frames)
drift_global = 1                       # polynomial drift order
lam_global   = 0.5                     # ridge penalty
lag_global   = 0                       # causal lag (frames)

INCLUDED_BASELINE = 25.0   # minutes of baseline before drug onset in the GLM fit (= baseline length - EPOCH_DISCARD_MIN)
CLIP_ABS_DZ       = 50.0               # |ΔZ| clipping threshold
NULL_TAG          = "iaaft"            # "iaaft" or "shift"
RESPONDER_NULL_THRESH = 99             # null percentile for responder cutoff
PHASIC_RESPONDER_NULL_THRESH = 99          # IAAFT null percentile for the phasic d' responder validation
L_MIN = 10.0                           # plateau duration (minutes) for plateau-ΔZ


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

# Derived: subfolder name for all GLM outputs (defined once, used everywhere)
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
_HCRT_REF_DIR = "/resnick/home/ychiu/yun/lightsheet/hcrt-trpv1_hcrt-h2b-g8m_120min"
regressor_files = [
    f"{_HCRT_REF_DIR}/260426_hcrt-trpv1_hcrt-h2b-g8m_csn_10uM_fish{i}/hcrt_all.csv"
    for i in (1, 2, 3)
]
regressor_column         = "Mean"
regressor_preprocess     = "tonic_zscore"   # or "none"
regressor_tonic_window   = 600              # volumes
regressor_tonic_quantile = 0.20


# ============================================================
# GROUP TAGS + PLOT META
# ============================================================
EXPT_TAG = "YNT185"
CTRL_TAG = "CTRL"     # refers to ctrl_fish_csn from config_hcrt_trpv1_csn_120min

PLOT_META = {
    EXPT_TAG: {"label": EXPT_TAG, "color": "mediumpurple", "alpha": 0.7},
    CTRL_TAG: {"label": CTRL_TAG, "color": "grey",         "alpha": 0.7},
}

GROUP_OF_FISH = {fish: EXPT_TAG for fish in ynt_fish}


def get_group_meta(fish):
    tag = GROUP_OF_FISH.get(fish, CTRL_TAG)
    return tag, PLOT_META[tag]


# ============================================================
# COMPARISON TAG + FIGURE PATH
# ============================================================
COMPARISON_TAG = f"{EXPT_TAG}_vs_{CTRL_TAG}"   # "YNT185_vs_CTRL"


def comparison_fig_dir(dir_analysis, comparison_tag=None):
    from pathlib import Path
    tag = comparison_tag if comparison_tag is not None else COMPARISON_TAG
    p = Path(dir_analysis) / "comparisons" / tag / "figures"
    p.mkdir(parents=True, exist_ok=True)
    return p
