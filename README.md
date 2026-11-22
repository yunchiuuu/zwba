# zwba — zebrafish whole-brain analysis

Pipeline for whole-brain light-sheet GCaMP imaging in larval zebrafish with a
chemogenetic / pharmacological perturbation (baseline → drug → washout). It takes
raw volumes to per-cell effect sizes, responder calls, and registered brain maps,
and compares an experimental group against controls.

```
raw TIFF ──► HDF5 volumes ──► voluseg (cells + traces) ──► ANTs SyN to shared template
                                                              │
   brain maps ◄── responders (GLM + IAAFT null, d′) ◄── F_tonic / F_phasic ◄─┘
```

Written for the Prober lab Resnick HPC (SLURM). Everything below runs **from the repo root**.

---

## 0. Before you start

**You need**
- a Resnick HPC account with access to `/resnick/groups/Proberlab/`
- the `voluseg` conda env (contains ANTsPy, voluseg, scipy, joblib, h5py, tifffile)
- raw volumes for your fish, acquired at the settings in section 3

**Paths and names — what to change for your own account**

Everything lab-specific lives in the config; the group data folder is
`/resnick/groups/Proberlab/<your_name>/lightsheet/` (Yun's is `yun`, her HPC login is `ychiu`,
and the repo lives at `~/zwba`). Outputs for each person go under their own group folder.

| What | Where | Set to |
|---|---|---|
| `dir_voluseg`, `dir_registration`, `dir_analysis` | every `chemogenetic/config/config_*.py` | your `/resnick/groups/Proberlab/<name>/lightsheet/...` |
| `PYTHON_BIN` | every config (default `~/miniconda3/envs/voluseg/bin/python`, expands to your home) | your conda env, if it lives elsewhere |
| `regressor_files` | every config | the regressor CSVs you define (section 6) |
| `dir_voluseg`, `dir_registration`, `PYTHON_BIN` | `registration/config/config_registration.py` | same as above |
| `BASE` | `run_voluseg.py` | your `lightsheet/` folder (also holds the shared `voluseg_setup/install/bin` ANTs install) |
| ANTs scratch dir | the three registration scripts | automatic: `/resnick/scratch/<your login>/tmp_ants` (`getpass.getuser()`) |

The submit scripts read `PYTHON_BIN` from the config you pass with `--config`; they stop with a
clear message if the config cannot be imported or the interpreter does not exist.

---

## 1. Concepts you need

- **A fish is a `(proj_ID, expt_ID)` tuple.** All paths resolve through
  `utils.data_io.fish_dir(base_dir, fish)` = `base_dir / proj_ID / expt_ID`.
- **Three path trees** (set in the config):
  - `dir_voluseg` — raw HDF5 in, voluseg out (`input/volume1000.h5…`, `output/cells0_clean.hdf5`, `output/volume0.hdf5`)
  - `dir_registration` — ANTs warps/affines and the shared template
  - `dir_analysis` — everything the analysis writes (`f_tonic.npy`, GLM folders, figures)
- **Controls are processed once and reused.** Control fish live under their own
  `proj_ID` (e.g. `huc-h2b-g8m_csn_120min`); any config that lists them reads the same
  outputs. If your drug protocol matches an existing control cohort, list those fish as
  your `ctrl_fish` instead of re-imaging.
- **Module vs run script.** `chemogenetic/*.py` are modules (operate on one fish).
  `chemogenetic/run/run_*.py` loop over the fish list in the config and are the sbatch
  targets. `chemogenetic/submit/submit_*.sh` wrap them.
- **`--config` is always a module name** (no folder, no `.py`):
  `--config config_hcrt_trpv1_csn_120min` (chemogenetic, under `chemogenetic/config/`) or
  `--config config_registration` (registration, under `registration/config/`).

---

## 2. Run order for a new dataset

### Step 1 — register your fish in two configs
1. `cp chemogenetic/config/config_template.py chemogenetic/config/config_<exp>.py`, then fill every `TODO`
   (project folders, fish tuples, `drug_uM`, epoch windows, group tags, `PYTHON_BIN`).
2. Add the same fish to `registration/config/config_registration.py` (`all_fish`).
   `submit_registration_syn.sh` submits **only uncommented** fish.

### Step 2 — raw TIFF → HDF5 (analysis PC, not HPC)
Not in the repo yet (`run_preprocess_batch.py`, `tif_to_h5.py`). Output per fish:
`{dir_voluseg}/{proj}/{expt}/input/volume1000.h5, volume1001.h5, …` with dataset key
`'default'`. Run conversions **serially** (parallel conversion produced truncated files), then sync to the HPC.

### Step 3 — voluseg (HPC, one job per fish)
```bash
sbatch --job-name=voluseg_<fish> --cpus-per-task=24 --mem=200G --time=8:00:00 \
  --output=logs/voluseg/<expt_ID>.log \
  --wrap="<PYTHON_BIN> run_voluseg.py --proj_ID <proj_ID> --expt_ID <expt_ID>"
```
Writes `output/cells0_clean.hdf5` and `output/volume0.hdf5`. Never run more than one fish per process (Spark).
Imaging constants (`FRAMES_PER_VOLUME=40`, 250 µm, 1.52 µm/px) are set at the top of `run_voluseg.py`.

### Step 4 — registration to the shared template (HPC)
```bash
bash registration/submit_registration_syn.sh --config config_registration
```
One job per fish in `all_fish` (32 CPUs; the script requests `MEM="24G"` — raise it in the
script if jobs are OOM-killed). **Fish that already have a complete registration made with the current
template are skipped**: the three output files exist and `registration_meta.json` records the same template
sha256. Preview first with `--dry-run`, which lists every fish as SUBMIT or SKIP with the reason:
```bash
bash registration/submit_registration_syn.sh --config config_registration --dry-run
bash registration/submit_registration_syn.sh --config config_registration --overwrite        # redo every fish
bash registration/submit_registration_syn.sh --config config_registration --trust-existing   # also skip old fish with no meta file
```
A fish is re-registered when files are missing, when the template changed, or with `--overwrite`.
Per fish it writes `expt_to_mean_warp.nii.gz`,
`expt_to_mean_affine.mat`, a QC overlay PDF, and `registration_meta.json` (template sha256,
SyN parameters, git commit).

- **Do not regenerate the template.** `template_mean_brain.nii.gz` (version `zstack_v1`) is
  read-only and shared by all fish. `generate_template.py` and
  `register_mapzebrain_to_template.py` were one-time steps. If the template ever changes,
  re-register **all** fish and re-run medoids with overwrite.
- Open `qc_registration_overlay.pdf` for every fish before moving on.

### Step 5 — chemogenetic preprocessing (HPC)
```bash
bash chemogenetic/submit/submit_decompose.sh --config config_<exp>
#   -> f_dark_scalar.npy, f_tonic.npy, f_phasic.npy, QC_figures/decompose_QC.png
bash chemogenetic/submit/submit_medoids.sh --config config_<exp>      # needs step 4
#   -> medoids_all.npy, medoids_rot_all.npy, medoids_template_vox.npy
```
Both skip fish whose outputs exist. To recompute, flip `OVERWRITE_*` at the top of the run script.
`run_decompose.py` checks the recording length against the epochs, and warns when many samples are
clipped at `F_dark` (>5 %) or `F_phasic` is non-finite (>2 %). Look at `decompose_QC.png` for every fish.

### Step 6 — responder statistics (HPC)
```bash
bash chemogenetic/submit/submit_glm.sh    --config config_<exp>   # GLM + IAAFT null + tonic responders
bash chemogenetic/submit/submit_dprime.sh --config config_<exp>   # phasic d′ + IAAFT + group figures
bash chemogenetic/submit/submit_figures_glm.sh --config config_<exp>   # after the GLM finishes
```
`run_glm.py` and `run_dprime.py` take `--stages` and `--overwrite`; the submit scripts pass them through:

```bash
bash chemogenetic/submit/submit_glm.sh --config config_<exp> --stages responders             # redo only the responder call
bash chemogenetic/submit/submit_glm.sh --config config_<exp> --stages refit,ablation,null,responders --overwrite refit,ablation,null,responders
bash chemogenetic/submit/submit_dprime.sh --config config_<exp> --stages responders,figures
```
| Script | Valid stages | Default `--stages` | Default `--overwrite` |
|---|---|---|---|
| `run_glm.py` | `cv, refit, ablation, null, responders` | `refit,ablation,null,responders` | `responders` |
| `run_dprime.py` | `dprime, null, responders, figures` | all four | `responders` |

Existing outputs are skipped unless listed in `--overwrite`. The responder call uses
`RESPONDER_NULL_THRESH` (GLM) and `PHASIC_RESPONDER_NULL_THRESH` (d′) from the config — both are IAAFT null
percentiles you can change (default 99); if only another percentile's threshold file exists, it is
derived from the saved IAAFT null distribution (no surrogates are re-run).
GLM outputs land in `{dir_analysis}/{proj}/{expt}/glm/in{input_tag}_K{K}_drift{d}_lam{λ}_lag{l}/`, plus
`tonic_pos/neg_glm_iaaft_nullp{N}_idxs.npy` and `phasic_pos/neg_dprime_iaaft_nullp{N}_idxs.npy`.

**After changing `INCLUDED_BASELINE`** (or the regressor), the fits must be rebuilt. The GLM folder name does not
encode the fit window, so `run_glm.py` reads `RUN_META.json` and stops if existing fits used a different window.
Rebuild everything in one go:
```bash
bash chemogenetic/submit/submit_glm.sh --config config_<exp> \
  --stages cv,refit,ablation,null,responders --overwrite refit,ablation,null,responders
```
(`cv` re-checks K/λ for the new window; the globals in the config are still what the refit uses.)
Then rerun `submit_figures_glm.sh` and the maps. Use a different `param_folder_name` if you want to keep the old results.

### Step 7 — maps (all need `medoids_template_vox.npy`)
```bash
bash chemogenetic/submit/submit_temporal_intensity_map.sh --config config_<exp>   # epoch_dz.npy, epoch_dprime.npy, per-window maps
bash chemogenetic/submit/submit_glm_dprime_map.sh        --config config_<exp>   # N-fish presence maps
bash chemogenetic/submit/submit_kernel_map.sh            --config config_<exp>   # kernel maps
bash chemogenetic/submit/submit_combined_fish_map.sh     --config config_hcrt_trpv1_pooled_csn_120min
bash chemogenetic/submit/submit_survey_by_regions.sh [--config config_<exp>] [--overwrite] [--top_n N]
```
`submit_survey_by_regions.sh` without `--config` submits the three configs listed at the top of the script;
any other arguments are passed to `run_survey_by_regions.py`. It needs `epoch_*.npy` and the MapZebrain masks
in `{dir_registration}/300um_template/mapzebrain_in_template/`.
`run_glm_brainmap.py` has no submit script. `submit_individual_fish_map.sh` points to a run script that is not in
the repo (the script now stops with a message).

Figures go to `{dir_analysis}/{proj}/{expt}/figures/` (per fish) and
`{dir_analysis}/comparisons/{COMPARISON_TAG}/figures/` (group level).

### When something fails
Every run script uses the shared guard in `utils/run_guard.py`:
- **Config / input problems** are checked before heavy work and reported all at once (every missing
  config key, every missing file). They print `❌ <message>` and exit with code **2**.
- **Per-fish failures** are caught, printed with a traceback, and the batch continues. At the end of
  the log a summary lists every failed item, and the job exits with code **1** so `sacct` shows it as `FAILED`.
- A clean run ends with `✅ complete: N ok, 0 failed`. `grep -n "❌" logs/<stage>/<config>.log` shows all failures.

---

## 3. What your data must look like

The timing is read from the config (`volume_per_sec`, `baseline_start/end`, `drug_start/end`,
`wash_end`, `EPOCH_DISCARD_MIN`); nothing in the analysis assumes 1 volume/s or 120 min any more.
The 120-min CSN design is 1 vol/s, 45 min baseline / 45 min drug / 30 min washout (drug onset at vol 2700).
The scripts refuse to run when the recording is shorter than the windows need (`wash_end`), or when the
discard leaves an empty window.

Imaging constants in the config (`n_slices=40`, `depth=250`, `res_x/y=1.52`, `rotation_k=2`) are
still per-setup values; `run_voluseg.py` has its own copy at the top.

---

## 4. How the quantities are defined (from the code)

1. **`F_dark`** — one scalar per fish: pooled median of 10×10 corner patches (top-left and
   bottom-left, all z-planes) over 100 evenly spaced volumes.
2. **`F_corrected`** = `clip(F_raw − F_dark, min=0)`.
3. **`F_tonic`** — per cell, `scipy.ndimage.percentile_filter` of `F_corrected`: centered,
   `mode="reflect"`, window `f_tonic_window_size` (config; 600 s) and percentile
   `f_tonic_percentile` (config; 20), in counts.
4. **`F_phasic`** = `(F_corrected − F_tonic) / F_tonic`. Default `DENOM_MODE="legacy"`
   (NaN where `F_tonic == 0`, no floor); `"fixed_floor"` is available for new experiments.
5. **Nothing z-scores raw traces.** Normalisation is per cell, per fish, inside the effect-size step.

**One epoch rule everywhere** (defined once in `chemogenetic/epochs.py` from the config):
the first `EPOCH_DISCARD_MIN` (default 20) minutes of each phase are discarded.

| Epoch | Frames (120-min design) | Minutes |
|---|---|---|
| Baseline summary window | `[baseline_start + 20·60·volume_per_sec, baseline_end)` = 1200–2700 | 20–45 |
| Drug summary window | `[drug_start + 20·60·volume_per_sec, drug_end)` = 3900–5400 | 65–90 |
| Temporal windows (`epoch_dz.npy`, `epoch_dprime.npy`) | seven 15-min windows, 10-min step, **starting at drug onset (not discarded)**, each compared with the baseline window above | 0–15, 10–25, … 60–75 after onset |

Used by: `run_decompose.py` QC, `temporal_windows.py`, the temporal / survey maps, `run_dprime.py`,
`run_figures_glm.py` and the GLM responder sign. All effect sizes are
`(μ_window − μ_baseline) / σ_baseline` per cell. For the epoch path, **"d′" is not a true d′**: it
uses the baseline SD only (no pooled SD), i.e. it is ΔZ computed on `F_phasic`.

The **GLM fit window** follows the same rule. `INCLUDED_BASELINE` is the number of baseline minutes
*before drug onset* that enter the fit, so the fit window is

```
[drug_start − INCLUDED_BASELINE·60·volume_per_sec,  end of recording)      # includes drug AND washout
```
Set `INCLUDED_BASELINE = baseline length − EPOCH_DISCARD_MIN` (45 − 20 = **25**), so the fit starts at min 20,
exactly where the baseline summary window starts (min 20–120 for the 120-min design). The responder sign
(pos/neg) is `mean(ŷ over the drug window) − mean(ŷ over the baseline window)` using the epoch windows
above, expressed relative to the fit start (baseline min 0–25, drug min 45–70 of the fit). The run stops
with the minimum value needed if `INCLUDED_BASELINE` is too short, and warns if it reaches back into the
discarded minutes.

Epoch arrays are cached per fish with an `epoch_meta.json` sidecar recording the windows they were
computed with. If you change `EPOCH_DISCARD_MIN` (or the timing), the maps recompute automatically;
caches without a sidecar (made before this change) are recomputed once.

**Responder calls**
- Tonic: GLM with IAAFT null at the 99th percentile (`RESPONDER_NULL_THRESH = 99`).
- Phasic: fixed `|d′| ≥ 0.5`, with the ctrl-pooled IAAFT p99 as validation.
- Statistics must use **fish as the unit**, never cells.

---

## 5. Interpreting results — current caveats

- **Baseline contamination.** The centered `F_tonic` filter can pull drug-period
  fluorescence into the last part of the baseline (within half a filter window of onset).
  This inflates `μ_baseline` / `σ_baseline`.
- **Edge effects.** Reflect padding makes `F_tonic` near the end of the recording
  non-independent; the last windows (W5, W6) are most affected.
- **Tonic σ_baseline is mostly slow drift**, so tonic ΔZ is not a pure noise-normalised
  effect size and can differ between fish for non-biological reasons.
- **Signed averaging cancels** at the voxel mean, the cross-fish mean and the Z/X projection.
  Strong responders of opposite sign can look flat in a mean map.
- **No N-support masking in difference maps.** A voxel backed by 2 fish looks like one
  backed by 9, and the EXPT − CTRL map is NaN wherever either group is empty.
- **`clip_min=0` in background subtraction** truncates the lower tail of dim cells near `F_dark`.

---

## 6. The GLM regressor (you define it in the config)

There is no fallback regressor. Set in your config:

```python
regressor_files = ["/resnick/.../<fish>/hcrt_all.csv", ...]  # one CSV per reference recording
regressor_column         = "Mean"           # raw trace column, one value per volume
regressor_preprocess     = "tonic_zscore"   # or "none" to use the column as is
regressor_tonic_window   = 600              # volumes
regressor_tonic_quantile = 0.20
```
With `tonic_zscore` each trace becomes a centered rolling 20th-percentile (F_tonic), is z-scored on
`[0, drug_start)` (ddof=1) and averaged across files; this reproduces the original `run_glm.py`
regressor exactly. Each file must exist, have the column, contain no NaN and cover the whole recording
(`wash_end` volumes); all problems are listed together and the run stops before any fitting.
`input_tag` must stay `"HCRT"`: it is only the label that names the `glm/inHCRT_K…` folder.
The same files draw the reference row of the temporal maps (z-scored on the baseline window, ddof=0).

---

## 7. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `could not read PYTHON_BIN from chemogenetic/config/...` | run from the repo root and pass the config module name; your config needs `PYTHON_BIN` |
| `PYTHON_BIN does not exist` | point `PYTHON_BIN` in the config at your conda env's `bin/python` |
| `❌ <config> is missing required setting(s)` | copy the listed names from `config_template.py` (e.g. `EPOCH_DISCARD_MIN`, `regressor_files`) |
| `Cannot load the GLM regressor` | the message lists each bad file (missing / no column / NaN / too short) |
| `recording has N volumes but epochs need M` | wrong fish in the list or wrong `wash_end`/timing in the config |
| Job `FAILED` in `sacct` | exit code 1 = some fish failed (see the summary at the end of the log); exit code 2 = bad config/input |
| Registration submit says `--config is required` | pass a module name: `--config config_registration` |
| Registration submits 0 jobs | all fish are already registered with the current template (see `--dry-run`), or commented out in `config_registration.py` |
| Registered fish look misaligned | check `registration_meta.json` template sha256 vs. the current template; check the affine z-scale in the QC overlay (values far from ~1.0 are suspect) |
| `Existing GLM outputs were fit with a different window` | `INCLUDED_BASELINE` changed since those fits; rebuild as in Step 6 |
| GLM results differ from an older run | the fit window, regressor or threshold changed — rerun with `--stages refit,ablation,null,responders --overwrite refit,ablation,null,responders` |

---

## 8. Repo layout

```
run_voluseg.py                   voluseg steps, one fish per job
registration/
  config/config_registration.py  fish lists, template path/version, imaging specs
  run_registration_syn.py        SyN per fish (+ QC, meta sidecar)
  submit_registration_syn.sh     one sbatch per uncommented fish
  generate_template.py           one-time template build (do not re-run)
  register_mapzebrain_to_template.py, submit_mapzebrain_registration.sh   one-time atlas → template
  registration.py, coordinates.py  library code
chemogenetic/
  config/config_<exp>.py         per-experiment config (copy config_template.py)
  run/  submit/                  run scripts and their sbatch wrappers
  epochs.py regressor.py        the single epoch definition; the config-defined GLM regressor
  glm.py dprime.py temporal_windows.py drug_regressor.py brainmap.py kernel_map.py population.py scatter.py traces.py
utils/  data_io.py  preprocess.py  stats.py  run_guard.py   (run_guard: failure railings used by every run script)
```

Legacy Spearman / permutation / effect-size pipeline was removed (superseded by GLM + IAAFT);
recover from git tag `pre-spearman-removal`.

## 9. Maintainer TODO

- Commit the missing preprocessing scripts (`run_preprocess_batch.py`, `tif_to_h5.py`, `apply_voluseg_patch.py`) and `run_individual_fish_map.py`.
- Read the remaining hardcoded timing in `run_combined_fish_map.py` / `run_kernel_map.py` / `run_glm_brainmap.py` from the config (they only load existing outputs, but check before using a different design).
- Add `--overwrite` CLI flags to `run_decompose.py` and `run_medoids.py` (still `OVERWRITE_*` constants).
- Add N-support masks to the difference maps.
- Stop tracking `__pycache__/`, `.ipynb_checkpoints/` and `logs/` (`git rm -r --cached`; `.gitignore` is updated).
