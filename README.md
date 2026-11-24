# zwba — zebrafish whole-brain analysis

Pipeline for whole-brain light-sheet GCaMP imaging in larval zebrafish with a
chemogenetic / pharmacological perturbation (baseline → drug → washout). It takes
raw volumes to per-cell effect sizes, responder calls, and registered brain maps,
and compares an experimental group against controls.

```
raw TIFF ──► HDF5 volumes ──► voluseg (cells + traces) ──► ANTs SyN to shared template
                                                              │
   brain maps ◄── responders (tonic: top-p% |ΔZ|; phasic: d′ vs per-fish IAAFT null) ◄── F_tonic / F_phasic ◄─┘
```

Written for the Prober lab Resnick HPC (SLURM). Everything below runs **from the repo root**.

---

## Quick start

**Contents:** Quick start · Part 1 Before you start · Part 2 Concepts · Part 3 Data and quantities · **Part 4 Responder methods (F1, F2, G: nature, null, pros/cons)** ·
Part 5 Stages A–I (commands) · Part 6 Interpreting results · Part 7 GLM regressor · Part 8 Troubleshooting · Part 9 Repo layout · Part 10 Maintainer TODO


**Where am I?** Run this any time. It is read-only and takes seconds, and it prints a table of every fish × stage
(✔ done, ✗ missing, · blocked by an earlier stage, ⚠ stale) plus the next command to run:

```bash
cd ~/zwba
bash chemogenetic/submit/status.sh --config config_<exp>
```

**Run order** (each stage is explained in Part 5; the F1 / F2 / G methods in Part 4; `<exp>` = your config module name):

|  | Stage | Command | Needs | Check before moving on |
|---|---|---|---|---|
| A | TIFF → HDF5 (analysis PC) | `tif_to_h5.py` (not in repo yet) | raw volumes | files `input/volume1000.h5…` complete |
| B | voluseg (one fish per job) | `sbatch … run_voluseg.py --proj_ID P --expt_ID E` | A | `output/cells0_clean.hdf5`, `volume0.hdf5` exist |
| C | registration to template | `bash registration/submit_registration_syn.sh --config config_registration` | B | **open every `qc_registration_overlay.pdf`** |
| D | decompose (F_tonic / F_phasic) | `bash chemogenetic/submit/submit_decompose.sh --config config_<exp>` | B | **open every `decompose_QC.png`** |
| E | medoids → template coordinates | `bash chemogenetic/submit/submit_medoids.sh --config config_<exp>` | C | status shows ✔ in `med` |
| F1 | tonic responders (top-p% \|max ΔZ\| of `F_tonic`) | `bash chemogenetic/submit/submit_tonic_dz.sh --config config_<exp>` | D | log ends `✅ complete … 0 failed`; **argmax-window shares** printed per fish look sensible (see Stage F1) |
| F2 (optional) | GLM + IAAFT null + GLM tonic responders | `bash chemogenetic/submit/submit_glm.sh --config config_<exp>` | D | not in the default run; needs the regressor in the config (Part 7); log ends `✅ complete … 0 failed` |
| G | phasic d′ + per-fish IAAFT null + responders + group figures | `bash chemogenetic/submit/submit_dprime.sh --config config_<exp>` | D | log ends `✅ complete … 0 failed`; the per-fish thresholds T are listed (all within a similar range) |
| G2 | d′ null diagnostics | `bash chemogenetic/submit/submit_figures_dprime_null.sh --config config_<exp>` | G | real \|d′\| tail rises above the null in the expt fish and not in ctrl; **null SD spread** (p95/p5) is below ~3 |
| H1 | tonic ΔZ figures | `bash chemogenetic/submit/submit_figures_tonic_dz.sh --config config_<exp>` | F1 | **argmax-window histogram** is not dominated by washout windows; σ-vs-\|ΔZ\| shows no flat-baseline cloud |
| H2 (optional) | GLM figures | `bash chemogenetic/submit/submit_figures_glm.sh --config config_<exp>` | F2 | responder fraction in controls is near the false-positive rate (≈ 1 − p/100) |
| I | maps (the kernel map needs F2) | `submit_temporal_intensity_map.sh`, `submit_glm_dprime_map.sh` (+ optional `submit_tonic_overlap.sh`), `submit_kernel_map.sh`, `submit_combined_fish_map.sh`, `submit_survey_by_regions.sh` | E, F1, G | per-fish N shown on the maps |

**The three responder methods in one line each** (full notes, null definitions and pros/cons: Part 4):
- **F1** tonic ΔZ: top 5 % of |max ΔZ| per baseline-σ stratum, **per fish, rank-based, no null**; the responder fraction is fixed by p and carries no information.
- **F2** tonic GLM (optional): supervised by the HCRT regressor, **per-fish IAAFT null on the regressor**, one scalar threshold per fish.
- **G** phasic d′: **per-fish IAAFT null on the baseline/drug label**, one scalar threshold per fish (not control-pooled); in our 16-fish data the control fish still sit above the nominal 1 % (Part 4.4).

C and D are independent of each other and can run at the same time; F1 and G can too. The GLM (optional stage F2, `submit_glm.sh`, with its figures H2 `submit_figures_glm.sh` and the kernel map) is **not** part of the default run; see Stage F2.
Registration (C), decompose (D), medoids (E), tonic ΔZ (F1) and d′ (G) skip fish whose outputs already exist; use `--overwrite`, `--stages` or the `OVERWRITE_*` constants to redo. Maps cache their epoch arrays and recompute them when the epochs change.
A failed fish never stops the batch, but the job exits with code 1 and the end of the log lists every failure.

### Reading the status table

```bash
bash chemogenetic/submit/status.sh --config config_<exp>               # all fish in the config
bash chemogenetic/submit/status.sh --config config_<exp> --fish fish3   # only expt_IDs containing "fish3"
bash chemogenetic/submit/status.sh --config config_<exp> --strict       # exit code 1 if anything is missing or stale
bash chemogenetic/submit/status.sh --config config_<exp> --glm          # also show the optional GLM columns (fit/null/resp)
```
It only reads files (nothing is submitted or created) and takes a few seconds. One row per fish, one column per stage:

| Column | Stage | Done when |
|---|---|---|
| `vol` | voluseg | `cells0_clean.hdf5` and `volume0.hdf5` exist |
| `reg` | registration | warp, affine and registered image exist, and `registration_meta.json` records the current template |
| `dec` | decompose | `f_dark_scalar.npy`, `f_tonic.npy`, `f_phasic.npy` |
| `cell` | cell filter | `cell_valid.npy` made from the current `f_tonic.npy`, epochs and `CELL_*` thresholds (⚠ = any of them changed; ✗ = run decompose, it only builds the filter). Does not block other columns |
| `med` | medoids | `medoids_all.npy`, `medoids_rot_all.npy`, `medoids_template_vox.npy` |
| `tdz` | tonic ΔZ responders (F1) | `tonic_dz_*.npy` and `tonic_{pos,neg}_dz_top<p>_idxs.npy` made with the current epochs, cell filter and `TONIC_DZ_*` settings (⚠ = any of them changed; · = blocked until `dec` is done) |
| `fit`, `null`, `resp` | optional GLM (F2); shown only with `--glm` | files in `glm/<param_folder_name>/` and the `tonic_pos/neg_glm_…_nullp<N>_idxs.npy` files, with the p<N> from the config |
| `dpr`, `dnul`, `dres` | d′, its per-fish null, responder calls | `phasic_dprime_cells_*.npy`; `phasic_dprime_null__iaaft.npy` + sidecar made with the current cell filter (**every fish, ctrl and expt**); `phasic_pos/neg_dprime_…_idxs.npy` + `phasic_dprime_responders_meta.json` recording per-fish thresholds (⚠ = called with the old ctrl-pooled threshold or an older filter) |
| `map` | epoch maps | `epoch_dz.npy` and `epoch_dprime.npy`, made with the current epochs |

| Symbol | Meaning |
|---|---|
| ✔ | done and current |
| ✗ | missing, and everything before it is done: this is what to run next |
| · | blocked: an earlier stage in that row is not done yet |
| ⚠ | exists but stale or unverifiable; the "Details" list below the table says why |

Below the table, **Next** lists the command for each stage that has ✗ or ⚠ and how many fish it affects, and **Details**
explains each ⚠ (different template, GLM fit window changed, epoch cache made with other epochs). A stage with a ⚠ in
`map` just needs the map script re-run; a ⚠ in `tdz` means rerun `submit_tonic_dz.sh` (it recomputes by itself); a ⚠ in `fit` (`--glm` only) means the GLM must be rebuilt (see Stage F1).
A ⚠ on `tdz`, `fit`, `null`, `dpr`, `dnul` or `map` can also mean that output was made without (or with a different) cell filter; `tdz`, `dpr`, `dnul` and `map` are recomputed automatically by their scripts, the GLM needs the command shown under Stage F2.
`Next` does not know about fish you are deliberately leaving out; use `--exclude` (Stage F2). The GLM columns are hidden by default and then ignored by `Next` and `--strict`.

---
## Part 1 — Before you start

**Get the code** (once, on the HPC login node):
```bash
cd ~
git clone git@github.com:yunchiuuu/zwba.git zwba
cd ~/zwba            # every command in this README runs from here
```
To update later: `cd ~/zwba && git pull`.

**You need**
- a Resnick HPC account with access to `/resnick/groups/Proberlab/`
- the `voluseg` conda env (contains ANTsPy, voluseg, scipy, joblib, h5py, tifffile)
- raw volumes for your fish, acquired at the settings in Part 3

**Paths and names — what to change for your own account**

Everything lab-specific lives in the config; the group data folder is
`/resnick/groups/Proberlab/<your_name>/lightsheet/` (Yun's is `yun`, her HPC login is `ychiu`,
and the repo lives at `~/zwba`). Outputs for each person go under their own group folder.

| What | Where | Set to |
|---|---|---|
| `dir_voluseg`, `dir_registration`, `dir_analysis` | every `chemogenetic/config/config_*.py` | your `/resnick/groups/Proberlab/<name>/lightsheet/...` |
| `PYTHON_BIN` | every config (default `~/miniconda3/envs/voluseg/bin/python`, expands to your home) | your conda env, if it lives elsewhere |
| `regressor_files` | every config | the regressor CSVs you define (Part 7) |
| `dir_voluseg`, `dir_registration`, `PYTHON_BIN` | `registration/config/config_registration.py` | same as above |
| `BASE` | `run_voluseg.py` | your `lightsheet/` folder (also holds the shared `voluseg_setup/install/bin` ANTs install) |
| ANTs scratch dir | the three registration scripts | automatic: `/resnick/scratch/<your login>/tmp_ants` (`getpass.getuser()`) |

The submit scripts read `PYTHON_BIN` from the config you pass with `--config`; they stop with a
clear message if the config cannot be imported or the interpreter does not exist.

---

---

## Part 2 — Concepts you need

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

---

## Part 3 — Data and quantity definitions

### What your data must look like

The timing is read from the config (`volume_per_sec`, `baseline_start/end`, `drug_start/end`,
`wash_end`, `EPOCH_DISCARD_MIN`); nothing in the analysis assumes 1 volume/s or 120 min any more.
The 120-min CSN design is 1 vol/s, 45 min baseline / 45 min drug / 30 min washout (drug onset at vol 2700).
The scripts refuse to run when the recording is shorter than the windows need (`wash_end`), or when the
discard leaves an empty window.

Imaging constants in the config (`n_slices=40`, `depth=250`, `res_x/y=1.52`, `rotation_k=2`) are
still per-setup values; `run_voluseg.py` has its own copy at the top.

---

### How the quantities are defined (from the code)

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
`run_tonic_dz.py`, `run_figures_glm.py` and the GLM responder sign. All effect sizes are
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

The three responder methods built on these quantities are described in Part 4.


---

## Part 4 — Responder methods: F1, F2 and G (read this before running stages F–I)

Three methods call "responders". They answer **different questions**, use **different references**, and their responder
counts mean different things. Read the table and the notes of each method once; the notes are the place to come back to
when a result looks odd.

### 4.1 At a glance

| | **F1 — tonic ΔZ** (default) | **F2 — tonic GLM** (optional) | **G — phasic d′** |
|---|---|---|---|
| Question | Which cells shift their slow baseline level the most, in either direction? | Which cells' slow level follows the reference (HCRT) time course better than a drift-only model? | Which cells' fast activity (events) changes between baseline and drug? |
| Trace | `F_tonic` | `F_tonic` | `F_phasic` |
| Statistic | ΔZ per window; **max \|ΔZ\|** over the 7 windows, sign from the argmax window | ΔR² = R²(regressor + drift) − R²(drift only) | d′ = (drug mean − baseline mean) / pooled SD (baseline and drug epochs) |
| Supervised? | No | **Yes** (needs your regressor, Part 7) | No |
| Null / reference | **None**: rank within a σ stratum | IAAFT surrogates of the **regressor**, cells refit | IAAFT shuffle of the **baseline/drug label vector**, d′ recomputed |
| Unit the cutoff is calibrated on | **per fish** (strata of that fish's baseline σ) | **per fish**: one scalar T per fish | **per fish**: one scalar T per fish |
| Pooled over control fish? | No | No (not ctrl-pooled, not per cell) | No (the ctrl-pooled T is only printed and plotted as a reference) |
| Cutoff | top `TONIC_DZ_TOP_PERCENT` % (5) of each σ stratum | T = p`RESPONDER_NULL_THRESH` (99) of the pooled null ΔR² | T = p`PHASIC_RESPONDER_NULL_THRESH` (99) of the pooled \|null d′\| |
| Responder fraction | **fixed by p**, carries no information | ≈ (100 − p) % if there is no effect, so it can carry information | ≈ (100 − p) % if the null is right, so it can carry information; depends on T |
| Cost | seconds per fish | hours (24 h default walltime) | minutes per fish |
| Stage / figures | F1 / H1 | F2 / H2 | G / G2 |

The two nulls are **per fish and pooled over cells**: the surrogate × sampled-cell values of one fish are pooled into one
distribution and one scalar T per fish is read from it. Neither is estimated per cell (200 surrogates are far too few to estimate a
99th percentile per cell) and neither borrows from the control group.

### 4.2 F1 — tonic ΔZ, rank-based top-p % (default)

**What it does.** For each valid cell, `ΔZ_w = (mean F_tonic in window w − baseline mean) / baseline σ` for the seven temporal
windows (Part 3). `max ΔZ` is the ΔZ of the window with the largest |ΔZ|, keeping its sign. Valid cells are split into
`TONIC_DZ_N_STRATA` (5) quantile strata of baseline σ; within each stratum the top `TONIC_DZ_TOP_PERCENT` (5 %) of |max ΔZ| are the responders
(positive or negative by sign; one pooled ranking). Sensitivity sets at 2.5 % and 10 % are saved in the same run.

**Nature.** Per fish, rank-based, no null, no use of the control group. The responder *fraction* is p % of valid cells in every fish, so it is a
design constant. Controls therefore always have "responders".

**Pros**
- Unsupervised: makes no assumption about the time course of the effect, so it can find responses that are not HCRT-like.
- Cheap, deterministic, no regressor, no surrogates, nothing to tune except p and the number of strata.
- Equal set size in every fish, so set-based measures (region composition, overlap, maps) are not driven by one fish having more calls.
- The σ strata stop low-noise cells from filling the top set only by having a small denominator.

**Cons / what to watch**
- **Cannot say whether anything responded.** The count is fixed; compare groups with the per-fish ECDFs of |max ΔZ| (H1), the region composition of the top sets, or epoch maps. Fish are the unit.
- `max` is taken over all seven windows, including washout (windows 6–7) and window 5 (straddles drug end). A late whole-brain shift can set `max ΔZ`; check the argmax-window shares (printed on every run) and the H1 histogram.
- Baseline σ of `F_tonic` is mostly slow drift, so ΔZ is not a pure noise-normalised effect size and its scale can differ between fish for non-biological reasons (the strata are per fish, which absorbs the scale but not the shape).
- Tonic ΔZ is asymmetric (a 20th-percentile filter follows decreases and increases differently): positive and negative tails are not mirror images.
- In the first 16-fish run the maps were diffuse without regional specificity. That is a property of an unsupervised whole-distribution top-5 %, not necessarily a bug.
- *Not built yet:* a control-calibrated cutoff (a common σ-binned bar set by the control fish) that would make the fraction informative.

### 4.3 F2 — tonic GLM with a per-fish IAAFT null on the regressor (optional)

**What it does.** Each cell's `F_tonic` is fit with the reference regressor (HCRT trace; ridge, with a drift basis) over the fit window
`[drug_start − INCLUDED_BASELINE, end of recording)`. ΔR² = full model − drift-only model. For the null, the **regressor** is replaced by IAAFT
surrogates (same amplitude distribution and spectrum, timing scrambled), the sampled valid cells are refit, and ΔR² values are pooled over surrogates × cells.
T = p99 of that pooled distribution, **one scalar per fish**. Responders: ΔR² > T; the sign comes from the predicted response (drug window minus baseline window).

**Nature.** Supervised (needs the regressor), per fish, null pooled over cells, not control-pooled, not per cell.

**Pros**
- Looks for a *specific* time course, so it can give a sharper, more regionally specific pattern than an unsupervised ΔZ. In the first comparison the old GLM maps looked more distinct in HCRT-TRPV1 fish; that run predates the cell filter and the current fit window, so re-check with a fresh F2 run and the overlap figure.
- Has a false-positive control: under no effect about (100 − p) % of cells pass, so a control-group fraction near 1 % is a sanity check (H2).
- Drift is modelled explicitly instead of inflating σ.

**Cons / what to watch**
- **Circular for "does the drug drive HCRT-like activity"**: the regressor defines what a responder is, so the output cannot reveal responses with another time course.
- Depends on regressor choice and preprocessing, K, λ, lag, drift order, and the fit window (`INCLUDED_BASELINE`). The fit includes washout; changing the window requires a rebuild.
- Slow (24 h default) and the most parameters of the three; the null uses a sample of valid cells (`N_CELLS_NULL`), not all.
- One T per fish for cells of very different noise: lenient for noisy cells, strict for quiet ones.
- Because T is per fish, a fish with a strong real effect or different noise gets a different T and a different number of calls; the fraction is only comparable between groups if the T values are.
- Not directly comparable with the old (pre-filter) GLM maps; rebuild before comparing.

### 4.4 G — phasic d′ with a per-fish IAAFT null on the epoch label (default for phasic)

**What it does.** Per cell on `F_phasic`: `d′ = (drug mean − baseline mean) / pooled SD` (baseline epoch min 20–45, drug epoch min 65–90; variance floor 1e-4).
Every fish, ctrl and expt, gets its own null: the 0/1 baseline/drug label vector is IAAFT-shuffled (same counts, timing scrambled), d′ is recomputed on those labels
for `N_CELLS_NULL` (20000) sampled valid cells × `N_SURROGATES` (200) surrogates, and T = p`PHASIC_RESPONDER_NULL_THRESH` (99) of |null d′|, **one scalar per fish**.
A cell is a responder if its real d′ is above +T (pos) or below −T (neg) of **its own fish**. The ctrl-pooled threshold is printed and drawn only as a reference.

**Nature.** Unsupervised, per fish, null pooled over cells, not control-pooled, not per cell.

**Pros**
- No regressor and no assumption on the response shape beyond "drug epoch differs from baseline epoch".
- Calibrated by construction for a stationary trace: each fish's own noise sets its own bar, so differences in imaging noise do not create group differences in the false-positive rate.
- Two-sided and fast; null reused unless the cell filter or settings change; another percentile is derived from the saved null without re-running surrogates.

**Cons / what to watch** (all seen in the first 16-fish run)
- **The label-shuffle null assumes stationarity.** A whole-brain, onset-locked component (positive early, negative late in controls) is not reproduced by shuffling the labels, so control fish still get one-sided "responders" (1.4–3.8 % instead of ~1 %). Read control fish as the empirical floor, not as 1 %.
- **T can be inflated by real effects**: strong responders widen the d′ distribution the null is compared with; T looked higher in HCRT-TRPV1 fish (median ≈ 0.37 vs 0.30; Mann-Whitney p = 0.31, n = 6 vs 10), which would make the per-fish rule more conservative exactly where the effect is largest. Not significant, but check the thresholds figure (G2) before reading a fraction difference as biology.
- One T per fish while the per-cell null SD varies 4–7× (p95/p5): too lenient for noisy cells, too strict for quiet ones. A per-cell standardisation or a circular-shift null would address this; not built.
- The presence maps showed a bilateral anterolateral band of POS and NEG calls in both groups (possible edge / high-variance artefact); an edge-distance check is still open.
- Responder fractions are comparable between groups only as far as T is comparable.

### 4.5 Which to use, and how to compare groups

| You want to… | Use | Why |
|---|---|---|
| Describe where the biggest slow-level shifts are, without assuming a time course | F1 (maps, region composition, ECDFs) | unsupervised, equal set size |
| Test an HCRT-like tonic response | F2 | supervised; null gives a false-positive control |
| Describe fast-activity changes | G | the only method on `F_phasic` |
| Compare group means of a responder *fraction* | G (or F2) only, and only with the T values checked | F1's fraction is fixed |
| Check whether two tonic methods agree | `submit_tonic_overlap.sh` (Stage I) | enrichment over chance per fish and sign |

Rules for every method: **fish is the unit, never cells**; controls will have responders under all three (the question is the size and location of the difference);
tonic and phasic maps are separate files and are never merged; each method's sets are saved with its parameters in the file name so they can sit side by side.

*Discussed but not built:* a control-calibrated cutoff for ΔZ and a control-calibrated d′ threshold (leave-one-out for control fish, bootstrap CI on the bar,
a config switch between `ctrl_pooled`, `per_fish` and `fixed`); a global / local split of the whole-brain onset-locked component; a median-centred F1 sensitivity set.

---

## Part 5 — The stages, A to I

The letters match the Quick start table. Every command runs from the repo root; a stage that says *skips fish whose outputs exist* is safe to re-run.


### Before stage A — register your fish in two configs
1. `cp chemogenetic/config/config_template.py chemogenetic/config/config_<exp>.py`, then fill every `TODO`
   (project folders, fish tuples, `drug_uM`, epoch windows, group tags, `PYTHON_BIN`).
2. Add the same fish to `registration/config/config_registration.py` (`all_fish`).
   `submit_registration_syn.sh` submits **only uncommented** fish.

### Stage A — raw TIFF → HDF5 (analysis PC, not HPC)
Not in the repo yet (`run_preprocess_batch.py`, `tif_to_h5.py`). Output per fish:
`{dir_voluseg}/{proj}/{expt}/input/volume1000.h5, volume1001.h5, …` with dataset key
`'default'`. Run conversions **serially** (parallel conversion produced truncated files), then sync to the HPC.

### Stage B — voluseg (HPC, one job per fish)
```bash
sbatch --job-name=voluseg_<fish> --cpus-per-task=24 --mem=200G --time=8:00:00 \
  --output=logs/voluseg/<expt_ID>.log \
  --wrap="<PYTHON_BIN> run_voluseg.py --proj_ID <proj_ID> --expt_ID <expt_ID>"
```
Writes `output/cells0_clean.hdf5` and `output/volume0.hdf5`. Never run more than one fish per process (Spark).
Imaging constants (`FRAMES_PER_VOLUME=40`, 250 µm, 1.52 µm/px) are set at the top of `run_voluseg.py`.

### Stage C — registration to the shared template (HPC)
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

**Registrations made before `registration_meta.json` existed** (old fish show ⚠ in `status.sh`) do not need to be
redone if they used the current template. Verify and stamp them once:
```bash
# 1. when was the CURRENT template built? (mtime of the ORIGINAL file, e.g. zstack_template/template_mean_brain.nii.gz)
ls -l --time-style=long-iso <original template file>
# 2. dry run: shows the decision per fish and what it would write
python3 registration/backfill_registration_meta.py --config config_registration --registered-after "2026-07-11 19:01"
#    (the fish list is all_fish in config_registration; if you keep old fish commented out there, take it from
#     the analysis config instead: add  --fish-config config_hcrt_trpv1_csn_120min)
# 3. write the meta files
python3 registration/backfill_registration_meta.py --config config_registration --registered-after "2026-07-11 19:01" --write   # same options as the dry run
```
A fish is stamped only if its warp was written at or after that time and its job log (if there is one) loaded a
file with the same sha256 as the current template. Fish that fail either test are listed as BLOCKED and are left
alone, so `submit_registration_syn.sh` will re-register them. Stamped files carry `"backfilled": true` and the
evidence used. The check cannot prove the template identity, only make it very likely, so treat the stamp as a record
of your inference.
Per fish it writes `expt_to_mean_warp.nii.gz`,
`expt_to_mean_affine.mat`, a QC overlay PDF, and `registration_meta.json` (template sha256,
SyN parameters, git commit).

- **Do not regenerate the template.** `template_mean_brain.nii.gz` (version `zstack_v1`) is
  read-only and shared by all fish. `generate_template.py` and
  `register_mapzebrain_to_template.py` were one-time steps. If the template ever changes,
  re-register **all** fish and re-run medoids with overwrite.
- Open `qc_registration_overlay.pdf` for every fish before moving on.

### Stage D — decompose: F_tonic / F_phasic and the cell filter (HPC)
```bash
bash chemogenetic/submit/submit_decompose.sh --config config_<exp>
#   -> f_dark_scalar.npy, f_tonic.npy, f_phasic.npy, cell_valid.npy, QC_figures/decompose_QC.png
```
Skips fish whose outputs exist. To recompute, flip `OVERWRITE_*` at the top of the run script.
`run_decompose.py` checks the recording length against the epochs, and warns when many samples are
clipped at `F_dark` (>5 %) or `F_phasic` is non-finite (>2 %). **Open `decompose_QC.png` for every fish.**

**Cell filter (`cell_valid.npy`).** After the decomposition `run_decompose.py` writes a boolean mask per fish
(True = keep). `f_tonic.npy` / `f_phasic.npy` are never modified. If they already exist the decomposition is
skipped and only the mask is built, so re-running decompose on old fish is cheap. Every later stage applies the mask by
setting the excluded cells' results to NaN (rows are kept, so indices still match the medoids): ΔZ / d′ maps and
region surveys ignore them, they can never be responders, and the IAAFT nulls are sampled from valid cells only.
The criteria use the **baseline epoch only** (min 20-45 by default), so they cannot depend on the drug effect:

| Config setting (default) | Excludes a cell when |
|---|---|
| `CELL_MAX_ZERO_FRAC = 0.10` | more than 10 % of baseline samples have `F_tonic <= 0` (dim cells clipped by the `F_dark` subtraction) |
| `CELL_MIN_BASELINE_SD = 1e-3` | the baseline SD of `F_tonic` is below this (ΔZ divides by it, a flat baseline gives ΔZ ~ 1e7) |
| `CELL_MIN_BASELINE_MEAN = 1.0` | the baseline mean of `F_tonic` is below this many counts |
| non-finite values | any NaN / inf in the baseline epoch (always on) |
| `CELL_MAX_ABS_DZ = 200` | backstop for cells that passed the four checks above: max over the temporal windows of \|tonic ΔZ\| exceeds this (`None` = off) |

The log, `cell_valid_meta.json` and the status table report how many cells each criterion removed, and the ids of
cells caught only by the ΔZ backstop. `run_decompose.py` warns when more than 1 % of a fish's cells are excluded.
The thresholds are deliberately conservative (they target cells at the camera floor, not dim cells); check the counts
before tightening them. The mask is rebuilt automatically when the thresholds, the epochs or `f_tonic.npy` change, and
its fingerprint is stored with every downstream output, so results made with another mask are detected:
`d′`, its null and the maps recompute themselves; the GLM stops and prints
`--stages ablation,null,responders --overwrite ablation,null,responders` (the refit does not depend on the mask).
Not covered: the CV stage (`--stages cv`) still samples its 2000 hyperparameter cells from all cells.

### Stage E — medoids → template coordinates (HPC)
```bash
bash chemogenetic/submit/submit_medoids.sh --config config_<exp>      # needs stage C (registration)
#   -> medoids_all.npy, medoids_rot_all.npy, medoids_template_vox.npy
```
Skips fish whose outputs exist (`OVERWRITE_*` at the top of the run script to redo). Every map in stage I needs `medoids_template_vox.npy`.

### Stages F1, F2, G and their figures (HPC)

What each method is, its null and its pros/cons: **Part 4**. This section is the operational part (commands, outputs, caching, options).

#### Stage F1 — tonic ΔZ responders (default)
```bash
bash chemogenetic/submit/submit_tonic_dz.sh --config config_<exp>
```
Outputs per fish (`.npy`, in the fish analysis folder): `tonic_dz_windows.npy` (cells × 7), `tonic_dz_max.npy`,
`tonic_dz_argmax_window.npy` (−1 = excluded cell), `tonic_dz_sigma_baseline.npy`, `tonic_dz_stratum.npy`,
`tonic_{pos,neg}_dz_top{p}_idxs.npy` (`2p5` for 2.5) and `tonic_dz_meta.json` (cutoffs per stratum, inputs, counts).
A fish is skipped when its meta matches the current epochs, cell filter and `TONIC_DZ_*` settings, and recomputed otherwise;
`--overwrite` forces it. Also takes `--fish TEXT`, `--exclude TEXT` (repeatable) and `--n-jobs`.

**Check the argmax window.** Every run prints a per-fish table of which window each cell peaks in (read from the saved files, so skipped fish are included;
a re-run overwrites the sbatch log, so this table is what to read) and warns for each fish with more than half of its cells peaking in the last two (washout) windows.
If a washout or late whole-brain shift dominates, the top set is "late-shifted cells", not "drug responders".

#### Stage F2 — tonic GLM + IAAFT null (optional)
Not part of the default run; needed for the GLM figures (H2), the GLM presence map and the kernel map. Needs the regressor in the config (Part 7).
```bash
bash chemogenetic/submit/submit_glm.sh --config config_<exp>      # refit, ablation, null, responders; default walltime 24 h
```
`run_glm.py` takes `--stages`, `--overwrite`, `--exclude TEXT` and `--fish TEXT`; the submit script passes them through:
```bash
bash chemogenetic/submit/submit_glm.sh --config config_<exp> --stages responders             # redo only the responder call
bash chemogenetic/submit/submit_glm.sh --config config_<exp> --stages refit,ablation,null,responders --overwrite refit,ablation,null,responders
bash chemogenetic/submit/submit_glm.sh --config config_<exp> --exclude 260824 --exclude 260825   # leave fish out (text from their expt_ID)
```
The GLM needs `f_tonic.npy` for every fish in the config and stops (exit 2, nothing computed) if some are missing; fish are fitted independently.
Outputs land in `{dir_analysis}/{proj}/{expt}/glm/in{input_tag}_K{K}_drift{d}_lam{λ}_lag{l}/`, plus `tonic_{pos,neg}_glm_iaaft_nullp{N}_idxs.npy`.

**After changing `INCLUDED_BASELINE`** (or the regressor), the fits must be rebuilt. The GLM folder name does not
encode the fit window, so `run_glm.py` reads `RUN_META.json` and stops if existing fits used a different window.
```bash
bash chemogenetic/submit/submit_glm.sh --config config_<exp> \
  --stages cv,refit,ablation,null,responders --overwrite refit,ablation,null,responders
```
(`cv` re-checks K/λ for the new window; the globals in the config are still what the refit uses.)
Then rerun `submit_figures_glm.sh` and the maps. Use a different `param_folder_name` if you want to keep the old results.

#### Stage G — phasic d′ + per-fish IAAFT null (default)
```bash
bash chemogenetic/submit/submit_dprime.sh --config config_<exp>    # stages: dprime, null, responders, figures
bash chemogenetic/submit/submit_dprime.sh --config config_<exp> --stages responders,figures
```
Outputs per fish: `phasic_dprime_cells_*.npy`, `phasic_dprime_null__iaaft.npy` (surrogates × cells),
`phasic_dprime_iaaft_nullp{N}_thresh.npy` (T), `phasic_dprime_null_meta.json` (cell filter sha + null settings),
`phasic_{pos,neg}_dprime_iaaft_nullp{N}_idxs.npy` and `phasic_dprime_responders_meta.json` (T, percentile, counts).
A null is reused only if its cell filter and settings match; a new percentile is derived from the saved null without re-running surrogates;
the responders stop with a message if the null is stale. The log lists every fish's T: check they are in a similar range (Part 4.4).

#### Stage G2 — d′ null diagnostics
```bash
bash chemogenetic/submit/submit_figures_dprime_null.sh --config config_<exp>    # after G; reads the saved files
```
`phasic_dprime_null_thresholds_p{N}.png` (T per fish, T vs the width of real d′, fraction vs T), `phasic_dprime_null_vs_real_p{N}.png`
(per fish: real vs null |d′| survival curves) and `phasic_dprime_null_heterogeneity_p{N}.png` (spread of the per-cell null SD with a split-half
reliability; p95/p5 above ~3 means one T per fish is too lenient for noisy cells and too strict for quiet ones).

#### Stage H1 — tonic ΔZ figures
```bash
bash chemogenetic/submit/submit_figures_tonic_dz.sh --config config_<exp>      # after F1; open all three PNGs
```
Argmax-window histogram (is the top set dominated by washout?), σ-vs-|max ΔZ| density (flat-baseline cloud?), and per-fish ECDFs of |max ΔZ| and signed max ΔZ,
CTRL vs EXPT (this, not the responder count, is the group comparison).

#### Stage H2 — GLM figures (optional)
```bash
bash chemogenetic/submit/submit_figures_glm.sh --config config_<exp>           # after F2
```
Mean traces, responder-fraction boxplot, ΔZ amplitude boxplot. The control responder fraction should be near 1 − p/100.

#### Run options shared by the stage scripts

| Script | Valid stages | Default `--stages` | Default `--overwrite` |
|---|---|---|---|
| `run_glm.py` (F2) | `cv, refit, ablation, null, responders` | `refit,ablation,null,responders` | `responders` |
| `run_dprime.py` (G) | `dprime, null, responders, figures` | all four | `responders` |
| `run_tonic_dz.py` (F1) | (single stage) | — | skips current fish; `--overwrite` redoes all |

Existing outputs are skipped unless listed in `--overwrite`. The responder call uses `RESPONDER_NULL_THRESH` (F2) and `PHASIC_RESPONDER_NULL_THRESH` (G)
from the config: both are IAAFT null percentiles you can change (default 99); if only another percentile's threshold file exists, it is derived from the saved null (no surrogates re-run).

### Stage I — maps (all need `medoids_template_vox.npy`)
```bash
bash chemogenetic/submit/submit_temporal_intensity_map.sh --config config_<exp>   # epoch_dz.npy, epoch_dprime.npy, per-window maps
bash chemogenetic/submit/submit_glm_dprime_map.sh        --config config_<exp>   # N-fish presence maps: ΔZ (F1), GLM (F2, if run), phasic d′ (G), one file each
bash chemogenetic/submit/submit_tonic_overlap.sh         --config config_<exp>   # optional: per-fish overlap of the ΔZ and GLM sets (needs F1 and F2)
bash chemogenetic/submit/submit_kernel_map.sh            --config config_<exp>   # kernel maps (needs F2)
bash chemogenetic/submit/submit_combined_fish_map.sh     --config config_hcrt_trpv1_pooled_csn_120min
bash chemogenetic/submit/submit_survey_by_regions.sh [--config config_<exp>] [--overwrite] [--top_n N]
```
**Presence maps** (`run_glm_dprime_map.py`): one file per method, parameters in the name, never overwriting each other:

| File | Method |
|---|---|
| `tonic_dz_presence_map_top{p}.png` | F1 |
| `tonic_glm_presence_map_{NULL_TAG}_nullp{N}.png` | F2 |
| `phasic_dprime_presence_map_{NULL_TAG}_nullp{N}.png` | G |

Each is made in two steps: the per-voxel fish counts are saved as `<same name>_counts.npz`, then the PNG is drawn from it.
Options (the submit script passes them through): `--tonic-source dz|glm|both` (default `both`; a method with no per-fish files is skipped with a message that names the stage to run)
and `--plot-only` (re-draw from the saved counts). The old `tonic_presence_map.png` / `phasic_presence_map.png` on disk are stale; delete them.

**Overlap** (`run_tonic_overlap.py`) writes `tonic_method_overlap_dz_top{p}_vs_glm_…` (`.npz`, `.csv`, `.png`): per fish and per sign the set sizes, Jaccard and the enrichment over
chance (1 = chance; Jaccard alone is capped because the ΔZ set size is fixed by p), plus the number of sign-conflicting cells (POS in one method, NEG in the other).

`submit_survey_by_regions.sh` without `--config` submits the three configs listed at the top of the script;
any other arguments are passed to `run_survey_by_regions.py`. It needs `epoch_*.npy` and the MapZebrain masks
in `{dir_registration}/300um_template/mapzebrain_in_template/`.
`run_glm_brainmap.py` has no submit script. `submit_individual_fish_map.sh` points to a run script that is not in
the repo (the script now stops with a message).

Figures go to `{dir_analysis}/{proj}/{expt}/figures/` (per fish) and
`{dir_analysis}/comparisons/{COMPARISON_TAG}/figures/` (group level).

### When something fails (all stages)
Every run script uses the shared guard in `utils/run_guard.py`:
- **Config / input problems** are checked before heavy work and reported all at once (every missing
  config key, every missing file). They print `❌ <message>` and exit with code **2**.
- **Per-fish failures** are caught, printed with a traceback, and the batch continues. At the end of
  the log a summary lists every failed item, and the job exits with code **1** so `sacct` shows it as `FAILED`.
- `bash chemogenetic/submit/status.sh --config config_<exp>` shows which fish still have missing or stale outputs.
- A clean run ends with `✅ complete: N ok, 0 failed`. `grep -n "❌" logs/<stage>/<config>.log` shows all failures.

---

---

## Part 6 — Interpreting results: current caveats

- **Baseline contamination.** The centered `F_tonic` filter can pull drug-period
  fluorescence into the last part of the baseline (within half a filter window of onset).
  This inflates `μ_baseline` / `σ_baseline`.
- **Edge effects.** Reflect padding makes `F_tonic` near the end of the recording
  non-independent; the last windows (W5, W6) are most affected.
- **Tonic σ_baseline is mostly slow drift**, so tonic ΔZ is not a pure noise-normalised
  effect size and can differ between fish for non-biological reasons.
- **Washout windows are inside the |ΔZ| max.** Windows 6–7 (95–120 min) are washout and window 5 straddles drug end; a late shift can set `max ΔZ`. Check the argmax-window histogram (H1) per fish.
- **Tonic ΔZ is asymmetric.** A 20th-percentile `F_tonic` follows decreases and increases differently, so positive and negative tails are not mirror images; compare signed `max ΔZ` ECDFs, not only |ΔZ|.
- **Flat baselines inflate ΔZ.** The σ strata stop low-σ cells from filling the top set only if σ is not itself ~0; the cell filter removes the degenerate ones, the σ-vs-|ΔZ| plot shows the rest. Dim, edge and surface cells have the lowest σ.
- **Top-p % always returns p %.** Controls will have responders too; judge a group by the ECDF shift and by which regions the top sets sit in, never by the count.
- **Signed averaging cancels** at the voxel mean, the cross-fish mean and the Z/X projection.
  Strong responders of opposite sign can look flat in a mean map.
- **No N-support masking in difference maps.** A voxel backed by 2 fish looks like one
  backed by 9, and the EXPT − CTRL map is NaN wherever either group is empty.
- **`clip_min=0` in background subtraction** truncates the lower tail of dim cells near `F_dark`.

---

---

## Part 7 — The GLM regressor (optional stage F2 only; you define it in the config)

The default pipeline (F1) does not use it. For F2 there is no fallback regressor. Set in your config:

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

---

## Part 8 — Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `could not read PYTHON_BIN from chemogenetic/config/...` | run from the repo root and pass the config module name; your config needs `PYTHON_BIN` |
| `PYTHON_BIN does not exist` | point `PYTHON_BIN` in the config at your conda env's `bin/python` |
| `❌ <config> is missing required setting(s)` | copy the listed names from `config_template.py` (e.g. `EPOCH_DISCARD_MIN`, `regressor_files`) |
| `TONIC_DZ_N_STRATA must be >= 1` (or other `TONIC_DZ_*` error) | fix the setting in the config; valid: `0 < TONIC_DZ_TOP_PERCENT < 50`, `TONIC_DZ_N_STRATA >= 1` |
| `no usable cells` for a fish (F1) | the cell filter removed everything: check `cell_valid_meta.json` and the `CELL_*` thresholds |
| `⚠ more than 50 % of cells peak in the last two windows` (F1) | washout / late shift dominates `max ΔZ`; inspect the argmax histogram and the fish's traces before using the sets |
| `the d′ null was sampled without / with a different cell filter` (d′ responders) | run `submit_dprime.sh` (stages `null` and `responders`); the stale null is recomputed |
| Status `dres` ⚠: responders called with the old ctrl-pooled threshold | rerun `submit_dprime.sh`; responders are redone with each fish's own T |
| `Cannot load the GLM regressor` (F2) | the message lists each bad file (missing / no column / NaN / too short) |
| `recording has N volumes but epochs need M` | wrong fish in the list or wrong `wash_end`/timing in the config |
| Job `FAILED` in `sacct` | exit code 1 = some fish failed (see the summary at the end of the log); exit code 2 = bad config/input |
| Registration submit says `--config is required` | pass a module name: `--config config_registration` |
| Registration submits 0 jobs | all fish are already registered with the current template (see `--dry-run`), or commented out in `config_registration.py` |
| Registered fish look misaligned | check `registration_meta.json` template sha256 vs. the current template; check the affine z-scale in the QC overlay (values far from ~1.0 are suspect) |
| `Existing GLM outputs were fit with a different window` (F2) | `INCLUDED_BASELINE` changed since those fits; rebuild as in Stage F2 |
| `cell_valid.npy is missing` | run `run_decompose.py` for that fish (existing `f_tonic` / `f_phasic` are skipped; only the filter is built) |
| `Existing GLM outputs were made without (or with a different) cell filter` | the mask changed or the GLM predates it: `run_glm.py --stages ablation,null,responders --overwrite ablation,null,responders` |
| Cell filter excludes > 1 % of cells | check `decompose_QC.png` and the `F_dark` clipping warning; thresholds are the `CELL_*` settings in the config |
| GLM results differ from an older run | the fit window, regressor or threshold changed — rerun with `--stages refit,ablation,null,responders --overwrite refit,ablation,null,responders` |

---

---

## Part 9 — Repo layout

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
  cell_filter.py                valid-cell mask (cell_valid.npy) used by every stage after decompose
  tonic_dz.py                   tonic ΔZ responders (stage F1): library, per-fish outputs, status check
  dprime.py                     phasic d′, the per-fish IAAFT null, responder calls (stage G)
  glm.py dprime.py temporal_windows.py drug_regressor.py brainmap.py kernel_map.py population.py scatter.py traces.py
utils/  data_io.py  preprocess.py  stats.py  run_guard.py   (run_guard: failure railings used by every run script)
```

Legacy Spearman / permutation / effect-size pipeline was removed (superseded by GLM + IAAFT);
recover from git tag `pre-spearman-removal`.

## Part 10 — Maintainer TODO

- Commit the missing preprocessing scripts (`run_preprocess_batch.py`, `tif_to_h5.py`, `apply_voluseg_patch.py`) and `run_individual_fish_map.py`.
- Read the remaining hardcoded timing in `run_combined_fish_map.py` / `run_kernel_map.py` / `run_glm_brainmap.py` from the config (they only load existing outputs, but check before using a different design).
- Add `--overwrite` CLI flags to `run_decompose.py` and `run_medoids.py` (still `OVERWRITE_*` constants).
- Add N-support masks to the difference maps.
- Stop tracking `__pycache__/`, `.ipynb_checkpoints/` and `logs/` (`git rm -r --cached`; `.gitignore` is updated).