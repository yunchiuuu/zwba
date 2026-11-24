"""
run_glm.py
==========
SLURM entry point: run the full GLM pipeline for all fish.

Stages:
    A. CV            — per-fish hyperparameter selection
    B. Choose global — aggregate CV results
    C. Refit         — refit all cells with global params
    D. Ablation      — drift-only ΔR²
    E. IAAFT null    — surrogate null distribution (now parallel across fish;
                        BLAS threads capped per worker to avoid oversubscription)
    F. Responders    — threshold + save pos/neg idx

Regressor: defined in the config (regressor_files / regressor_column / ...), loaded
           once and passed as u_ext to all stages. There is NO fallback: if it cannot be
           read the run stops with a message (see chemogenetic/regressor.py).

Stages are selected on the command line (default: refit,ablation,null,responders;
existing outputs are skipped unless listed in --overwrite):
    --stages cv,refit,ablation,null,responders
    --overwrite responders            (default)

Usage:
    bash chemogenetic/submit/submit_glm.sh --config config_hcrt_trpv1_csn_120min
    python chemogenetic/run/run_glm.py --config config_x --stages responders --overwrite responders

Failures: every fish/stage is reported with a traceback; the job exits 1 if any failed.

Location:
    ~/zwba/chemogenetic/run/run_glm.py
"""

import argparse
import json
import os
import sys
from pathlib import Path

# ── cap BLAS threads per worker BEFORE numpy import ─────────────────────────
# With N_JOBS_STAGE_E parallel fish workers each doing ridge regression,
# uncapped BLAS threading causes massive core oversubscription and can be
# SLOWER than serial. Set so N_JOBS_STAGE_E * BLAS_THREADS_PER_WORKER ≈
# total allocated CPUs (see submit_glm.sh --cpus-per-task).
BLAS_THREADS_PER_WORKER = "4"
os.environ["OMP_NUM_THREADS"]        = BLAS_THREADS_PER_WORKER
os.environ["MKL_NUM_THREADS"]        = BLAS_THREADS_PER_WORKER
os.environ["OPENBLAS_NUM_THREADS"]   = BLAS_THREADS_PER_WORKER
os.environ["NUMEXPR_NUM_THREADS"]    = BLAS_THREADS_PER_WORKER
os.environ["VECLIB_MAXIMUM_THREADS"] = BLAS_THREADS_PER_WORKER

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

# ── CLI ───────────────────────────────────────────────────────────────────────
VALID_STAGES = ("cv", "refit", "ablation", "null", "responders")
parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True,
                    help="Config module under chemogenetic/config/")
parser.add_argument("--stages", default="refit,ablation,null,responders",
                    help=f"comma list from {VALID_STAGES}; cv is only needed to pick new "
                         f"global hyperparameters (otherwise the config globals are used)")
parser.add_argument("--overwrite", default="responders",
                    help="comma list of stages to recompute even if outputs exist")
args = parser.parse_args()


def _stage_set(text, flag):
    out = {t.strip() for t in text.split(",") if t.strip()}
    bad = out - set(VALID_STAGES)
    if bad:
        parser.error(f"{flag}: unknown stage(s) {sorted(bad)}; valid: {VALID_STAGES}")
    return out


STAGES = _stage_set(args.stages, "--stages")
OVERWRITE_STAGES = _stage_set(args.overwrite, "--overwrite")

import importlib
cfg = importlib.import_module(f"chemogenetic.config.{args.config}")

all_fish         = cfg.all_fish
dir_analysis     = cfg.dir_analysis
drug_end         = cfg.drug_end
drug_start       = cfg.drug_start
drug_uM          = cfg.drug_uM
INCLUDED_BASELINE = cfg.INCLUDED_BASELINE
input_tag        = cfg.input_tag
K_global         = cfg.K_global
drift_global     = cfg.drift_global
lam_global       = cfg.lam_global
lag_global       = cfg.lag_global
NULL_TAG         = cfg.NULL_TAG
param_folder_name = cfg.param_folder_name
Q_ml_min         = cfg.Q_ml_min
RESPONDER_NULL_THRESH = cfg.RESPONDER_NULL_THRESH
sampling_rate_hz = cfg.sampling_rate_hz
V_ml             = cfg.V_ml

from chemogenetic.glm import (
    choose_global_params,
    cv_one_fish,
    refit_one_fish,
    ablation_one_fish,
    iaaft_null_one_fish,
    save_responder_idx,
    get_fit_window_indices,
    get_run_dir,
)
from chemogenetic.cell_filter import ensure_cell_valid, filter_params
from chemogenetic.epochs import get_epochs
from chemogenetic.regressor import build_regressor, check_regressor_vs_recording
from utils.data_io import fish_dir
from utils.run_guard import (
    PipelineError, RunGuard, guarded, require_paths, run_main, setup_or_exit,
)

# ============================================================
# STAGES / SETTINGS (stages and overwrite come from the command line)
# ============================================================
RUN_CV         = "cv" in STAGES
RUN_REFIT      = "refit" in STAGES
RUN_ABLATION   = "ablation" in STAGES
RUN_NULL       = "null" in STAGES
RUN_RESPONDERS = "responders" in STAGES

OVERWRITE_CV         = "cv" in OVERWRITE_STAGES
OVERWRITE_REFIT      = "refit" in OVERWRITE_STAGES
OVERWRITE_ABLATION   = "ablation" in OVERWRITE_STAGES
OVERWRITE_NULL       = "null" in OVERWRITE_STAGES
OVERWRITE_RESPONDERS = "responders" in OVERWRITE_STAGES

K_LIST        = (60, 120, 300, 600, 900, 1200)
DRIFT_ORDERS  = (1, 2, 3)
LAM_LIST      = (1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3)
BLOCK_LEN_SEC = 10 * 60
N_CELLS_CV    = 2000
CHUNK_CELLS   = 2000

N_SURROGATES  = 200
N_ITER_IAAFT  = 50
N_CELLS_NULL  = 20000

N_JOBS_FISH = 16

# Stage E (IAAFT null) is parallel across fish. Keep this modest since each worker also
# spawns BLAS_THREADS_PER_WORKER threads internally —
# N_JOBS_STAGE_E * int(BLAS_THREADS_PER_WORKER) should stay close to your allocated
# --cpus-per-task in submit_glm.sh (e.g. 8 * 4 = 32).
N_JOBS_STAGE_E = 8

FIT_BASELINE_SEC = INCLUDED_BASELINE * 60

# ============================================================
# EPOCHS for the responder sign (Stage F)
# ============================================================
# Sign = mean(yhat in drug window) - mean(yhat in baseline window). Both windows come
# from the config epochs (first EPOCH_DISCARD_MIN of each phase discarded) and are
# converted to minutes RELATIVE TO THE START OF THE FIT WINDOW, which is what
# save_responder_idx expects.
EPOCHS = setup_or_exit(get_epochs, cfg)


def _fit_relative_windows(epochs):
    fit_start, _, _ = get_fit_window_indices(
        int(cfg.wash_end), sampling_rate_hz, drug_start, FIT_BASELINE_SEC)
    b0, b1 = epochs.baseline
    d0, d1 = epochs.drug
    if b0 < fit_start:
        raise PipelineError(
            f"The GLM fit window starts at vol {fit_start} (INCLUDED_BASELINE="
            f"{INCLUDED_BASELINE:g} min before drug onset) but the baseline summary "
            f"window starts at vol {b0}. Increase INCLUDED_BASELINE to at least "
            f"{(drug_start - b0) / (60 * sampling_rate_hz):g} min.")
    if fit_start < b0:
        print(f"\n⚠️  INCLUDED_BASELINE={INCLUDED_BASELINE:g} min starts the GLM fit at vol {fit_start}, "
              f"before the baseline summary window (vol {b0}): the discarded first "
              f"{epochs.minutes(epochs.discard_frames):g} min enter the fit. Use "
              f"{(drug_start - b0) / (60 * sampling_rate_hz):g} to follow the epoch rule.\n",
              file=sys.stderr)
    to_min = lambda f: (f - fit_start) / (60.0 * sampling_rate_hz)   # noqa: E731
    return (to_min(b0), to_min(b1)), (to_min(d0), to_min(d1))


BASELINE_WIN_MIN, DRUG_WIN_MIN = setup_or_exit(_fit_relative_windows, EPOCHS)

# ============================================================
# REGRESSOR (user-defined in the config; no fallback)
# ============================================================
if cfg.input_tag != "HCRT":
    print(f"\n❌ config input_tag is {cfg.input_tag!r}; it must be 'HCRT' (the label of the "
          f"user-defined regressor, it names the glm/inHCRT_... output folder).\n",
          file=sys.stderr)
    sys.exit(2)

u_ext = setup_or_exit(build_regressor, cfg)
print(f"✅ regressor loaded from {len(cfg.regressor_files)} file(s): {u_ext.size} samples, "
      f"peak={u_ext.max():.2f} at vol {u_ext.argmax()}")


# ============================================================
# HELPERS
# ============================================================
def _common_kwargs(fish):
    return dict(
        fish=fish,
        dir_analysis=dir_analysis,
        sampling_rate_hz=sampling_rate_hz,
        drug_start_frame_full=drug_start,
        drug_end_frame_full=drug_end,
        drug_uM=drug_uM,
        V_ml=V_ml,
        Q_ml_min=Q_ml_min,
        input_tag=input_tag,
        lag_global=lag_global,
        fit_baseline_sec=FIT_BASELINE_SEC,
        chunk_cells_fit=CHUNK_CELLS,
        u_ext=u_ext,
    )


# ============================================================
# MAIN
# ============================================================
def _preflight():
    """Fail BEFORE any compute, listing every problem at once."""
    require_paths(
        [(fish_dir(dir_analysis, f) / "f_tonic.npy", "run run_decompose.py first")
         for f in all_fish],
        "GLM input F_tonic",
    )
    for f in all_fish:
        T = np.load(str(fish_dir(dir_analysis, f) / "f_tonic.npy"), mmap_mode="r").shape[1]
        EPOCHS.check_recording(T, f"{f[1]} f_tonic")
        check_regressor_vs_recording(u_ext, T, f[1])
    # Stale-fit guard: the GLM folder name does not encode the fit window, so outputs made with a
    # different INCLUDED_BASELINE would otherwise be reused silently.
    exp_fit_start, _, _ = get_fit_window_indices(
        int(cfg.wash_end), sampling_rate_hz, drug_start, FIT_BASELINE_SEC)
    rebuilt = {"refit", "ablation", "null"} <= OVERWRITE_STAGES
    stale = []
    for f in all_fish:
        _, rd = get_run_dir(dir_analysis, f, K_global, drift_global, lam_global,
                            input_tag, lag_global)
        mp = rd / "RUN_META.json"
        if mp.exists():
            try:
                old = int(json.loads(mp.read_text()).get("fit_start", exp_fit_start))
            except Exception:    # noqa: BLE001 - unreadable meta: treat as unknown, not stale
                continue
            if old != exp_fit_start:
                stale.append(f"    {f[1]}: existing fit starts at vol {old}, config gives vol {exp_fit_start}")
    # Cell filter: cell_valid.npy is normally written by run_decompose.py; built here if missing/stale.
    # Fit outputs made with a different mask are detected (the GLM folder name does not encode it).
    mask_sha = {}
    for f in all_fish:
        _, sha, info = ensure_cell_valid(fish_dir(dir_analysis, f), EPOCHS, cfg, verbose=False)
        mask_sha[f] = sha
        print(f"  cell filter {f[1]}: excluded {info['n_excluded']:,}/{info['n_cells']:,}")
    mask_stale = []
    for f in all_fish:
        _, rd = get_run_dir(dir_analysis, f, K_global, drift_global, lam_global,
                            input_tag, lag_global)
        for tag, name in (("ablation", "RUN_META_ABLATION.json"), ("null", f"RUN_META_NULL__{NULL_TAG}.json")):
            mp = rd / name
            if mp.exists():
                try:
                    made_with = json.loads(mp.read_text()).get("cell_valid_sha")
                except Exception:    # noqa: BLE001
                    made_with = None
                if made_with != mask_sha[f]:
                    mask_stale.append(f"    {f[1]}: {tag} made without / with a different cell filter")
    if mask_stale and not ({"ablation", "null"} <= OVERWRITE_STAGES and RUN_ABLATION and RUN_NULL):
        raise PipelineError(
            "Existing GLM outputs were made without (or with a different) cell filter (cell_valid.npy):\n"
            + "\n".join(mask_stale) + "\n  Recompute the mask-dependent stages (refit is not affected): "
            "--stages ablation,null,responders --overwrite ablation,null,responders")
    if stale and not rebuilt:
        raise PipelineError(
            "Existing GLM outputs were fit with a different window (INCLUDED_BASELINE changed):\n"
            + "\n".join(stale) + "\n  Rebuild them: --stages cv,refit,ablation,null,responders "
            "--overwrite refit,ablation,null,responders  (or point param_folder_name elsewhere "
            "to keep the old results).")
    if RUN_RESPONDERS and not (RUN_REFIT and RUN_ABLATION):
        # fit outputs are not being (re)built in this run: they must already exist
        need = []
        for f in all_fish:
            _, rd = get_run_dir(dir_analysis, f, K_global, drift_global, lam_global,
                                input_tag, lag_global)
            need.append((rd / "kernel_delta_r2_fit.npy",
                         "run: run_glm.py --stages refit,ablation first"))
        require_paths(need, "GLM fit outputs needed by --stages responders")


def main():
    guard = RunGuard("glm")
    print(f"GLM pipeline | {len(all_fish)} fish")
    print(f"  config={args.config}  input_tag={input_tag}  stages={sorted(STAGES)}  "
          f"overwrite={sorted(OVERWRITE_STAGES)}")
    print(f"  K={K_global}  drift={drift_global}  lam={lam_global}  "
          f"lag={lag_global}  null={NULL_TAG} p{RESPONDER_NULL_THRESH}")
    print(f"  epochs: {EPOCHS.describe()}")
    print(f"  cell filter: {filter_params(cfg)}")
    print(f"  responder sign windows (min from fit start): "
          f"baseline {BASELINE_WIN_MIN[0]:g}-{BASELINE_WIN_MIN[1]:g}, "
          f"drug {DRUG_WIN_MIN[0]:g}-{DRUG_WIN_MIN[1]:g}\n")

    _preflight()

    # ── A: CV ──────────────────────────────────────────────────
    if RUN_CV:
        print("── Stage A: CV ─────────────────────────────────────")

        @guarded
        def _cv(fish):
            return cv_one_fish(
                **_common_kwargs(fish),
                K_list=K_LIST, drift_orders=DRIFT_ORDERS, lam_list=LAM_LIST,
                block_len_sec=BLOCK_LEN_SEC, n_cells_cv=N_CELLS_CV,
                overwrite=OVERWRITE_CV, show_progress=False, verbose_cv=False,
            )

        results = Parallel(n_jobs=N_JOBS_FISH, backend="loky")(
            delayed(_cv)(fish) for fish in tqdm(all_fish, desc="CV")
        )
        for r in results:
            print(f"  {r['fish'][1]:50s}  {r['status']}")
        guard.record_results("cv", results)

    # ── B: Choose globals ──────────────────────────────────────
    print("\n── Stage B: Using config globals ───────────────────────")
    print(f"  K={K_global}  drift={drift_global}  lam={lam_global}")

    # ── C: Refit ───────────────────────────────────────────────
    if RUN_REFIT:
        print("\n── Stage C: Global refit ───────────────────────────────")

        @guarded
        def _refit(fish):
            return refit_one_fish(
                **_common_kwargs(fish),
                K_global=K_global, drift_global=drift_global, lam_global=lam_global,
                overwrite=OVERWRITE_REFIT, show_progress=False,
            )

        results = Parallel(n_jobs=N_JOBS_FISH, backend="loky")(
            delayed(_refit)(fish) for fish in tqdm(all_fish, desc="Refit")
        )
        for r in results:
            print(f"  {r['fish'][1]:50s}  {r['status']}")
        guard.record_results("refit", results)

    # ── D: Ablation ────────────────────────────────────────────
    if RUN_ABLATION:
        print("\n── Stage D: Drift-only ablation ────────────────────────")

        @guarded
        def _ablation(fish):
            return ablation_one_fish(
                **_common_kwargs(fish),
                K_global=K_global, drift_global=drift_global, lam_global=lam_global,
                overwrite=OVERWRITE_ABLATION, show_progress=False,
            )

        results = Parallel(n_jobs=N_JOBS_FISH, backend="loky")(
            delayed(_ablation)(fish) for fish in tqdm(all_fish, desc="Ablation")
        )
        for r in results:
            print(f"  {r['fish'][1]:50s}  {r['status']}")
        guard.record_results("ablation", results)

    # ── E: IAAFT null ──────────────────────────────────────────
    if RUN_NULL:
        print(f"\n── Stage E: IAAFT null (parallel, n_jobs={N_JOBS_STAGE_E}, "
              f"{BLAS_THREADS_PER_WORKER} BLAS threads/worker) ──────────")

        @guarded
        def _null(fish):
            result = iaaft_null_one_fish(
                **_common_kwargs(fish),
                K_global=K_global, drift_global=drift_global, lam_global=lam_global,
                n_surrogates=N_SURROGATES, n_iter_iaaft=N_ITER_IAAFT,
                n_cells_null=N_CELLS_NULL,
                null_percentile=RESPONDER_NULL_THRESH,
                overwrite=OVERWRITE_NULL, show_progress=False,
            )
            return {"fish": fish, "status": result["status"]}

        results = Parallel(n_jobs=N_JOBS_STAGE_E, backend="loky")(
            delayed(_null)(fish) for fish in tqdm(all_fish, desc="IAAFT null")
        )
        for r in results:
            print(f"  {r['fish'][1]:50s}  {r['status']}")
        guard.record_results("null", results)

    # ── F: Responder indices ───────────────────────────────────
    if RUN_RESPONDERS:
        print("\n── Stage F: Save responder indices ─────────────────────")

        @guarded
        def _responders(fish):
            return save_responder_idx(
                fish=fish,
                dir_analysis=dir_analysis,
                sampling_rate_hz=sampling_rate_hz,
                K_global=K_global, drift_global=drift_global, lam_global=lam_global,
                input_tag=input_tag, lag_global=lag_global,
                null_tag=NULL_TAG,
                null_percentile=RESPONDER_NULL_THRESH,
                baseline_win_min=BASELINE_WIN_MIN,
                drug_win_min=DRUG_WIN_MIN,
                u_ext=u_ext,
                overwrite=OVERWRITE_RESPONDERS,
            )

        results = Parallel(n_jobs=N_JOBS_FISH, backend="loky")(
            delayed(_responders)(fish) for fish in tqdm(all_fish, desc="Responders")
        )
        for r in results:
            print(f"  {r['fish'][1]:50s}  {r['status']}")
        guard.record_results("responders", results)
        for r in results:
            if r.get("n_resp") == 0 or r.get("n_pos", 1) + r.get("n_neg", 1) == 0:
                guard.warn(r["fish"][1], f"0 responders at {NULL_TAG} p{RESPONDER_NULL_THRESH} "
                                         f"— check the regressor and the null threshold")

    guard.finish()


if __name__ == "__main__":
    run_main(main)
