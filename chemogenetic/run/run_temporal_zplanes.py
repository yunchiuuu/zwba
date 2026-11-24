"""
run_temporal_zplanes.py   (stage I, optional)
=============================================
The temporal intensity maps (run_temporal_intensity_map.py) average over all z planes in the XY panel, so a
signed effect that has opposite sign in different planes cancels. This script draws the SAME figure (HCRT reference
row, CTRL | EXPT | EXPT − CTRL, seven windows, same colour scale) for single z planes, one PDF page per plane.

Two steps (both run by default):
    1. calculate: per-fish epoch ΔZ / d′ are voxelized, averaged per group and saved as the full 3-D coarse volumes
           temporal_group_volumes_DS{dx}x{dy}x{dz}_clip{C}.npz      (comparisons/<COMPARISON_TAG>/figures/)
    2. plot: pages are drawn from that file, so other planes can be shown without recomputing (--plot-only):
           temporal_dz_zplanes_DS…_slab{S}_z{bins}.pdf       page 1 = cells per z bin (+ the chosen planes), then one page per plane
           temporal_dprime_zplanes_DS…_slab{S}_z{bins}.pdf

Needs the epoch arrays made by run_temporal_intensity_map.py with the CURRENT epochs and cell filter (stops with a
message for fish whose cache is missing or stale) and medoids_template_vox.npy.

z planes: coarse-z bin k covers template planes dz·k … dz·k+dz−1 (dz = 2 with the default DS). Choose them with
--z-bins 3,6,9,12, or --z-bins all (every bin that holds at least one cell; default: --n-planes 6 spread evenly over
the bins that hold cells). --slab N averages N neighbouring bins (default 1; with --z-bins all and --slab N the pages step by N). Each page title gives the bin and the template planes.

Usage:
    python chemogenetic/run/run_temporal_zplanes.py --config config_hcrt_trpv1_csn_120min
    python chemogenetic/run/run_temporal_zplanes.py --config … --plot-only --z-bins all
"""
import argparse
import gc
import importlib
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from matplotlib.backends.backend_pdf import PdfPages

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="Config module under chemogenetic/config/")
parser.add_argument("--z-bins", default=None, help="comma-separated coarse-z bins, e.g. 3,6,9,12, or 'all' = every bin with cells (default: evenly spaced)")
parser.add_argument("--n-planes", type=int, default=6, help="number of planes when --z-bins is not given (default 6)")
parser.add_argument("--slab", type=int, default=1, help="coarse-z bins averaged per page (default 1)")
parser.add_argument("--plot-only", action="store_true", help="draw from the saved npz; do not touch the per-fish files")
parser.add_argument("--n-jobs", type=int, default=8)
args = parser.parse_args()
cfg = importlib.import_module(f"chemogenetic.config.{args.config}")

from utils.data_io import fish_dir                                            # noqa: E402
from utils.run_guard import RunGuard, setup_or_exit                           # noqa: E402
from chemogenetic.cell_filter import load_cell_valid                          # noqa: E402
from chemogenetic.epochs import get_epochs                                    # noqa: E402
from chemogenetic.regressor import reference_traces                           # noqa: E402
from chemogenetic.temporal_windows import (                                   # noqa: E402
    N_WINDOWS, WINDOW_LABELS, epoch_cache_is_current, load_epoch_dz, load_epoch_dprime)
from chemogenetic.brainmap import _coarse_shape, _voxelize_cell_values, plot_group_comparison  # noqa: E402

EPOCHS = setup_or_exit(get_epochs, cfg)
GUARD  = RunGuard("temporal_zplanes")

dir_analysis = cfg.dir_analysis
CTRL_FISH, EXPT_FISH = cfg.ctrl_fish, cfg.expt_fish
CTRL_TAG, EXPT_TAG   = cfg.CTRL_TAG, cfg.EXPT_TAG
CLIP_ABS_DZ = cfg.CLIP_ABS_DZ

# keep identical to run_temporal_intensity_map.py so the pages are comparable with the projection figure
DS          = (5, 5, 2)
BRAIN_SHAPE = (288, 568, 40)
VMIN = {"dz": -3.0, "dprime": -0.6}
VMAX = {"dz":  3.0, "dprime":  0.6}

FIG_DIR = Path(dir_analysis) / "comparisons" / cfg.COMPARISON_TAG / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DS_TAG  = "x".join(str(v) for v in DS)
NPZ     = FIG_DIR / f"temporal_group_volumes_DS{DS_TAG}_clip{CLIP_ABS_DZ:g}.npz"
GROUPS  = {CTRL_TAG: CTRL_FISH, EXPT_TAG: EXPT_FISH}


# ── step 1: voxelize + average ────────────────────────────────────────────────

def _one_fish(fish, metric):
    """-> (list of N_WINDOWS coarse volumes or None, cells per coarse-z bin) ; (None, None) if the fish cannot be used."""
    f_dir = fish_dir(dir_analysis, fish)
    cpath = f_dir / "medoids_template_vox.npy"
    if not cpath.exists():
        return None, f"{fish[1]}: medoids_template_vox.npy is missing"
    try:
        _, sha = load_cell_valid(f_dir)
    except Exception as e:                                       # noqa: BLE001
        return None, f"{fish[1]}: {e}"
    if not epoch_cache_is_current(EPOCHS, f_dir, sha):
        return None, (f"{fish[1]}: epoch_dz / epoch_dprime missing or made with other epochs / cell filter "
                      f"(run submit_temporal_intensity_map.sh first)")
    coords = np.load(str(cpath))
    arr = load_epoch_dz(f_dir) if metric == "dz" else load_epoch_dprime(f_dir)
    if arr.shape[0] != coords.shape[0]:
        return None, f"{fish[1]}: {arr.shape[0]} epoch rows vs {coords.shape[0]} medoids"
    vols = [_voxelize_cell_values(coords, np.clip(arr[:, w], -CLIP_ABS_DZ, CLIP_ABS_DZ), BRAIN_SHAPE, DS)
            for w in range(N_WINDOWS)]
    # cells per coarse-z bin (finite value in window 1, finite coordinates, inside the brain)
    ok = np.isfinite(arr[:, 0]) & ~np.isnan(coords).any(axis=1)
    z = np.floor(coords[ok, 2] / DS[2]).astype(int)
    zc = _coarse_shape(BRAIN_SHAPE, DS)[2]
    cnt = np.bincount(np.clip(z[(z >= 0) & (z < BRAIN_SHAPE[2] // DS[2] + 1)], 0, zc - 1), minlength=zc)
    return vols, cnt


def _group(fish_list, metric, tag):
    print(f"  {tag} {metric} ({len(fish_list)} fish)")
    res = Parallel(n_jobs=min(args.n_jobs, len(fish_list)), backend="loky")(
        delayed(_one_fish)(f, metric) for f in fish_list)
    used, counts = [], 0
    for f, (vols, info) in zip(fish_list, res):
        if vols is None:
            GUARD.fail(f"{metric}/{tag}: {f[1]}", info + " — fish left out")
            continue
        used.append(vols); counts = counts + info
    out = np.full((N_WINDOWS,) + _coarse_shape(BRAIN_SHAPE, DS), np.nan, dtype=np.float32)
    for w in range(N_WINDOWS):
        v = [u[w] for u in used if u[w] is not None]
        if v:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out[w] = np.nanmean(np.stack(v), axis=0)
    return out, len(used), counts


def _reference():
    d = {"has_ref": np.array(False)}
    if getattr(cfg, "regressor_files", None):
        _, tonic, phasic = setup_or_exit(reference_traces, cfg, EPOCHS.baseline, 0)
        d.update(has_ref=np.array(True), tonic_mat=tonic, phasic_mat=phasic,
                 time_min=np.arange(tonic.shape[1]) / (60.0 * cfg.sampling_rate_hz))
    return d


def compute():
    print("Step 1: voxelize and average per group")
    data = {"window_vols": np.array([(s, e) for _, s, e in EPOCHS.windows]),
            "window_labels": np.array(WINDOW_LABELS), "epochs": np.array(EPOCHS.describe()),
            "ctrl_tag": np.array(CTRL_TAG), "expt_tag": np.array(EXPT_TAG),
            "ds": np.array(DS), "brain_shape": np.array(BRAIN_SHAPE), "clip": np.array(CLIP_ABS_DZ)}
    for metric in ("dz", "dprime"):
        for tag, key in ((CTRL_TAG, "ctrl"), (EXPT_TAG, "expt")):
            vol, n, cnt = _group(GROUPS[tag], metric, tag)
            data[f"{metric}_{key}"] = vol
            data[f"n_{metric}_{key}"] = np.array(n)
            data[f"zcount_{metric}_{key}"] = np.asarray(cnt)
            gc.collect()
    data.update({f"ref_{k}": v for k, v in _reference().items()})
    np.savez_compressed(str(NPZ), **data)
    print(f"  → saved: {NPZ.name}")


# ── step 2: plot ──────────────────────────────────────────────────────────────

def _ref_df(z, metric):
    if not bool(z["ref_has_ref"]):
        return None, None
    mat = z["ref_tonic_mat"] if metric == "dz" else z["ref_phasic_mat"]
    n = mat.shape[0]
    sem = mat.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros(mat.shape[1])
    df = pd.DataFrame({"time_min": z["ref_time_min"], "m": mat.mean(axis=0), "sem": sem})
    cols = []
    for i in range(n):
        df[f"f{i}"] = mat[i]; cols.append(f"f{i}")
    return df, cols


def _choose_bins(z):
    cnt = z["zcount_dz_ctrl"] + z["zcount_dz_expt"]
    if args.z_bins and args.z_bins.strip().lower() == "all":
        occ = np.where(cnt > 0)[0]
        bins = list(range(int(occ.min()), int(occ.max()) + 1, args.slab))
    elif args.z_bins:
        bins = [int(b) for b in args.z_bins.split(",")]
    else:
        occ = np.where(cnt >= 0.10 * cnt.max())[0]
        bins = sorted(set(np.linspace(occ.min(), occ.max(), args.n_planes).round().astype(int).tolist()))
    zc = len(cnt)
    bins = [b for b in bins if b + args.slab <= zc] if (args.z_bins or "").strip().lower() == "all" else bins
    bad = [b for b in bins if b < 0 or b + args.slab > zc]
    if bad:
        sys.exit(f"❌ z bins {bad} with --slab {args.slab} fall outside 0..{zc - 1}")
    return bins, cnt


def _summary_page(z, bins, cnt):
    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(cnt)); w = 0.4
    ax.bar(x - w / 2, z["zcount_dz_ctrl"], w, color="#4477aa", label=f"{z['ctrl_tag']} (N={int(z['n_dz_ctrl'])})")
    ax.bar(x + w / 2, z["zcount_dz_expt"], w, color="#cc6677", label=f"{z['expt_tag']} (N={int(z['n_dz_expt'])})")
    for b in bins:
        ax.axvspan(b - 0.5, b + args.slab - 0.5, color="k", alpha=0.12)
    ax.set_xlabel(f"coarse-z bin (bin k = template planes {DS[2]}k … {DS[2]}k+{DS[2] - 1}); grey = pages in this PDF")
    ax.set_ylabel("valid cells (all fish pooled)")
    ax.set_xticks(x); ax.legend(frameon=False)
    ax.set_title(f"Cells per z bin | DS={DS} | slab={args.slab} | pages: z bins {bins}")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def plot():
    z = np.load(str(NPZ), allow_pickle=False)
    bins, cnt = _choose_bins(z)
    labels = [str(s) for s in z["window_labels"]]
    wvols = [tuple(v) for v in z["window_vols"].tolist()]
    btag = "all" if (args.z_bins and args.z_bins.strip().lower() == "all") else "-".join(str(b) for b in bins)
    for metric, mlabel, title in (("dz", "Tonic ΔZ", "Tonic ΔZ"), ("dprime", "Phasic d′", "Phasic d′")):
        df, cols = _ref_df(z, metric)
        out = FIG_DIR / f"temporal_{metric}_zplanes_DS{DS_TAG}_slab{args.slab}_z{btag}.pdf"
        with PdfPages(str(out)) as pdf:
            fig = _summary_page(z, bins, cnt); pdf.savefig(fig); plt.close(fig)
            for b in bins:
                lo, hi = DS[2] * b, DS[2] * (b + args.slab) - 1
                ctrl = [z[f"{metric}_ctrl"][w] for w in range(N_WINDOWS)]
                expt = [z[f"{metric}_expt"][w] for w in range(N_WINDOWS)]
                diff = [e - c for e, c in zip(expt, ctrl)]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    fig = plot_group_comparison(
                        ctrl_maps=ctrl, expt_maps=expt, diff_maps=diff, window_labels=labels,
                        brain_shape=tuple(int(v) for v in z["brain_shape"]), DS=DS,
                        ctrl_tag=f"{z['ctrl_tag']}\n(N={int(z[f'n_{metric}_ctrl'])})",
                        expt_tag=f"{z['expt_tag']}\n(N={int(z[f'n_{metric}_expt'])})",
                        diff_tag="EXPT − CTRL", metric_label=mlabel,
                        title=f"{title} | z bin {b}" + (f"–{b + args.slab - 1}" if args.slab > 1 else "")
                              + f" (template planes {lo}–{hi}) | DS={DS} | 15min window 10min step | clip={float(z['clip']):g}",
                        vmin=VMIN[metric], vmax=VMAX[metric], z_slice=b, z_slab=args.slab,
                        hcrt_df=df, hcrt_trace_col="m", hcrt_sem_col="sem", hcrt_fish_cols=cols,
                        hcrt_trace_label="Hcrt F_tonic" if metric == "dz" else "Hcrt F_phasic",
                        hcrt_window_vols=wvols, fig_dir=None)
                fig.set_size_inches(fig.get_size_inches())
                pdf.savefig(fig)
                print(f"    page: z bin {b} (template planes {lo}–{hi})")
        print(f"  → saved: {out.name}")


def main():
    print(f"Config: {args.config}\nepochs: {EPOCHS.describe()}\nOutput → {FIG_DIR}\n")
    if not args.plot_only:
        compute()
    elif not NPZ.exists():
        sys.exit(f"❌ {NPZ.name} not found — run without --plot-only first")
    print("\nStep 2: plot")
    plot()
    GUARD.finish()


if __name__ == "__main__":
    main()
