"""
tonic_dz.py   (stage F1)
=======================
Tonic responders defined by ΔZ of F_tonic. The GLM ΔR² is no longer part of the definition
(the GLM is the optional stage F2).

Per cell
--------
    ΔZ_w   = (mean F_tonic in window w − baseline mean) / baseline σ        w = 1..7 (the temporal windows)
    max ΔZ = ΔZ at the window with the largest |ΔZ|   (sign taken from that window; "argmax window")

    baseline = min 20-45 (epoch rule), windows = the seven 15-min windows used by the temporal maps
    (they start at drug onset, step 10 min, and the last three reach into the washout).
    Cells excluded by the cell filter (cell_valid.npy) are NaN and never selected.

Responders (rank-based, within each fish)
-----------------------------------------
    Cells are split into N_STRATA strata by baseline σ (quantiles of σ among valid cells, so high-variance
    cells do not dominate). Inside each stratum the top p % of |max ΔZ| are responders; positive ones are
    "pos", negative ones "neg" (pooled ranking, so the pos:neg ratio is a result, not a constraint).
    p = TONIC_DZ_TOP_PERCENT (default 5); TONIC_DZ_TOP_PERCENT_SENS (default 2.5, 10) are saved alongside.

    There is NO null and NO calibration against ctrl: the fraction of responders is fixed by p, so it
    cannot be compared between groups. The set answers "which cells changed most", not "how many".
    Compare groups with the per-fish distribution of |max ΔZ| (run_figures_tonic_dz.py) and with the
    composition of the top sets, using fish as the unit.

Files per fish (fish_dir(dir_analysis, fish)):
    tonic_dz_windows.npy          (n_cells, n_windows) float32  ΔZ per window (NaN = excluded cell)
    tonic_dz_max.npy              (n_cells,) float32            signed ΔZ at the argmax window
    tonic_dz_argmax_window.npy    (n_cells,) int8               0-based window of max |ΔZ| (-1 = excluded)
    tonic_dz_sigma_baseline.npy   (n_cells,) float32            baseline σ of F_tonic
    tonic_dz_stratum.npy          (n_cells,) int8               σ stratum (-1 = excluded)
    tonic_pos_dz_top{ptag}_idxs.npy / tonic_neg_dz_top{ptag}_idxs.npy   int64 cell indices
    tonic_dz_meta.json            inputs (epochs, cell_valid fingerprint, parameters), per-stratum cutoffs
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from chemogenetic.cell_filter import load_cell_valid
from chemogenetic.temporal_windows import compute_epoch_dz
from utils.data_io import fish_dir
from utils.run_guard import PipelineError

DEFAULT_TOP_PERCENT = 5.0
DEFAULT_TOP_PERCENT_SENS = (2.5, 10.0)
DEFAULT_N_STRATA = 5

META_NAME = "tonic_dz_meta.json"
_ARRAYS = {
    "windows": "tonic_dz_windows.npy",
    "max": "tonic_dz_max.npy",
    "argmax": "tonic_dz_argmax_window.npy",
    "sigma": "tonic_dz_sigma_baseline.npy",
    "stratum": "tonic_dz_stratum.npy",
}


def ptag(p: float) -> str:
    """5.0 -> '5', 2.5 -> '2p5', 10 -> '10' (file-name tag)."""
    return f"{float(p):g}".replace(".", "p")


def responder_paths(f_dir, p: float):
    f_dir = Path(f_dir)
    t = ptag(p)
    return (f_dir / f"tonic_pos_dz_top{t}_idxs.npy", f_dir / f"tonic_neg_dz_top{t}_idxs.npy")


def tonic_dz_params(cfg) -> dict:
    """Settings from the config, falling back to the defaults for keys an older config lacks."""
    main = float(getattr(cfg, "TONIC_DZ_TOP_PERCENT", DEFAULT_TOP_PERCENT))
    sens = [float(x) for x in getattr(cfg, "TONIC_DZ_TOP_PERCENT_SENS", DEFAULT_TOP_PERCENT_SENS)]
    percents = [main] + [x for x in sens if x != main]
    n_strata = int(getattr(cfg, "TONIC_DZ_N_STRATA", DEFAULT_N_STRATA))
    bad = []
    if any(not (0 < x < 50) for x in percents):
        bad.append(f"TONIC_DZ_TOP_PERCENT / _SENS must be between 0 and 50 (got {percents})")
    if n_strata < 1:
        bad.append(f"TONIC_DZ_N_STRATA must be >= 1 (got {n_strata})")
    if bad:
        raise PipelineError("Invalid tonic-ΔZ setting(s) in the config:\n  " + "\n  ".join(bad))
    return {"top_percent": main, "percents": percents, "n_strata": n_strata}


# ──────────────────────────────────────────────────────────────────────────────
# computation
# ──────────────────────────────────────────────────────────────────────────────

def max_abs_window(dz: np.ndarray):
    """
    dz (n_cells, n_windows) with NaN for unusable entries.
    Returns (signed ΔZ at the argmax window, argmax window int8, usable bool).
    A cell is usable when at least one window is finite.
    """
    fin = np.isfinite(dz)
    usable = fin.any(axis=1)
    absdz = np.where(fin, np.abs(dz), -1.0)
    arg = absdz.argmax(axis=1)
    signed = dz[np.arange(dz.shape[0]), arg].astype(np.float32)
    signed[~usable] = np.nan
    arg = np.where(usable, arg, -1).astype(np.int8)
    return signed, arg, usable


def assign_strata(sigma: np.ndarray, usable: np.ndarray, n_strata: int):
    """Quantile strata of baseline σ among usable cells (0..n_strata-1; -1 = not usable)."""
    stratum = np.full(sigma.shape[0], -1, dtype=np.int8)
    if n_strata <= 1 or usable.sum() == 0:
        stratum[usable] = 0
        return stratum, np.array([])
    edges = np.quantile(sigma[usable], np.linspace(0, 1, n_strata + 1)[1:-1])
    stratum[usable] = np.digitize(sigma[usable], edges).astype(np.int8)
    return stratum, edges


def select_top_percent(dz_max: np.ndarray, stratum: np.ndarray, p: float, n_strata: int):
    """Top p % of |ΔZ| inside each stratum. Returns (pos_idx, neg_idx, cutoffs per stratum)."""
    absmax = np.abs(dz_max)
    pos, neg, cutoffs = [], [], []
    for k in range(max(n_strata, 1)):
        m = np.flatnonzero(stratum == k)
        if m.size == 0:
            cutoffs.append(None)
            continue
        cut = float(np.percentile(absmax[m], 100.0 - p))
        cutoffs.append(cut)
        sel = m[absmax[m] >= cut]
        pos.append(sel[dz_max[sel] > 0])
        neg.append(sel[dz_max[sel] < 0])
    cat = lambda xs: (np.sort(np.concatenate(xs)).astype(np.int64) if xs else np.array([], dtype=np.int64))
    return cat(pos), cat(neg), cutoffs


def compute_tonic_dz(f_tonic, epochs, valid: np.ndarray, params: dict, chunk_cells: int = 20000):
    """
    f_tonic may be a memory-mapped (n_cells, T) array; it is processed in chunks of cells
    (the per-cell statistics do not depend on other cells, so the result is chunk-independent).
    """
    n_cells = f_tonic.shape[0]
    n_w = epochs.n_windows
    dz = np.empty((n_cells, n_w), dtype=np.float32)
    sigma = np.empty(n_cells, dtype=np.float32)
    for s in range(0, n_cells, chunk_cells):
        e = min(s + chunk_cells, n_cells)
        d, _, sd = compute_epoch_dz(np.asarray(f_tonic[s:e], dtype=np.float32), epochs, valid[s:e])
        dz[s:e], sigma[s:e] = d, sd

    dz_max, argmax, usable = max_abs_window(dz)
    usable &= valid
    dz_max[~usable] = np.nan
    argmax[~usable] = -1
    sigma_u = np.where(usable, sigma, np.nan).astype(np.float32)
    stratum, edges = assign_strata(sigma, usable, params["n_strata"])

    sets, cutoffs = {}, {}
    for p in params["percents"]:
        pos, neg, cut = select_top_percent(dz_max, stratum, p, params["n_strata"])
        sets[p], cutoffs[p] = (pos, neg), cut
    return {
        "windows": dz, "max": dz_max, "argmax": argmax, "sigma": sigma_u, "stratum": stratum,
        "sets": sets, "cutoffs": cutoffs, "sigma_edges": edges, "n_usable": int(usable.sum()),
    }


# ──────────────────────────────────────────────────────────────────────────────
# per-fish driver
# ──────────────────────────────────────────────────────────────────────────────

def _inputs(epochs, valid_sha: str, params: dict) -> dict:
    return {
        "baseline": [int(x) for x in epochs.baseline],
        "windows": [[lab, int(s), int(e)] for lab, s, e in epochs.windows],
        "cell_valid_sha": valid_sha,
        "percents": [float(x) for x in params["percents"]],
        "n_strata": int(params["n_strata"]),
        "statistic": "max|dZ| over windows, sign from argmax window; pooled |dZ| top-p% within baseline-sigma strata",
    }


def tonic_dz_is_current(f_dir, epochs, valid_sha: str, params: dict) -> bool:
    """True when every output exists and was made from exactly these inputs (read-only; used by run_status)."""
    f_dir = Path(f_dir)
    mp = f_dir / META_NAME
    if not mp.exists():
        return False
    for name in _ARRAYS.values():
        if not (f_dir / name).exists():
            return False
    for p in params["percents"]:
        if not all(x.exists() for x in responder_paths(f_dir, p)):
            return False
    try:
        return json.loads(mp.read_text()).get("inputs") == _inputs(epochs, valid_sha, params)
    except Exception:           # noqa: BLE001
        return False


def tonic_dz_one_fish(fish, dir_analysis, epochs, params: dict, overwrite: bool = False) -> dict:
    proj_ID, expt_ID = fish
    f_dir = fish_dir(dir_analysis, fish)
    ft_path = f_dir / "f_tonic.npy"
    if not ft_path.exists():
        raise PipelineError(f"{expt_ID}: f_tonic.npy missing ({ft_path}) — run run_decompose.py first")

    ft = np.load(str(ft_path), mmap_mode="r")
    epochs.check_recording(ft.shape[1], "f_tonic")
    valid, valid_sha = load_cell_valid(f_dir, expected_n_cells=ft.shape[0])

    if not overwrite and tonic_dz_is_current(f_dir, epochs, valid_sha, params):
        print(f"⏩ {expt_ID}: tonic ΔZ outputs exist and match the current epochs / cell filter / settings")
        return {"fish": fish, "status": "skipped"}

    out = compute_tonic_dz(ft, epochs, valid, params)
    if out["n_usable"] == 0:
        raise PipelineError(f"{expt_ID}: no usable cells (cell filter or non-finite data removed all of them)")

    np.save(str(f_dir / _ARRAYS["windows"]), out["windows"])
    np.save(str(f_dir / _ARRAYS["max"]), out["max"])
    np.save(str(f_dir / _ARRAYS["argmax"]), out["argmax"])
    np.save(str(f_dir / _ARRAYS["sigma"]), out["sigma"])
    np.save(str(f_dir / _ARRAYS["stratum"]), out["stratum"])
    counts = {}
    for p, (pos, neg) in out["sets"].items():
        pp, pn = responder_paths(f_dir, p)
        np.save(str(pp), pos)
        np.save(str(pn), neg)
        counts[ptag(p)] = {"n_pos": int(pos.size), "n_neg": int(neg.size)}
    meta = {
        "inputs": _inputs(epochs, valid_sha, params),
        "n_cells": int(ft.shape[0]),
        "n_usable": out["n_usable"],
        "sigma_stratum_edges": [float(x) for x in out["sigma_edges"]],
        "cutoff_abs_dz_per_stratum": {ptag(p): c for p, c in out["cutoffs"].items()},
        "counts": counts,
    }
    (f_dir / META_NAME).write_text(json.dumps(meta, indent=2))   # written last: marks the outputs as complete

    main = ptag(params["top_percent"])
    c = counts[main]
    am = out["argmax"][out["argmax"] >= 0]
    hist = np.bincount(am, minlength=epochs.n_windows)
    print(f"✅ {expt_ID}: tonic ΔZ | usable {out['n_usable']:,}/{ft.shape[0]:,} | top {params['top_percent']:g}%: "
          f"{c['n_pos']:,} pos + {c['n_neg']:,} neg | argmax window counts "
          f"{[f'{100 * h / max(hist.sum(), 1):.0f}%' for h in hist]}")
    return {"fish": fish, "status": "ok", "n_pos": c["n_pos"], "n_neg": c["n_neg"],
            "argmax_hist": hist.tolist()}
