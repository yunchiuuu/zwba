"""
cell_filter.py
==============
Per-fish valid-cell mask.  F_tonic / F_phasic are NEVER modified: the mask is a separate file
(`cell_valid.npy`, bool, True = keep) that every downstream stage reads and applies by setting
the excluded cells' results to NaN (rows are kept so cell indices still match the medoid
coordinates).

Why a mask and not a threshold on the result (|ΔZ|, d′)?
    The exclusion has to be independent of the drug effect. A cut on the readout also removes the
    strongest genuine responders, and it removes more of them in the group that responds more.
    So the criteria look at the BASELINE epoch only (min 20-45 by default):

        nonfinite   NaN / inf anywhere in the baseline epoch
        zero_frac   fraction of baseline samples with F_tonic <= 0 is above CELL_MAX_ZERO_FRAC
                    (dim cells that the F_dark subtraction clipped to 0)
        flat        baseline SD of F_tonic is below CELL_MIN_BASELINE_SD
                    (ΔZ divides by this SD: a flat baseline gives ΔZ ~ 1e7)
        dim         baseline mean of F_tonic is below CELL_MIN_BASELINE_MEAN counts

    plus one backstop on the result, applied only to cells that passed the four checks above:

        backstop    max over the temporal windows of |tonic ΔZ| is above CELL_MAX_ABS_DZ
                    (set CELL_MAX_ABS_DZ = None to switch it off). Far above the real tail
                    (genuine cells reach |ΔZ| of a few tens), so it only catches degenerate
                    cells that slipped through. The cell ids are listed in cell_valid_meta.json.

All thresholds are config parameters (defaults below, so older configs keep working).

Files (in fish_dir(dir_analysis, fish)):
    cell_valid.npy          bool (n_cells,)  True = keep
    cell_valid_meta.json    parameters, per-criterion counts, example cell ids, the mask's
                            fingerprint (`sha`) and the inputs it was computed from.

The fingerprint (`sha`) is stamped into every downstream cache (epoch_meta.json, d′ and GLM
sidecars) so that results made with a different mask are detected instead of reused.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

from utils.run_guard import PipelineError

MASK_NAME = "cell_valid.npy"
META_NAME = "cell_valid_meta.json"

DEFAULTS = {
    "CELL_MAX_ZERO_FRAC":     0.10,    # baseline samples with F_tonic <= 0
    "CELL_MIN_BASELINE_SD":   1e-3,    # counts
    "CELL_MIN_BASELINE_MEAN": 1.0,     # counts
    "CELL_MAX_ABS_DZ":        200.0,   # backstop; None disables
}
_SD_FLOOR_DZ = 1e-6                    # same floor temporal_windows uses for ΔZ


def filter_params(cfg) -> dict:
    """Thresholds from the config, falling back to DEFAULTS for keys the config lacks."""
    return {k: getattr(cfg, k, v) for k, v in DEFAULTS.items()}


def _check_params(p: dict) -> None:
    bad = []
    if not (p["CELL_MAX_ZERO_FRAC"] is not None and 0 <= float(p["CELL_MAX_ZERO_FRAC"]) <= 1):
        bad.append(f"CELL_MAX_ZERO_FRAC={p['CELL_MAX_ZERO_FRAC']!r} (must be between 0 and 1)")
    for k in ("CELL_MIN_BASELINE_SD", "CELL_MIN_BASELINE_MEAN"):
        if p[k] is None or float(p[k]) < 0:
            bad.append(f"{k}={p[k]!r} (must be >= 0)")
    if p["CELL_MAX_ABS_DZ"] is not None and float(p["CELL_MAX_ABS_DZ"]) <= 0:
        bad.append(f"CELL_MAX_ABS_DZ={p['CELL_MAX_ABS_DZ']!r} (must be > 0, or None to disable)")
    if bad:
        raise PipelineError("Invalid cell-filter setting(s) in the config:\n  " + "\n  ".join(bad))


def mask_sha(valid: np.ndarray) -> str:
    """Short fingerprint of a boolean mask (identical masks give identical strings)."""
    v = np.asarray(valid, dtype=bool).ravel()
    h = hashlib.sha256()
    h.update(str(v.size).encode())
    h.update(np.packbits(v).tobytes())
    return h.hexdigest()[:16]


def mask_rows(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Copy of `arr` (cells on axis 0) with the rows of excluded cells set to NaN."""
    valid = np.asarray(valid, dtype=bool).ravel()
    if arr.shape[0] != valid.size:
        raise PipelineError(f"cell_valid has {valid.size} cells but the array has {arr.shape[0]} rows "
                            f"— the mask belongs to a different decomposition")
    out = np.array(arr, dtype=np.float32, copy=True)
    out[~valid] = np.nan
    return out


# ──────────────────────────────────────────────────────────────────────────────
# computation
# ──────────────────────────────────────────────────────────────────────────────

def compute_cell_valid(f_tonic, epochs, params: dict, chunk_cells: int = 20000):
    """
    Returns (valid, info).  `f_tonic` may be a memory-mapped (n_cells, T) array; it is read in
    chunks of cells and only the columns the criteria need (baseline epoch .. last window).
    """
    _check_params(params)
    n_cells, T = f_tonic.shape
    epochs.check_recording(T, "f_tonic")
    b0, b1 = (int(x) for x in epochs.baseline)
    wins = [(lab, int(s), int(e)) for lab, s, e in epochs.windows]
    hi = max([b1] + [e for _, _, e in wins])
    lo = min([b0] + [s for _, s, _ in wins])
    max_dz = params["CELL_MAX_ABS_DZ"]

    c_nonfinite = np.zeros(n_cells, bool)
    c_zero      = np.zeros(n_cells, bool)
    c_flat      = np.zeros(n_cells, bool)
    c_dim       = np.zeros(n_cells, bool)
    c_backstop  = np.zeros(n_cells, bool)

    for s in range(0, n_cells, chunk_cells):
        e = min(s + chunk_cells, n_cells)
        X = np.asarray(f_tonic[s:e, lo:hi], dtype=np.float32)
        base = X[:, b0 - lo:b1 - lo]

        nonfin = ~np.isfinite(base).all(axis=1)
        base_c = np.where(np.isfinite(base), base, 0.0)           # keeps the maths warning-free
        zero   = (base_c <= 0).mean(axis=1) > float(params["CELL_MAX_ZERO_FRAC"])
        mu     = base_c.mean(axis=1)
        sd     = base_c.std(axis=1)
        flat   = sd < float(params["CELL_MIN_BASELINE_SD"])
        dim    = mu < float(params["CELL_MIN_BASELINE_MEAN"])

        c_nonfinite[s:e], c_zero[s:e], c_flat[s:e], c_dim[s:e] = nonfin, zero, flat, dim

        if max_dz is not None:
            ok = ~(nonfin | zero | flat | dim)
            if ok.any():
                sd_safe = np.where(sd < _SD_FLOOR_DZ, _SD_FLOOR_DZ, sd)
                big = np.zeros(e - s, bool)
                for _, ws, we in wins:
                    wm = X[:, ws - lo:we - lo].mean(axis=1)
                    with np.errstate(invalid="ignore"):
                        dz = np.abs((wm - mu) / sd_safe)
                    big |= np.isfinite(dz) & (dz > float(max_dz))
                c_backstop[s:e] = ok & big

    excluded = c_nonfinite | c_zero | c_flat | c_dim | c_backstop
    valid = ~excluded

    def _ids(m, k=20):
        return [int(i) for i in np.flatnonzero(m)[:k]]

    info = {
        "n_cells": int(n_cells),
        "n_valid": int(valid.sum()),
        "n_excluded": int(excluded.sum()),
        "counts": {                                   # a cell can fail several checks
            "nonfinite": int(c_nonfinite.sum()),
            "zero_frac": int(c_zero.sum()),
            "flat":      int(c_flat.sum()),
            "dim":       int(c_dim.sum()),
            "backstop_only": int(c_backstop.sum()),   # passed the baseline checks, failed |ΔZ|
        },
        "example_ids": {
            "nonfinite": _ids(c_nonfinite), "zero_frac": _ids(c_zero),
            "flat": _ids(c_flat), "dim": _ids(c_dim),
            "backstop_only": _ids(c_backstop, 50),
        },
    }
    return valid, info


def describe(info: dict) -> str:
    c = info["counts"]
    n, x = info["n_cells"], info["n_excluded"]
    s = (f"excluded {x:,}/{n:,} cells ({100 * x / max(n, 1):.3f}%) | baseline checks: "
         f"nonfinite {c['nonfinite']:,}, zero-frac {c['zero_frac']:,}, flat {c['flat']:,}, "
         f"dim {c['dim']:,} (a cell can fail several) | |ΔZ| backstop only: {c['backstop_only']:,}")
    if c["backstop_only"]:
        s += f"  ids: {info['example_ids']['backstop_only'][:10]}"
    return s


# ──────────────────────────────────────────────────────────────────────────────
# files
# ──────────────────────────────────────────────────────────────────────────────

def _inputs(f_dir: Path, epochs, params: dict) -> dict:
    st = (Path(f_dir) / "f_tonic.npy").stat()
    return {
        "params": {k: (None if v is None else float(v)) for k, v in params.items()},
        "baseline": [int(x) for x in epochs.baseline],
        "windows": [[lab, int(s), int(e)] for lab, s, e in epochs.windows],
        "f_tonic_size": int(st.st_size),
        "f_tonic_mtime_ns": int(st.st_mtime_ns),
    }


def ensure_cell_valid(f_dir, epochs, cfg, overwrite: bool = False, verbose: bool = True):
    """
    Return (valid, sha, info) for one fish, computing and saving the mask only when it is missing
    or was made from different inputs (other thresholds, other epochs, or a re-made f_tonic.npy).
    """
    f_dir = Path(f_dir)
    ft_path = f_dir / "f_tonic.npy"
    if not ft_path.exists():
        raise PipelineError(f"Cannot build the cell filter: {ft_path} is missing — run run_decompose.py first.")

    params = filter_params(cfg)
    _check_params(params)
    inputs = _inputs(f_dir, epochs, params)
    mpath, jpath = f_dir / MASK_NAME, f_dir / META_NAME

    if not overwrite and mpath.exists() and jpath.exists():
        try:
            meta = json.loads(jpath.read_text())
            if meta.get("inputs") == inputs:
                valid = np.load(str(mpath)).astype(bool)
                if valid.size == meta["info"]["n_cells"] and mask_sha(valid) == meta["sha"]:
                    if verbose:
                        print(f"  ⏩ {f_dir.name}: cell filter up to date | {describe(meta['info'])}")
                    return valid, meta["sha"], meta["info"]
        except Exception:                      # noqa: BLE001 - unreadable cache == recompute
            pass
        if verbose:
            print(f"  ♻️  {f_dir.name}: cell filter inputs changed (thresholds, epochs or f_tonic) — recomputing")

    ft = np.load(str(ft_path), mmap_mode="r")
    valid, info = compute_cell_valid(ft, epochs, params)
    sha = mask_sha(valid)
    # write via temp files + atomic rename: two stages (e.g. d′ and GLM jobs) may build the same
    # mask at the same time, and a reader must never see a half-written file
    tmp_m = f_dir / f"{MASK_NAME[:-4]}.tmp{os.getpid()}.npy"
    tmp_j = f_dir / f"{META_NAME}.tmp{os.getpid()}"
    np.save(str(tmp_m), valid)
    tmp_j.write_text(json.dumps({"sha": sha, "inputs": inputs, "info": info}, indent=2))
    os.replace(tmp_m, mpath)
    os.replace(tmp_j, jpath)           # meta last: it is what marks the mask as complete
    if verbose:
        print(f"  ✅ {f_dir.name}: cell filter | {describe(info)}\n     saved → {mpath}")
    return valid, sha, info


def mask_status(f_dir, epochs, cfg):
    """
    Read-only check used by run_status.py: ("ok" | "stale" | "missing", sha_or_None, info_or_None).
    Never computes or writes anything.
    """
    f_dir = Path(f_dir)
    jpath = f_dir / META_NAME
    if not (f_dir / MASK_NAME).exists() or not jpath.exists():
        return "missing", None, None
    try:
        meta = json.loads(jpath.read_text())
        if (f_dir / "f_tonic.npy").exists():
            cur = _inputs(f_dir, epochs, filter_params(cfg))
            if meta.get("inputs") != cur:
                return "stale", meta.get("sha"), meta.get("info")
        return "ok", meta.get("sha"), meta.get("info")
    except Exception:                          # noqa: BLE001 - unreadable == treat as stale
        return "stale", None, None


def load_cell_valid(f_dir, expected_n_cells: int | None = None):
    """
    Load an existing mask: (valid, sha).  Raises PipelineError when it does not exist, because every
    downstream stage must use the same mask (run_decompose.py creates it, also for fish whose
    F_tonic / F_phasic already exist).
    """
    f_dir = Path(f_dir)
    mpath = f_dir / MASK_NAME
    if not mpath.exists():
        raise PipelineError(
            f"{f_dir.name}: {MASK_NAME} is missing in {f_dir}\n"
            f"  Run run_decompose.py (it only builds the missing cell filter when F_tonic / F_phasic "
            f"already exist).")
    valid = np.load(str(mpath)).astype(bool)
    if expected_n_cells is not None and valid.size != int(expected_n_cells):
        raise PipelineError(
            f"{f_dir.name}: {MASK_NAME} has {valid.size} cells but the data has {expected_n_cells} "
            f"— re-run run_decompose.py to rebuild it.")
    return valid, mask_sha(valid)
