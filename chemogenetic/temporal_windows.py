"""
temporal_windows.py
===================
Compute per-cell metrics for each temporal window.

Tonic ΔZ:
    ΔZ_i(w) = (μ_window_w - μ_baseline) / σ_baseline        on F_tonic

Phasic d' (same formula on F_phasic; NOT a pooled-SD d′):
    d'_i(w) = (μ_phasic_window_w - μ_phasic_baseline) / σ_phasic_baseline

Baseline and windows come from chemogenetic/epochs.py (config-driven):
    baseline = [baseline_start + EPOCH_DISCARD_MIN, baseline_end)   (vols 1200-2700 for the 120-min CSN design)
    windows  = seven 15-min windows, 10-min step, starting at drug onset:
        W0 0-15, W1 10-25, W2 20-35, W3 30-45, W4 40-55, W5 50-65, W6 60-75 min

I/O conventions (must match run_decompose.py):
    Input  reads from:  fish_dir(dir_analysis, fish) / "f_tonic.npy"
                        fish_dir(dir_analysis, fish) / "f_phasic.npy"
    Output writes to:   fish_dir(dir_analysis, fish) / "epoch_dz.npy"
                        fish_dir(dir_analysis, fish) / "epoch_dprime.npy"

Location:
    ~/zwba/chemogenetic/temporal_windows.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from chemogenetic.cell_filter import mask_rows
from chemogenetic.epochs import (   # noqa: F401  (re-exported for the run scripts)
    Epochs,
    N_WINDOWS,
    WINDOW_LABELS,
)


def _baseline_stats(arr: np.ndarray, epochs: Epochs):
    b0, b1 = epochs.baseline
    epochs.check_recording(arr.shape[1], "f_tonic/f_phasic")
    seg = arr[:, b0:b1]
    mu, sd = seg.mean(axis=1), seg.std(axis=1)
    sd_safe = np.where(sd < 1e-6, 1e-6, sd)     # guard against silent / flat cells
    return mu, sd, sd_safe


def _windowed_effect(arr: np.ndarray, epochs: Epochs):
    mu_bl, sd_bl, sd_safe = _baseline_stats(arr, epochs)
    out = np.empty((arr.shape[0], epochs.n_windows), dtype=np.float32)
    for w, (_, start, end) in enumerate(epochs.windows):
        out[:, w] = (arr[:, start:end].mean(axis=1) - mu_bl) / sd_safe
    return out, mu_bl, sd_bl


# ── tonic ΔZ ──────────────────────────────────────────────────────────────

def compute_epoch_dz(
    f_tonic: np.ndarray, epochs: Epochs, valid: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Per-cell tonic ΔZ for each temporal window.

    valid : optional bool mask (chemogenetic/cell_filter.py). Cells with valid == False get
            NaN in `dz` (rows are kept, so indices still match the medoid coordinates).

    Returns
    -------
    dz : (n_cells, n_windows) float32
    mu_baseline, sigma_baseline : (n_cells,) float32   (for QC)
    """
    dz, mu_bl, sd_bl = _windowed_effect(f_tonic, epochs)
    if valid is not None:
        dz = mask_rows(dz, valid)
    return dz, mu_bl.astype(np.float32), sd_bl.astype(np.float32)


# ── phasic d' per window ───────────────────────────────────────────────────

def compute_epoch_dprime(f_phasic: np.ndarray, epochs: Epochs,
                         valid: np.ndarray | None = None) -> np.ndarray:
    """Per-cell phasic effect size for each temporal window, (n_cells, n_windows) float32.
    Cells with valid == False (cell_filter) get NaN."""
    dprime, _, _ = _windowed_effect(f_phasic, epochs)
    if valid is not None:
        dprime = mask_rows(dprime, valid)
    return dprime.astype(np.float32)


# ── I/O helpers ────────────────────────────────────────────────────────────
# All paths follow run_decompose.py conventions:
#   fish_dir(dir_analysis, fish)  →  the per-fish output directory
#   filenames are lowercase (f_tonic.npy, not F_tonic.npy)

def load_ftonic(fish_analysis_dir: Path) -> np.ndarray:
    """Load f_tonic.npy from the fish analysis directory."""
    p = Path(fish_analysis_dir) / "f_tonic.npy"
    if not p.exists():
        raise FileNotFoundError(f"f_tonic.npy not found: {p}\n"
                                "Run run_decompose.py first.")
    return np.load(str(p))


def load_fphasic(fish_analysis_dir: Path) -> np.ndarray:
    """Load f_phasic.npy from the fish analysis directory."""
    p = Path(fish_analysis_dir) / "f_phasic.npy"
    if not p.exists():
        raise FileNotFoundError(f"f_phasic.npy not found: {p}\n"
                                "Run run_decompose.py first.")
    return np.load(str(p))


def save_epoch_dz(dz: np.ndarray, fish_analysis_dir: Path) -> Path:
    """Save dz (n_cells × n_windows) → epoch_dz.npy"""
    out = Path(fish_analysis_dir) / "epoch_dz.npy"
    np.save(str(out), dz)
    print(f"  Saved epoch_dz:     {out}  shape={dz.shape}")
    return out


def save_epoch_dprime(dprime: np.ndarray, fish_analysis_dir: Path) -> Path:
    """Save dprime (n_cells × n_windows) → epoch_dprime.npy"""
    out = Path(fish_analysis_dir) / "epoch_dprime.npy"
    np.save(str(out), dprime)
    print(f"  Saved epoch_dprime: {out}  shape={dprime.shape}")
    return out


def load_epoch_dz(fish_analysis_dir: Path) -> np.ndarray:
    """Load epoch_dz.npy, shape (n_cells, n_windows)."""
    p = Path(fish_analysis_dir) / "epoch_dz.npy"
    if not p.exists():
        raise FileNotFoundError(f"epoch_dz.npy not found: {p}\n"
                                "Run run_temporal_intensity_map.py first.")
    return np.load(str(p))


# ── cache provenance: epoch arrays are only valid for the epochs that made them ──

def epochs_meta(epochs: Epochs, cell_valid_sha: str | None = None) -> dict:
    return {
        "baseline": [int(x) for x in epochs.baseline],
        "windows": [[lab, int(s), int(e)] for lab, s, e in epochs.windows],
        "discard_frames": int(epochs.discard_frames),
        "cell_valid_sha": cell_valid_sha,     # fingerprint of cell_valid.npy used for the NaN rows
    }


def save_epoch_meta(epochs: Epochs, fish_analysis_dir: Path,
                    cell_valid_sha: str | None = None) -> Path:
    import json
    out = Path(fish_analysis_dir) / "epoch_meta.json"
    out.write_text(json.dumps(epochs_meta(epochs, cell_valid_sha), indent=2))
    return out


def epoch_cache_is_current(epochs: Epochs, fish_analysis_dir: Path,
                           cell_valid_sha: str | None = None) -> bool:
    """True only if epoch_dz/epoch_dprime exist AND were made with exactly these epochs and
    exactly this cell_valid mask (arrays made before the cell filter existed are recomputed once)."""
    import json
    d = Path(fish_analysis_dir)
    meta = d / "epoch_meta.json"
    if not (d / "epoch_dz.npy").exists() or not (d / "epoch_dprime.npy").exists():
        return False
    if not meta.exists():
        return False        # legacy arrays of unknown provenance -> recompute once
    try:
        return json.loads(meta.read_text()) == epochs_meta(epochs, cell_valid_sha)
    except Exception:       # noqa: BLE001 — unreadable sidecar == not current
        return False


def load_epoch_dprime(fish_analysis_dir: Path) -> np.ndarray:
    """Load epoch_dprime.npy, shape (n_cells, n_windows)."""
    p = Path(fish_analysis_dir) / "epoch_dprime.npy"
    if not p.exists():
        raise FileNotFoundError(f"epoch_dprime.npy not found: {p}\n"
                                "Run run_temporal_intensity_map.py first.")
    return np.load(str(p))
