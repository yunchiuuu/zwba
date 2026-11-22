"""
regressor.py
============
User-defined tonic-GLM regressor, read from the files named in the experiment config.
There is NO fallback: if the regressor is not defined or cannot be read, the run stops
with a message saying what to fix (the old silent CSTR-capsaicin fallback is gone).

Config settings (see config_template.py)
----------------------------------------
    regressor_files           list of CSV paths, one per reference recording
                              (e.g. the hcrt_all.csv of each HCRT-imaged fish)
    regressor_column          CSV column holding the raw trace            (default "Mean")
    regressor_preprocess      "tonic_zscore" — per file: centered rolling percentile
                                (= F_tonic), z-score on [0, drug_start), then mean
                                across files. Same as hcrt_regressor.ipynb.
                              "none" — the column is already the final regressor
                                (e.g. a measured concentration trace); mean across files
    regressor_tonic_window    rolling window, volumes                     (default 600)
    regressor_tonic_quantile  rolling quantile, 0-1                       (default 0.20)

Every trace must be sampled once per volume (same time base as the imaging) and be at
least `wash_end` volumes long.

Location:
    ~/zwba/chemogenetic/regressor.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from utils.run_guard import PipelineError, require_cfg

_DEFAULTS = dict(regressor_column="Mean", regressor_preprocess="tonic_zscore",
                 regressor_tonic_window=600, regressor_tonic_quantile=0.20)


def _setting(cfg, name):
    return getattr(cfg, name, _DEFAULTS[name])


def rolling_tonic(raw: np.ndarray, window: int, quantile: float) -> np.ndarray:
    """Centered rolling quantile (pandas; min_periods=1) — the notebook's F_tonic."""
    return (pd.Series(raw).rolling(int(window), center=True, min_periods=1)
            .quantile(float(quantile)).values)


def zscore(sig: np.ndarray, baseline: tuple, ddof: int) -> np.ndarray:
    seg = sig[baseline[0]:baseline[1]]
    return (sig - seg.mean()) / max(seg.std(ddof=ddof), 1e-6)


def load_raw_traces(cfg):
    """
    Read the raw trace of every regressor file.
    Returns (labels, [np.ndarray, ...]). Raises PipelineError with ALL problems listed.
    """
    require_cfg(cfg, ["regressor_files", "wash_end"], getattr(cfg, "__name__", "config"))
    files = list(cfg.regressor_files or [])
    if not files:
        raise PipelineError(
            "regressor_files is empty. The tonic GLM needs a regressor trace you define "
            "in the config, e.g.\n"
            "    regressor_files = ['/resnick/.../<fish>/hcrt_all.csv', ...]\n"
            "(there is no CSTR-capsaicin fallback any more)."
        )
    col = _setting(cfg, "regressor_column")
    n_needed = int(cfg.wash_end)

    problems, labels, raws = [], [], []
    for f in files:
        p = Path(f)
        if not p.exists():
            problems.append(f"{p} — file not found")
            continue
        try:
            df = pd.read_csv(p)
        except Exception as e:   # noqa: BLE001
            problems.append(f"{p} — cannot read as CSV ({type(e).__name__}: {e})")
            continue
        if col not in df.columns:
            problems.append(f"{p} — no column {col!r} (columns: {list(df.columns)[:8]})")
            continue
        raw = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        n_bad = int((~np.isfinite(raw)).sum())
        if n_bad:
            problems.append(f"{p} — {n_bad} non-numeric/NaN value(s) in {col!r}")
            continue
        if raw.size < n_needed:
            problems.append(f"{p} — {raw.size} rows but the recording has {n_needed} "
                            f"volumes (wash_end); trace must cover the whole recording")
            continue
        labels.append(p.parent.name if p.stem == "hcrt_all" else p.stem)
        raws.append(raw)
    if problems:
        raise PipelineError("Cannot load the GLM regressor:\n  - " + "\n  - ".join(problems))
    return labels, raws


def build_regressor(cfg) -> np.ndarray:
    """
    The regressor u(t) passed to every GLM stage (float32, length >= wash_end).

    "tonic_zscore" reproduces the legacy run_glm.py computation exactly:
    pandas centered rolling 20th percentile, z-scored on [0, drug_start) with ddof=1,
    averaged across files. (This normalisation is deliberately NOT tied to
    EPOCH_DISCARD_MIN: changing it would change every GLM fit.)
    """
    mode = _setting(cfg, "regressor_preprocess")
    if mode not in ("tonic_zscore", "none"):
        raise PipelineError(f"regressor_preprocess must be 'tonic_zscore' or 'none', "
                            f"got {mode!r}")
    _, raws = load_raw_traces(cfg)
    n_min = min(r.size for r in raws)

    if mode == "tonic_zscore":
        win, q = _setting(cfg, "regressor_tonic_window"), _setting(cfg, "regressor_tonic_quantile")
        base = (0, int(cfg.drug_start))
        traces = [zscore(rolling_tonic(r, win, q), base, ddof=1) for r in raws]
    else:
        traces = [r for r in raws]
    u = np.mean([t[:n_min] for t in traces], axis=0).astype(np.float32)
    return u


def check_regressor_vs_recording(u: np.ndarray, n_vols: int, name: str) -> None:
    if u.size < n_vols:
        raise PipelineError(f"regressor has {u.size} samples but {name} has {n_vols} volumes")


def reference_traces(cfg, baseline: tuple, ddof: int = 0):
    """
    Per-file F_tonic and F_phasic reference traces, z-scored on `baseline`
    (frames) — used only to draw the reference row of the temporal maps next to the
    cells' ΔZ, so it takes the same baseline summary window and ddof as the cells.

    Returns (labels, tonic_z [n_files, T], phasic_z [n_files, T]).
    """
    labels, raws = load_raw_traces(cfg)
    win, q = _setting(cfg, "regressor_tonic_window"), _setting(cfg, "regressor_tonic_quantile")
    n_min = min(r.size for r in raws)
    tonic, phasic = [], []
    for raw in raws:
        ft = rolling_tonic(raw, win, q)
        fp = (raw - ft) / np.maximum(ft, 1e-6)
        tonic.append(zscore(ft, baseline, ddof)[:n_min])
        phasic.append(zscore(fp, baseline, ddof)[:n_min])
    return labels, np.stack(tonic), np.stack(phasic)
