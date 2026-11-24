"""
run_figures_dprime_null.py   (phasic d′ per-fish null diagnostics; step 2: reads the saved .npy files)
=====================================================================================================
Everything is read from the per-fish files written by run_dprime.py (nothing is recomputed from F_phasic):
    phasic_dprime_cells_raw.npy, phasic_dprime_null__iaaft.npy, phasic_dprime_iaaft_nullp{N}_thresh.npy,
    phasic_dprime_null_meta.json, phasic_pos/neg_dprime_iaaft_nullp{N}_idxs.npy

Figures (dir_analysis/comparisons/<COMPARISON_TAG>/figures/, parameters in names and titles):
    phasic_dprime_null_thresholds_p{N}.png     per-fish threshold T by group; T vs the width of the real
                                               d′ distribution; responder fraction vs T. Tells you whether a
                                               group difference in the fraction is just a difference in noise floor.
    phasic_dprime_null_vs_real_p{N}.png        one panel per fish: survival curve of |real d′| (valid cells)
                                               against |null d′|, T marked. Excess over the null = responders.
    phasic_dprime_null_heterogeneity_p{N}.png  per-cell null SD (SD over surrogates of each sampled cell):
                                               its spread per fish. A wide spread means one pooled T per fish is too
                                               lenient for noisy cells and too strict for quiet ones.

Usage:
    python chemogenetic/run/run_figures_dprime_null.py --config config_hcrt_trpv1_csn_120min
"""
import argparse
import importlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="Config module under chemogenetic/config/")
args = parser.parse_args()
cfg = importlib.import_module(f"chemogenetic.config.{args.config}")

from utils.data_io import fish_dir                               # noqa: E402
from utils.run_guard import RunGuard, run_main                   # noqa: E402

P       = int(getattr(cfg, "PHASIC_RESPONDER_NULL_THRESH", 99))
MODE    = "raw"
FIG_DIR = cfg.comparison_fig_dir(cfg.dir_analysis, cfg.COMPARISON_TAG)
GROUPS  = [(cfg.CTRL_TAG, cfg.ctrl_fish), (cfg.EXPT_TAG, cfg.expt_fish)]
META    = cfg.PLOT_META
DPI     = 150
GUARD   = RunGuard("figures_dprime_null")


def _load(fish):
    d = fish_dir(cfg.dir_analysis, fish)
    need = [d / f"phasic_dprime_cells_{MODE}.npy", d / "phasic_dprime_null__iaaft.npy",
            d / f"phasic_dprime_iaaft_nullp{P}_thresh.npy", d / "phasic_dprime_null_meta.json",
            d / f"phasic_pos_dprime_iaaft_nullp{P}_idxs.npy", d / f"phasic_neg_dprime_iaaft_nullp{P}_idxs.npy"]
    missing = [p.name for p in need if not p.exists()]
    if missing:
        GUARD.fail(fish[1], f"missing {missing} — run run_dprime.py (stages null, responders) first")
        return None
    dp = np.load(str(need[0]))
    null = np.load(str(need[1]))
    thr = float(np.load(str(need[2])))
    meta = json.loads(need[3].read_text())
    n_pos, n_neg = (np.load(str(p)).size for p in need[4:6])
    fin = dp[np.isfinite(dp)]
    nf = null[np.isfinite(null)]
    # per-cell null SD: SD over the surrogate axis of each sampled cell
    with np.errstate(all="ignore"):
        cell_sd = np.nanstd(null, axis=0)
    cell_sd = cell_sd[np.isfinite(cell_sd) & (cell_sd > 0)]
    # is the spread real? split-half reliability of the per-cell SD (surrogates split in two)
    half = null.shape[0] // 2
    r_half = np.nan
    if half >= 20:
        with np.errstate(all="ignore"):
            a, b = np.nanstd(null[:half], axis=0), np.nanstd(null[half:2 * half], axis=0)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() > 20:
            r_half = float(np.corrcoef(a[ok], b[ok])[0, 1])
    return {"r_half": r_half, "dp": fin, "null_abs": np.abs(nf), "thr": thr, "meta": meta, "cell_sd": cell_sd,
            "n_pos": n_pos, "n_neg": n_neg, "n_valid": fin.size}


def _mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.size < 2 or b.size < 2:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def _plabel(p):
    return "n/a" if not np.isfinite(p) else f"MWU p={p:.3f}"


def _survival(x, n_pts=400):
    """P(|x| >= v) on a log-spaced grid (for log-y plotting)."""
    x = np.sort(np.abs(x))
    v = np.linspace(0, float(x[-1]) if x.size else 1.0, n_pts)
    s = 1.0 - np.searchsorted(x, v, side="left") / max(x.size, 1)
    return v, np.clip(s, 1.0 / max(x.size, 1), 1)


def main():
    data = {tag: [(f, _load(f)) for f in fish] for tag, fish in GROUPS}
    data = {tag: [(f, d) for f, d in lst if d is not None] for tag, lst in data.items()}
    if not any(data.values()):
        GUARD.fail("all", "no fish with a d′ null and responders")
        GUARD.finish()
        return
    any_meta = next(d["meta"] for lst in data.values() for _, d in lst)
    tmeta = (f"per-fish IAAFT null p{P} | {any_meta.get('n_surrogates', '?')} surrogates x "
             f"≤{any_meta.get('n_cells_null', '?')} cells | {args.config}")

    # ── 1. thresholds ────────────────────────────────────────────────────
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axs[0]
    xs = {}
    for k, (tag, lst) in enumerate(data.items()):
        v = np.array([d["thr"] for _, d in lst])
        xs[tag] = v
        if v.size:
            ax.boxplot([v], positions=[k], widths=0.45, showfliers=False)
            ax.scatter(k + np.random.default_rng(0).uniform(-0.1, 0.1, v.size), v, s=34,
                       color=META[tag]["color"], edgecolor="k", lw=0.5, zorder=3)
    tags = list(data.keys())
    ax.set_xticks(range(len(tags)))
    ax.set_xticklabels([f"{t}\n(n={len(data[t])})" for t in tags])
    ax.set_ylabel(f"per-fish threshold T (p{P} of |null d′|)")
    ax.set_title(f"T per fish | {_plabel(_mwu(*[xs[t] for t in tags[:2]]))}" if len(tags) >= 2 else "T per fish",
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axs[1]
    for tag, lst in data.items():
        for f, d in lst:
            q = np.percentile(d["dp"], [25, 75]) if d["dp"].size else [np.nan, np.nan]
            ax.scatter((q[1] - q[0]) / 1.349, d["thr"], s=40, color=META[tag]["color"], edgecolor="k", lw=0.5)
    ax.set_xlabel("robust SD of real d′ over valid cells (IQR/1.349)")
    ax.set_ylabel("T")
    ax.set_title("noise floor: T should track the width of the bulk of real d′", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axs[2]
    for tag, lst in data.items():
        for f, d in lst:
            ax.scatter(d["thr"], 100 * (d["n_pos"] + d["n_neg"]) / max(d["n_valid"], 1), s=40,
                       color=META[tag]["color"], edgecolor="k", lw=0.5, label=tag)
    h, l = ax.get_legend_handles_labels()
    ax.legend(dict(zip(l, h)).values(), dict(zip(l, h)).keys(), fontsize=8, frameon=False)
    ax.set_xlabel("T")
    ax.set_ylabel("responders, % of valid cells (pos + neg)")
    ax.set_title("responder fraction vs T (a fraction that only tracks T is a noise-floor effect)", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Phasic d′ thresholds per fish | {tmeta}", fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / f"phasic_dprime_null_thresholds_p{P}.png"
    fig.savefig(str(out), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved: {out}")

    # ── 2. real vs null survival, one panel per fish ──────────────────────
    allf = [(tag, f, d) for tag, lst in data.items() for f, d in lst]
    ncol = 4
    nrow = int(np.ceil(len(allf) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax in axs.ravel():
        ax.axis("off")
    for ax, (tag, f, d) in zip(axs.ravel(), allf):
        ax.axis("on")
        v, s_real = _survival(d["dp"])
        vn, s_null = _survival(d["null_abs"])
        ax.semilogy(v, s_real, color=META[tag]["color"], lw=1.6, label="real |d′|")
        ax.semilogy(vn, s_null, color="k", lw=1.2, ls="--", label="null |d′|")
        ax.axvline(d["thr"], color="r", lw=0.9, ls=":")
        ax.axhline(1 - P / 100, color="gray", lw=0.6, ls=":")
        ax.set_title(f"{f[1][:34]}\n{tag} | T=±{d['thr']:.3f} | pos {d['n_pos']} neg {d['n_neg']} "
                     f"/ {d['n_valid']:,}", fontsize=7)
        ax.set_xlim(0, max(float(np.percentile(np.abs(d["dp"]), 99.9)) if d["dp"].size else 1, d["thr"] * 1.5))
        ax.set_xlabel("|d′|", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    axs.ravel()[0].legend(fontsize=6, frameon=False)
    fig.suptitle(f"Real vs null |d′| (survival; red dotted = T, grey dotted = {100 - P}% null tail) | {tmeta}",
                 fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / f"phasic_dprime_null_vs_real_p{P}.png"
    fig.savefig(str(out), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved: {out}")

    # ── 3. null heterogeneity ─────────────────────────────────────────────
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axs[0]
    pos = 0
    ticks, tl = [], []
    for tag, lst in data.items():
        for f, d in lst:
            if d["cell_sd"].size:
                ax.boxplot([d["cell_sd"]], positions=[pos], widths=0.6, showfliers=False,
                           patch_artist=True, boxprops=dict(facecolor=META[tag]["color"], alpha=0.5))
            ticks.append(pos)
            tl.append(f[1][:18])
            pos += 1
        pos += 1
    ax.set_xticks(ticks)
    ax.set_xticklabels(tl, rotation=90, fontsize=6)
    ax.set_ylabel("per-cell null SD of d′ (over surrogates)")
    ax.set_title("null SD across sampled cells, per fish (box = IQR, whiskers = 1.5 IQR)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axs[1]
    rat = {}
    for k, (tag, lst) in enumerate(data.items()):
        v = np.array([np.percentile(d["cell_sd"], 95) / np.percentile(d["cell_sd"], 5)
                      for _, d in lst if d["cell_sd"].size > 20])
        rat[tag] = v
        if v.size:
            ax.scatter(k + np.random.default_rng(1).uniform(-0.1, 0.1, v.size), v, s=38,
                       color=META[tag]["color"], edgecolor="k", lw=0.5)
    ax.axhline(1.5, color="gray", ls=":", lw=0.8)
    ax.axhline(3.0, color="r", ls=":", lw=0.8)
    ax.set_xticks(range(len(tags)))
    ax.set_xticklabels(tags)
    ax.set_ylabel("p95 / p5 of per-cell null SD")
    _rh = [d["r_half"] for lst in data.values() for _, d in lst if np.isfinite(d["r_half"])]
    ax.set_title("spread ratio per fish (grey 1.5 = fine, red 3 = calibrate per cell)\n"
                 f"split-half r of per-cell SD = {np.median(_rh) if _rh else float('nan'):.2f} "
                 "(near 1 = the spread is real, near 0 = sampling noise)", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Heterogeneity of the d′ null across cells | {tmeta}", fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / f"phasic_dprime_null_heterogeneity_p{P}.png"
    fig.savefig(str(out), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → saved: {out}")
    for tag, v in rat.items():
        if v.size:
            rh = [d["r_half"] for _, d in data[tag] if np.isfinite(d["r_half"])]
            print(f"  spread ratio p95/p5 ({tag}): median {np.median(v):.2f}, max {v.max():.2f}"
                  f" | split-half r of per-cell null SD: {np.median(rh) if rh else float('nan'):.2f}")
    GUARD.finish()


if __name__ == "__main__":
    run_main(main)
