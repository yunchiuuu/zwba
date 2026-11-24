"""
run_figures_tonic_dz.py   (stage F1 figures; step 2: loads the .npy files made by run_tonic_dz.py)
==============================================================================================
Diagnostics for the tonic ΔZ responder sets and the group comparison. Everything is read from the
per-fish files written by run_tonic_dz.py (nothing is recomputed from F_tonic).

Figures (saved in dir_analysis/comparisons/<COMPARISON_TAG>/figures/, parameters in names and titles):
    tonic_dz_argmax_window_top{p}.png      which window each cell peaks in (all valid cells, and the top-p% set),
                                           one line per fish, washout windows shaded. If most cells peak in the
                                           late windows the max-over-windows ranking is dominated by washout.
    tonic_dz_sigma_vs_absdz_top{p}.png     baseline σ vs |max ΔZ| per group (pooled cells, density); a flat or
                                           tiny σ inflates ΔZ. Dashed lines = the per-stratum top-p% cutoffs.
    tonic_dz_ecdf_top{p}.png               per-fish ECDFs of |max ΔZ| and signed max ΔZ, CTRL vs EXPT overlaid,
                                           and the per-fish 95th percentile of |max ΔZ| (fish = the unit, MWU).

The responder FRACTION is fixed by p, so it is not plotted as a group comparison; use the ECDFs.

Usage:
    python chemogenetic/run/run_figures_tonic_dz.py --config config_hcrt_trpv1_csn_120min
"""
import argparse
import importlib
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

from chemogenetic.epochs import get_epochs                                   # noqa: E402
from chemogenetic.tonic_dz import ptag, responder_paths, tonic_dz_params     # noqa: E402
from utils.data_io import fish_dir                                           # noqa: E402
from utils.run_guard import RunGuard, run_main, setup_or_exit                # noqa: E402

EPOCHS = setup_or_exit(get_epochs, cfg)
PARAMS = setup_or_exit(tonic_dz_params, cfg)
P      = PARAMS["top_percent"]
PT     = ptag(P)
FIG_DIR = cfg.comparison_fig_dir(cfg.dir_analysis, cfg.COMPARISON_TAG)
GROUPS = [(cfg.CTRL_TAG, cfg.ctrl_fish), (cfg.EXPT_TAG, cfg.expt_fish)]
META   = cfg.PLOT_META
DPI    = 150
GUARD  = RunGuard("figures_tonic_dz")
SR     = float(cfg.sampling_rate_hz)
TITLE_META = f"top {P:g}% |max ΔZ| | {PARAMS['n_strata']} σ strata | {args.config}"


def _load(fish):
    d = fish_dir(cfg.dir_analysis, fish)
    need = [d / "tonic_dz_max.npy", d / "tonic_dz_argmax_window.npy", d / "tonic_dz_sigma_baseline.npy",
            d / "tonic_dz_stratum.npy", d / "tonic_dz_meta.json", *responder_paths(d, P)]
    missing = [p.name for p in need if not p.exists()]
    if missing:
        GUARD.fail(fish[1], f"tonic ΔZ outputs missing {missing} — run run_tonic_dz.py first")
        return None
    import json
    pos, neg = (np.load(str(p)) for p in responder_paths(d, P))
    return {
        "max": np.load(str(d / "tonic_dz_max.npy")),
        "arg": np.load(str(d / "tonic_dz_argmax_window.npy")),
        "sig": np.load(str(d / "tonic_dz_sigma_baseline.npy")),
        "strat": np.load(str(d / "tonic_dz_stratum.npy")),
        "pos": pos, "neg": neg,
        "meta": json.loads((d / "tonic_dz_meta.json").read_text()),
    }


def _ecdf(x):
    x = np.sort(x[np.isfinite(x)])
    return x, np.arange(1, x.size + 1) / max(x.size, 1)


def main():
    data = {tag: [(f, _load(f)) for f in fish] for tag, fish in GROUPS}
    data = {tag: [(f, d) for f, d in lst if d is not None] for tag, lst in data.items()}
    if not any(data.values()):
        GUARD.fail("all", "no fish with tonic ΔZ outputs")
        GUARD.finish()
        return
    wins = EPOCHS.windows
    drug_end = EPOCHS.drug[1]
    labels = [f"{lab}\n(abs {s / SR / 60:g}-{e / SR / 60:g})" for lab, s, e in wins]
    washout = [s >= drug_end for _, s, _ in wins]
    straddle = [(s < drug_end <= e - 1) for _, s, e in wins]

    # ── 1. argmax window ────────────────────────────────────────────────
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), sharey="row")
    for j, (tag, lst) in enumerate(data.items()):
        for i, (rowname, sel) in enumerate([("all valid cells", None), (f"top {P:g}% set", "set")]):
            ax = axs[i, j]
            for w, (wo, st) in enumerate(zip(washout, straddle)):
                if wo:
                    ax.axvspan(w - 0.5, w + 0.5, color="#dddddd", zorder=0)
                elif st:
                    ax.axvspan(w - 0.5, w + 0.5, color="#f1f1f1", zorder=0)
            fr = []
            for f, d in lst:
                idx = np.where(d["arg"] >= 0)[0] if sel is None else np.union1d(d["pos"], d["neg"]).astype(int)
                a = d["arg"][idx]
                a = a[a >= 0]
                h = np.bincount(a, minlength=len(wins)) / max(a.size, 1)
                fr.append(h)
                ax.plot(range(len(wins)), h, "-o", ms=3, lw=0.8, color=META[tag]["color"], alpha=0.5)
            if fr:
                ax.plot(range(len(wins)), np.mean(fr, 0), "-o", color="k", lw=2, ms=5, label="mean of fish")
            ax.set_title(f"{tag} (n={len(lst)} fish) | {rowname}", fontsize=9)
            ax.set_xticks(range(len(wins)))
            ax.set_xticklabels(labels, fontsize=6)
            if j == 0:
                ax.set_ylabel("fraction of cells peaking in window")
            ax.legend(fontsize=7, loc="upper left")
    fig.suptitle(f"Argmax window of |ΔZ| (grey = washout, light grey = straddles drug end) | {TITLE_META}", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIG_DIR / f"tonic_dz_argmax_window_top{PT}.png"
    fig.savefig(str(out), dpi=DPI, facecolor="white"); plt.close(fig); print(f"  → saved: {out}")

    # ── 2. sigma vs |max dZ| ────────────────────────────────────────────
    rng = np.random.default_rng(0)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    for ax, (tag, lst) in zip(axs, data.items()):
        sx, sy, cut_lines = [], [], []
        for f, d in lst:
            ok = np.flatnonzero(np.isfinite(d["max"]) & np.isfinite(d["sig"]) & (d["sig"] > 0))
            if ok.size:
                k = rng.choice(ok, size=min(50_000, ok.size), replace=False)
                sx.append(d["sig"][k]); sy.append(np.abs(d["max"][k]))
            c = d["meta"]["cutoff_abs_dz_per_stratum"].get(PT, [])
            e = d["meta"].get("sigma_stratum_edges", [])
            cut_lines.append((e, c))
        if sx:
            x, y = np.concatenate(sx), np.concatenate(sy)
            hb = ax.hexbin(x, np.clip(y, 1e-3, None), gridsize=70, xscale="log", yscale="log", bins="log",
                           cmap="viridis", mincnt=1)
            fig.colorbar(hb, ax=ax, label="cells (log)")
        ax.set_title(f"{tag} (n={len(lst)} fish; ≤50k cells/fish)", fontsize=9)
        ax.set_xlabel("baseline σ of F_tonic (counts)")
        ax.axhline(cfg.__dict__.get("CELL_MAX_ABS_DZ", 200.0) or 1e9, color="r", lw=0.8, ls=":")
    axs[0].set_ylabel("|max ΔZ|")
    fig.suptitle(f"Baseline σ vs |max ΔZ| (flat baselines inflate ΔZ; red dotted = cell-filter backstop) | {TITLE_META}",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / f"tonic_dz_sigma_vs_absdz_top{PT}.png"
    fig.savefig(str(out), dpi=DPI, facecolor="white"); plt.close(fig); print(f"  → saved: {out}")

    # ── 3. ECDFs + per-fish p95 ─────────────────────────────────────────
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw={"width_ratios": [1, 1, 0.6]})
    p95 = {}
    for tag, lst in data.items():
        col, al = META[tag]["color"], META[tag]["alpha"]
        vals = []
        for f, d in lst:
            m = d["max"]
            m = m[np.isfinite(m)]
            if m.size == 0:
                continue
            x, y = _ecdf(np.abs(m)); axs[0].plot(x, y, color=col, alpha=0.5, lw=1)
            x, y = _ecdf(m);         axs[1].plot(x, y, color=col, alpha=0.5, lw=1)
            vals.append(np.percentile(np.abs(m), 95))
        p95[tag] = np.array(vals)
        axs[0].plot([], [], color=col, label=f"{tag} (n={len(vals)})")
    axs[0].set_xlabel("|max ΔZ|"); axs[0].set_ylabel("ECDF (per fish)"); axs[0].legend(fontsize=8)
    axs[0].set_xlim(0, max(1.0, np.nanpercentile(np.concatenate(
        [np.abs(d["max"][np.isfinite(d["max"])]) for lst in data.values() for _, d in lst]), 99.9)))
    axs[0].set_title("|max ΔZ|", fontsize=9)
    axs[1].set_xlabel("max ΔZ (signed)"); axs[1].set_title("signed max ΔZ (direction asymmetry)", fontsize=9)
    lo, hi = np.nanpercentile(np.concatenate([d["max"][np.isfinite(d["max"])] for lst in data.values() for _, d in lst]),
                              [0.1, 99.9])
    axs[1].set_xlim(lo, hi)
    tags = list(p95)
    for i, tag in enumerate(tags):
        v = p95[tag]
        axs[2].scatter(np.full(v.size, i) + rng.uniform(-0.08, 0.08, v.size), v, color=META[tag]["color"], s=28)
        if v.size:
            axs[2].hlines(np.median(v), i - 0.25, i + 0.25, color="k")
    axs[2].set_xticks(range(len(tags))); axs[2].set_xticklabels(tags)
    axs[2].set_ylabel("per-fish 95th percentile of |max ΔZ|")
    if len(tags) == 2 and all(p95[t].size >= 2 for t in tags):
        pv = mannwhitneyu(p95[tags[0]], p95[tags[1]], alternative="two-sided").pvalue
        axs[2].set_title(f"fish = unit | MWU p={pv:.3g}", fontsize=9)
        print(f"  per-fish p95 |max ΔZ|: {tags[0]} median {np.median(p95[tags[0]]):.2f}, "
              f"{tags[1]} median {np.median(p95[tags[1]]):.2f}, MWU p={pv:.3g}")
    fig.suptitle(f"Tonic ΔZ: group comparison without a responder fraction | {TITLE_META}", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = FIG_DIR / f"tonic_dz_ecdf_top{PT}.png"
    fig.savefig(str(out), dpi=DPI, facecolor="white"); plt.close(fig); print(f"  → saved: {out}")
    GUARD.finish()


if __name__ == "__main__":
    run_main(main)
