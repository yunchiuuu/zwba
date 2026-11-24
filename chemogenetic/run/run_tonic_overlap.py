"""
run_tonic_overlap.py   (optional; needs stage F1 and stage F2 outputs)
======================================================================
How much do the two tonic methods pick the same cells?  Per fish, for POS and NEG sets separately:

    n_dz, n_glm, n_both        set sizes and intersection (indices into the full cell array)
    jaccard                    n_both / n_union
    enrichment                 n_both / expected,  expected = n_dz * n_glm / N_valid  (hypergeometric mean; 1 = chance)
    sign_conflict              cells that are POS in one method and NEG in the other (should be ~0)

Why enrichment and not only Jaccard: the ΔZ set size is fixed by p (top-5% per σ stratum) while the GLM set size
is set by a null, so Jaccard is capped by the size mismatch; enrichment says whether the overlap beats chance.

Files (saved in dir_analysis/comparisons/<COMPARISON_TAG>/figures/, parameters in names and titles):
    tonic_method_overlap_dz_top{p}_vs_glm_{NULL_TAG}_nullp{N}.npz   the table (step 1)
    tonic_method_overlap_dz_top{p}_vs_glm_{NULL_TAG}_nullp{N}.csv
    tonic_method_overlap_dz_top{p}_vs_glm_{NULL_TAG}_nullp{N}.png   drawn from the npz (step 2)

Usage:
    python chemogenetic/run/run_tonic_overlap.py --config config_hcrt_trpv1_csn_120min
"""
import argparse
import importlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="Config module under chemogenetic/config/")
args = parser.parse_args()
cfg = importlib.import_module(f"chemogenetic.config.{args.config}")

from chemogenetic.tonic_dz import ptag, responder_paths, tonic_dz_params     # noqa: E402
from utils.data_io import fish_dir                                           # noqa: E402

P       = tonic_dz_params(cfg)["top_percent"]
GLM_N   = int(cfg.RESPONDER_NULL_THRESH)
NULL    = cfg.NULL_TAG
FIG_DIR = cfg.comparison_fig_dir(cfg.dir_analysis, cfg.COMPARISON_TAG)
STEM    = f"tonic_method_overlap_dz_top{ptag(P)}_vs_glm_{NULL}_nullp{GLM_N}"
GROUPS  = [(cfg.CTRL_TAG, cfg.ctrl_fish, "#4477aa"), (cfg.EXPT_TAG, cfg.expt_fish, "#cc6677")]
TITLE   = f"tonic ΔZ top {P:g}%  vs  GLM ΔR² {NULL} null p{GLM_N}  |  {args.config}"


def _n_cells(d):
    """Cell count from the mask (always present after stage C), else from the medoids."""
    for name in ("cell_valid.npy", "medoids_template_vox.npy"):
        if (d / name).exists():
            return np.load(str(d / name), mmap_mode="r").shape[0]
    return None


def _n_valid(d):
    m = d / "cell_valid.npy"
    return int(np.load(str(m)).sum()) if m.exists() else _n_cells(d)


def _one(fish):
    d = fish_dir(cfg.dir_analysis, fish)
    dz_p, dz_n = responder_paths(d, P)
    gl_p = d / f"tonic_pos_glm_{NULL}_nullp{GLM_N}_idxs.npy"
    gl_n = d / f"tonic_neg_glm_{NULL}_nullp{GLM_N}_idxs.npy"
    miss = [p.name for p in (dz_p, dz_n, gl_p, gl_n) if not p.exists()]
    if miss:
        print(f"  ⚠️  {fish[1]}: missing {miss} — skipped")
        return None
    S = {k: set(np.load(str(p)).astype(np.int64).tolist())
         for k, p in (("dzP", dz_p), ("dzN", dz_n), ("glP", gl_p), ("glN", gl_n))}
    N = _n_valid(d)
    row = {"N_valid": N}
    for sign, a, b in (("pos", "dzP", "glP"), ("neg", "dzN", "glN")):
        nd, ng, nb = len(S[a]), len(S[b]), len(S[a] & S[b])
        union = nd + ng - nb
        exp = nd * ng / N if N else np.nan
        row[f"n_dz_{sign}"], row[f"n_glm_{sign}"], row[f"n_both_{sign}"] = nd, ng, nb
        row[f"jaccard_{sign}"] = nb / union if union else np.nan
        row[f"enrich_{sign}"] = nb / exp if exp and exp > 0 else np.nan
    row["conflict_dzPos_glmNeg"] = len(S["dzP"] & S["glN"])
    row["conflict_dzNeg_glmPos"] = len(S["dzN"] & S["glP"])
    return row


def compute():
    rows = []
    for tag, fish_list, _ in GROUPS:
        for fish in fish_list:
            r = _one(fish)
            if r is not None:
                r.update(group=tag, fish=fish[1])
                rows.append(r)
    return rows


def save(rows):
    out = Path(FIG_DIR); out.mkdir(parents=True, exist_ok=True)
    keys = ["group", "fish"] + [k for k in rows[0] if k not in ("group", "fish")]
    np.savez_compressed(str(out / f"{STEM}.npz"), **{k: np.array([r[k] for r in rows]) for k in keys})
    with open(out / f"{STEM}.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    print(f"  → saved: {STEM}.npz / .csv")


def plot():
    z = np.load(str(Path(FIG_DIR) / f"{STEM}.npz"))
    grp = z["group"]
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.5))
    for i, sign in enumerate(("pos", "neg")):
        for j, (key, lab, log) in enumerate((
                (f"n_dz_{sign}", "set size (dz = solid, glm = open)", True),
                (f"jaccard_{sign}", "Jaccard (dz ∩ glm / dz ∪ glm)", False),
                (f"enrich_{sign}", "enrichment over chance (1 = chance)", True))):
            ax = axs[i, j]
            for k, (tag, _, col) in enumerate(GROUPS):
                m = grp == tag
                x = k + (np.random.default_rng(0).random(m.sum()) - 0.5) * 0.25
                if j == 0:
                    ax.scatter(x - 0.0, z[f"n_dz_{sign}"][m], color=col, s=22)
                    ax.scatter(x + 0.0, z[f"n_glm_{sign}"][m], facecolor="none", edgecolor=col, s=22)
                else:
                    ax.scatter(x, z[key][m], color=col, s=26)
                    if np.isfinite(z[key][m]).any():
                        ax.hlines(np.nanmedian(z[key][m]), k - 0.3, k + 0.3, color="k", lw=1.5)
            if log:
                ax.set_yscale("log")
            if j == 2:
                ax.axhline(1, color="grey", ls=":")
            ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels([g[0] for g in GROUPS], fontsize=8)
            ax.set_ylabel(lab, fontsize=8); ax.tick_params(labelsize=7)
            ax.set_title(f"{'▲ POS' if sign == 'pos' else '▼ NEG'}", fontsize=9,
                         color="#cc2222" if sign == "pos" else "#1155cc", loc="left")
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
    conf = z["conflict_dzPos_glmNeg"] + z["conflict_dzNeg_glmPos"]
    fig.suptitle(f"{TITLE}\nfish = unit (n={len(grp)}); sign-conflicting cells: median {np.median(conf):g}, max {conf.max():g}",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = Path(FIG_DIR) / f"{STEM}.png"
    fig.savefig(str(out), dpi=150, facecolor="white"); plt.close(fig)
    print(f"  → saved: {out.name}")


def main():
    print(f"Config: {args.config}\n{TITLE}\nOutput → {FIG_DIR}\n")
    rows = compute()
    if not rows:
        print("⚠️  no fish with both ΔZ (stage F1) and GLM (stage F2) responder files — nothing to compare")
        sys.exit(1)
    print(f"  {len(rows)} fish compared")
    save(rows)
    plot()
    print("\nAll done.")


if __name__ == "__main__":
    main()
