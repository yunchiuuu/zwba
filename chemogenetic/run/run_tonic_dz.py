"""
run_tonic_dz.py   (stage F1)
============================
Tonic responders from ΔZ of F_tonic (replaces the GLM ΔR² as the definition; the GLM is the optional
stage F2). Definition, outputs and what the sets can and cannot be used for: see chemogenetic/tonic_dz.py.

    ΔZ_w = (mean F_tonic in window w − baseline mean) / baseline σ          (seven temporal windows)
    max ΔZ = ΔZ at the window of largest |ΔZ|; responders = top p % of |max ΔZ| within baseline-σ strata.

Prerequisites (per fish):
    run_decompose.py  -> f_tonic.npy and cell_valid.npy   (the cell filter is built here if missing)

Settings (config; defaults apply when a config lacks them):
    TONIC_DZ_TOP_PERCENT = 5.0           main set
    TONIC_DZ_TOP_PERCENT_SENS = [2.5, 10.0]   sensitivity sets, saved next to the main one
    TONIC_DZ_N_STRATA = 5                baseline-σ strata per fish (1 = no stratification)

Usage:
    bash chemogenetic/submit/submit_tonic_dz.sh --config config_hcrt_trpv1_csn_120min
    python chemogenetic/run/run_tonic_dz.py --config config_hcrt_trpv1_csn_120min [--overwrite] [--fish TEXT] [--exclude TEXT]

Fish whose outputs exist AND were made with the same epochs, cell filter and settings are skipped.
Failures are reported with tracebacks and the job exits 1 if any fish failed (exit 2 = bad config / missing input).
"""
import argparse
import importlib
import sys
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="Config module under chemogenetic/config/")
parser.add_argument("--overwrite", action="store_true", help="recompute even if current outputs exist")
parser.add_argument("--fish", action="append", default=[],
                    help="only fish whose expt_ID contains this text (repeatable)")
parser.add_argument("--exclude", action="append", default=[],
                    help="leave out fish whose expt_ID contains this text (repeatable)")
parser.add_argument("--n-jobs", type=int, default=4, help="fish processed in parallel (default 4)")
args = parser.parse_args()

cfg = importlib.import_module(f"chemogenetic.config.{args.config}")

from chemogenetic.cell_filter import ensure_cell_valid, filter_params       # noqa: E402
from chemogenetic.epochs import get_epochs                                    # noqa: E402
from chemogenetic.tonic_dz import tonic_dz_one_fish, tonic_dz_params          # noqa: E402
from utils.data_io import fish_dir                                            # noqa: E402
from utils.run_guard import (RunGuard, guarded, require_paths, run_main,      # noqa: E402
                             setup_or_exit)

EPOCHS = setup_or_exit(get_epochs, cfg)
PARAMS = setup_or_exit(tonic_dz_params, cfg)
dir_analysis = cfg.dir_analysis


@guarded
def _run_one(fish, epochs, params, overwrite):
    return tonic_dz_one_fish(fish, dir_analysis, epochs, params, overwrite=overwrite)


def main():
    guard = RunGuard("tonic_dz")
    all_fish = [f for f in cfg.all_fish
                if (not args.fish or any(t in f[1] for t in args.fish))
                and not any(t in f[1] for t in args.exclude)]
    if not all_fish:
        print("❌ no fish left after --fish / --exclude")
        sys.exit(2)

    print(f"Tonic ΔZ responders for {len(all_fish)} fish")
    print(f"  dir_analysis : {dir_analysis}")
    print(f"  epochs       : {EPOCHS.describe()}")
    print(f"  statistic    : max |ΔZ| over the {EPOCHS.n_windows} windows "
          f"{[lab for lab, _, _ in EPOCHS.windows]}")
    print(f"  responders   : top {PARAMS['top_percent']:g}% of |max ΔZ| per fish, within "
          f"{PARAMS['n_strata']} baseline-σ strata; sensitivity {PARAMS['percents'][1:]}")
    print(f"  cell filter  : {filter_params(cfg)}")
    print("  NOTE: rank-based, no null: the responder FRACTION is fixed by p and not comparable between groups.\n")

    require_paths(
        [(fish_dir(dir_analysis, f) / "f_tonic.npy", "run run_decompose.py first") for f in all_fish],
        "tonic ΔZ input F_tonic",
    )
    for fish in all_fish:                      # builds cell_valid.npy when missing / stale
        with guard.attempt(f"cell filter: {fish[1]}"):
            ensure_cell_valid(fish_dir(dir_analysis, fish), EPOCHS, cfg, verbose=False)

    results = Parallel(n_jobs=max(1, args.n_jobs), backend="loky")(
        delayed(_run_one)(f, EPOCHS, PARAMS, args.overwrite) for f in all_fish)
    guard.record_results("tonic ΔZ", results)

    ok = [r for r in results if r.get("status") == "ok"]
    if ok:
        tot = np.sum([r["argmax_hist"] for r in ok], axis=0)
        share = 100 * tot / max(tot.sum(), 1)
        print("\nArgmax-window share of ALL cells (0-based window: share):")
        for (lab, s, e), sh in zip(EPOCHS.windows, share):
            print(f"   {lab:>10s}  (abs {s / 60:g}-{e / 60:g}): {sh:5.1f}%")
        late = share[-2:].sum()
        if late > 50:
            guard.warn("argmax windows", f"{late:.0f}% of cells peak in the last two windows (washout). The "
                                         f"max-over-windows ranking is then dominated by washout / late "
                                         f"whole-brain shifts rather than the drug phase; look at "
                                         f"run_figures_tonic_dz.py before using the sets.")
    guard.finish()


if __name__ == "__main__":
    run_main(main)
