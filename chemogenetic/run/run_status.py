"""
run_status.py
=============
Where am I in the pipeline? One table: a row per fish, a column per stage.

    bash chemogenetic/submit/status.sh --config config_<exp>
    (or: <PYTHON_BIN> chemogenetic/run/run_status.py --config config_<exp>)

Read-only. Runs on the login node in a few seconds; submits nothing and creates no folders.

Symbols
-------
    ✔  stage output exists and is current
    ✗  missing (its prerequisites are done, so this is what to run next)
    ·  blocked: a prerequisite stage is not done yet
    ⚠  output exists but is stale / unverifiable (reason listed under the table)

Stages (column) -> what is checked
----------------------------------
    vol   voluseg      {dir_voluseg}/{proj}/{expt}/output/cells0_clean.hdf5 + volume0.hdf5
    reg   registration warp + affine + registered_final + registration_meta.json whose template
                       sha256 equals the current template                      (⚠ = other template / no meta)
    dec   decompose    f_dark_scalar.npy, f_tonic.npy, f_phasic.npy
    med   medoids      medoids_all.npy, medoids_rot_all.npy, medoids_template_vox.npy
    fit   GLM fit      glm/<param_folder_name>/kernel_delta_r2_fit.npy; ⚠ if RUN_META.json says
                       the fit window differs from INCLUDED_BASELINE in the config
    null  GLM null     kernel_delta_r2_null__<NULL_TAG>.npy
    resp  GLM resp.    tonic_pos/neg_glm_<NULL_TAG>_nullp<RESPONDER_NULL_THRESH>_idxs.npy
    dpr   d′           phasic_dprime_cells_*.npy
    dnul  d′ null      phasic_dprime_null__iaaft.npy
    dres  d′ resp.     phasic_pos/neg_dprime_iaaft_nullp<PHASIC_RESPONDER_NULL_THRESH>_idxs.npy
    map   epoch maps   epoch_dz.npy + epoch_dprime.npy; ⚠ if made with different epochs

Options
-------
    --fish TEXT   only fish whose expt_ID contains TEXT
    --strict      exit 1 if any ✗ or ⚠ is shown (for scripts)
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="Module name under chemogenetic/config/")
parser.add_argument("--fish", default=None, help="only fish whose expt_ID contains this text")
parser.add_argument("--strict", action="store_true", help="exit 1 if anything is missing or stale")
args = parser.parse_args()

if args.config.endswith(".py") or "/" in args.config:
    sys.exit(f"--config must be a module name (e.g. config_hcrt_trpv1_csn_120min), not a path: {args.config}")
try:
    cfg = importlib.import_module(f"chemogenetic.config.{args.config}")
except Exception as e:      # noqa: BLE001
    sys.exit(f"Cannot import chemogenetic/config/{args.config}.py: {type(e).__name__}: {e}")

from chemogenetic.epochs import get_epochs                    # noqa: E402
from utils.run_guard import PipelineError                     # noqa: E402

try:
    EPOCHS = get_epochs(cfg)
except PipelineError as e:
    sys.exit(f"\n❌ {e}\n")

OK, MISSING, BLOCKED, STALE = "✔", "✗", "·", "⚠"

dir_voluseg = Path(cfg.dir_voluseg)
dir_registration = Path(cfg.dir_registration)
dir_analysis = Path(cfg.dir_analysis)
NULL_TAG = getattr(cfg, "NULL_TAG", "iaaft")
GLM_P = int(getattr(cfg, "RESPONDER_NULL_THRESH", 99))
DPR_P = int(getattr(cfg, "PHASIC_RESPONDER_NULL_THRESH", 99))
PARAM_FOLDER = cfg.param_folder_name

# ── template identity (registration column) ─────────────────────────────────
template_fname = "template_mean_brain.nii.gz"
try:
    sys.path.insert(0, str(REPO / "registration" / "config"))
    import config_registration as _regcfg                      # noqa: E402
    template_fname = getattr(_regcfg, "MEAN_BRAIN_FNAME", template_fname)
except Exception:           # noqa: BLE001 - fall back to the conventional filename
    pass
template_path = dir_registration / template_fname              # same folder the registration writes to
TEMPLATE_SHA = None
if template_path.exists():
    h = hashlib.sha256()
    with open(template_path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    TEMPLATE_SHA = h.hexdigest()

# expected GLM fit start (vol) for the stale-window check
_fit_start = int(round(cfg.drug_start - float(cfg.INCLUDED_BASELINE) * 60 * cfg.sampling_rate_hz))
EXPECTED_FIT_START = max(0, _fit_start)

notes: list[str] = []


def _all(d: Path, names) -> bool:
    return all((d / n).exists() for n in names)


def _json(p: Path):
    try:
        return json.loads(p.read_text())
    except Exception:       # noqa: BLE001
        return None


def check_fish(fish):
    proj, expt = fish
    vol_d = dir_voluseg / proj / expt / "output"
    reg_d = dir_registration / proj / expt
    an_d = dir_analysis / proj / expt
    glm_d = an_d / "glm" / PARAM_FOLDER
    r: dict[str, str] = {}

    # raw prerequisites (no prerequisites of their own)
    r["vol"] = OK if _all(vol_d, ["cells0_clean.hdf5", "volume0.hdf5"]) else MISSING

    # registration: needs vol
    if r["vol"] != OK:
        r["reg"] = BLOCKED
    elif not _all(reg_d, ["expt_to_mean_warp.nii.gz", "expt_to_mean_affine.mat",
                          "expt_to_mean_registered_final.nii.gz"]):
        r["reg"] = MISSING
    else:
        meta = _json(reg_d / "registration_meta.json")
        if meta is None:
            r["reg"] = STALE
            notes.append(f"{expt}: registration has no readable registration_meta.json (template unverified)")
        elif TEMPLATE_SHA is None:
            r["reg"] = STALE
            notes.append(f"template not found at {template_path}; cannot verify registrations")
        elif meta.get("template", {}).get("sha256") != TEMPLATE_SHA:
            r["reg"] = STALE
            notes.append(f"{expt}: registered to a different template than {template_path.name}")
        else:
            r["reg"] = OK

    # decompose: needs vol
    if r["vol"] != OK:
        r["dec"] = BLOCKED
    else:
        r["dec"] = OK if _all(an_d, ["f_dark_scalar.npy", "f_tonic.npy", "f_phasic.npy"]) else MISSING

    # medoids: needs reg
    if r["reg"] not in (OK, STALE):
        r["med"] = BLOCKED
    else:
        r["med"] = OK if _all(an_d, ["medoids_all.npy", "medoids_rot_all.npy",
                                     "medoids_template_vox.npy"]) else MISSING

    # GLM chain: needs dec
    if r["dec"] != OK:
        r["fit"] = r["null"] = r["resp"] = BLOCKED
    else:
        if (glm_d / "kernel_delta_r2_fit.npy").exists():
            meta = _json(glm_d / "RUN_META.json")
            old = meta.get("fit_start") if isinstance(meta, dict) else None
            if old is not None and int(old) != EXPECTED_FIT_START:
                r["fit"] = STALE
                notes.append(f"{expt}: GLM fit starts at vol {old} but INCLUDED_BASELINE="
                             f"{cfg.INCLUDED_BASELINE:g} gives vol {EXPECTED_FIT_START} (rebuild the GLM)")
            else:
                r["fit"] = OK
        else:
            r["fit"] = MISSING
        if r["fit"] not in (OK, STALE):
            r["null"] = r["resp"] = BLOCKED
        else:
            r["null"] = OK if (glm_d / f"kernel_delta_r2_null__{NULL_TAG}.npy").exists() else MISSING
            if r["null"] != OK:
                r["resp"] = BLOCKED
            else:
                r["resp"] = OK if _all(an_d, [f"tonic_pos_glm_{NULL_TAG}_nullp{GLM_P}_idxs.npy",
                                              f"tonic_neg_glm_{NULL_TAG}_nullp{GLM_P}_idxs.npy"]) else MISSING

    # d' chain: needs dec
    if r["dec"] != OK:
        r["dpr"] = r["dnul"] = r["dres"] = BLOCKED
    else:
        r["dpr"] = OK if any(an_d.glob("phasic_dprime_cells_*.npy")) else MISSING
        if r["dpr"] != OK:
            r["dnul"] = r["dres"] = BLOCKED
        else:
            r["dnul"] = OK if (an_d / "phasic_dprime_null__iaaft.npy").exists() else MISSING
            if r["dnul"] != OK:
                r["dres"] = BLOCKED
            else:
                r["dres"] = OK if _all(an_d, [f"phasic_pos_dprime_iaaft_nullp{DPR_P}_idxs.npy",
                                              f"phasic_neg_dprime_iaaft_nullp{DPR_P}_idxs.npy"]) else MISSING

    # epoch maps: need dec + med
    if r["dec"] != OK or r["med"] != OK:
        r["map"] = BLOCKED
    elif not _all(an_d, ["epoch_dz.npy", "epoch_dprime.npy"]):
        r["map"] = MISSING
    else:
        from chemogenetic.temporal_windows import epochs_meta      # numpy -> needs the conda env
        if _json(an_d / "epoch_meta.json") == epochs_meta(EPOCHS):
            r["map"] = OK
        else:
            r["map"] = STALE
            notes.append(f"{expt}: epoch_dz/epoch_dprime were made with different (or unrecorded) "
                         f"epochs; the map script recomputes them")
    return r


COLS = ["vol", "reg", "dec", "med", "fit", "null", "resp", "dpr", "dnul", "dres", "map"]
NEXT = {
    "vol":  "voluseg: sbatch run_voluseg.py per fish (README step 3)",
    "reg":  "bash registration/submit_registration_syn.sh --config config_registration",
    "dec":  f"bash chemogenetic/submit/submit_decompose.sh --config {args.config}",
    "med":  f"bash chemogenetic/submit/submit_medoids.sh --config {args.config}",
    "fit":  f"bash chemogenetic/submit/submit_glm.sh --config {args.config}",
    "null": f"bash chemogenetic/submit/submit_glm.sh --config {args.config}",
    "resp": f"bash chemogenetic/submit/submit_glm.sh --config {args.config}",
    "dpr":  f"bash chemogenetic/submit/submit_dprime.sh --config {args.config}",
    "dnul": f"bash chemogenetic/submit/submit_dprime.sh --config {args.config}",
    "dres": f"bash chemogenetic/submit/submit_dprime.sh --config {args.config}",
    "map":  f"bash chemogenetic/submit/submit_temporal_intensity_map.sh --config {args.config}",
}

fish_list = [f for f in cfg.all_fish if args.fish is None or args.fish in f[1]]
if not fish_list:
    sys.exit(f"No fish match {args.fish!r} in {args.config}.")
group_of = getattr(cfg, "GROUP_OF_FISH", {})

rows = [(f, check_fish(f)) for f in fish_list]

w_name = max(len(f[1]) for f in fish_list)
print(f"\nconfig: {args.config}   {len(fish_list)} fish   {EPOCHS.describe()}")
print(f"GLM: {PARAM_FOLDER}  fit window starts at vol {EXPECTED_FIT_START}  null={NULL_TAG} p{GLM_P}   d′ null p{DPR_P}")
print(f"template: {template_path}" + ("" if TEMPLATE_SHA else "   (NOT FOUND)") + "\n")
print(f"{'':{w_name}}  {'group':<10}" + " ".join(f"{c:>4}" for c in COLS))
for fish, r in rows:
    g = str(group_of.get(fish, ""))[:10]
    print(f"{fish[1]:{w_name}}  {g:<10}" + " ".join(f"{r[c]:>4}" for c in COLS))
print("\n✔ done   ✗ missing (next to run)   · blocked by an earlier stage   ⚠ stale or unverifiable")

# ── what to do next ─────────────────────────────────────────────────────────
todo = []
for c in COLS:
    n_missing = sum(1 for _, r in rows if r[c] == MISSING)
    n_stale = sum(1 for _, r in rows if r[c] == STALE)
    if n_missing or n_stale:
        todo.append((c, n_missing, n_stale))
print()
if not todo:
    print("Everything is done and current for these fish.")
else:
    print("Next:")
    seen = set()
    for c, n_missing, n_stale in todo:
        cmd = NEXT[c]
        if cmd in seen:        # fit/null/resp share one command
            continue
        seen.add(cmd)
        parts = []
        if n_missing:
            parts.append(f"{n_missing} fish missing")
        if n_stale:
            parts.append(f"{n_stale} stale")
        print(f"  {c:<5} ({', '.join(parts)})\n         {cmd}")
    n_blocked = sum(1 for _, r in rows for c in COLS if r[c] == BLOCKED)
    if n_blocked:
        print(f"\n  ({n_blocked} cells are blocked (·): they unlock once the stage before them is done.)")
if notes:
    print("\nDetails:")
    for n in dict.fromkeys(notes):
        print(f"  ⚠ {n}")
print()

if args.strict and any(r[c] in (MISSING, STALE) for _, r in rows for c in COLS):
    sys.exit(1)
