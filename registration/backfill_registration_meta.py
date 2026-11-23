"""
backfill_registration_meta.py
=============================
Write a registration_meta.json for fish that were registered BEFORE the provenance
sidecar existed, so that submit_registration_syn.sh skips them and status.sh shows them
as done instead of "template unverified".

This asserts something the files alone cannot prove (that the old registration used the
current template), so it only stamps a fish when the evidence you give holds:

  1. the three registration outputs exist and there is no registration_meta.json yet;
  2. the warp file was written AT OR AFTER --registered-after (the time the current template
     was built; for zstack_v1 that is the mtime of the ORIGINAL template file);
  3. if a job log (logs/registration/<expt_ID>.log) says which template it loaded, that
     file must exist and have the same sha256 as the current template (a log that loaded a
     different file, or a file that no longer exists, blocks the fish).

Everything it writes is marked  "backfilled": true  with the evidence used.
Default is a DRY RUN; nothing is written without --write.

Usage (from the repo root, login node, stdlib only):
    python3 registration/backfill_registration_meta.py --config config_registration \\
        --registered-after "2026-07-11 19:01"                 # dry run: show the decision per fish
    python3 registration/backfill_registration_meta.py --config config_registration \\
        --registered-after "2026-07-11 19:01" --write
Options:
    --fish TEXT         only fish whose expt_ID contains TEXT (repeatable)
    --fish-config NAME  take the fish list from chemogenetic/config/NAME.py (its all_fish) instead of
                        all_fish in the registration config. Use this when the registration config only
                        lists the fish you are registering now and the old ones are commented out.
    --log-dir DIR       job logs to cross-check (default: logs/registration)
    --ignore-logs       do not cross-check logs
    --registered-after  "YYYY-MM-DD HH:MM" in the cluster's local time (same as `ls -l`)
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True, help="module name under registration/config/")
parser.add_argument("--registered-after", required=True,
                    help='"YYYY-MM-DD HH:MM": the time the current template was built')
parser.add_argument("--write", action="store_true", help="actually write the files (default: dry run)")
parser.add_argument("--fish", action="append", default=[], help="only fish whose expt_ID contains this text")
parser.add_argument("--fish-config", default=None,
                    help="module name under chemogenetic/config/ whose all_fish is used as the fish list")
parser.add_argument("--log-dir", default=str(REPO / "logs" / "registration"))
parser.add_argument("--ignore-logs", action="store_true")
args = parser.parse_args()

if args.config.endswith(".py") or "/" in args.config:
    sys.exit(f"--config must be a module name, not a path: {args.config}")
try:
    after = datetime.datetime.strptime(args.registered_after, "%Y-%m-%d %H:%M")
except ValueError:
    sys.exit('--registered-after must look like "2026-07-11 19:01"')

sys.path.insert(0, str(REPO / "registration" / "config"))
try:
    cfg = __import__(args.config)
except Exception as e:      # noqa: BLE001
    sys.exit(f"Cannot import registration/config/{args.config}.py: {type(e).__name__}: {e}")

template = Path(cfg.MEAN_BRAIN_PATH)
if not template.exists():
    sys.exit(f"Template not found: {template}")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


TEMPLATE_SHA = sha256(template)
tstat = template.stat()
_sha_cache: dict[str, str | None] = {str(template): TEMPLATE_SHA}


def sha_of(path: str):
    if path not in _sha_cache:
        _sha_cache[path] = sha256(path) if os.path.exists(path) else None
    return _sha_cache[path]


def template_in_log(expt_id: str):
    """(status, path): 'none' (no log / no line), or 'found' with the LAST template path the log printed."""
    lp = Path(args.log_dir) / f"{expt_id}.log"
    if args.ignore_logs or not lp.exists():
        return "none", None
    text = lp.read_bytes().decode("utf-8", errors="ignore")      # some logs contain binary bytes
    hits = re.findall(r"Loading template brain:\s*(\S+)", text)
    return ("found", hits[-1]) if hits else ("none", None)


NEED = ["expt_to_mean_warp.nii.gz", "expt_to_mean_affine.mat", "expt_to_mean_registered_final.nii.gz"]
all_fish = cfg.all_fish
fish_source = f"all_fish in {args.config}"
if args.fish_config:
    if args.fish_config.endswith(".py") or "/" in args.fish_config:
        sys.exit(f"--fish-config must be a module name, not a path: {args.fish_config}")
    sys.path.insert(0, str(REPO))
    try:
        import importlib
        all_fish = importlib.import_module(f"chemogenetic.config.{args.fish_config}").all_fish
    except Exception as e:      # noqa: BLE001
        sys.exit(f"Cannot import chemogenetic/config/{args.fish_config}.py: {type(e).__name__}: {e}")
    fish_source = f"all_fish in {args.fish_config}"
fish_list = [f for f in all_fish if not args.fish or any(t in f[1] for t in args.fish)]
if not fish_list:
    sys.exit("No fish match --fish.")

print(f"fish     : {len(fish_list)} from {fish_source}")
print(f"template : {template}  sha256 {TEMPLATE_SHA[:12]}…")
print(f"evidence : warp written at/after {after:%Y-%m-%d %H:%M}"
      + ("" if args.ignore_logs else f"; job log in {args.log_dir} (if any) loaded a file with this sha256"))
print(f"mode     : {'WRITE' if args.write else 'dry run (add --write to write)'}\n")

n_stamp = n_skip = n_block = 0
for proj, expt in fish_list:
    d = Path(str(cfg.dir_registration)) / proj / expt
    meta_path = d / "registration_meta.json"
    missing = [n for n in NEED if not (d / n).exists()]
    if missing:
        print(f"  not registered   {expt}  (missing {', '.join(missing)})"); n_block += 1; continue
    if meta_path.exists():
        print(f"  already has meta {expt}"); n_skip += 1; continue

    warp_mtime = datetime.datetime.fromtimestamp((d / "expt_to_mean_warp.nii.gz").stat().st_mtime)
    if warp_mtime < after:
        print(f"  BLOCKED          {expt}  warp {warp_mtime:%Y-%m-%d %H:%M} is older than the template build "
              f"({after:%Y-%m-%d %H:%M}): registered to an earlier template, re-register it")
        n_block += 1; continue

    log_status, log_path = template_in_log(expt)
    if log_status == "found":
        s = sha_of(log_path)
        if s is None:
            print(f"  BLOCKED          {expt}  its log loaded {log_path}, which no longer exists (cannot verify)")
            n_block += 1; continue
        if s != TEMPLATE_SHA:
            print(f"  BLOCKED          {expt}  its log loaded {log_path}, a DIFFERENT file than the current "
                  f"template (sha256 {s[:12]}…): re-register it")
            n_block += 1; continue
        log_note = f"log loaded {log_path} (same sha256)"
    else:
        log_note = "no log line to cross-check"

    meta = {
        "expt_ID": expt,
        "proj_ID": proj,
        "backfilled": True,
        "backfilled_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "backfill_evidence": {
            "warp_mtime": warp_mtime.isoformat(timespec="minutes"),
            "registered_after": after.isoformat(timespec="minutes"),
            "log": log_note,
            "note": "Registration predates registration_meta.json. Template identity was inferred from "
                    "timestamps and the job log, not recorded by the job.",
        },
        "config_module": args.config,
        "template": {
            "version": getattr(cfg, "TEMPLATE_VERSION", None),
            "path": str(template),
            "sha256": TEMPLATE_SHA,
            "size_bytes": tstat.st_size,
            "mtime_utc": datetime.datetime.fromtimestamp(
                tstat.st_mtime, datetime.timezone.utc).isoformat(timespec="seconds"),
        },
    }
    if args.write:
        meta_path.write_text(json.dumps(meta, indent=2))
        print(f"  WROTE            {expt}  (warp {warp_mtime:%Y-%m-%d %H:%M}; {log_note})")
    else:
        print(f"  would write      {expt}  (warp {warp_mtime:%Y-%m-%d %H:%M}; {log_note})")
    n_stamp += 1

verb = "written" if args.write else "would be written"
print(f"\n{n_stamp} {verb}, {n_skip} already had a meta file, {n_block} blocked or not registered.")
if n_block:
    print("Blocked fish keep showing ⚠ / get re-registered by submit_registration_syn.sh; nothing was changed for them.")
sys.exit(0)
