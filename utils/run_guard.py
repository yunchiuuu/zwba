"""
run_guard.py
============
Failure railings shared by every run script.

Why this exists
---------------
Run scripts loop over many fish. Historically a per-fish exception was printed as
one line, the job still exited 0 and printed "complete", so a broken fish looked
like a finished one in `sacct`. This module makes failures loud and keeps them
visible until the very end of the log:

    guard = RunGuard("decompose")
    with guard.attempt(expt_ID):          # catches, prints traceback, records
        ...
    guard.record_results("glm refit", results)   # joblib workers returning dicts
    guard.finish()                        # summary + exit code 1 if anything failed

    @guarded                              # wraps a joblib worker taking `fish` first;
    def _worker(fish): ...                # exceptions become an ERROR result dict

    require_cfg(cfg, ["EPOCH_DISCARD_MIN", ...], "config_x")   # all missing names at once
    require_paths([(path, "hint"), ...], "what")               # all missing files at once
    run_main(main)                        # PipelineError -> clean message, exit code 2
    setup_or_exit(get_epochs, cfg)        # same, for module-level setup before main()

PipelineError means "the inputs/config are wrong, here is what to do" (no traceback
needed). Anything else is a bug or an environment problem and gets a traceback.
"""

from __future__ import annotations

import functools
import sys
import time
import traceback
from contextlib import contextmanager
from pathlib import Path


class PipelineError(RuntimeError):
    """Input or config problem with an actionable message."""


def _say(msg: str = "") -> None:
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# up-front validation (fail BEFORE heavy work, list ALL problems at once)
# ---------------------------------------------------------------------------

def require_cfg(cfg, names, config_name: str | None = None) -> None:
    """Raise PipelineError listing every missing config attribute."""
    missing = [n for n in names if not hasattr(cfg, n)]
    if missing:
        where = config_name or getattr(cfg, "__name__", "config")
        raise PipelineError(
            f"{where} is missing required setting(s): {', '.join(missing)}.\n"
            f"  Copy them from chemogenetic/config/config_template.py."
        )


def require_paths(items, what: str) -> None:
    """
    items: iterable of (path, hint) or just path.
    Raise PipelineError listing every path that does not exist.
    """
    missing = []
    for it in items:
        path, hint = (it if isinstance(it, (tuple, list)) else (it, ""))
        if not Path(path).exists():
            missing.append(f"    {path}" + (f"   <- {hint}" if hint else ""))
    if missing:
        raise PipelineError(f"{what}: {len(missing)} required file(s) not found:\n"
                            + "\n".join(missing))


# ---------------------------------------------------------------------------
# worker wrapper for joblib (loky) pools
# ---------------------------------------------------------------------------

def guarded(fn):
    """
    Wrap a worker whose first argument is `fish`. Any exception is returned as
    {"fish", "status": "ERROR: <Type>: <msg>", "traceback": str} instead of being
    raised, so one bad fish cannot kill the pool and cannot be silently lost.
    Pass the returned list to RunGuard.record_results().
    """
    @functools.wraps(fn)
    def wrapper(fish, *args, **kwargs):
        try:
            return fn(fish, *args, **kwargs)
        except PipelineError as e:
            return {"fish": fish, "status": f"ERROR: {e}", "traceback": None}
        except Exception as e:   # noqa: BLE001 — deliberate: report, don't swallow
            return {"fish": fish,
                    "status": f"ERROR: {type(e).__name__}: {e}",
                    "traceback": traceback.format_exc()}
    return wrapper


# ---------------------------------------------------------------------------
# failure tracker
# ---------------------------------------------------------------------------

class RunGuard:
    def __init__(self, name: str):
        self.name = name
        self.n_ok = 0
        self.failures: list[tuple[str, str]] = []
        self.warnings: list[tuple[str, str]] = []
        self._t0 = time.time()

    # -- recording --------------------------------------------------------
    def ok(self, n: int = 1) -> None:
        self.n_ok += n

    def fail(self, label: str, message: str, tb: str | None = None) -> None:
        self.failures.append((label, message))
        _say(f"\n❌ [{self.name}] {label}: {message}")
        if tb:
            _say(tb.rstrip())

    def warn(self, label: str, message: str) -> None:
        self.warnings.append((label, message))
        _say(f"⚠️  [{self.name}] {label}: {message}")

    @contextmanager
    def attempt(self, label: str):
        """Catch, report with traceback, record — and let the caller continue."""
        try:
            yield
        except PipelineError as e:
            self.fail(label, str(e))
        except Exception as e:   # noqa: BLE001
            self.fail(label, f"{type(e).__name__}: {e}", traceback.format_exc())
        else:
            self.n_ok += 1

    def record_results(self, stage: str, results, label=lambda r: r["fish"][1]) -> None:
        """Register joblib worker result dicts; status starting 'ERROR' = failure."""
        for r in results:
            status = str(r.get("status", ""))
            if status.startswith("ERROR"):
                self.fail(f"{stage}: {label(r)}", status[len("ERROR: "):], r.get("traceback"))
            else:
                self.n_ok += 1

    # -- end of run -------------------------------------------------------
    def finish(self) -> None:
        """Print a summary; exit non-zero if any failure was recorded."""
        dt = time.time() - self._t0
        _say("\n" + "=" * 64)
        if self.warnings:
            _say(f"[{self.name}] {len(self.warnings)} warning(s):")
            for label, msg in self.warnings:
                _say(f"   ⚠️  {label}: {msg}")
        if self.failures:
            _say(f"[{self.name}] FINISHED WITH FAILURES: "
                 f"{self.n_ok} ok, {len(self.failures)} failed  ({dt/60:.1f} min)")
            for label, msg in self.failures:
                _say(f"   ❌ {label}: {msg.splitlines()[0] if msg else ''}")
            _say("   Outputs for the failed items are missing or stale. "
                 "Exiting with code 1 so SLURM marks this job FAILED.")
            _say("=" * 64)
            sys.exit(1)
        _say(f"[{self.name}] ✅ complete: {self.n_ok} ok, 0 failed  ({dt/60:.1f} min)")
        _say("=" * 64)


def setup_or_exit(fn, *args, **kwargs):
    """Run module-level setup (get_epochs, build_regressor, ...); PipelineError -> exit 2."""
    try:
        return fn(*args, **kwargs)
    except PipelineError as e:
        _say(f"\n❌ {e}\n")
        sys.exit(2)


def run_main(main) -> None:
    """Entry-point wrapper: PipelineError -> readable message + exit code 2."""
    try:
        main()
    except PipelineError as e:
        _say(f"\n❌ {e}\n")
        sys.exit(2)
