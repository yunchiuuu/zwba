"""
epochs.py
=========
ONE definition of the analysis epochs, derived from the experiment config.

Rule (all effect-size summaries: ΔZ, d′, responder sign, figures)
-----------------------------------------------------------------
Discard the first EPOCH_DISCARD_MIN of each phase:

    baseline summary window = [baseline_start + discard, baseline_end)
    drug     summary window = [drug_start     + discard, drug_end)

with discard = EPOCH_DISCARD_MIN * 60 * volume_per_sec frames. For the 120-min CSN
design (45 min baseline, 45 min drug) and EPOCH_DISCARD_MIN = 20 that is
baseline min 20-45 (vols 1200-2700) and drug min 65-90 (vols 3900-5400).

Temporal-window analysis (run_temporal_intensity_map / epoch_dz / epoch_dprime)
-------------------------------------------------------------------------------
Seven 15-min windows with a 10-min step, starting at drug onset (NOT discarded: the
point of the windows is to follow the response over time), each compared with the
SAME baseline summary window above:

    W0 0-15, W1 10-25, W2 20-35, W3 30-45, W4 40-55, W5 50-65, W6 60-75 min
    after drug onset.

The GLM fit window is set by INCLUDED_BASELINE (minutes of baseline before drug onset) and
should equal baseline length - EPOCH_DISCARD_MIN, so the fit starts where the baseline
summary window starts; it still contains the drug onset and runs to the end of the recording.

Location:
    ~/zwba/chemogenetic/epochs.py
"""

from __future__ import annotations

from dataclasses import dataclass

from utils.run_guard import PipelineError, require_cfg

# ── temporal-window design (minutes; window starts are relative to drug onset) ──
WINDOW_WIDTH_MIN: int = 15
WINDOW_STEP_MIN:  int = 10
N_WINDOWS:        int = 7
WINDOW_LABELS: list[str] = [
    f"{k * WINDOW_STEP_MIN}–{k * WINDOW_STEP_MIN + WINDOW_WIDTH_MIN} min"
    for k in range(N_WINDOWS)
]

_REQUIRED_CFG = [
    "volume_per_sec", "baseline_start", "baseline_end",
    "drug_start", "drug_end", "wash_end", "EPOCH_DISCARD_MIN",
]


@dataclass(frozen=True)
class Epochs:
    sampling_rate_hz: float
    discard_frames: int
    baseline: tuple          # (start, end) frames, end exclusive — summary window
    drug: tuple              # (start, end) frames, end exclusive — summary window
    windows: tuple           # ((label, start, end), ...) temporal windows
    n_vols_required: int     # recording must have at least this many volumes

    @property
    def n_windows(self) -> int:
        return len(self.windows)

    def minutes(self, frames: int) -> float:
        return frames / (60.0 * self.sampling_rate_hz)

    def describe(self) -> str:
        b, d = self.baseline, self.drug
        return (f"baseline vols [{b[0]}, {b[1]}) = min {self.minutes(b[0]):.0f}-"
                f"{self.minutes(b[1]):.0f} | drug vols [{d[0]}, {d[1]}) = min "
                f"{self.minutes(d[0]):.0f}-{self.minutes(d[1]):.0f} | "
                f"{self.n_windows} temporal windows | "
                f"discard first {self.minutes(self.discard_frames):.0f} min of each phase")

    def check_recording(self, n_vols: int, name: str = "recording") -> None:
        """Raise if the recording is shorter than the epochs/windows need."""
        if n_vols < self.n_vols_required:
            raise PipelineError(
                f"{name} has {n_vols} volumes but the configured epochs need "
                f"{self.n_vols_required} (config wash_end / temporal windows). "
                f"Wrong config for this fish, truncated recording, or a failed "
                f"voluseg run?"
            )


def get_epochs(cfg) -> Epochs:
    """Build and validate the epochs for a config module."""
    require_cfg(cfg, _REQUIRED_CFG, getattr(cfg, "__name__", "config"))

    sr = float(cfg.volume_per_sec)
    discard_min = float(cfg.EPOCH_DISCARD_MIN)
    if discard_min < 0:
        raise PipelineError(f"EPOCH_DISCARD_MIN must be >= 0, got {discard_min}")
    discard = int(round(discard_min * 60 * sr))

    b = (int(cfg.baseline_start) + discard, int(cfg.baseline_end))
    d = (int(cfg.drug_start) + discard, int(cfg.drug_end))
    problems = []
    if b[1] <= b[0]:
        problems.append(f"baseline window is empty after discarding {discard_min:g} min "
                        f"(vols [{b[0]}, {b[1]})) — baseline phase is too short")
    if d[1] <= d[0]:
        problems.append(f"drug window is empty after discarding {discard_min:g} min "
                        f"(vols [{d[0]}, {d[1]})) — drug phase is too short")
    if int(cfg.baseline_end) > int(cfg.drug_start):
        problems.append(f"baseline_end ({cfg.baseline_end}) is after drug_start "
                        f"({cfg.drug_start})")

    w_step = int(round(WINDOW_STEP_MIN * 60 * sr))
    w_width = int(round(WINDOW_WIDTH_MIN * 60 * sr))
    windows = []
    for k in range(N_WINDOWS):
        s = int(cfg.drug_start) + k * w_step
        windows.append((WINDOW_LABELS[k], s, s + w_width))
    n_required = max(int(cfg.wash_end), windows[-1][2])
    if windows[-1][2] > int(cfg.wash_end):
        problems.append(
            f"the {N_WINDOWS} temporal windows end at vol {windows[-1][2]} but wash_end "
            f"is {cfg.wash_end}. The windows (7 x 15 min, 10 min step from drug onset) "
            f"are fixed in chemogenetic/epochs.py; change N_WINDOWS there for a "
            f"shorter recording."
        )
    if problems:
        raise PipelineError("Invalid epoch configuration:\n  - " + "\n  - ".join(problems))

    return Epochs(sampling_rate_hz=sr, discard_frames=discard, baseline=b, drug=d,
                  windows=tuple(windows), n_vols_required=n_required)
