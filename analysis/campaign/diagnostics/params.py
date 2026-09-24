"""Named, documented detector thresholds.

These are conservative *screening* defaults, not physical constants.  Every
value is recorded in each diagnostic's provenance; later profiles/policy may
override them.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class DiagnosticParams:
    # ── timeline ───────────────────────────────────────────────────────────
    #: a positive frame interval >= this × the median interval is reported as a gap
    gap_ratio_threshold: float = 5.0
    #: intervals within this relative tolerance of the median count as "regular"
    sampling_rel_tol: float = 1e-3

    # ── box ────────────────────────────────────────────────────────────────
    #: a frame-to-frame change must exceed BOTH an absolute fraction ...
    box_volume_change_fraction: float = 0.05
    box_edge_change_fraction: float = 0.10
    box_offdiag_change_nm: float = 0.10
    #: ... AND this multiple of the trajectory's own median change, so steady
    #: barostat deformation (e.g. semi-isotropic membrane boxes) is not flagged
    box_change_vs_typical_ratio: float = 10.0
    #: cumulative first→last edge change reported (INFO) as box-shape drift
    box_drift_fraction: float = 0.10

    # ── component motion ───────────────────────────────────────────────────
    #: a displacement is "large" when some fractional (box) coordinate moves by >= this
    jump_box_fraction: float = 0.30
    #: a large displacement is a *likely periodic wrap* only if the minimum-image
    #: displacement is <= wrap_residual_nm AND <= wrap_residual_ratio × raw displacement
    wrap_residual_nm: float = 1.0
    wrap_residual_ratio: float = 0.25
    #: raw vs minimum-image evidence "strong" if residual <= this × raw displacement
    strong_wrap_ratio: float = 0.10

    # ── partner ↔ receptor ─────────────────────────────────────────────────
    #: |Δ minimum-image separation| <= this while the raw separation jumps => periodic evidence
    separation_stable_nm: float = 0.5
    #: raw − minimum-image separation >= this fraction of the shortest box edge
    #: => partner stored in a non-nearest periodic image
    image_discrepancy_box_fraction: float = 0.30

    # ── reporting ──────────────────────────────────────────────────────────
    max_events_reported: int = 25

    def to_dict(self) -> dict:
        return asdict(self)
