"""Level 1 — box diagnostics from the per-frame boxes of the Phase 1 time index."""
from __future__ import annotations

import numpy as np

from analysis.campaign.diagnostics.context import DiagnosticContext, box_matrix
from analysis.campaign.diagnostics.registry import detector
from analysis.campaign.models import (
    Diagnostic, EvidenceConfidence, SamplingStrategy, Severity, ValidationState,
)

_EDGE_NAMES = ("a", "b", "c")


def _has_timeline(ctx: DiagnosticContext):
    ti = ctx.get_time_index()
    if ti.state == ValidationState.INVALID or ti.n_frames == 0:
        return False, "no valid time index"
    if not ti.boxes:
        return False, "time index carries no box records"
    return True, ""


def _sampling(n: int) -> dict:
    return {"strategy": SamplingStrategy.ALL_FRAMES, "frames_examined": n, "total_frames": n}


@detector("box_behaviour", "1", 1,
          "missing boxes, abrupt box discontinuities and cumulative box-shape drift",
          applicable=_has_timeline, sampling=SamplingStrategy.ALL_FRAMES)
def box_behaviour(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    raw = ctx.boxes()
    n = len(raw)
    p = ctx.params
    mats = np.stack([box_matrix(b) for b in raw])
    vol = np.abs(np.linalg.det(mats))
    edges = np.linalg.norm(mats, axis=2)                    # (n, 3)
    offdiag = raw[:, 3:]
    out: list[Diagnostic] = []

    missing = np.where(vol <= 1e-9)[0]
    if missing.size:
        out.append(Diagnostic(
            code="box_missing", severity=Severity.ERROR, scope="box",
            message=f"{missing.size} frame(s) have no (zero-volume) box — periodic "
                    f"reasoning is impossible there",
            frames=[int(missing[0]), int(missing[-1])],
            times_ps=[ti.times_ps[int(missing[0])], ti.times_ps[int(missing[-1])]],
            measured={"n_frames_without_box": int(missing.size)},
            sampling=_sampling(n)))
        return out
    if n < 2:
        return out

    dvol = np.abs(np.diff(vol)) / vol[:-1]
    dedge = np.abs(np.diff(edges, axis=0)) / edges[:-1]     # (n-1, 3)
    doff = np.abs(np.diff(offdiag, axis=0)).max(axis=1)
    med_vol = float(np.median(dvol))
    med_edge = np.median(dedge, axis=0)
    events = []
    for i in range(n - 1):
        reasons = []
        if dvol[i] >= p.box_volume_change_fraction and \
                dvol[i] >= p.box_change_vs_typical_ratio * max(med_vol, 1e-12):
            reasons.append("volume")
        for k in range(3):
            if dedge[i, k] >= p.box_edge_change_fraction and \
                    dedge[i, k] >= p.box_change_vs_typical_ratio * max(med_edge[k], 1e-12):
                reasons.append(f"edge_{_EDGE_NAMES[k]}")
        if doff[i] >= p.box_offdiag_change_nm:
            reasons.append("shape")
        if reasons:
            events.append({
                "frames": [i, i + 1], "times_ps": [ti.times_ps[i], ti.times_ps[i + 1]],
                "volume_nm3": [float(vol[i]), float(vol[i + 1])],
                "relative_volume_change": float((vol[i + 1] - vol[i]) / vol[i]),
                "edges_nm": [edges[i].round(5).tolist(), edges[i + 1].round(5).tolist()],
                "max_offdiagonal_change_nm": float(doff[i]),
                "triggered_by": reasons})
    if events:
        lim = p.max_events_reported
        e0 = events[0]
        out.append(Diagnostic(
            code="box_discontinuity", severity=Severity.REVIEW, scope="box",
            message=(f"{len(events)} abrupt box change(s); first frame {e0['frames'][0]}→"
                     f"{e0['frames'][1]}: volume {e0['volume_nm3'][0]:.1f}→"
                     f"{e0['volume_nm3'][1]:.1f} nm³ "
                     f"({100 * e0['relative_volume_change']:+.1f}%)"),
            frames=[events[0]["frames"][0], events[-1]["frames"][1]],
            times_ps=[events[0]["times_ps"][0], events[-1]["times_ps"][1]],
            measured={"n_events": len(events), "events": events[:lim],
                      "events_truncated": len(events) > lim,
                      "typical_relative_volume_change": med_vol,
                      "typical_relative_edge_change": med_edge.tolist()},
            interpretation="abrupt_box_change", confidence=EvidenceConfidence.MODERATE,
            uncertainty=["may be a legitimate change (e.g. restart with a different box / "
                         "ensemble) — not assumed to be corruption"],
            candidate_remediations=["segment_selection"],
            sampling=_sampling(n)))

    drift = np.abs(edges[-1] - edges[0]) / edges[0]
    if (drift >= p.box_drift_fraction).any():
        out.append(Diagnostic(
            code="box_shape_drift", severity=Severity.INFO, scope="box",
            message=("box edges change cumulatively over the trajectory: "
                     + ", ".join(f"{_EDGE_NAMES[k]} {edges[0, k]:.3f}→{edges[-1, k]:.3f} nm"
                                 for k in range(3))),
            frames=[0, n - 1], times_ps=[ti.times_ps[0], ti.times_ps[-1]],
            measured={"edges_first_nm": edges[0].tolist(), "edges_last_nm": edges[-1].tolist(),
                      "relative_edge_change": drift.tolist(),
                      "volume_first_nm3": float(vol[0]), "volume_last_nm3": float(vol[-1])},
            interpretation="progressive_box_deformation", confidence=EvidenceConfidence.STRONG,
            uncertainty=["typical of anisotropic/semi-isotropic pressure coupling; "
                         "coordinates are rescaled with the box"],
            sampling=_sampling(n)))
    return out
