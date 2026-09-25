"""Level 2 — component centre-motion diagnostics (raw vs minimum-image).

Positions are centres of the *stored* coordinates (``gmx trajectory -nopbc``),
so a component re-imaged by a box vector appears as a displacement ≈ one box
translation while its minimum-image displacement stays small.  That pattern is
reported as *evidence* of periodic wrapping; a displacement that stays large
under the minimum image is reported as uncertain large motion — never as
"nonphysical".
"""
from __future__ import annotations

import numpy as np

from analysis.campaign.diagnostics.context import (
    MOTION_ROLES, DiagnosticContext, box_matrix, box_valid, min_edge, minimum_image,
)
from analysis.campaign.diagnostics.registry import detector
from analysis.campaign.models import (
    Diagnostic, EvidenceConfidence, SamplingStrategy, Severity, ValidationState,
)

PARTNER_ROLES = ("ligand", "peptide", "cofactor")

_SPLIT_CAVEAT = ("centres are computed from stored (possibly wrapped) coordinates; a "
                 "component split across the boundary has a displaced centre; "
                 "molecule_periodic_image_change examines the molecules individually")


def _present_roles(ctx: DiagnosticContext, candidates) -> list[str]:
    return [r for r in candidates if ctx.has_role(r)]


def _timeline_ok(ctx: DiagnosticContext) -> tuple[bool, str]:
    ti = ctx.get_time_index()
    if ti.state == ValidationState.INVALID or ti.n_frames < 2:
        return False, "no valid time index with at least two frames"
    return True, ""


def _sampling(ctx: DiagnosticContext) -> dict:
    n = ctx.get_time_index().n_frames
    return {"strategy": SamplingStrategy.ALL_FRAMES, "frames_examined": n, "total_frames": n,
            "pairs": "consecutive frames in file order"}


def _boxes(ctx: DiagnosticContext, n: int):
    raw = ctx.boxes()
    if raw is None or len(raw) != n:
        return None
    return [box_matrix(b) for b in raw]


# ═══════════════════════════════════════════════════════════════════════════════
# Pure classification (unit-testable without gmx)
# ═══════════════════════════════════════════════════════════════════════════════

def classify_displacements(pos: np.ndarray, boxes, params) -> list[dict]:
    """Consecutive-frame centre displacements that are large in box units."""
    events = []
    for i in range(len(pos) - 1):
        B = boxes[i + 1] if boxes is not None else None
        d = pos[i + 1] - pos[i]
        raw = float(np.linalg.norm(d))
        if B is None or not box_valid(B):
            continue
        d_mi, f = minimum_image(d, B)
        frac = float(np.abs(f).max())
        if frac < params.jump_box_fraction:
            continue
        mi = float(np.linalg.norm(d_mi))
        wrap = mi <= params.wrap_residual_nm and mi <= params.wrap_residual_ratio * raw
        events.append({
            "frames": [i, i + 1], "raw_displacement_nm": d.round(4).tolist(),
            "raw_magnitude_nm": round(raw, 4),
            "min_image_displacement_nm": d_mi.round(4).tolist(),
            "min_image_magnitude_nm": round(mi, 4),
            "max_fractional_displacement": round(frac, 4),
            "box_translation": np.round(f).astype(int).tolist(),
            "box_edges_nm": np.linalg.norm(B, axis=1).round(4).tolist(),
            "classification": "likely_periodic_wrap" if wrap else "large_motion_uncertain",
            "strong": wrap and mi <= params.strong_wrap_ratio * raw,
        })
    return events


def separation_series(p_partner: np.ndarray, p_receptor: np.ndarray, boxes):
    """Per-frame raw and minimum-image partner–receptor centre distances."""
    raw = np.linalg.norm(p_partner - p_receptor, axis=1)
    mi = np.full(len(raw), np.nan)
    edge = np.full(len(raw), np.nan)
    if boxes is not None:
        for i, B in enumerate(boxes):
            if box_valid(B):
                d_mi, _f = minimum_image(p_partner[i] - p_receptor[i], B)
                mi[i] = np.linalg.norm(d_mi)
                edge[i] = min_edge(B)
    return raw, mi, edge


def classify_separation(raw, mi, edge, params) -> tuple[list[dict], list[int]]:
    """``(events, discrepancy_frames)`` for a partner–receptor pair."""
    events = []
    for i in range(len(raw) - 1):
        L = edge[i + 1]
        if not np.isfinite(L):
            continue
        d_raw = float(raw[i + 1] - raw[i])
        if abs(d_raw) < params.jump_box_fraction * L:
            continue
        d_mi = float(mi[i + 1] - mi[i])
        periodic = abs(d_mi) <= params.separation_stable_nm
        events.append({
            "frames": [i, i + 1],
            "raw_distance_nm": [round(float(raw[i]), 4), round(float(raw[i + 1]), 4)],
            "min_image_distance_nm": [round(float(mi[i]), 4), round(float(mi[i + 1]), 4)],
            "raw_change_nm": round(d_raw, 4), "min_image_change_nm": round(d_mi, 4),
            "shortest_box_edge_nm": round(float(L), 4),
            "classification": "periodic_separation" if periodic else "separation_change",
        })
    disc = [i for i in range(len(raw)) if np.isfinite(edge[i])
            and raw[i] - mi[i] >= params.image_discrepancy_box_fraction * edge[i]]
    return events, disc


def _ranges(frames: list[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for f in frames:
        if out and f == out[-1][1] + 1:
            out[-1][1] = f
        else:
            out.append([f, f])
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Detectors
# ═══════════════════════════════════════════════════════════════════════════════

def _motion_applicable(ctx):
    ok, why = _timeline_ok(ctx)
    if not ok:
        return ok, why
    roles = _present_roles(ctx, MOTION_ROLES)
    return (True, "") if roles else (False, "no semantic component groups available")


@detector("component_periodic_jump", "1", 2,
          "centre jumps of whole components / the System: periodic wrap vs large motion",
          applicable=_motion_applicable,
          roles=lambda ctx: _present_roles(ctx, MOTION_ROLES),
          sampling=SamplingStrategy.ALL_FRAMES)
def component_periodic_jump(ctx: DiagnosticContext) -> list[Diagnostic]:
    roles = _present_roles(ctx, MOTION_ROLES)
    series = ctx.motion(roles)
    ti = ctx.get_time_index()
    boxes = _boxes(ctx, ti.n_frames)
    p = ctx.params
    out: list[Diagnostic] = []
    for role in roles:
        s = series[role]
        events = classify_displacements(s.positions, boxes, p)
        for e in events:
            a, b = e["frames"]
            e["times_ps"] = [ti.times_ps[a], ti.times_ps[b]]
        base_prov = {"group": s.group, "group_evidence": s.group_evidence,
                     "centre_convention": s.convention, "mass_source": s.mass_source,
                     "coordinates": s.coordinates, "probe_cache_key": s.cache_key}
        wraps = [e for e in events if e["classification"] == "likely_periodic_wrap"]
        large = [e for e in events if e["classification"] == "large_motion_uncertain"]
        scope = "system" if role == "system" else f"component:{role}"
        if wraps:
            all_strong = all(e["strong"] for e in wraps)
            system = role == "system"
            out.append(Diagnostic(
                code="system_periodic_translation" if system else "component_periodic_wrap",
                severity=Severity.WARN, scope=scope,
                message=(f"{role} centre jumps by ≈ a box translation {len(wraps)} time(s) "
                         f"while its minimum-image displacement stays small "
                         f"(first: frames {wraps[0]['frames']}, raw {wraps[0]['raw_magnitude_nm']} nm, "
                         f"minimum-image {wraps[0]['min_image_magnitude_nm']} nm)"),
                frames=[wraps[0]["frames"][0], wraps[-1]["frames"][1]],
                times_ps=[wraps[0]["times_ps"][0], wraps[-1]["times_ps"][1]],
                measured={"n_events": len(wraps), "events": wraps[:p.max_events_reported],
                          "events_truncated": len(wraps) > p.max_events_reported},
                interpretation="likely_periodic_wrap",
                confidence=EvidenceConfidence.STRONG if all_strong else EvidenceConfidence.MODERATE,
                uncertainty=[] if system else [_SPLIT_CAVEAT],
                candidate_remediations=(["center", "nojump"] if system
                                        else ["nojump", "make_whole", "minimum_image_analysis"]),
                provenance=dict(base_prov), sampling=_sampling(ctx)))
        if large:
            out.append(Diagnostic(
                code="component_large_displacement", severity=Severity.REVIEW, scope=scope,
                message=(f"{role} centre moves ≥ {p.jump_box_fraction} of a box edge between "
                         f"frames {len(large)} time(s) and the minimum image does NOT explain it "
                         f"(first: frames {large[0]['frames']}, minimum-image "
                         f"{large[0]['min_image_magnitude_nm']} nm) — cause undetermined"),
                frames=[large[0]["frames"][0], large[-1]["frames"][1]],
                times_ps=[large[0]["times_ps"][0], large[-1]["times_ps"][1]],
                measured={"n_events": len(large), "events": large[:p.max_events_reported],
                          "events_truncated": len(large) > p.max_events_reported},
                interpretation="large_motion_uncertain", confidence=EvidenceConfidence.WEAK,
                uncertainty=[_SPLIT_CAVEAT, "may be genuine motion"],
                candidate_remediations=["make_whole", "minimum_image_analysis"],
                provenance=dict(base_prov), sampling=_sampling(ctx)))
    return out


def _partner_applicable(ctx):
    ok, why = _timeline_ok(ctx)
    if not ok:
        return ok, why
    if not ctx.has_role("receptor"):
        return False, "no resolved receptor group"
    if not _present_roles(ctx, PARTNER_ROLES):
        return False, "no resolved ligand / peptide / cofactor group"
    return True, ""


@detector("partner_receptor_separation", "1", 2,
          "partner (ligand/peptide/cofactor) vs receptor: periodic-image vs genuine separation",
          applicable=_partner_applicable,
          roles=lambda ctx: ["receptor"] + _present_roles(ctx, PARTNER_ROLES),
          sampling=SamplingStrategy.ALL_FRAMES)
def partner_receptor_separation(ctx: DiagnosticContext) -> list[Diagnostic]:
    partners = _present_roles(ctx, PARTNER_ROLES)
    series = ctx.motion(["receptor"] + partners)
    ti = ctx.get_time_index()
    boxes = _boxes(ctx, ti.n_frames)
    p = ctx.params
    rec = series["receptor"]
    out: list[Diagnostic] = []
    for role in partners:
        s = series[role]
        raw, mi, edge = separation_series(s.positions, rec.positions, boxes)
        events, disc = classify_separation(raw, mi, edge, p)
        for e in events:
            a, b = e["frames"]
            e["times_ps"] = [ti.times_ps[a], ti.times_ps[b]]
        prov = {"partner_group": s.group, "receptor_group": rec.group,
                "partner_group_evidence": s.group_evidence,
                "receptor_group_evidence": rec.group_evidence,
                "centre_convention": s.convention, "coordinates": s.coordinates}
        periodic = [e for e in events if e["classification"] == "periodic_separation"]
        real = [e for e in events if e["classification"] == "separation_change"]
        lim = p.max_events_reported
        if periodic:
            e0 = periodic[0]
            out.append(Diagnostic(
                code=f"{role}_receptor_periodic_separation", severity=Severity.WARN,
                scope=f"pair:{role}-receptor",
                message=(f"raw {role}–receptor centre distance jumps "
                         f"{e0['raw_distance_nm'][0]}→{e0['raw_distance_nm'][1]} nm at frames "
                         f"{e0['frames']} while the minimum-image distance stays "
                         f"{e0['min_image_distance_nm'][0]}→{e0['min_image_distance_nm'][1]} nm "
                         f"({len(periodic)} such event(s)) — periodic-image evidence"),
                frames=[periodic[0]["frames"][0], periodic[-1]["frames"][1]],
                times_ps=[periodic[0]["times_ps"][0], periodic[-1]["times_ps"][1]],
                measured={"n_events": len(periodic), "events": periodic[:lim],
                          "events_truncated": len(periodic) > lim},
                interpretation="periodic_image_separation", confidence=EvidenceConfidence.STRONG,
                uncertainty=[_SPLIT_CAVEAT],
                candidate_remediations=["complex_reconstruction", "minimum_image_analysis"],
                provenance=dict(prov), sampling=_sampling(ctx)))
        if real:
            e0 = real[0]
            out.append(Diagnostic(
                code=f"{role}_receptor_separation_change", severity=Severity.INFO,
                scope=f"pair:{role}-receptor",
                message=(f"{role}–receptor distance changes abruptly in BOTH raw and "
                         f"minimum-image geometry ({len(real)} event(s); first frames "
                         f"{e0['frames']}: minimum-image {e0['min_image_distance_nm'][0]}→"
                         f"{e0['min_image_distance_nm'][1]} nm) — possibly genuine separation, "
                         f"not diagnosed as periodic"),
                frames=[real[0]["frames"][0], real[-1]["frames"][1]],
                times_ps=[real[0]["times_ps"][0], real[-1]["times_ps"][1]],
                measured={"n_events": len(real), "events": real[:lim],
                          "events_truncated": len(real) > lim},
                interpretation="possible_physical_separation", confidence=EvidenceConfidence.WEAK,
                uncertainty=[_SPLIT_CAVEAT],
                provenance=dict(prov), sampling=_sampling(ctx)))
        if disc:
            rng = _ranges(disc)
            out.append(Diagnostic(
                code=f"{role}_receptor_image_discrepancy", severity=Severity.WARN,
                scope=f"pair:{role}-receptor",
                message=(f"in {len(disc)} frame(s) the stored {role} lies in a non-nearest "
                         f"periodic image of the receptor (raw − minimum-image distance ≥ "
                         f"{p.image_discrepancy_box_fraction} of the shortest box edge)"),
                frames=[disc[0], disc[-1]], times_ps=[ti.times_ps[disc[0]], ti.times_ps[disc[-1]]],
                measured={"n_frames": len(disc), "frame_ranges": rng[:lim],
                          "ranges_truncated": len(rng) > lim,
                          "first_frame_raw_nm": round(float(raw[disc[0]]), 4),
                          "first_frame_min_image_nm": round(float(mi[disc[0]]), 4)},
                interpretation="stored_in_non_nearest_image", confidence=EvidenceConfidence.MODERATE,
                uncertainty=[_SPLIT_CAVEAT],
                candidate_remediations=["complex_reconstruction", "minimum_image_analysis"],
                provenance=dict(prov), sampling=_sampling(ctx)))
    return out
