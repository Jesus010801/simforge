"""Level 0 — timeline diagnostics, derived from the Phase 1 FrameTimeIndex.

Times are never re-parsed here; Phase 1 validation results are translated into
diagnostics with their quantitative evidence.
"""
from __future__ import annotations

import statistics

from analysis.campaign.diagnostics.context import DiagnosticContext
from analysis.campaign.diagnostics.registry import detector
from analysis.campaign.models import (
    Diagnostic, EvidenceConfidence, SamplingStrategy, Severity, ValidationState,
)


def _sampling(n: int) -> dict:
    return {"strategy": SamplingStrategy.METADATA_ONLY, "frames_examined": n, "total_frames": n}


@detector("timeline_integrity", "1", 0,
          "empty / invalid / truncated timelines (from the Phase 1 time index)",
          sampling=SamplingStrategy.METADATA_ONLY)
def timeline_integrity(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    n = ti.n_frames
    out: list[Diagnostic] = []
    codes = {w.code: w for w in ti.warnings}
    if ti.state == ValidationState.INVALID:
        code = "empty_trajectory" if "empty_trajectory" in codes or n == 0 else "time_index_invalid"
        out.append(Diagnostic(
            code=code, severity=Severity.ERROR, scope="timeline",
            message="the trajectory timeline cannot be established: "
                    + "; ".join(w.message for w in ti.warnings),
            measured={"n_frames": n, "index_warnings": sorted(codes),
                      "backend_reported_frames": ti.backend_reported_frames},
            sampling=_sampling(n)))
    if ti.truncated:
        out.append(Diagnostic(
            code="truncated_trajectory", severity=Severity.WARN, scope="timeline",
            message=f"trajectory ends with an incomplete frame; {n} complete frame(s) usable",
            frames=[n - 1, n - 1] if n else None,
            times_ps=[ti.end_time_ps, ti.end_time_ps] if n else None,
            measured={"n_complete_frames": n, "evidence": ti.truncation_evidence},
            sampling=_sampling(n)))
    return out


@detector("duplicate_timestamps", "1", 0,
          "distinct frames sharing one simulation time", sampling=SamplingStrategy.METADATA_ONLY)
def duplicate_timestamps(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    out: list[Diagnostic] = []
    lim = ctx.params.max_events_reported
    for group in ti.duplicate_groups:
        t = ti.times_ps[group[0]]
        steps = [ti.steps[i] for i in group] if ti.steps else []
        if steps and len(set(steps)) == 1:
            interp, sev, conf = "duplicated_state_same_md_step", Severity.WARN, EvidenceConfidence.STRONG
            msg = (f"frames {group} repeat time {t} ps with the same MD step "
                   f"({steps[0]}) — duplicated frames (e.g. overlapping concatenation/restart)")
        elif steps:
            interp, sev, conf = "time_precision_collapse_or_retimed", Severity.REVIEW, EvidenceConfidence.MODERATE
            msg = (f"frames {group} share time {t} ps but carry different MD steps {steps} — "
                   f"stored time precision may be insufficient, or the timeline was re-timed")
        else:
            interp, sev, conf = "duplicate_time_unknown_cause", Severity.REVIEW, EvidenceConfidence.WEAK
            msg = f"frames {group} share time {t} ps (no MD-step evidence)"
        out.append(Diagnostic(
            code="duplicate_timestamps", severity=sev, scope="timeline", message=msg,
            frames=[group[0], group[-1]], times_ps=[t, t],
            measured={"time_ps": t, "frames": group[:lim], "md_steps": steps[:lim],
                      "n_frames_sharing": len(group)},
            interpretation=interp, confidence=conf,
            candidate_remediations=["segment_selection", "deduplicate_frames"],
            sampling=_sampling(ti.n_frames)))
    return out


@detector("non_monotonic_time", "1", 0,
          "simulation time decreasing between consecutive frames",
          sampling=SamplingStrategy.METADATA_ONLY)
def non_monotonic_time(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    pairs = ti.non_monotonic_steps
    if not pairs:
        return []
    lim = ctx.params.max_events_reported
    trans = [{"frames": [a, b], "times_ps": [ti.times_ps[a], ti.times_ps[b]],
              "md_steps": [ti.steps[a], ti.steps[b]] if ti.steps else None}
             for a, b in pairs[:lim]]
    return [Diagnostic(
        code="non_monotonic_time", severity=Severity.REVIEW, scope="timeline",
        message=f"time decreases at {len(pairs)} transition(s) (first: frame {pairs[0][0]} "
                f"{ti.times_ps[pairs[0][0]]} ps → frame {pairs[0][1]} {ti.times_ps[pairs[0][1]]} ps); "
                f"time→frame mapping is ambiguous",
        frames=[pairs[0][0], pairs[-1][1]],
        times_ps=[ti.times_ps[pairs[0][0]], ti.times_ps[pairs[-1][1]]],
        measured={"n_transitions": len(pairs), "transitions": trans,
                  "transitions_truncated": len(pairs) > lim},
        interpretation="concatenated_or_reordered_segments", confidence=EvidenceConfidence.MODERATE,
        candidate_remediations=["segment_selection", "segment_ordering"],
        sampling=_sampling(ti.n_frames))]


@detector("time_gaps", "1", 0,
          "unusually large gaps and irregular output intervals",
          sampling=SamplingStrategy.METADATA_ONLY)
def time_gaps(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    t = ti.times_ps
    dts = [(i, t[i + 1] - t[i]) for i in range(len(t) - 1) if t[i + 1] > t[i]]
    if len(dts) < 2:
        return []
    p = ctx.params
    typical = statistics.median(dt for _i, dt in dts)
    out: list[Diagnostic] = []
    gaps = [(i, dt) for i, dt in dts if dt >= p.gap_ratio_threshold * typical]
    if gaps:
        ev = [{"frames": [i, i + 1], "times_ps": [t[i], t[i + 1]], "gap_ps": dt,
               "ratio_to_typical": dt / typical} for i, dt in gaps[:p.max_events_reported]]
        out.append(Diagnostic(
            code="time_gap", severity=Severity.WARN, scope="timeline",
            message=f"{len(gaps)} interval(s) ≥ {p.gap_ratio_threshold}× the typical "
                    f"{typical} ps spacing (largest {max(dt for _i, dt in gaps)} ps); "
                    f"frames may be missing — or output frequency changed",
            frames=[gaps[0][0], gaps[-1][0] + 1],
            times_ps=[t[gaps[0][0]], t[gaps[-1][0] + 1]],
            measured={"typical_dt_ps": typical, "n_gaps": len(gaps), "gaps": ev,
                      "gaps_truncated": len(gaps) > p.max_events_reported},
            interpretation="missing_frames_or_output_change", confidence=EvidenceConfidence.WEAK,
            uncertainty=["a gap may be a legitimate output-frequency change"],
            sampling=_sampling(ti.n_frames)))
    tol = p.sampling_rel_tol * typical
    irregular = [dt for _i, dt in dts if abs(dt - typical) > tol]
    if irregular:
        hist: dict[float, int] = {}
        for _i, dt in dts:
            key = round(dt, 6)
            hist[key] = hist.get(key, 0) + 1
        top = sorted(hist.items(), key=lambda kv: -kv[1])[:5]
        out.append(Diagnostic(
            code="irregular_sampling", severity=Severity.INFO, scope="timeline",
            message=f"{len(irregular)} of {len(dts)} forward intervals differ from the "
                    f"typical {typical} ps",
            measured={"typical_dt_ps": typical, "n_irregular": len(irregular),
                      "n_intervals": len(dts),
                      "interval_histogram_ps": [{"dt_ps": k, "count": v} for k, v in top]},
            interpretation="variable_output_interval", confidence=EvidenceConfidence.STRONG,
            sampling=_sampling(ti.n_frames)))
    return out
