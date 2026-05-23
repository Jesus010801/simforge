"""StudyAnalyzer — comparative multi-system MD study analysis.

Entry point: parse_study(path) → Study

Discovers all XVG files under ``path``, extracts system/replica/observable
from filenames, computes per-series statistics, aggregates across replicas,
detects outliers, and generates comparative findings.

Supported naming convention (auto-detected):
  SYSTEM-REPLICAobservable.xvg       AA-A1rmsd_protein.xvg
  SYSTEM-REPLICA_observable.xvg      LP-A4_rmsd-ligand.xvg
  SYSTEM-REPLICA-observable.xvg      HMG-A5-mindist_lig.xvg

Fallback: if no files match the study pattern, the entire directory is treated
as a single system (directory name) with one replica ("R1").
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from runtime.study_models import (
    Study, SystemGroup, Replica, ObservableSeries, AggregateMetrics,
    ComparativeSummary, ComparativeFinding,
)
from runtime.observable_resolver import ObservableResolver, ResolvedObservable
from runtime.xvg_parser import parse_xvg, XVGData


# ─── Pure-Python statistics helpers ──────────────────────────────────────────

def _mean(vals: list[float]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def _std(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5


def _median(vals: list[float]) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def _linear_slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    sx = sum(xs); sy = sum(ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    sxx = sum(x * x for x in xs)
    denom = n * sxx - sx * sx
    return (n * sxy - sx * sy) / denom if abs(denom) > 1e-12 else 0.0


# ─── Filename pattern recognition ─────────────────────────────────────────────
# Matches: SYSTEM sep REPLICA [sep?] observable_hint
#   where SYSTEM  = 1–6 alpha chars
#         sep     = _ or -
#         REPLICA = one letter + one or more digits  (A1, A10, B3, …)
#         rest    = anything

_STUDY_RE = re.compile(r'^([A-Za-z]{1,6})[_-]([A-Za-z]\d+)[_\-]?(.+)$')


def _parse_filename(stem: str) -> Optional[tuple[str, str, str]]:
    """Return (system, replica, obs_hint) or None if pattern doesn't match."""
    m = _STUDY_RE.match(stem)
    if not m:
        return None
    system   = m.group(1).upper()
    replica  = m.group(2).upper()
    obs_hint = m.group(3).lstrip('_-')
    return system, replica, obs_hint


# ─── Per-series statistics ────────────────────────────────────────────────────

def _compute_series_stats(series: ObservableSeries) -> None:
    """Fill in mean, std, median, drift, plateau_std, convergence_score."""
    vals  = series.values
    times = series.time_ns
    n = len(vals)

    if not vals:
        return

    series.mean   = _mean(vals)
    series.std    = _std(vals)
    series.median = _median(vals)

    if len(times) >= 2:
        series.drift = _linear_slope(times, vals)

    # Last-20% window
    start      = max(1, int(n * 0.8))
    last_vals  = vals[start:]
    last_times = times[start:]

    series.plateau_std = _std(last_vals)

    if len(last_vals) >= 3:
        last_mean  = _mean(last_vals)
        last_std   = _std(last_vals)
        last_drift = _linear_slope(last_times, last_vals) if len(last_times) >= 2 else 0.0

        ref        = abs(last_mean) if abs(last_mean) > 1e-6 else max(last_std, 1e-6)
        rel_std    = last_std / ref
        time_span  = (last_times[-1] - last_times[0]) if len(last_times) >= 2 else 1.0
        rel_drift  = abs(last_drift * time_span) / ref

        std_penalty   = min(1.0, rel_std   / 0.20)
        drift_penalty = min(1.0, rel_drift / 0.10)
        series.convergence_score = max(0.0, 1.0 - 0.6 * std_penalty - 0.4 * drift_penalty)
    else:
        series.convergence_score = 0.5


# ─── Aggregate statistics ────────────────────────────────────────────────────

def _compute_aggregate(sys_group: SystemGroup, observable_units: dict[str, str]) -> None:
    for obs_name in sys_group.observables:
        replica_means:  list[float] = []
        conv_scores:    list[float] = []

        for replica in sys_group.replicas.values():
            if obs_name in replica.observables:
                s = replica.observables[obs_name]
                replica_means.append(s.mean)
                conv_scores.append(s.convergence_score)

        if not replica_means:
            continue

        sys_group.aggregate[obs_name] = AggregateMetrics(
            observable       = obs_name,
            system           = sys_group.name,
            n_replicas       = len(replica_means),
            mean             = _mean(replica_means),
            std              = _std(replica_means),
            median           = _median(replica_means),
            min_val          = min(replica_means),
            max_val          = max(replica_means),
            mean_convergence = _mean(conv_scores) if conv_scores else 0.5,
            units            = observable_units.get(obs_name, ""),
        )


# ─── Outlier detection ────────────────────────────────────────────────────────

# Grubbs test critical values (α≈0.05) for small n.
# For n ≥ 9 we fall back to the classical 2σ rule.
_GRUBBS_CRITICAL: dict[int, float] = {
    3: 1.15, 4: 1.48, 5: 1.71, 6: 1.89, 7: 2.02, 8: 2.13,
}


def _outlier_threshold(n: int) -> float:
    return _GRUBBS_CRITICAL.get(n, 2.0)


def _detect_outliers(study: Study) -> None:
    """Tag replicas whose mean deviates beyond the Grubbs critical value from
    the system group mean (n-dependent threshold, more sensitive at small n)."""
    for sys_group in study.systems.values():
        for obs_name, agg in sys_group.aggregate.items():
            n = agg.n_replicas
            if n < 3:
                continue

            group_std = agg.std
            if group_std < 1e-10:
                continue

            threshold = _outlier_threshold(n)

            for replica in sys_group.replicas.values():
                if obs_name not in replica.observables:
                    continue
                replica_mean = replica.observables[obs_name].mean
                z = (replica_mean - agg.mean) / group_std
                rel_dev = (abs(replica_mean - agg.mean) / abs(agg.mean)
                           if abs(agg.mean) > 1e-6 else 0.0)
                # Require both statistical significance and a meaningful effect size
                if abs(z) > threshold and rel_dev > 0.10:
                    direction = "above" if z > 0 else "below"
                    display  = study.observable_display.get(obs_name, obs_name)
                    units    = study.observable_units.get(obs_name, "")
                    val_str  = f"{replica_mean:.3g}{' ' + units if units else ''}"
                    mean_str = f"{agg.mean:.3g}{' ' + units if units else ''}"
                    reason = (
                        f"{display}: {val_str} is {abs(z):.1f}σ {direction} "
                        f"system mean ({mean_str})"
                    )
                    replica.outlier_reasons.append(reason)
                    replica.is_outlier = True


# ─── Comparative findings ────────────────────────────────────────────────────

def _generate_findings(study: Study) -> ComparativeSummary:
    findings:         list[ComparativeFinding]   = []
    outlier_replicas: list[tuple[str, str, str]] = []

    # Collect outlier entries
    for sys_group in study.systems.values():
        for replica in sys_group.replicas.values():
            if replica.is_outlier:
                for reason in replica.outlier_reasons:
                    outlier_replicas.append((sys_group.name, replica.label, reason))

    # System stability ranking: mean of convergence scores across observables
    system_ranking: dict[str, float] = {}
    for sys_group in study.systems.values():
        scores = [agg.mean_convergence for agg in sys_group.aggregate.values()]
        if scores:
            system_ranking[sys_group.name] = _mean(scores)

    # Per-observable comparative findings
    for obs_name in study.observables_detected:
        sys_metrics: dict[str, AggregateMetrics] = {
            sn: sg.aggregate[obs_name]
            for sn, sg in study.systems.items()
            if obs_name in sg.aggregate
        }
        if len(sys_metrics) < 2:
            continue

        display  = study.observable_display.get(obs_name, obs_name)
        units    = study.observable_units.get(obs_name, "")
        unit_sfx = f" {units}" if units else ""

        all_means = [m.mean for m in sys_metrics.values()]
        all_stds  = [m.std  for m in sys_metrics.values()]
        mean_of_means = _mean(all_means)
        range_mean    = max(all_means) - min(all_means)

        # Only generate findings when there's meaningful between-system variation
        rel_range = range_mean / mean_of_means if abs(mean_of_means) > 1e-6 else 0.0

        if rel_range > 0.15:
            max_sys = max(sys_metrics.items(), key=lambda x: x[1].mean)
            min_sys = min(sys_metrics.items(), key=lambda x: x[1].mean)

            findings.append(ComparativeFinding(
                message=f"{max_sys[0]} shows the highest {display} mean "
                        f"({max_sys[1].mean:.3g}{unit_sfx} across {max_sys[1].n_replicas} replicas)",
                level="info",
                system=max_sys[0],
                observable=obs_name,
            ))
            findings.append(ComparativeFinding(
                message=f"{min_sys[0]} shows the lowest {display} mean "
                        f"({min_sys[1].mean:.3g}{unit_sfx} across {min_sys[1].n_replicas} replicas)",
                level="highlight",
                system=min_sys[0],
                observable=obs_name,
            ))

        # Most consistent replicas (lowest inter-replica std)
        if max(all_stds, default=0.0) > 1e-6:
            min_std_sys = min(sys_metrics.items(), key=lambda x: x[1].std)
            findings.append(ComparativeFinding(
                message=f"{min_std_sys[0]} shows the most consistent {display} across replicas "
                        f"(inter-replica std={min_std_sys[1].std:.3g}{unit_sfx})",
                level="highlight",
                system=min_std_sys[0],
                observable=obs_name,
            ))

    # Outlier findings
    for system, replica, reason in outlier_replicas:
        findings.append(ComparativeFinding(
            message=f"{system} replica {replica} may be a statistical outlier — {reason}",
            level="warning",
            system=system,
            replica=replica,
        ))

    return ComparativeSummary(
        findings=findings,
        outlier_replicas=outlier_replicas,
        system_ranking=system_ranking,
    )


# ─── Main entry point ────────────────────────────────────────────────────────

def parse_study(path: Path) -> Study:
    """Discover, parse, and analyze a multi-system MD study directory.

    Returns a populated Study object. Never raises — errors are captured in
    study.parse_errors and study.n_xvg_ungrouped.
    """
    resolver = ObservableResolver()
    study    = Study(path=path)

    xvg_files = sorted(path.rglob("*.xvg"))
    study.n_xvg_discovered = len(xvg_files)

    # (system, replica, resolved, data, xvg_path) tuples
    matched:  list[tuple[str, str, ResolvedObservable, XVGData, Path]] = []
    ungrouped: list[tuple[Path, XVGData]] = []

    for xvg_path in xvg_files:
        try:
            data = parse_xvg(xvg_path)
        except Exception as exc:
            study.parse_errors.append(f"{xvg_path.name}: {exc}")
            continue

        if not data.series or not data.time_ps:
            continue

        parsed = _parse_filename(xvg_path.stem)
        if parsed is None:
            ungrouped.append((xvg_path, data))
            continue

        system, replica, obs_hint = parsed
        resolved = resolver.resolve(obs_hint, data.title, data.ylabel)
        matched.append((system, replica, resolved, data, xvg_path))

    # Fallback: single-system mode when no files matched the study pattern
    if not matched and ungrouped:
        dir_name = re.sub(r'[^A-Z0-9]', '', path.name.upper())[:6] or "SIM"
        for xvg_path, data in ungrouped:
            resolved = resolver.resolve_from_path(xvg_path, data.title, data.ylabel)
            matched.append((dir_name, "R1", resolved, data, xvg_path))
    else:
        study.n_xvg_ungrouped = len(ungrouped)

    study.n_xvg_parsed = len(matched)

    # ── Build Study structure ─────────────────────────────────────────────────
    obs_all: set[str] = set()

    for system_name, replica_label, resolved, data, xvg_path in matched:
        times_ns = [t / 1000.0 for t in data.time_ps]
        values   = data.series[0].values[: len(times_ns)]

        series = ObservableSeries(
            observable = resolved.canonical,
            replica    = replica_label,
            system     = system_name,
            xvg_path   = xvg_path,
            time_ns    = times_ns,
            values     = values,
        )
        _compute_series_stats(series)

        if system_name not in study.systems:
            study.systems[system_name] = SystemGroup(name=system_name)
        sys_group = study.systems[system_name]

        if replica_label not in sys_group.replicas:
            sys_group.replicas[replica_label] = Replica(label=replica_label, system=system_name)
        replica_obj = sys_group.replicas[replica_label]

        # Keep first occurrence per (replica, observable)
        if resolved.canonical not in replica_obj.observables:
            replica_obj.observables[resolved.canonical] = series

        obs_all.add(resolved.canonical)
        study.observable_display[resolved.canonical] = resolved.display
        study.observable_units[resolved.canonical]   = resolved.units

    study.observables_detected = sorted(obs_all)

    # ── Aggregate + outliers + findings ───────────────────────────────────────
    for sys_group in study.systems.values():
        _compute_aggregate(sys_group, study.observable_units)

    _detect_outliers(study)
    study.summary = _generate_findings(study)

    return study
