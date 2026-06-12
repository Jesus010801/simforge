"""Sprint 2 — Scientific summary generator.

Scans a SimForge workspace for XVG files produced by analysis steps,
runs convergence/stability analysis, and writes:
  metadata/scientific_summary.json
  metadata/scientific_summary.md

Best-effort: if no XVG files are found, returns an empty summary without error.

Phase A addition: analyze_trajectory(path) → (ScientificSummary, QualityReport)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, TYPE_CHECKING

from runtime.xvg_parser import parse_xvg, XVGData
from runtime.convergence_analyzer import (
    analyze_rmsd_convergence,
    analyze_energy_stability,
    ConvergenceResult,
    StabilityResult,
)

if TYPE_CHECKING:
    from runtime.quality_classifier import QualityReport


# ─────────────────────────────────────────────────────────────────────────────
# Result dataclass
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ScientificSummary:
    workspace:      Path
    converged:      bool
    rmsd_verdict:   str
    energy_verdict: str
    runtime_ns:     Optional[float]   # total simulation time (from XVG time axis)
    wall_time_s:    Optional[float]   # not yet extracted — placeholder
    analyses:       list[dict]        # per-analysis results
    warnings:       list[str]         = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "workspace":      str(self.workspace),
            "converged":      self.converged,
            "rmsd_verdict":   self.rmsd_verdict,
            "energy_verdict": self.energy_verdict,
            "runtime_ns":     self.runtime_ns,
            "wall_time_s":    self.wall_time_s,
            "n_analyses":     len(self.analyses),
            "analyses":       self.analyses,
            "warnings":       self.warnings,
        }

    def as_markdown(self) -> str:
        lines: list[str] = [
            "# Scientific Summary",
            "",
            f"**Workspace:** `{self.workspace}`",
            f"**Converged:** {'Yes' if self.converged else 'No'}",
        ]
        if self.runtime_ns is not None:
            lines.append(f"**Runtime:** {self.runtime_ns:.3f} ns")

        lines += ["", "## RMSD Convergence", "", self.rmsd_verdict or "_No RMSD data._"]
        lines += ["", "## Energy Stability",  "", self.energy_verdict or "_No energy data._"]

        if self.analyses:
            lines += ["", "## Per-Analysis Results", ""]
            for a in self.analyses:
                xvg  = a.get("xvg_file", "?")
                kind = a.get("kind", "?")
                verd = a.get("verdict", "")
                lines.append(f"### `{xvg}` ({kind})")
                lines.append(verd)
                lines.append("")

        if self.warnings:
            lines += ["", "## Warnings", ""]
            for w in self.warnings:
                lines.append(f"- {w}")

        return "\n".join(lines) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# XVG classification heuristics
# ─────────────────────────────────────────────────────────────────────────────

def _classify_xvg(path: Path, data: XVGData) -> str:
    """Return a coarse kind label: 'rmsd', 'energy', 'rmsf', or 'other'."""
    name  = path.name.lower()
    title = data.title.lower()
    ylabel = data.ylabel.lower()

    if "rmsd" in name or "rmsd" in title:
        return "rmsd"
    if "rmsf" in name or "rmsf" in title:
        return "rmsf"
    if any(k in name or k in title or k in ylabel
           for k in ("energy", "potential", "epot", "temperature", "pressure")):
        return "energy"
    return "other"


# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────

def generate_summary(workspace_path: Path) -> ScientificSummary:
    """Scan workspace for XVG files, analyze, return ScientificSummary.

    Best-effort: missing files / parse errors are captured as warnings.
    """
    analyses:  list[dict] = []
    warnings:  list[str]  = []

    rmsd_result:   Optional[ConvergenceResult] = None
    energy_result: Optional[StabilityResult]   = None
    runtime_ns:    Optional[float]             = None

    # Scan steps/ and analysis/ directories for .xvg files
    xvg_files: list[Path] = []
    for sub in ("steps", "analysis"):
        d = workspace_path / sub
        if d.exists():
            xvg_files.extend(sorted(d.rglob("*.xvg")))

    for xvg_path in xvg_files:
        try:
            data = parse_xvg(xvg_path)
        except Exception as exc:
            warnings.append(f"Failed to parse {xvg_path.name}: {exc}")
            continue

        kind = _classify_xvg(xvg_path, data)

        entry: dict = {
            "xvg_file":   xvg_path.name,
            "xvg_path":   str(xvg_path),
            "kind":       kind,
            "n_points":   len(data.time_ps),
            "n_series":   len(data.series),
            "title":      data.title,
            "verdict":    "",
        }

        # Track total simulation time from the last RMSD or any XVG time axis
        if data.time_ps and (runtime_ns is None or kind == "rmsd"):
            candidate_ns = data.times_ns()[-1]
            if runtime_ns is None or candidate_ns > runtime_ns:
                runtime_ns = candidate_ns

        if kind == "rmsd" and rmsd_result is None:
            try:
                res = analyze_rmsd_convergence(data)
                rmsd_result    = res
                entry["verdict"] = res.verdict
                entry["converged"] = res.converged
                entry["plateau_ns"] = res.plateau_ns
                entry["std_last20pct"] = res.std_last20pct
            except Exception as exc:
                warnings.append(f"RMSD analysis failed for {xvg_path.name}: {exc}")

        elif kind == "energy" and energy_result is None:
            try:
                res = analyze_energy_stability(data)
                energy_result  = res
                entry["verdict"]      = res.verdict
                entry["stable"]       = res.stable
                entry["drift_per_ns"] = res.drift_per_ns
                entry["mean"]         = res.mean
            except Exception as exc:
                warnings.append(f"Energy analysis failed for {xvg_path.name}: {exc}")

        analyses.append(entry)

    # Aggregate convergence decision
    converged = False
    if rmsd_result is not None:
        converged = rmsd_result.converged
    elif energy_result is not None:
        converged = energy_result.stable
    # If no analyses at all → converged stays False (best-effort)

    return ScientificSummary(
        workspace      = workspace_path,
        converged      = converged,
        rmsd_verdict   = rmsd_result.verdict   if rmsd_result   else "",
        energy_verdict = energy_result.verdict if energy_result else "",
        runtime_ns     = runtime_ns,
        wall_time_s    = None,
        analyses       = analyses,
        warnings       = warnings,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Phase A: analyze_trajectory
# ─────────────────────────────────────────────────────────────────────────────

def analyze_trajectory(
    path: Path,
    context: "str | None" = None,
) -> "tuple[ScientificSummary, QualityReport]":
    """Discover, parse, and classify simulation quality from any directory.

    1. Calls discover_trajectory(path) to find all MD output files.
    2. Calls load_xvg_files(manifest) to parse them.
    3. Calls classify_run() to produce a QualityReport (with optional context).
    4. Calls generate_summary() for a ScientificSummary (workspace-aware).

    Returns (ScientificSummary, QualityReport). Never raises — degrades
    gracefully on missing or unreadable data.
    """
    from runtime.trajectory_ingestor import discover_trajectory, load_xvg_files
    from runtime.quality_classifier import classify_run

    manifest = discover_trajectory(path)
    xvg_data = load_xvg_files(manifest)

    rmsd_data   = xvg_data.get("rmsd")
    energy_data = xvg_data.get("potential_energy")
    pressure    = xvg_data.get("pressure")
    temperature = xvg_data.get("temperature")

    quality_report = classify_run(
        rmsd_data        = rmsd_data,
        energy_data      = energy_data,
        pressure_data    = pressure,
        temperature_data = temperature,
        context          = context,
    )

    # Build ScientificSummary using the workspace-aware generator
    try:
        summary = generate_summary(path)
    except Exception:
        summary = ScientificSummary(
            workspace      = path,
            converged      = False,
            rmsd_verdict   = "",
            energy_verdict = "",
            runtime_ns     = None,
            wall_time_s    = None,
            analyses       = [],
            warnings       = ["Summary generation failed; using empty summary."],
        )

    return summary, quality_report
