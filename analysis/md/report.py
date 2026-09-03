"""Human-readable Markdown report for dry-run analysis planning."""
from __future__ import annotations

from pathlib import Path

from analysis.md.models import DiscoveredMDRun
from analysis.md.observables import AnalysisPlan
from analysis.md.provenance import Provenance

_STATUS_ICONS = {
    "planned":          "✓ planned",
    "unavailable":      "✗ unavailable",
    "needs_selection":  "⚠ needs selection",
    "already_available":"✓ already available",
}


def _display_path(run: DiscoveredMDRun, role: str, path: object) -> str | None:
    """Return a display-friendly path string for the report."""
    if path is None:
        return None
    if role in run.user_overrides:
        return run.user_overrides[role]
    p = path if isinstance(path, Path) else Path(str(path))
    try:
        return str(p.relative_to(run.run_dir))
    except ValueError:
        return str(p)


def write_report(
    run: DiscoveredMDRun,
    plan: AnalysisPlan,
    provenance: Provenance,
    output_dir: Path,
) -> None:
    lines: list[str] = [
        "# SimForge MD Analysis — Dry Run Report",
        "",
        f"**Run directory:** `{run.run_dir}`  ",
        f"**Generated:** {provenance.timestamp}  ",
        f"**SimForge version:** {provenance.simforge_version}",
        "",
        "---",
        "",
        "## Detected Files",
        "",
    ]

    def row(label: str, role: str, path: object) -> str:
        display = _display_path(run, role, path)
        v = f"`{display}`" if display else "_not found_"
        suffix = " _(user-specified)_" if role in run.user_overrides else ""
        return f"- **{label}:** {v}{suffix}"

    lines += [
        row("Topology input (.tpr)", "tpr",        run.tpr),
        row("Primary trajectory",    "trajectory", run.trajectory),
        row("Reference structure",   "structure",  run.structure),
        row("Energy file (.edr)",    "edr",        run.edr),
        row("Topology text (.top)",  "topology",   run.topology),
        row("Index file (.ndx)",     "index",      run.index),
        row("Primary log",           "log",        run.log),
        f"- **Existing XVG files:** {len(run.xvg_files)}",
    ]

    if len(run.trajectory_candidates) > 1:
        lines += ["", "**Trajectory candidates:**"]
        for p in run.trajectory_candidates:
            marker = " ← selected" if p == run.trajectory else ""
            disp   = _display_path(run, "_candidate", p)
            lines.append(f"  - `{disp}`{marker}")

    if run.warnings:
        lines += ["", "## Warnings", ""]
        for w in run.warnings:
            icon = "⚠" if w.severity == "warn" else "ℹ"
            lines.append(f"- {icon} `{w.code}`: {w.message}")

    lines += [
        "",
        "## Planned Observables",
        "",
        "| Observable | Status | Notes |",
        "|---|---|---|",
    ]
    for obs in plan.observables:
        status = _STATUS_ICONS.get(obs.status, obs.status)
        notes  = "; ".join(obs.warnings) if obs.warnings else ""
        lines.append(f"| {obs.name} | {status} | {notes} |")

    active = [o for o in plan.observables if o.gmx_command and o.status in ("planned", "needs_selection")]
    if active:
        lines += ["", "## Proposed GROMACS Commands", "", "```bash"]
        for obs in active:
            lines.append(f"# {obs.name}")
            lines.append(obs.gmx_command)  # type: ignore[arg-type]
            lines.append("")
        lines.append("```")

    lines += ["", "## Next Steps", ""]
    if not run.tpr:
        lines.append("1. **Locate or regenerate `.tpr`** — required for most GROMACS analysis.")
    if not run.trajectory:
        lines.append("1. **Locate trajectory** (`.xtc` or `.trr`) before running analysis.")
    if active:
        lines.append("1. Run the proposed GROMACS commands above to generate observables.")
        lines.append("2. Re-run `simforge analyze-md` (without `--dry-run`) to analyze the results.")
    else:
        lines.append("1. Ensure all required files are present, then re-run without `--dry-run`.")

    (output_dir / "report.md").write_text("\n".join(lines) + "\n")
