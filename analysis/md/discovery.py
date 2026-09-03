"""MD run discovery: scan a directory for GROMACS artifacts and infer their roles."""
from __future__ import annotations

from pathlib import Path

from analysis.md.models import AnalysisWarning, DiscoveredFile, DiscoveredMDRun

_GROMACS_SUFFIXES = {
    ".tpr", ".xtc", ".trr", ".gro", ".pdb",
    ".edr", ".xvg", ".top", ".itp", ".ndx", ".log",
}
_JOB_LOG_NAMES = {"slurm.out", "nohup.out"}


def discover_md_run(
    run_dir: Path,
    overrides: dict[str, str] | None = None,
) -> DiscoveredMDRun:
    """Scan run_dir recursively for GROMACS artifacts and return a DiscoveredMDRun.

    overrides maps role names ("tpr", "trajectory", "structure", "edr", "index",
    "log") to path strings exactly as the user supplied them.  Non-absolute paths
    are resolved relative to run_dir.  Issues AnalysisWarnings instead of raising
    for missing files.
    """
    run_dir_abs = run_dir.resolve()
    run = DiscoveredMDRun(run_dir=run_dir_abs)

    candidates: list[Path] = []
    for p in sorted(run_dir_abs.rglob("*")):
        if p.is_file() and (p.suffix in _GROMACS_SUFFIXES or p.name in _JOB_LOG_NAMES):
            candidates.append(p)

    by_suffix: dict[str, list[Path]] = {}
    for p in candidates:
        key = p.name if p.name in _JOB_LOG_NAMES else p.suffix
        by_suffix.setdefault(key, []).append(p)

    def get(suffix: str) -> list[Path]:
        return by_suffix.get(suffix, [])

    # ── TPR ──────────────────────────────────────────────────────────────────
    tpr_files = get(".tpr")
    if not tpr_files:
        run.warnings.append(AnalysisWarning("no_tpr", "No .tpr file found — topology input missing"))
    else:
        run.tpr = _pick_primary(tpr_files, preferred_stem="md")

    # ── Trajectory ───────────────────────────────────────────────────────────
    xtc = get(".xtc")
    trr = get(".trr")
    traj_candidates = xtc + trr
    if not traj_candidates:
        run.warnings.append(AnalysisWarning("no_trajectory", "No trajectory file (.xtc or .trr) found"))
    else:
        run.trajectory = _pick_primary(xtc if xtc else trr, preferred_stem="md")
        run.trajectory_candidates = traj_candidates
        if len(traj_candidates) > 1:
            run.warnings.append(AnalysisWarning(
                "multiple_trajectories",
                f"Multiple trajectories found — selected {run.trajectory.name} as primary",
            ))

    # ── Structure ────────────────────────────────────────────────────────────
    gro = get(".gro")
    pdb = get(".pdb")
    if gro:
        run.structure = _pick_primary(gro, preferred_stem="md")
    elif pdb:
        run.structure = _pick_primary(pdb, preferred_stem="md")
    else:
        run.warnings.append(AnalysisWarning("no_structure", "No structure file (.gro or .pdb) found"))

    # ── Energy ───────────────────────────────────────────────────────────────
    edr = get(".edr")
    if not edr:
        run.warnings.append(AnalysisWarning("no_edr", "No .edr file found — energy QC unavailable"))
    else:
        run.edr = _pick_primary(edr, preferred_stem="md")

    # ── Topology text ────────────────────────────────────────────────────────
    top = get(".top")
    if top:
        run.topology = _pick_primary(top, preferred_stem="topol")

    # ── Index ────────────────────────────────────────────────────────────────
    ndx = get(".ndx")
    if ndx:
        run.index = _pick_primary(ndx, preferred_stem="index")
    else:
        run.warnings.append(AnalysisWarning(
            "no_ndx",
            "No .ndx file found — selections will use default GROMACS groups",
            severity="info",
        ))

    # ── Logs ─────────────────────────────────────────────────────────────────
    log_files = get(".log") + get("slurm.out") + get("nohup.out")
    if log_files:
        run.log = _pick_primary(log_files, preferred_stem="md")

    # ── Existing XVG ─────────────────────────────────────────────────────────
    run.xvg_files = get(".xvg")

    # ── Apply user overrides (before building all_files so roles are correct) ─
    if overrides:
        _apply_overrides(run, overrides, run_dir_abs)

    # ── Build annotated file list ─────────────────────────────────────────────
    for p in candidates:
        run.all_files.append(DiscoveredFile(
            path=p,
            role=_infer_role(p, run),
            size_bytes=p.stat().st_size,
        ))

    # Add user-specified files that were not in the discovered candidates
    discovered_paths = {f.path for f in run.all_files}
    for role, path_str in run.user_overrides.items():
        resolved = _resolve_override_path(path_str, run_dir_abs)
        if resolved is not None and resolved not in discovered_paths:
            run.all_files.append(DiscoveredFile(
                path=resolved,
                role=_infer_role(resolved, run),
                size_bytes=resolved.stat().st_size,
            ))

    return run


def _resolve_override_path(path_str: str, run_dir_abs: Path) -> Path | None:
    p = Path(path_str)
    if not p.is_absolute():
        p = run_dir_abs / p
    return p if p.exists() else None


def _apply_overrides(
    run: DiscoveredMDRun,
    overrides: dict[str, str],
    run_dir_abs: Path,
) -> None:
    """Apply user-specified artifact overrides in-place, updating warnings."""
    _role_attr = {
        "tpr":        "tpr",
        "trajectory": "trajectory",
        "structure":  "structure",
        "edr":        "edr",
        "index":      "index",
        "log":        "log",
    }

    for role, path_str in overrides.items():
        if role not in _role_attr:
            continue
        resolved = _resolve_override_path(path_str, run_dir_abs)
        if resolved is None:
            run.warnings.append(AnalysisWarning(
                "override_not_found",
                f"User-specified {role} not found: {path_str}",
            ))
            continue

        setattr(run, _role_attr[role], resolved)
        run.user_overrides[role] = path_str

        if role == "trajectory":
            _update_trajectory_warning(run, path_str)


def _update_trajectory_warning(run: DiscoveredMDRun, user_path_str: str) -> None:
    """Replace the auto-selection multiple_trajectories warning with a user-specified notice."""
    had_multiple = any(w.code == "multiple_trajectories" for w in run.warnings)
    run.warnings = [w for w in run.warnings if w.code != "multiple_trajectories"]
    if had_multiple:
        run.warnings.append(AnalysisWarning(
            "multiple_trajectories",
            f"Multiple trajectories found, but user-specified trajectory was used: {user_path_str}",
            severity="info",
        ))


def _pick_primary(files: list[Path], preferred_stem: str) -> Path:
    """Return the file matching preferred_stem, falling back to the largest file."""
    for p in files:
        if p.stem == preferred_stem:
            return p
    return max(files, key=lambda p: p.stat().st_size)


def _infer_role(p: Path, run: DiscoveredMDRun) -> str:
    if p == run.tpr:
        return "topology_input"
    if p == run.trajectory:
        return "primary_trajectory"
    if p in run.trajectory_candidates:
        return "trajectory_candidate"
    if p == run.structure:
        return "structure"
    if p == run.edr:
        return "energy"
    if p == run.topology:
        return "topology_text"
    if p == run.index:
        return "index"
    if p == run.log:
        return "log"
    if p.suffix == ".xvg":
        return "xvg_analysis"
    if p.suffix == ".itp":
        return "include_topology"
    if p.name in _JOB_LOG_NAMES:
        return "job_log"
    return "other"
