"""Dry-run observable planning: given discovered files, plan possible GROMACS analyses."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from analysis.md.models import AnalysisWarning, DiscoveredMDRun


@dataclass
class AnalysisObservablePlan:
    name: str
    status: str                          # planned | unavailable | needs_selection | already_available
    required_inputs: list[str]
    expected_outputs: list[str]
    gmx_command: str | None = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "status": self.status,
            "required_inputs": self.required_inputs,
            "expected_outputs": self.expected_outputs,
            "gmx_command": self.gmx_command,
            "warnings": self.warnings,
        }


@dataclass
class AnalysisPlan:
    observables: list[AnalysisObservablePlan] = field(default_factory=list)
    warnings: list[AnalysisWarning] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "observables": [o.to_dict() for o in self.observables],
            "warnings": [w.to_dict() for w in self.warnings],
        }


def _cmd_path(run: DiscoveredMDRun, role: str, path: Path | None) -> str:
    """Return the path string to embed in a GROMACS command.

    Preference order:
    1. User-supplied path string (preserves exactly what the user typed).
    2. Path relative to run_dir (keeps subdirectory context, e.g. R-1/md.xtc).
    3. Fallback basename placeholder when the artifact was not found.
    """
    if role in run.user_overrides:
        return run.user_overrides[role]
    if path is None:
        _fallbacks = {"tpr": "md.tpr", "trajectory": "md.xtc", "edr": "md.edr"}
        return _fallbacks.get(role, role)
    try:
        return str(path.relative_to(run.run_dir))
    except ValueError:
        return str(path)


def plan_observables(run: DiscoveredMDRun) -> AnalysisPlan:
    """Return an AnalysisPlan with planned/unavailable observables. No GROMACS is run."""
    plan = AnalysisPlan()

    has_tpr  = run.tpr is not None
    has_traj = run.trajectory is not None
    has_edr  = run.edr is not None
    has_ndx  = run.index is not None
    has_xvg  = bool(run.xvg_files)

    tpr_and_traj = has_tpr and has_traj
    tpr_s  = _cmd_path(run, "tpr",        run.tpr)
    traj_s = _cmd_path(run, "trajectory", run.trajectory)

    def _traj_obs(name: str, cmd_template: str, outputs: list[str]) -> AnalysisObservablePlan:
        if tpr_and_traj:
            return AnalysisObservablePlan(
                name=name,
                status="planned",
                required_inputs=["tpr", "trajectory"],
                expected_outputs=outputs,
                gmx_command=cmd_template.format(tpr=tpr_s, traj=traj_s),
            )
        return AnalysisObservablePlan(
            name=name,
            status="unavailable",
            required_inputs=["tpr", "trajectory"],
            expected_outputs=outputs,
            warnings=["requires .tpr and trajectory"],
        )

    # RMSD
    plan.observables.append(_traj_obs(
        "rmsd",
        "echo Backbone Backbone | gmx rms -s {tpr} -f {traj} -o rmsd.xvg",
        ["rmsd.xvg"],
    ))

    # RMSF
    plan.observables.append(_traj_obs(
        "rmsf",
        "echo Backbone | gmx rmsf -s {tpr} -f {traj} -o rmsf.xvg -res",
        ["rmsf.xvg"],
    ))

    # Radius of gyration
    plan.observables.append(_traj_obs(
        "radius_of_gyration",
        "echo Protein | gmx gyrate -s {tpr} -f {traj} -o gyrate.xvg",
        ["gyrate.xvg"],
    ))

    # SASA
    plan.observables.append(_traj_obs(
        "sasa",
        "echo Protein | gmx sasa -s {tpr} -f {traj} -o sasa.xvg",
        ["sasa.xvg"],
    ))

    # H-bonds — needs selection warning when no .ndx
    if tpr_and_traj:
        hbond_status = "needs_selection" if not has_ndx else "planned"
        hbond_warn   = ["no .ndx — will use default GROMACS groups (Protein-Protein)"] if not has_ndx else []
        plan.observables.append(AnalysisObservablePlan(
            name="hydrogen_bonds",
            status=hbond_status,
            required_inputs=["tpr", "trajectory"],
            expected_outputs=["hbnum.xvg"],
            gmx_command=f"echo Protein Protein | gmx hbond -s {tpr_s} -f {traj_s} -num hbnum.xvg",
            warnings=hbond_warn,
        ))
    else:
        plan.observables.append(AnalysisObservablePlan(
            name="hydrogen_bonds",
            status="unavailable",
            required_inputs=["tpr", "trajectory"],
            expected_outputs=["hbnum.xvg"],
            warnings=["requires .tpr and trajectory"],
        ))

    # Energy / temperature / pressure QC
    if has_edr:
        edr_s = _cmd_path(run, "edr", run.edr)
        plan.observables.append(AnalysisObservablePlan(
            name="energy_temperature_pressure_qc",
            status="planned",
            required_inputs=["edr"],
            expected_outputs=["energy.xvg", "temperature.xvg", "pressure.xvg"],
            gmx_command=(
                f"echo Potential Temperature Pressure | gmx energy -f {edr_s} -o energy.xvg"
            ),
        ))
    else:
        plan.observables.append(AnalysisObservablePlan(
            name="energy_temperature_pressure_qc",
            status="unavailable",
            required_inputs=["edr"],
            expected_outputs=["energy.xvg", "temperature.xvg", "pressure.xvg"],
            warnings=["requires .edr file"],
        ))

    # Existing XVG summary
    plan.observables.append(AnalysisObservablePlan(
        name="existing_xvg_summary",
        status="already_available" if has_xvg else "unavailable",
        required_inputs=["xvg_files"],
        expected_outputs=["xvg_summary (in report)"],
        warnings=[] if has_xvg else ["no .xvg files found in run directory"],
    ))

    return plan
