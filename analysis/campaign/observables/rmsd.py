"""RMSD observables.

Four scientifically distinct analyses (spec §48) — the reference frame is part
of the definition and is recorded in every result:

Reference coordinates for the least-squares fit and the measurement are taken
from the production topology (``.tpr`` — the start-of-production frame written
by ``grompp`` from the equilibrated system), falling back to a ``.gro``/``.pdb``
only when no ``.tpr`` is available.

``rmsd-receptor``
    Fit receptor backbone to the reference; measure receptor backbone RMSD.
    Intrinsic receptor conformational change relative to the production start.

``rmsd-complex``
    Fit the whole complex backbone to the reference; measure complex backbone RMSD.

``rmsd-peptide-intrinsic``
    Fit peptide backbone to the reference; measure peptide backbone RMSD.
    The peptide's *own* conformational drift, receptor motion removed.

``rmsd-peptide-receptor-frame``
    Fit the trajectory to the *receptor* backbone, then measure peptide backbone
    RMSD with **no further fitting**.  Captures how the peptide moves *relative
    to the receptor* (binding-pose stability).

The first three use a ``whole`` trajectory view and let ``gmx rms`` do the
least-squares fit to the measurement group.  The last uses a receptor-fitted
view and ``gmx rms -nofit``.
"""
from __future__ import annotations

from pathlib import Path

from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, CampaignWarning, ComponentType, Severity, ViewKind,
)
from analysis.campaign.observables.base import AnalysisContext, ObservableSpec
from analysis.campaign.observables.registry import register
from analysis.campaign.trajectory import requirements as reqs


def _summarise_xvg(path: Path) -> dict:
    try:
        from runtime.xvg_parser import parse_xvg
        data = parse_xvg(path)
    except Exception as exc:  # noqa: BLE001
        return {"parse_error": str(exc)}
    if not data.series:
        return {"parse_error": "no data series in xvg"}
    ys = data.series[0].values if hasattr(data.series[0], "values") else data.series[0]
    ys = [float(v) for v in ys]
    ts = [float(t) for t in getattr(data, "time_ps", [])]
    if not ys:
        return {"parse_error": "empty rmsd series"}
    n = len(ys)
    mean = sum(ys) / n
    tail = ys[max(1, int(0.8 * n)):]
    return {
        "n_frames": n,
        "rmsd_unit": "nm",
        "time_unit": "ps",
        "mean_nm": round(mean, 5),
        "min_nm": round(min(ys), 5),
        "max_nm": round(max(ys), 5),
        "final_nm": round(ys[-1], 5),
        "plateau_mean_nm": round(sum(tail) / len(tail), 5) if tail else None,
        "start_time_ps": ts[0] if ts else None,
        "end_time_ps": ts[-1] if ts else None,
    }


class _RmsdBase(ObservableSpec):
    fit_group: str = "Receptor_Backbone"
    measure_group: str = "Receptor_Backbone"
    reference_frame_desc: str = ""
    use_prefit_view: bool = False       # True => view fitted, gmx rms -nofit

    def trajectory_requirements(self, ctx_params: dict):
        if self.use_prefit_view:
            return reqs.fit_to(
                self.fit_group,
                f"{self.id}: measure '{self.measure_group}' in the "
                f"'{self.fit_group}'-aligned frame",
            )
        return reqs.whole_only(
            f"{self.id}: whole molecules; gmx rms performs the least-squares fit "
            f"to '{self.measure_group}'"
        )

    def execute(self, ctx: AnalysisContext) -> AnalysisResult:
        result = AnalysisResult(
            analysis_id=self.id, system_id=ctx.system.system_id,
            status=AnalysisStatus.PLANNED,
            fit_selection=self.fit_group, measure_selection=self.measure_group,
            reference_frame=self.reference_frame_desc,
            trajectory_view_kind=ctx.trajectory_view.kind,
            parameters=dict(ctx.parameters),
        )

        index = ctx.semantic_index
        groups = set(index.group_names()) if index else set()
        needed = {self.measure_group}
        if not self.use_prefit_view:
            needed.add(self.fit_group)
        missing = sorted(needed - groups)
        if missing or not index or not index.path:
            result.status = AnalysisStatus.REVIEW_REQUIRED
            result.message = (
                f"semantic index is missing group(s) {missing or list(needed)}; "
                f"cannot construct an explicit RMSD selection"
            )
            return result

        view = ctx.trajectory_view
        if not view.safe or not view.path:
            result.status = AnalysisStatus.SKIPPED
            result.message = (
                "required trajectory view could not be built safely: "
                + "; ".join(w.message for w in view.warnings)
            )
            return result

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        out_xvg = ctx.output_dir / "rmsd.xvg"
        # Reference = the production topology (.tpr): its coordinates are the
        # start-of-production frame that `grompp` wrote from the equilibrated
        # system.  This is the standard "RMSD relative to the production start".
        # Only fall back to a .gro/.pdb structure if there is no .tpr.
        ref = ctx.topology_path
        if not ref or not ref.lower().endswith(".tpr"):
            ref = ctx.structure_path or ctx.topology_path
        result.reference_frame = (
            f"{self.reference_frame_desc}; reference coordinates = {Path(ref).name}"
        )
        args = [
            "rms", "-s", ref,
            "-f", view.path, "-n", index.path, "-o", str(out_xvg),
            "-tu", "ps",
        ]
        if self.use_prefit_view:
            args.append("-nofit")
            stdin = f"{self.measure_group}\n"
        else:
            stdin = f"{self.fit_group}\n{self.measure_group}\n"

        if ctx.dry_run:
            result.status = AnalysisStatus.PLANNED
            result.message = f"planned: gmx {' '.join(args)}  (stdin: {stdin!r})"
            result.data_summary = {"planned_command": ["gmx", *args], "stdin": stdin}
            return result

        res = run_gmx(args, stdin=stdin, gmx=ctx.gmx, timeout=3600)
        (ctx.output_dir / "gmx_rms.log").write_text(
            f"$ {' '.join(res.argv)}\nSTDIN:\n{stdin}\n\nSTDERR:\n{res.stderr}\n\nSTDOUT:\n{res.stdout}\n"
        )
        if not res.ok or not out_xvg.is_file():
            result.status = AnalysisStatus.FAILED
            result.message = f"gmx rms failed (rc={res.returncode}): {res.stderr.strip()[-400:]}"
            result.warnings.append(CampaignWarning("gmx_rms_failed", result.message, Severity.ERROR))
            return result

        result.status = AnalysisStatus.SUCCESS
        result.output_files = [str(out_xvg.resolve()), str((ctx.output_dir / "gmx_rms.log").resolve())]
        result.data_summary = _summarise_xvg(out_xvg)
        result.data_summary["gmx_version"] = res.gmx_version
        result.message = f"receptor-frame RMSD computed" if self.use_prefit_view else "RMSD computed"
        return result


class RmsdReceptor(_RmsdBase):
    id = "rmsd-receptor"
    display_name = "RMSD — receptor"
    description = ("Intrinsic receptor backbone RMSD relative to the production start "
                   "(reference = .tpr coordinates), after least-squares backbone fit.")
    category = "structural"
    required_components = (ComponentType.RECEPTOR,)
    fit_group = "Receptor_Backbone"
    measure_group = "Receptor_Backbone"
    reference_frame_desc = "receptor backbone least-squares fitted to the reference structure, then receptor backbone RMSD measured"


class RmsdComplex(_RmsdBase):
    id = "rmsd-complex"
    display_name = "RMSD — complex"
    description = ("Whole-complex backbone RMSD relative to the production start "
                   "(reference = .tpr coordinates), after least-squares backbone fit.")
    category = "structural"
    required_components = (ComponentType.COMPLEX,)
    fit_group = "Complex_Backbone"
    measure_group = "Complex_Backbone"
    reference_frame_desc = "complex backbone least-squares fitted to the reference structure, then complex backbone RMSD measured"


class RmsdPeptideIntrinsic(_RmsdBase):
    id = "rmsd-peptide-intrinsic"
    display_name = "RMSD — peptide (intrinsic)"
    description = ("Peptide backbone RMSD relative to the production start (reference = .tpr "
                   "coordinates), after least-squares peptide-backbone fit. Receptor motion "
                   "is removed — this is the peptide's own conformational drift.")
    category = "structural"
    required_components = (ComponentType.PEPTIDE,)
    fit_group = "Peptide_Backbone"
    measure_group = "Peptide_Backbone"
    reference_frame_desc = "peptide backbone least-squares fitted to the reference structure, then peptide backbone RMSD measured"


class RmsdPeptideReceptorFrame(_RmsdBase):
    id = "rmsd-peptide-receptor-frame"
    display_name = "RMSD — peptide (receptor frame)"
    description = ("Peptide backbone RMSD measured in the receptor-aligned frame "
                   "(trajectory fitted to the receptor backbone, no further fitting). "
                   "Captures peptide motion relative to the receptor / binding-pose stability.")
    category = "interaction"
    required_components = (ComponentType.RECEPTOR, ComponentType.PEPTIDE)
    fit_group = "Receptor_Backbone"
    measure_group = "Peptide_Backbone"
    reference_frame_desc = "trajectory least-squares fitted to the receptor backbone of the reference structure; peptide backbone RMSD then measured with -nofit"
    use_prefit_view = True


BUILTIN_OBSERVABLES = (
    RmsdReceptor, RmsdComplex, RmsdPeptideIntrinsic, RmsdPeptideReceptorFrame,
)


def register_builtins() -> None:
    """Idempotently register the built-in RMSD observables.

    Safe to call after ``registry._REGISTRY.clear()`` (test isolation): unlike a
    module-level import side-effect it re-registers every time.
    """
    from analysis.campaign.observables.registry import is_registered, register
    for cls in BUILTIN_OBSERVABLES:
        obs = cls()
        if not is_registered(obs.id):
            register(obs)


register_builtins()
