"""Generic, annotation-aware observables (Trajectory Review Phase 7).

Every observable here is an ordinary :class:`ObservableSpec` on top of the
existing contract: selections (``component:<type>`` / ``annotation:<id>``)
resolve through :mod:`selections`; coordinates come from the policy-planned
view; results are :class:`ResultArray` objects whose definition tokens carry
the exact atom evidence of every selection used (and nothing else).

PBC / coordinate semantics are part of each definition.  The GROMACS 2025
analysis tools used here (gyrate, sasa, distance, hbond) reassemble molecules
in memory from ``.tpr`` bonds (``-rmpbc``) and apply periodic minimum-image
geometry (``-pbc``), so they run on the raw trajectory — no derived
"whole" trajectory is generated.  A ``.tpr`` is therefore required.
``mindist`` measures minimum-image interatomic distances (``-pbc``).
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

from analysis.campaign.fingerprint import fingerprint_file
from analysis.campaign.gmx import gmx_version, run_gmx
from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, Axis, AxisKind, CampaignWarning, ClassificationState,
    FrameAlignment, FrameAlignmentMode, MissingSemantics, ObservablePurpose, ResultArray,
    Severity, StorageFormat, StorageRef, TrajectoryRequirements,
)
from analysis.campaign.observables.base import AnalysisContext, ObservableSpec
from analysis.campaign.observables.selections import (
    ResolvedSelection, SelectionRef, resolve_selection,
)

GENERIC_DEFINITION_VERSION = "generic/v1"


def _params(ctx_params: Optional[dict], obs_id: str) -> dict:
    """Flat parameters, overridden by a per-observable block ``{obs_id: {...}}``."""
    p = {k: v for k, v in (ctx_params or {}).items() if not isinstance(v, dict)}
    p.update((ctx_params or {}).get(obs_id) or {})
    return p


# ═══════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ═══════════════════════════════════════════════════════════════════════════════

def time_alignment(ctx: AnalysisContext, times: list[float]) -> tuple[FrameAlignment, dict]:
    """Verify one-output-row-per-frame against the view's own FrameTimeIndex.

    Times are compared within the output's printing resolution; a mismatch is
    recorded (mode UNKNOWN), never papered over.
    """
    from analysis.campaign.trajectory.preprocessor import load_view_time_index
    from analysis.campaign.trajectory.time_index import GromacsTimeIndexBackend, index_trajectory
    view = ctx.trajectory_view
    idx = load_view_time_index(view)
    source = "view time index (Phase 3)"
    if idx is None and view.path:
        idx, _ = index_trajectory(view.path, backend=GromacsTimeIndexBackend(ctx.gmx),
                                  cache_dir=ctx.time_index_cache_dir)
        source = "trajectory time index (Phase 1)"
    attrs: dict = {"alignment_source": source}
    if idx is None or not idx.times_ps:
        attrs["alignment_check"] = "no time index available"
        return FrameAlignment(mode=FrameAlignmentMode.UNKNOWN, view_ref=view.cache_key or None), attrs
    ok = len(times) == idx.n_frames and all(
        abs(a - b) <= max(1e-3, 6e-6 * abs(b)) for a, b in zip(times, idx.times_ps))
    attrs.update({"rows": len(times), "frames": idx.n_frames,
                  "alignment_check": "verified: one row per frame, times match" if ok else
                  "MISMATCH: rows/times differ from the trajectory timeline"})
    return FrameAlignment(
        mode=FrameAlignmentMode.ONE_ROW_PER_FRAME if ok else FrameAlignmentMode.UNKNOWN,
        trajectory_digest=idx.fingerprint.digest if idx.fingerprint else None,
        view_ref=view.cache_key or None), attrs


class _GenericObservable(ObservableSpec):
    """Selection-parameterised observable executed by one GROMACS tool."""
    tool: str = ""
    quantity: str = ""
    expected_unit: str = ""
    selection_params: tuple[str, ...] = ()
    value_column: int = 1
    #: rationale / semantics recorded in the definition
    coordinate_semantics: dict = {}

    # ── selections ─────────────────────────────────────────────────────────
    def default_selection(self, name: str, system) -> Optional[str]:
        return None

    def selection_refs(self, params: dict, system) -> tuple[dict, list[str]]:
        refs, errors = {}, []
        for name in self.selection_params:
            raw = params.get(name) or self.default_selection(name, system)
            if not raw:
                errors.append(f"parameter '{name}' is required (component:<type> or "
                              f"annotation:<id>); no unambiguous default exists")
                continue
            try:
                refs[name] = SelectionRef.parse(raw)
            except ValueError as exc:
                errors.append(str(exc))
        return refs, errors

    def resolve(self, ctx_or_system, params: dict, semantic_index=None
                ) -> tuple[dict[str, ResolvedSelection], list[str]]:
        system = getattr(ctx_or_system, "system", ctx_or_system)
        # an explicitly passed index wins (SystemRecord also has a semantic_index field)
        index = semantic_index if semantic_index is not None else getattr(
            ctx_or_system, "semantic_index", None)
        refs, errors = self.selection_refs(params, system)
        resolved = {}
        for name, ref in refs.items():
            r = resolve_selection(ref, system, index)
            resolved[name] = r
            if not r.ok:
                errors.append(f"{name} ({ref}): {r.reason}")
        return resolved, errors

    def applicability(self, system, params=None, semantic_index=None) -> dict:
        _, errors = self.resolve(system, _params(params, self.id), semantic_index)
        missing_a = [e.split("'")[1] for e in errors if "annotation '" in e]
        return {"applicable": not errors, "reasons": errors, "missing_components": [],
                "missing_annotations": missing_a}

    def parameters_schema(self) -> dict:
        return {n: "component:<type> | annotation:<id>" for n in self.selection_params}

    def extra_parameters(self, params: dict) -> dict:
        return {}

    # ── contract ───────────────────────────────────────────────────────────
    def output_schema(self, params=None) -> list[ResultArray]:
        return [ResultArray(name=self.quantity, quantity=self.quantity, unit=self.expected_unit,
                            axes=[Axis(name="time", kind=AxisKind.TIME, unit="ps")],
                            missing=MissingSemantics.NONE_EXPECTED)]

    def definition_evidence(self, ctx: AnalysisContext):
        params = _params(ctx.parameters, self.id)
        resolved, errors = self.resolve(ctx, params)
        tpr = ctx.topology_path
        if errors or not tpr or not tpr.lower().endswith(".tpr") or not Path(tpr).is_file():
            return None
        view = ctx.trajectory_view
        return {
            "observable": self.id, "definition_version": GENERIC_DEFINITION_VERSION,
            "quantity": self.quantity,
            "selections": {n: r.evidence for n, r in sorted(resolved.items())},
            "parameters": self.extra_parameters(params),
            "coordinate_semantics": self.coordinate_semantics,
            "topology": {"digest": fingerprint_file(Path(tpr)).digest,
                         "role": "masses / bonds for in-memory whole molecules"},
            "trajectory_view": {"kind": view.kind, "cache_key": view.cache_key},
            "backend": {"tool": f"gmx {self.tool}", "version": gmx_version(ctx.gmx)},
        }

    def build_command(self, ctx, resolved, params, out: Path) -> tuple[list[str], Optional[str]]:
        raise NotImplementedError

    def unit_fallback(self) -> tuple[Optional[str], str]:
        return None, ""

    # ── execution ──────────────────────────────────────────────────────────
    def execute(self, ctx: AnalysisContext) -> AnalysisResult:
        from analysis.campaign.results import definition_token, normalize_unit, xvg_units
        params = _params(ctx.parameters, self.id)
        result = AnalysisResult(analysis_id=self.id, system_id=ctx.system.system_id,
                                status=AnalysisStatus.PLANNED,
                                trajectory_view_kind=ctx.trajectory_view.kind,
                                parameters=dict(params))
        resolved, errors = self.resolve(ctx, params)
        if errors:
            result.status = AnalysisStatus.REVIEW_REQUIRED
            result.message = "; ".join(errors)
            return result
        tpr = ctx.topology_path
        if not tpr or not tpr.lower().endswith(".tpr") or not Path(tpr).is_file():
            result.status = AnalysisStatus.SKIPPED
            result.message = (f"gmx {self.tool} needs a .tpr (masses and bonds for in-memory "
                              f"whole molecules); none available")
            return result
        view = ctx.trajectory_view
        if not view.safe or not view.path:
            result.status = AnalysisStatus.SKIPPED
            result.message = ("trajectory view not available: "
                              + "; ".join(w.message for w in view.warnings))
            return result
        evidence = self.definition_evidence(ctx)
        result.definition_evidence = evidence
        result.definition_token = definition_token(evidence)
        result.measure_selection = " / ".join(r.group for r in resolved.values())

        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        out = ctx.output_dir / f"{self.quantity}.xvg"
        argv, stdin = self.build_command(ctx, resolved, params, out)
        result.data_summary = {"command": ["gmx", *argv], "stdin": stdin}
        if ctx.dry_run:
            result.message = f"planned: gmx {' '.join(argv)}"
            return result
        res = run_gmx(argv, stdin=stdin, gmx=ctx.gmx, timeout=4 * 3600)
        log = ctx.output_dir / f"gmx_{self.tool}.log"
        log.write_text(f"$ {' '.join(res.argv)}\nSTDIN:\n{stdin or ''}\n\nSTDERR:\n{res.stderr}\n"
                       f"\nSTDOUT:\n{res.stdout}\n")
        if not res.ok or not out.is_file():
            result.status = AnalysisStatus.FAILED
            result.message = f"gmx {self.tool} failed (rc={res.returncode}): {res.stderr.strip()[-400:]}"
            result.warnings.append(CampaignWarning(f"gmx_{self.tool}_failed", result.message,
                                                   Severity.ERROR))
            return result

        from runtime.xvg_parser import parse_xvg
        data = parse_xvg(out)
        x_unit, y_unit_raw = xvg_units(out)
        unit, unit_source = normalize_unit(y_unit_raw), "xvg header"
        if unit is None:
            unit, unit_source = self.unit_fallback()
        time_attrs = {"unit_source": "xvg header" if x_unit else "command flag -tu ps"}
        alignment, align_attrs = time_alignment(ctx, list(data.time_ps))
        time_attrs.update(align_attrs)
        if alignment.mode != FrameAlignmentMode.ONE_ROW_PER_FRAME:
            result.warnings.append(CampaignWarning(
                "time_alignment_unverified", align_attrs.get("alignment_check", ""), Severity.WARN))
        values = data.series[self.value_column - 1].values if data.series else []
        path = str(out.resolve())
        result.arrays = [ResultArray(
            name=self.quantity, quantity=self.quantity, unit=unit,
            axes=[Axis(name="time", kind=AxisKind.TIME, unit=x_unit or "ps",
                       values_ref=StorageRef(path=path, format=StorageFormat.XVG, column=0),
                       alignment=alignment, attrs=time_attrs)],
            storage=StorageRef(path=path, format=StorageFormat.XVG, column=self.value_column,
                               fingerprint=fingerprint_file(out)),
            view_ref=view.cache_key or None, definition_token=result.definition_token,
            missing=MissingSemantics.NONE_EXPECTED,
            attrs={"unit_source": unit_source, "selections": {n: str(r.ref) for n, r in resolved.items()},
                   "groups": {n: r.group for n, r in resolved.items()},
                   "coordinate_semantics": self.coordinate_semantics})]
        finite = [v for v in values if math.isfinite(v)]
        result.data_summary.update({
            "n_rows": len(values), "unit": unit,
            "mean": round(sum(finite) / len(finite), 6) if finite else None,
            "min": min(finite) if finite else None, "max": max(finite) if finite else None,
            "gmx_version": res.gmx_version})
        result.output_files = [path, str(log.resolve())]
        result.status = AnalysisStatus.SUCCESS
        result.message = f"{self.quantity} computed ({len(values)} rows)"
        return result


def _default_receptor(system) -> Optional[str]:
    comps = [c for c in system.components if c.component_type == "receptor"]
    if len(comps) == 1 and comps[0].classification_state == ClassificationState.RESOLVED:
        return "component:receptor"
    return None


def _g(r: ResolvedSelection) -> str:
    return f'group "{r.group}"'


# ═══════════════════════════════════════════════════════════════════════════════
# Observables
# ═══════════════════════════════════════════════════════════════════════════════

class RadiusOfGyration(_GenericObservable):
    id = "rg"
    display_name = "Radius of gyration"
    description = "Mass-weighted radius of gyration of one selection over time (gmx gyrate)."
    category = "structural"
    purpose = ObservablePurpose.INTRAMOLECULAR_SHAPE
    tool = "gyrate"
    quantity = "radius_of_gyration"
    expected_unit = "nm"
    selection_params = ("selection",)
    value_column = 1                                  # total Rg (then Rx, Ry, Rz)
    coordinate_semantics = {"molecules": "made whole in memory from .tpr bonds (-rmpbc)",
                            "weighting": "mass (-mode mass)"}

    def default_selection(self, name, system):
        return _default_receptor(system)

    def trajectory_requirements(self, ctx_params):
        return TrajectoryRequirements(rationale="gyrate makes molecules whole in memory (-rmpbc); "
                                                "raw coordinates suffice")

    def build_command(self, ctx, resolved, params, out):
        return (["gyrate", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
                 "-n", ctx.semantic_index.path, "-sel", _g(resolved["selection"]),
                 "-mode", "mass", "-o", str(out), "-tu", "ps"], None)


class SolventAccessibleSurface(_GenericObservable):
    id = "sasa"
    display_name = "SASA"
    description = ("Solvent-accessible surface area of an output selection computed within a "
                   "surface selection (gmx sasa).")
    category = "structural"
    purpose = ObservablePurpose.INTRAMOLECULAR_SHAPE
    tool = "sasa"
    quantity = "sasa"
    expected_unit = "nm^2"
    selection_params = ("selection", "surface")
    value_column = 2                                  # the output selection's area
    coordinate_semantics = {"molecules": "made whole in memory (-rmpbc)",
                            "pbc": "periodic neighbour search (-pbc)"}

    def default_selection(self, name, system):
        return _default_receptor(system)

    def selection_refs(self, params, system):
        # the surface defaults to the measured selection itself (explicitly recorded)
        if not params.get("surface") and params.get("selection"):
            params = {**params, "surface": params["selection"]}
        return super().selection_refs(params, system)

    def extra_parameters(self, params):
        return {"probe_nm": float(params.get("probe", 0.14)),
                "ndots": int(params.get("ndots", 24))}

    def trajectory_requirements(self, ctx_params):
        return TrajectoryRequirements(rationale="sasa handles PBC and whole molecules itself")

    def build_command(self, ctx, resolved, params, out):
        x = self.extra_parameters(params)
        return (["sasa", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
                 "-n", ctx.semantic_index.path, "-surface", _g(resolved["surface"]),
                 "-output", _g(resolved["selection"]), "-probe", str(x["probe_nm"]),
                 "-ndots", str(x["ndots"]), "-o", str(out), "-tu", "ps"], None)


class _PairObservable(_GenericObservable):
    purpose = ObservablePurpose.INTER_COMPONENT_GEOMETRY
    selection_params = ("selection_a", "selection_b")

    def trajectory_requirements(self, ctx_params):
        # minimum-image geometry on raw coordinates — never a reconstructed trajectory
        return TrajectoryRequirements(minimum_image_distances=True,
                                      rationale=f"gmx {self.tool} applies minimum-image geometry")


class ComDistance(_PairObservable):
    id = "com-distance"
    display_name = "COM distance"
    description = "Distance between the mass-weighted centres of two selections (gmx distance)."
    category = "interaction"
    tool = "distance"
    quantity = "com_distance"
    expected_unit = "nm"
    coordinate_semantics = {"metric": "centre of mass (.tpr masses)",
                            "molecules": "made whole in memory (-rmpbc)",
                            "pbc": "minimum-image distance vector (-pbc)",
                            "output_resolution_nm": 0.001}

    def build_command(self, ctx, resolved, params, out):
        a, b = resolved["selection_a"], resolved["selection_b"]
        return (["distance", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
                 "-n", ctx.semantic_index.path,
                 "-select", f'com of group "{a.group}" plus com of group "{b.group}"',
                 "-oall", str(out), "-tu", "ps"], None)


class MinimumDistance(_PairObservable):
    id = "min-distance"
    display_name = "Minimum distance"
    description = "Minimum interatomic distance between two selections (gmx mindist)."
    category = "interaction"
    tool = "mindist"
    quantity = "minimum_distance"
    expected_unit = "nm"
    coordinate_semantics = {"metric": "minimum interatomic distance",
                            "pbc": "per 'pbc' parameter: minimum image over periodic images "
                                   "(-pbc) or plain stored coordinates (-nopbc)"}

    @staticmethod
    def _pbc(params) -> bool:
        v = params.get("pbc", True)
        return v if isinstance(v, bool) else str(v).lower() not in ("0", "false", "no")

    def extra_parameters(self, params):
        return {"pbc_minimum_image": self._pbc(params)}

    def trajectory_requirements(self, ctx_params):
        if not self._pbc(_params(ctx_params, self.id)):
            return TrajectoryRequirements(rationale="mindist -nopbc: stored coordinates, no "
                                                    "minimum image")
        return super().trajectory_requirements(ctx_params)

    def build_command(self, ctx, resolved, params, out):
        a, b = resolved["selection_a"], resolved["selection_b"]
        return (["mindist", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
                 "-n", ctx.semantic_index.path, "-od", str(out), "-tu", "ps",
                 "-pbc" if self._pbc(params) else "-nopbc"],
                f"{a.group}\n{b.group}\n")


class HydrogenBondCount(_PairObservable):
    id = "hbond-count"
    display_name = "Hydrogen-bond count"
    description = ("Number of hydrogen bonds between a reference and a target selection "
                   "(gmx hbond, geometric criteria recorded explicitly).")
    category = "interaction"
    tool = "hbond"
    quantity = "hydrogen_bond_count"
    expected_unit = "count"
    coordinate_semantics = {"molecules": "made whole in memory (-rmpbc)",
                            "pbc": "minimum-image donor/acceptor geometry (-pbc)"}

    def extra_parameters(self, params):
        return {"distance_cutoff_nm": float(params.get("distance_cutoff", 0.35)),
                "angle_cutoff_deg": float(params.get("angle_cutoff", 30.0)),
                "donor_elements": str(params.get("donor_elements", "N O")),
                "acceptor_elements": str(params.get("acceptor_elements", "N O")),
                "merge_donors (-m)": False, "per_frame_file (-pf)": False}

    def unit_fallback(self):
        return "count", "fallback: hydrogen-bond count is dimensionless (header has no unit)"

    def build_command(self, ctx, resolved, params, out):
        a, b = resolved["selection_a"], resolved["selection_b"]
        x = self.extra_parameters(params)
        return (["hbond", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
                 "-n", ctx.semantic_index.path, "-r", _g(a), "-t", _g(b),
                 "-hbr", str(x["distance_cutoff_nm"]), "-hba", str(x["angle_cutoff_deg"]),
                 "-cutoff", str(x["distance_cutoff_nm"]),
                 # multi-value options: one argv token per element ("N O" as a single
                 # token would be read as one element literally named "N O")
                 "-de", *x["donor_elements"].split(), "-ae", *x["acceptor_elements"].split(),
                 "-num", str(out), "-o", str(out.with_name("hbond_pairs.ndx")),
                 "-tu", "ps"], None)


BUILTIN_GENERIC = (RadiusOfGyration, SolventAccessibleSurface, ComDistance, MinimumDistance,
                   HydrogenBondCount)


def register_builtins() -> None:
    from analysis.campaign.observables.registry import is_registered, register
    for cls in BUILTIN_GENERIC:
        obs = cls()
        if not is_registered(obs.id):
            register(obs)
