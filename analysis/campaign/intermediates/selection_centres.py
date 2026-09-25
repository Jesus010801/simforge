"""``selection-centre-series`` — the centre of one resolved selection per frame.

    selection   component:<type> | annotation:<id>   (Phase 7 resolver: only RESOLVED
                                                       components / ACTIVE annotations)
    convention  com | cog                             (required; COM needs .tpr masses —
                                                       never silently replaced by COG)
    pbc         whole (default) | stored              whole: molecules made whole in
                                                       memory from .tpr bonds (-rmpbc);
                                                       stored: stored coordinates (-normpbc)

Output ``[time, coordinate]`` in nm, box coordinates, one row per frame of the
raw view — verified against the view's FrameTimeIndex (duplicate timestamps
kept).  Stored as ``.npy`` (deterministic bytes, so a republish can be checked
for reproducibility).  Computed with ``gmx trajectory``.

Relation to the Phase 4 diagnostic probes: those run ``gmx trajectory -nopbc``
with the default ``-rmpbc`` (with a .tpr: whole molecules) and keep their own
cache keyed by role; with a .tpr and COM the numbers are the same definition
as ``pbc=whole``, but the probes are not rerouted through this store (their
cache/key scheme and failure handling belong to diagnostics).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from analysis.campaign.intermediates.base import (
    IntermediateContext, IntermediateScope, IntermediateSpec,
)
from analysis.campaign.intermediates.registry import register_intermediate
from analysis.campaign.models import (
    Axis, AxisKind, FrameAlignmentMode, MissingSemantics, ObservablePurpose, ResultArray,
    StorageRef, TrajectoryRequirements,
)

CONVENTIONS = ("com", "cog")
PBC_MODES = {"whole": ["-rmpbc", "-nopbc"], "stored": ["-normpbc", "-nopbc"]}


def _assembly(view):
    from analysis.campaign.trajectory.assembly import assembly_evidence
    return assembly_evidence(view)


def _tpr(path) -> bool:
    return bool(path) and str(path).lower().endswith(".tpr") and Path(path).is_file()


class SelectionCentreSeries(IntermediateSpec):
    id = "selection-centre-series"
    version = 1
    purpose = ObservablePurpose.INTER_COMPONENT_GEOMETRY
    scope = IntermediateScope.VIEW
    description = "Per-frame centre (COM or COG) of one resolved selection, in box coordinates."

    def parameters_schema(self) -> dict:
        return {"selection": "component:<type> | annotation:<id>",
                "convention": "com | cog (required)", "pbc": "whole (default) | stored"}

    def validate(self, params: dict) -> list[str]:
        from analysis.campaign.observables.selections import SelectionRef
        errors = super().validate(params)
        try:
            SelectionRef.parse(params.get("selection", ""))
        except ValueError as exc:
            errors.append(str(exc))
        if params.get("convention") not in CONVENTIONS:
            errors.append(f"convention must be one of {CONVENTIONS} (no default)")
        if params.get("pbc", "whole") not in PBC_MODES:
            errors.append(f"pbc must be one of {sorted(PBC_MODES)}")
        return errors

    def _resolve(self, system, params, semantic_index):
        from analysis.campaign.observables.selections import SelectionRef, resolve_selection
        return resolve_selection(SelectionRef.parse(params["selection"]), system, semantic_index)

    def assembly_groups(self, system, params, semantic_index) -> tuple:
        # a centre is a global property of the selection: all its molecules must
        # share one periodic image (single-molecule selections are unaffected)
        r = self._resolve(system, params, semantic_index)
        return (r.group,) if r.ok and r.group else ()

    def applicability(self, system, params, semantic_index) -> dict:
        reasons = []
        r = self._resolve(system, params, semantic_index)
        if not r.ok:
            reasons.append(f"selection {params['selection']}: {r.reason}")
        needs_tpr = params["convention"] == "com" or params.get("pbc", "whole") == "whole"
        if needs_tpr and not _tpr(system.topology_path):
            reasons.append("COM masses / whole-molecule connectivity need a .tpr; none available "
                           "(use convention=cog, pbc=stored explicitly if that is what you mean)")
        return {"applicable": not reasons, "reasons": reasons}

    def trajectory_requirements(self, params) -> TrajectoryRequirements:
        return TrajectoryRequirements(
            rationale="centres of the raw coordinates; whole molecules in memory (-rmpbc) "
                      "when pbc=whole — no trajectory transformation is written")

    def output_schema(self, params) -> list[ResultArray]:
        return [ResultArray("selection_center", "selection_center", "nm",
                            axes=[Axis("time", AxisKind.TIME, "ps"),
                                  Axis("coordinate", AxisKind.COORDINATE, "nm",
                                       values=["x", "y", "z"], frame="box")],
                            missing=MissingSemantics.NONE_EXPECTED)]

    def definition_evidence(self, ctx: IntermediateContext):
        from analysis.campaign.fingerprint import fingerprint_file
        from analysis.campaign.gmx import gmx_version
        p = ctx.parameters
        r = self._resolve(ctx.system, p, ctx.semantic_index)
        if not r.ok:
            return None
        pbc = p.get("pbc", "whole")
        tpr = ctx.topology_path if _tpr(ctx.topology_path) else None
        return {"selection": r.evidence, "convention": p["convention"], "pbc": pbc,
                "gmx_flags": PBC_MODES[pbc],
                "mass_source": "tpr masses" if p["convention"] == "com" else
                               "none (geometric centre)",
                "topology": fingerprint_file(Path(tpr)).digest if tpr else None,
                "coordinate_frame": "box",
                "backend": {"tool": "gmx trajectory", "version": gmx_version(ctx.gmx)}} | (
                    {"assembly": ev} if (ev := _assembly(ctx.trajectory_view)) else {})

    def compute(self, ctx: IntermediateContext):
        from analysis.campaign.gmx import run_gmx
        from analysis.campaign.observables.generic import time_alignment
        p = ctx.parameters
        r = self._resolve(ctx.system, p, ctx.semantic_index)
        kw = "com" if p["convention"] == "com" else "cog"
        out = ctx.work_dir / "centres.xvg"                   # gmx needs the .xvg extension
        view = ctx.trajectory_view
        argv = ["trajectory", "-f", str(view.path),
                "-s", str(ctx.topology_path if _tpr(ctx.topology_path) else ctx.structure_path),
                "-n", str(ctx.semantic_index.path), *PBC_MODES[p.get("pbc", "whole")],
                "-tu", "ps", "-ox", str(out), "-select", f'{kw} of group "{r.group}"']
        res = run_gmx(argv, gmx=ctx.gmx, timeout=4 * 3600)
        if not res.ok or not out.is_file():
            raise RuntimeError(f"gmx trajectory failed (rc={res.returncode}): "
                               f"{res.stderr.strip()[-400:]}")
        rows = [[float(v) for v in line.split()] for line in out.read_text().splitlines()
                if line.strip() and line.lstrip()[0] not in "#@"]
        out.unlink()                                   # only the canonical .npy files persist
        if not rows or any(len(x) != 4 for x in rows):
            raise RuntimeError("gmx trajectory returned an unexpected column layout")
        data = np.asarray(rows, dtype=np.float64)
        times = data[:, 0].tolist()
        alignment, attrs = time_alignment(ctx, times)
        if alignment.mode != FrameAlignmentMode.ONE_ROW_PER_FRAME:
            raise RuntimeError("centre series does not align one row per frame with the "
                               f"view's time index: {attrs.get('alignment_check')}")
        tfile, cfile = ctx.work_dir / "time.npy", ctx.work_dir / "center.npy"
        np.save(tfile, data[:, 0], allow_pickle=False)
        np.save(cfile, data[:, 1:4], allow_pickle=False)
        arr = ResultArray(
            "selection_center", "selection_center", "nm",
            axes=[Axis("time", AxisKind.TIME, "ps", values_ref=StorageRef(str(tfile), "npy"),
                       alignment=alignment, attrs=attrs),
                  Axis("coordinate", AxisKind.COORDINATE, "nm", values=["x", "y", "z"],
                       frame="box")],
            storage=StorageRef(str(cfile), "npy"), view_ref=view.cache_key,
            missing=MissingSemantics.NONE_EXPECTED,
            attrs={"convention": p["convention"], "pbc": p.get("pbc", "whole"),
                   "group": r.group, "selection": p["selection"], "shape": list(data[:, 1:4].shape)})
        return [arr], [], {"command": ["gmx", *argv[:-1], argv[-1]],
                           "rows": len(rows), "note": "time.npy / center.npy are the canonical "
                                                      "artifacts; the gmx .xvg is not kept"}


register_intermediate(SelectionCentreSeries())
