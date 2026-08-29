"""Resolve :class:`TrajectoryRequirements` into a concrete :class:`TrajectoryView`.

One derived trajectory per *distinct* requirement set, cached by a key built
from the source fingerprints + the requirement token + the semantic selections
used.  Different observables therefore share a view only when it is
scientifically identical.

Preprocessing pipeline (only the steps a requirement asks for):

1. ``-pbc whole``   (make_whole)   — molecules made whole
2. ``-pbc nojump``  (nojump)       — unwrapped, diffusion-safe
3. ``-center``      (center)       — a semantic group centred in the box
4. ``-fit rot+trans`` (fit)        — least-squares fit to a semantic group

Each step is a separate ``gmx trjconv`` call fed the relevant *named* group
from the semantic index (never a numeric id), and every call is recorded in a
:class:`PreprocessingOperation`.  If a required semantic group is missing the
view is marked ``safe=False`` and the caller must fail the analysis rather than
fall back to raw coordinates.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

from analysis.campaign.gmx import gmx_available, run_gmx
from analysis.campaign.models import (
    CampaignWarning, PreprocessingOperation, SemanticIndex, Severity,
    SourceFingerprint, TrajectoryRequirements, TrajectoryView, ViewKind,
)


def _cache_key(
    fingerprints: list[SourceFingerprint],
    topology_fp: Optional[SourceFingerprint],
    req: TrajectoryRequirements,
) -> str:
    parts = [fp.digest for fp in fingerprints]
    if topology_fp:
        parts.append(topology_fp.digest)
    parts.append(req.cache_token())
    return hashlib.sha256("\x00".join(parts).encode()).hexdigest()[:16]


def _resolve_group(index: Optional[SemanticIndex], name: str) -> Optional[str]:
    if not index:
        return None
    return name if name in index.group_names() else None


def build_view(
    *,
    requirements: TrajectoryRequirements,
    trajectory_paths: list[str],
    topology_path: str,
    structure_path: Optional[str],
    semantic_index: Optional[SemanticIndex],
    source_fingerprints: list[SourceFingerprint],
    topology_fingerprint: Optional[SourceFingerprint],
    work_dir: Path,
    gmx: str = "gmx",
    force: bool = False,
    dry_run: bool = False,
) -> TrajectoryView:
    kind = requirements.view_kind()
    view = TrajectoryView(
        kind=kind, requirements=requirements, topology_path=topology_path,
        source_fingerprints=source_fingerprints,
    )
    view.cache_key = _cache_key(source_fingerprints, topology_fingerprint, requirements)

    if kind == ViewKind.RAW:
        if len(trajectory_paths) == 1:
            view.path = trajectory_paths[0]
        view.safe = len(trajectory_paths) == 1
        if not view.safe:
            view.warnings.append(CampaignWarning(
                "multi_segment_raw",
                "multiple trajectory segments and no preprocessing requested; "
                "an observable must handle segments explicitly",
                Severity.WARN,
            ))
        return view

    if not gmx_available(gmx):
        view.safe = False
        view.warnings.append(CampaignWarning(
            "gmx_unavailable", f"'{gmx}' not available; cannot build the "
            f"'{kind}' trajectory view", Severity.ERROR))
        return view

    if len(trajectory_paths) != 1:
        view.safe = False
        view.warnings.append(CampaignWarning(
            "segments_not_concatenated",
            "trajectory has multiple segments; automatic concatenation is not "
            "performed. Provide a single trajectory or a validated concatenation.",
            Severity.ERROR,
        ))
        return view

    src = trajectory_paths[0]
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = work_dir / view.cache_key
    final_out = cache_dir / f"{kind}.xtc"

    if final_out.is_file() and not force and not dry_run:
        view.path = str(final_out.resolve())
        view.reused = True
        return view

    ref_for_center = structure_path or topology_path
    current = src
    steps: list[tuple[str, list[str], Optional[str], str, str]] = []

    if requirements.requires_whole_molecules:
        out = cache_dir / "whole.xtc"
        steps.append((
            "make_whole",
            ["trjconv", "-s", topology_path, "-f", current, "-o", str(out),
             "-pbc", "whole"],
            "System",
            "make molecules whole across periodic boundaries",
            "required before any fitting or geometric measurement",
        ))
        current = str(out)

    if requirements.requires_nojump:
        out = cache_dir / "nojump.xtc"
        steps.append((
            "nojump",
            ["trjconv", "-s", topology_path, "-f", current, "-o", str(out),
             "-pbc", "nojump"],
            "System",
            "remove periodic jumps (unwrap) — preserves diffusion",
            "diffusion/MSD-safe unwrapping; no spatial fitting applied",
        ))
        current = str(out)

    if requirements.centering_target:
        grp = _resolve_group(semantic_index, requirements.centering_target)
        if grp is None:
            view.safe = False
            view.warnings.append(CampaignWarning(
                "missing_semantic_group",
                f"centering target '{requirements.centering_target}' is not in the "
                f"semantic index; refusing to guess",
                Severity.ERROR,
            ))
            return view
        out = cache_dir / "centered.xtc"
        steps.append((
            "center",
            ["trjconv", "-s", ref_for_center, "-f", current, "-o", str(out),
             "-pbc", "mol", "-center"],
            f"{grp}\nSystem",
            f"centre '{grp}' in the box",
            "centering only; no rotational fitting",
        ))
        current = str(out)

    if requirements.fit_selection:
        grp = _resolve_group(semantic_index, requirements.fit_selection)
        if grp is None:
            view.safe = False
            view.warnings.append(CampaignWarning(
                "missing_semantic_group",
                f"fit selection '{requirements.fit_selection}' is not in the "
                f"semantic index; refusing to guess",
                Severity.ERROR,
            ))
            return view
        out = cache_dir / "fitted.xtc"
        steps.append((
            "fit",
            ["trjconv", "-s", ref_for_center, "-f", current, "-o", str(out),
             "-fit", "rot+trans"],
            f"{grp}\nSystem",
            f"least-squares rot+trans fit to '{grp}'",
            f"reference frame = frame 0 of '{grp}'",
        ))
        current = str(out)

    index_args: list[str] = []
    if semantic_index and semantic_index.path:
        index_args = ["-n", semantic_index.path]

    if dry_run:
        for op_name, args, stdin, reason, validation in steps:
            full_args = list(args)
            if index_args and op_name in ("center", "fit"):
                full_args += index_args
            view.operations.append(PreprocessingOperation(
                operation=op_name, tool="gmx trjconv", command=["gmx", *full_args],
                stdin=stdin, semantic_selection=stdin.split("\n")[0] if stdin else None,
                input_trajectory=args[args.index("-f") + 1],
                input_topology=args[args.index("-s") + 1],
                output=args[args.index("-o") + 1],
                reason=reason, validation_note=validation,
                returncode=None, stdout_tail="", stderr_tail="(dry run — not executed)",
            ))
        view.path = str(final_out.resolve())   # planned location; not created
        view.warnings.append(CampaignWarning(
            "dry_run", "trajectory preprocessing planned but not executed", Severity.INFO))
        return view

    cache_dir.mkdir(parents=True, exist_ok=True)
    for op_name, args, stdin, reason, validation in steps:
        full_args = list(args)
        if index_args and op_name in ("center", "fit"):
            full_args += index_args
        res = run_gmx(full_args, stdin=(stdin + "\n") if stdin else None,
                      gmx=gmx, timeout=3600)
        op = PreprocessingOperation(
            operation=op_name, tool="gmx trjconv", command=res.argv,
            stdin=stdin, semantic_selection=stdin.split("\n")[0] if stdin else None,
            input_trajectory=args[args.index("-f") + 1],
            input_topology=args[args.index("-s") + 1],
            output=args[args.index("-o") + 1],
            reason=reason, validation_note=validation,
            returncode=res.returncode,
            stdout_tail=res.tail("stdout", 800), stderr_tail=res.tail("stderr", 1500),
        )
        view.operations.append(op)
        if not res.ok:
            view.safe = False
            view.warnings.append(CampaignWarning(
                "preprocessing_failed",
                f"trajectory preprocessing step '{op_name}' failed "
                f"(rc={res.returncode}); analysis must not proceed on raw coordinates",
                Severity.ERROR,
            ))
            return view

    if not Path(current).is_file():
        view.safe = False
        view.warnings.append(CampaignWarning(
            "preprocessing_no_output", "expected derived trajectory was not produced",
            Severity.ERROR))
        return view

    Path(current).replace(final_out)
    view.path = str(final_out.resolve())
    return view
