"""Study discovery: a directory tree of MD artifacts -> candidate systems.

Association uses progressively weaker evidence (spec §11).  **File names never
define scientific identity** — they are, at most, the weakest tie-breaker.
Directory structure is *evidence*, recorded as such, not truth.

Output is a list of :class:`SystemCandidate` — the raw material the manifest
builder turns into :class:`SystemRecord` objects after component inference and
trajectory inspection.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.fingerprint import fingerprint_map
from analysis.campaign.models import (
    CampaignWarning, Severity, SourceFingerprint, TrajectoryArtifact, TrajectoryStage,
)
from analysis.campaign.trajectory.stage_classifier import classify_trajectory_stage

_TRAJ_SUFFIXES = {".xtc", ".trr"}
_TOPOLOGY_SUFFIXES = {".tpr", ".top"}
_STRUCTURE_SUFFIXES = {".gro", ".pdb"}
_ALL_SUFFIXES = _TRAJ_SUFFIXES | _TOPOLOGY_SUFFIXES | _STRUCTURE_SUFFIXES | {
    ".ndx", ".cpt", ".edr", ".log",
}
_SKIP_DIR_NAMES = {
    "simforge_analysis", ".git", "__pycache__", ".ipynb_checkpoints", "venv",
}
_DERIVED_TRAJ_HINTS = (
    "nojump", "nopbc", "no_pbc", "_pbc", "-pbc", "whole", "center", "centered",
    "_fit", "-fit", "fitted", "_mol", "dt", "skip", "aligned", "wrapped",
)
_NONPRODUCTION_WORD_HINTS = ("em", "min", "nvt", "npt", "eq", "heat", "ion")
_NONPRODUCTION_SUBSTR_HINTS = ("minimiz", "equilibrat", "embedding", "anneal",
                               "solvat", "prepar", "restraint", "posre", "warmup")
_STEP_PREFIX_RE = re.compile(r"^\d{1,3}[_\-](.+)$")   # simforge "11_production_md"
_WATER_MODELS = {
    "tip3p": "TIP3P", "tip3": "TIP3P", "spce": "SPC/E", "spc": "SPC",
    "tip4p": "TIP4P", "tip4pew": "TIP4P-Ew", "tip5p": "TIP5P",
    "opc": "OPC", "opc3": "OPC3", "tips3p": "TIPS3P",
}
_REPLICATE_RE = re.compile(r"^(?:rep|replica|r|run|sim|seed|traj)[\-_]?(\d+)$", re.I)
_BARE_INT_RE = re.compile(r"^(\d{1,3})$")
_SEGMENT_PART_RE = re.compile(r"(?:part|seg|chunk)[\-_]?(\d+)", re.I)
_SEGMENT_TAIL_RE = re.compile(r"[._-](\d{1,4})$")


@dataclass
class SystemCandidate:
    sim_dir: str
    production_trajectory_paths: list[str] = field(default_factory=list)
    trajectory_artifacts: list[TrajectoryArtifact] = field(default_factory=list)
    derived_trajectory_paths: list[str] = field(default_factory=list)
    topology_path: Optional[str] = None
    structure_path: Optional[str] = None
    index_path: Optional[str] = None
    fingerprints: dict[str, SourceFingerprint] = field(default_factory=dict)
    path_segments: list[str] = field(default_factory=list)
    replicate_id: Optional[str] = None
    water_model: Optional[str] = None
    partner_hint: Optional[str] = None
    condition_dimensions: dict[str, str] = field(default_factory=dict)
    discovery_evidence: list[str] = field(default_factory=list)
    association_evidence: list[str] = field(default_factory=list)
    warnings: list[CampaignWarning] = field(default_factory=list)
    simforge_run_metadata: Optional[dict] = None

    @property
    def workflow_trajectory_paths(self) -> list[str]:
        return [a.path for a in self.trajectory_artifacts
                if a.stage in TrajectoryStage.NON_PRODUCTION
                and a.stage != TrajectoryStage.DERIVED]


@dataclass
class DiscoveryResult:
    study_root: str
    candidates: list[SystemCandidate]
    unassigned_files: list[str]
    warnings: list[CampaignWarning]
    settings: dict


def _iter_artifact_files(root: Path):
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if any(part in _SKIP_DIR_NAMES for part in p.relative_to(root).parts[:-1]):
            continue
        if p.suffix.lower() in _ALL_SUFFIXES:
            yield p


def _stage_classify_dir(paths: list[Path]) -> list[TrajectoryArtifact]:
    """Classify every trajectory in one directory by MD-workflow stage.

    Two passes: first classify each file, then use directory-relative frame
    spacing / duration (from paired mdps) to refine borderline cases.
    """
    first: list[tuple[Path, object]] = []
    for p in sorted(paths, key=lambda q: q.name):
        first.append((p, classify_trajectory_stage(p)))

    dts = [r.mdp_dt_out_ps for _p, r in first if r.mdp_dt_out_ps]
    durs = [r.mdp_duration_ps for _p, r in first if r.mdp_duration_ps]

    artifacts: list[TrajectoryArtifact] = []
    for p, r0 in first:
        r = classify_trajectory_stage(
            p,
            dt_ps=r0.mdp_dt_out_ps,
            duration_ps=r0.mdp_duration_ps,
            sibling_dts=[d for d in dts if d != r0.mdp_dt_out_ps],
            sibling_durations=[d for d in durs if d != r0.mdp_duration_ps],
        )
        artifacts.append(TrajectoryArtifact(
            path=str(p), stage=r.stage, stage_confidence=r.confidence,
            stage_evidence=r.evidence, part_index=r.part_index,
            dt_ps=r.mdp_dt_out_ps,
            end_time_ps=r.mdp_duration_ps,
        ))
    return artifacts


def _pair_stage_file(
    production_stems: list[str], hits: list[tuple[Path, str]], want_production: bool,
) -> tuple[Optional[Path], str, Optional[CampaignWarning]]:
    """Pick a topology/structure for the production trajectory by stage stem match.

    Returns (path, evidence, optional warning).  ``review_required`` (via a
    REVIEW-severity warning) is emitted rather than an arbitrary pick when
    multiple production-stage candidates remain.
    """
    if not hits:
        return None, "", None
    same_dir = [p for p, e in hits if e == "same_directory"] or [p for p, _e in hits]

    # 1. exact stem match with a production trajectory stem
    for stem in production_stems:
        exact = [p for p in same_dir if p.stem.lower() == stem]
        if len(exact) == 1:
            return exact[0], f"stage-matched stem '{stem}'", None

    # 2. classify each candidate; keep production-stage ones
    from analysis.campaign.trajectory.stage_classifier import classify_trajectory_stage
    staged = [(p, classify_trajectory_stage(p)) for p in same_dir]
    prod = [p for p, r in staged if r.stage == TrajectoryStage.PRODUCTION]
    nonprod_stages = {TrajectoryStage.MINIMIZATION, TrajectoryStage.PREPARATION,
                      *TrajectoryStage.EQUILIBRATION}
    if want_production:
        if len(prod) == 1:
            return prod[0], "only production-stage candidate", None
        if len(prod) > 1:
            # tie-break by preferred stem, but flag for review
            pref = next((p for p in prod if p.stem.lower() in ("md", "prod", "production")), prod[0])
            return pref, f"multiple production-stage candidates {[p.name for p in prod]}", (
                CampaignWarning(
                    "ambiguous_production_topology",
                    f"{len(prod)} production-stage candidates: {[p.name for p in prod]}; "
                    f"picked '{pref.name}'. Set the correct one in study_manifest.yaml.",
                    "review_required",
                ))
        # no production-stage candidate: drop obviously-non-production, then prefer stem
        remaining = [p for p, r in staged if r.stage not in nonprod_stages] or same_dir
        pref = next((p for p in remaining if p.stem.lower() in ("md", "prod", "production")),
                    remaining[0])
        return pref, "no clearly production-stage candidate; stem preference", (
            CampaignWarning(
                "topology_stage_unverified",
                f"could not confirm a production-stage file among {[p.name for p in same_dir]}; "
                f"using '{pref.name}'",
                Severity.WARN,
            ))
    return same_dir[0], "first candidate", None


def _order_segments(paths: list[Path]) -> list[Path]:
    def key(p: Path):
        s = p.name
        m = _SEGMENT_PART_RE.search(s) or _SEGMENT_TAIL_RE.search(s)
        return (0, int(m.group(1))) if m else (1, s)
    return sorted(paths, key=key)


def _load_run_metadata(sim_dir: Path, root: Path) -> Optional[dict]:
    d = sim_dir
    while True:
        for name in ("metadata.json", "provenance.json", "simforge_run.json"):
            f = d / name
            if f.is_file():
                try:
                    return {"path": str(f), "data": json.loads(f.read_text())}
                except (json.JSONDecodeError, OSError):
                    pass
        if d == root or d.parent == d:
            return None
        d = d.parent


def _find_nearby(sim_dir: Path, root: Path, suffixes: set[str],
                 prefer_stems: tuple[str, ...] = ()) -> list[tuple[Path, str]]:
    """Return [(path, evidence)] for files of the given suffixes, nearest first."""
    hits: list[tuple[Path, str]] = []
    d = sim_dir
    depth = 0
    while True:
        local = sorted(p for p in d.iterdir()
                       if p.is_file() and p.suffix.lower() in suffixes)
        if local:
            ev = "same_directory" if depth == 0 else f"ancestor_directory(+{depth})"
            # stable ordering: preferred stems first, then by name
            local.sort(key=lambda p: (p.stem.lower() not in prefer_stems, p.name))
            hits.extend((p, ev) for p in local)
        if d == root or d.parent == d or depth > 6:
            break
        d = d.parent
        depth += 1
    return hits


def _condition_from_segments(cand: SystemCandidate, root: Path) -> None:
    rel = Path(cand.sim_dir).relative_to(root)
    segments = [s for s in rel.parts]
    cand.path_segments = segments
    norm_segments: list[str] = []
    for seg in segments:
        sm = _STEP_PREFIX_RE.match(seg)
        norm_segments.append(sm.group(1) if sm else seg)

    # Two-pass replicate detection: an *explicit* rep/replica/run/r<N> marker wins.
    # A bare integer directory (e.g. "300" for a temperature axis, "150" for a
    # salt concentration) is only treated as the replicate when NO explicit
    # marker exists anywhere in the path — otherwise it is a condition dimension.
    has_explicit_rep = any(_REPLICATE_RE.match(s.lower()) for s in norm_segments)

    leftover: list[str] = []
    for seg in norm_segments:
        low = seg.lower()
        rep_m = _REPLICATE_RE.match(low)
        bare_m = _BARE_INT_RE.match(low)
        if rep_m and cand.replicate_id is None:
            cand.replicate_id = f"rep{int(rep_m.group(1)):02d}"
            cand.discovery_evidence.append(f"replicate from directory segment '{seg}'")
            continue
        if bare_m and not has_explicit_rep and cand.replicate_id is None:
            cand.replicate_id = f"rep{int(bare_m.group(1)):02d}"
            cand.discovery_evidence.append(
                f"replicate assumed from bare-integer directory segment '{seg}' "
                f"(no explicit rep/replica/run marker in path)"
            )
            continue
        if low in _WATER_MODELS and cand.water_model is None:
            cand.water_model = _WATER_MODELS[low]
            cand.condition_dimensions["water_model"] = cand.water_model
            cand.discovery_evidence.append(f"water model from directory segment '{seg}'")
            continue
        leftover.append(seg)
    if leftover:
        # first leftover segment = most likely the primary comparison axis (partner)
        cand.partner_hint = leftover[0]
        cand.condition_dimensions.setdefault("partner", leftover[0])
        cand.discovery_evidence.append(
            f"partner/condition hint from directory segment(s) {leftover} "
            f"(directory semantics — not confirmed structurally)"
        )
        for i, seg in enumerate(leftover[1:], start=2):
            cand.condition_dimensions[f"dimension_{i}"] = seg
    if cand.replicate_id is None:
        cand.replicate_id = "rep01"
        cand.warnings.append(CampaignWarning(
            "no_replicate_marker",
            f"no replicate marker in path '{rel}'; assumed a single replicate (rep01) "
            f"— documented harmless assumption",
            Severity.INFO, scope=cand.sim_dir,
        ))


def _resolve_production(
    prod: list[TrajectoryArtifact],
    workflow: list[TrajectoryArtifact],
    cand: SystemCandidate,
    sim_dir: Path,
) -> list[str]:
    """Decide the production trajectory list for a directory.

    * 0 production, ≥1 workflow  → no production (system unusable for observables).
    * 1 production               → that trajectory.
    * ≥2 production              → ONE segment collection ONLY with positive
      continuation evidence (all have part indices); otherwise keep them all
      separate and emit ``ambiguous_production_trajectories`` so the downstream
      "multiple segments, no safe concatenation" fail-safe applies.
    """
    if not prod:
        # No clearly-production trajectory.  Trajectories with NO stage evidence
        # (UNKNOWN) are production candidates by elimination — an arbitrary name
        # like "final.xtc" with no mdp is not evidence of equilibration.
        unknown = [a for a in workflow if a.stage == TrajectoryStage.UNKNOWN]
        classified_nonprod = [a for a in workflow if a.stage != TrajectoryStage.UNKNOWN]
        if unknown:
            # single trajectory, no equilibration siblings -> the common real
            # case (arbitrary production filename); promote silently (INFO).
            benign = len(unknown) == 1 and not classified_nonprod
            for a in unknown:
                a.stage = TrajectoryStage.PRODUCTION
                a.stage_confidence = min(a.stage_confidence, 0.35 if benign else 0.25)
                a.stage_evidence.append(
                    "no equilibration/minimisation evidence and "
                    + ("it is the only trajectory in the directory"
                       if benign else "other trajectories are classified")
                    + "; treated as production")
            cand.warnings.append(CampaignWarning(
                "production_stage_unverified",
                f"{sim_dir.name}: production stage of "
                f"{[Path(a.path).name for a in unknown]} not positively confirmed "
                f"(no paired mdp / stage marker); proceeding as production.",
                Severity.INFO if benign else Severity.WARN, scope=str(sim_dir)))
            prod = unknown
        else:
            if classified_nonprod:
                cand.warnings.append(CampaignWarning(
                    "no_production_trajectory",
                    f"{sim_dir.name}: only equilibration/preparation trajectories found "
                    f"({[Path(a.path).name for a in classified_nonprod]}); no production "
                    f"trajectory to analyse. If one IS production, rename it or supply a manifest.",
                    Severity.WARN, scope=str(sim_dir)))
            return []

    if len(prod) == 1:
        return [prod[0].path]

    parts = [a for a in prod if a.part_index is not None]
    if len(parts) == len(prod):
        ordered = sorted(prod, key=lambda a: a.part_index)
        for i, a in enumerate(ordered):
            a.stage_evidence.append(f"continuation segment {i + 1}/{len(ordered)} "
                                    f"(part {a.part_index})")
        cand.discovery_evidence.append(
            f"{len(ordered)} production continuation segments "
            f"(part markers): {[Path(a.path).name for a in ordered]}"
        )
        return [a.path for a in ordered]

    # multiple production trajectories, no continuation markers -> do NOT merge
    names = [Path(a.path).name for a in prod]
    cand.warnings.append(CampaignWarning(
        "ambiguous_production_trajectories",
        f"{sim_dir.name}: {len(prod)} production-stage trajectories {names} with no "
        f"continuation (partNNNN) markers. Same-directory membership alone is not "
        f"evidence of one run. They will NOT be concatenated; select one in "
        f"study_manifest.yaml.",
        "review_required", scope=str(sim_dir)))
    return sorted(a.path for a in prod)


def _pair_workflow_artifacts(
    artifacts: list[TrajectoryArtifact], sim_dir: Path, root: Path,
) -> None:
    """Attach a stem-matched .tpr / structure to each artifact for provenance."""
    tprs = {p.stem.lower(): p for p in sim_dir.glob("*.tpr")}
    gros = {p.stem.lower(): p for p in sim_dir.glob("*.gro")}
    for a in artifacts:
        stem = Path(a.path).stem.lower()
        if stem in tprs:
            a.topology_path = str(tprs[stem])
        if stem in gros:
            a.structure_path = str(gros[stem])


def discover_study(
    study_root: str | Path, *, strong_fingerprint: bool = False,
) -> DiscoveryResult:
    root = Path(study_root).resolve()
    warnings: list[CampaignWarning] = []
    if not root.is_dir():
        return DiscoveryResult(str(root), [], [], [CampaignWarning(
            "study_root_missing", f"not a directory: {root}", Severity.ERROR)], {})

    all_files = list(_iter_artifact_files(root))
    traj_by_dir: dict[Path, list[Path]] = {}
    for f in all_files:
        if f.suffix.lower() in _TRAJ_SUFFIXES:
            traj_by_dir.setdefault(f.parent, []).append(f)

    if not traj_by_dir:
        hint = ""
        if list(root.rglob("*.xvg")):
            hint = (" — this directory contains .xvg files; for pre-computed XVG "
                    "comparative analysis use `simforge study <dir>` (no sub-command)")
        warnings.append(CampaignWarning(
            "no_trajectories",
            f"no .xtc/.trr trajectories found under {root}{hint}", Severity.ERROR))
        return DiscoveryResult(str(root), [], [str(f) for f in all_files], warnings,
                               {"strong_fingerprint": strong_fingerprint})

    assigned: set[Path] = set()
    candidates: list[SystemCandidate] = []

    for sim_dir in sorted(traj_by_dir):
        trajs = traj_by_dir[sim_dir]
        artifacts = _stage_classify_dir(trajs)
        for a in trajs:
            assigned.add(a)

        prod = [a for a in artifacts if a.stage == TrajectoryStage.PRODUCTION]
        non_prod = [a for a in artifacts if a.stage in TrajectoryStage.EQUILIBRATION
                    or a.stage in (TrajectoryStage.MINIMIZATION, TrajectoryStage.PREPARATION,
                                   TrajectoryStage.UNKNOWN)]

        cand = SystemCandidate(sim_dir=str(sim_dir), trajectory_artifacts=artifacts)

        # ── resolve the production trajectory / segment collection ─────────
        prod_paths = _resolve_production(prod, non_prod, cand, sim_dir)
        cand.production_trajectory_paths = prod_paths
        prod_stems = [Path(p).stem.lower() for p in prod_paths]

        # recompute after _resolve_production may have promoted UNKNOWN -> PRODUCTION
        derived = [a for a in artifacts if a.stage == TrajectoryStage.DERIVED]
        workflow = [a for a in artifacts
                    if a.stage in TrajectoryStage.NON_PRODUCTION
                    and a.stage != TrajectoryStage.DERIVED]
        cand.derived_trajectory_paths = [a.path for a in derived]

        for a in derived:
            a.derived_from = prod_paths[0] if prod_paths else None
        if derived:
            cand.discovery_evidence.append(
                f"excluded {len(derived)} pre-processed trajectory file(s) "
                f"(PBC/fit derivatives): {[Path(a.path).name for a in derived]}"
            )
        wf_by_stage = {}
        for a in workflow:
            wf_by_stage.setdefault(a.stage, []).append(Path(a.path).name)
        if wf_by_stage:
            cand.discovery_evidence.append(
                f"workflow (non-production) trajectories classified: {wf_by_stage}"
            )

        # ── SimForge run metadata (strongest) ──────────────────────────────
        meta = _load_run_metadata(sim_dir, root)
        if meta:
            cand.simforge_run_metadata = meta
            cand.association_evidence.append(f"simforge_run_metadata:{Path(meta['path']).name}")

        # ── stage-matched topology (production) ────────────────────────────
        tpr_hits = _find_nearby(sim_dir, root, {".tpr"}, prefer_stems=tuple(prod_stems) or ("md",))
        top_hits = _find_nearby(sim_dir, root, {".top"}, prefer_stems=("topol", "system"))
        if tpr_hits:
            best, ev, warn = _pair_stage_file(prod_stems, tpr_hits, want_production=bool(prod_paths))
            if best is not None:
                cand.topology_path = str(best)
                cand.association_evidence.append(f"topology(.tpr):{ev}")
            if warn is not None:
                warn.scope = str(sim_dir)
                cand.warnings.append(warn)
        elif top_hits:
            best, ev = top_hits[0]
            cand.topology_path = str(best)
            cand.association_evidence.append(f"topology(.top):{ev}")
        else:
            cand.warnings.append(CampaignWarning(
                "no_topology",
                f"no .tpr or .top for {sim_dir.name}; RMSD-type analyses need a topology",
                Severity.WARN, scope=str(sim_dir)))

        # ── stage-matched reference structure (production) ─────────────────
        struct_hits = _find_nearby(sim_dir, root, _STRUCTURE_SUFFIXES,
                                   prefer_stems=tuple(prod_stems) or ("md", "confout"))
        if struct_hits:
            best, ev, warn = _pair_stage_file(prod_stems, struct_hits, want_production=bool(prod_paths))
            if best is not None:
                cand.structure_path = str(best)
                cand.association_evidence.append(f"structure({Path(best).suffix}):{ev}")
        else:
            cand.warnings.append(CampaignWarning(
                "no_structure",
                f"no .gro/.pdb reference structure for {sim_dir.name}; semantic analysis limited",
                Severity.WARN, scope=str(sim_dir)))

        # ── per-artifact stage-matched topology/structure (for provenance) ─
        _pair_workflow_artifacts(cand.trajectory_artifacts, sim_dir, root)

        # ── index ─────────────────────────────────────────────────────────
        ndx_hits = _find_nearby(sim_dir, root, {".ndx"})
        if ndx_hits:
            cand.index_path = str(ndx_hits[0][0])
            cand.association_evidence.append(f"user_index:{ndx_hits[0][1]}")

        # ── condition dimensions from directory structure ─────────────────
        _condition_from_segments(cand, root)

        # ── fingerprints (production sources only) ────────────────────────
        fp_targets: dict[str, Path] = {}
        for i, p in enumerate(prod_paths):
            fp_targets[f"trajectory_{i}" if len(prod_paths) > 1 else "trajectory"] = Path(p)
        if cand.topology_path:
            fp_targets["topology"] = Path(cand.topology_path)
        if cand.structure_path:
            fp_targets["structure"] = Path(cand.structure_path)
        cand.fingerprints = fingerprint_map(fp_targets, strong=strong_fingerprint)
        for a in cand.trajectory_artifacts:
            if Path(a.path).is_file():
                from analysis.campaign.fingerprint import fingerprint_file
                a.fingerprint = fingerprint_file(Path(a.path), strong=strong_fingerprint)

        candidates.append(cand)

    # ── unassigned-file inventory (invariant: assigned ∩ unassigned == ∅) ──
    for cand in candidates:
        for p in (cand.topology_path, cand.structure_path, cand.index_path):
            if p:
                assigned.add(Path(p))
        for a in cand.trajectory_artifacts:
            if a.topology_path:
                assigned.add(Path(a.topology_path))
            if a.structure_path:
                assigned.add(Path(a.structure_path))
    tracked = _TRAJ_SUFFIXES | _TOPOLOGY_SUFFIXES | _STRUCTURE_SUFFIXES
    unassigned = sorted(str(f) for f in all_files
                        if f not in assigned and f.suffix.lower() in tracked)

    return DiscoveryResult(
        study_root=str(root),
        candidates=candidates,
        unassigned_files=unassigned,
        warnings=warnings,
        settings={"strong_fingerprint": strong_fingerprint},
    )
