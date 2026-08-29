"""Build / serialise / reload the :class:`StudyManifest`.

``build_manifest`` runs discovery, then per candidate: component inference,
trajectory inspection, canonical-id assignment and condition grouping.  The
manifest is written as both YAML and JSON and can be fed back in with
``load_manifest`` — including a hand-corrected one (ambiguities resolved by the
user are honoured, not rediscovered).
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from typing import Optional

import yaml

from analysis.campaign.canonical_ids import (
    canonical_system_id, deduplicate_ids, sanitize_token,
)
from analysis.campaign.discovery import DiscoveryResult, SystemCandidate, discover_study
from analysis.campaign.models import (
    Ambiguity, CampaignWarning, ClassificationState, ComponentType, ConditionRecord,
    MANIFEST_SCHEMA_VERSION, Severity, SourceFingerprint, StudyManifest, SystemRecord,
)
from analysis.campaign.structure.component_detector import detect_components


def _simforge_version() -> str:
    try:
        from simforge import __version__
        return __version__
    except Exception:  # noqa: BLE001
        return "unknown"


def _study_receptor_hint(root: Path) -> Optional[str]:
    name = root.name.strip()
    generic = {"", ".", "..", "study", "studies", "data", "runs", "md", "analysis",
               "simulations", "sims", "production", "results"}
    if name.lower() in generic or len(name) < 2:
        return None
    return name


def _candidate_to_record(
    cand: SystemCandidate,
    receptor_hint: Optional[str],
    *,
    inspect_trajectories: bool,
    gmx: str,
) -> tuple[SystemRecord, list[Ambiguity]]:
    ambiguities: list[Ambiguity] = []

    rec = SystemRecord(
        system_id="(pending)",
        condition_id="(pending)",
        replicate_id=cand.replicate_id or "rep01",
        receptor=receptor_hint,
        partner=cand.partner_hint,
        water_model=cand.water_model,
        condition_dimensions=dict(cand.condition_dimensions),
        production_trajectory_paths=list(cand.production_trajectory_paths),
        trajectory_artifacts=list(cand.trajectory_artifacts),
        topology_path=cand.topology_path,
        structure_path=cand.structure_path,
        index_path=cand.index_path,
        source_fingerprints=dict(cand.fingerprints),
        discovery_evidence=list(cand.discovery_evidence),
        association_evidence=list(cand.association_evidence),
        warnings=list(cand.warnings),
    )
    if receptor_hint:
        rec.discovery_evidence.append(
            f"study label '{receptor_hint}' taken from the study-root directory name "
            f"and used as the canonical-id receptor slot — this is a STUDY LABEL, not a "
            f"verified molecular/biological receptor identity"
        )

    # ── component inference ────────────────────────────────────────────────
    det = detect_components(
        structure_path=cand.structure_path,
        topology_path=cand.topology_path,
        membrane_present_hint=None,
    )
    rec.components = det.components
    rec.membrane_present = det.membrane_present
    for w in det.warnings:
        rec.warnings.append(CampaignWarning("component_inference", w, Severity.INFO,
                                            scope=cand.sim_dir))

    # partner_type from resolved components
    peptide = det.by_type(ComponentType.PEPTIDE)
    ligand = det.by_type(ComponentType.LIGAND)
    if peptide and peptide.classification_state == ClassificationState.RESOLVED:
        rec.partner_type = "peptide"
    elif ligand:
        rec.partner_type = "ligand"

    # ── classification state roll-up ──────────────────────────────────────
    states = [c.classification_state for c in det.components
              if c.component_type in (ComponentType.RECEPTOR, ComponentType.PEPTIDE,
                                      ComponentType.COMPLEX, ComponentType.PROTEIN_PARTNER)]
    if any(s == ClassificationState.REVIEW_REQUIRED for s in states):
        rec.classification_state = ClassificationState.REVIEW_REQUIRED
    elif any(s == ClassificationState.AMBIGUOUS for s in states):
        rec.classification_state = ClassificationState.AMBIGUOUS
    elif states and all(s == ClassificationState.RESOLVED for s in states):
        rec.classification_state = ClassificationState.RESOLVED
    else:
        rec.classification_state = ClassificationState.UNKNOWN

    for comp in det.components:
        if comp.classification_state in (ClassificationState.AMBIGUOUS,
                                         ClassificationState.REVIEW_REQUIRED):
            ambiguities.append(Ambiguity(
                system_id="(pending)", kind="component_classification",
                message=f"{comp.label}: {comp.classification_state}. "
                        + "; ".join(comp.warnings or [e.detail for e in comp.evidence]),
                options=[f"chain {c}" for c in comp.chain_ids],
            ))

    # ── trajectory inspection (PRODUCTION trajectories only) ─────────────
    if inspect_trajectories and cand.production_trajectory_paths:
        from analysis.campaign.trajectory.inspector import inspect_trajectory
        insp = inspect_trajectory(
            trajectory_paths=cand.production_trajectory_paths,
            topology_path=cand.topology_path,
            structure_path=cand.structure_path,
            membrane_present=rec.membrane_present,
            gmx=gmx,
        )
        rec.trajectory_inspection = insp
        rec.warnings.extend(insp.warnings)
        # fold measured timing back into the production artifacts
        for seg in insp.segments:
            art = rec.artifact(seg.path)
            if art is not None:
                art.n_frames = seg.n_frames if seg.n_frames is not None else art.n_frames
                art.start_time_ps = seg.start_time_ps
                art.end_time_ps = seg.end_time_ps if seg.end_time_ps is not None else art.end_time_ps
                art.dt_ps = seg.dt_ps if seg.dt_ps is not None else art.dt_ps

    return rec, ambiguities


def build_manifest(
    study_root: str | Path,
    *,
    inspect_trajectories: bool = True,
    strong_fingerprint: bool = False,
    gmx: str = "gmx",
    discovery: Optional[DiscoveryResult] = None,
) -> StudyManifest:
    root = Path(study_root).resolve()
    disc = discovery or discover_study(root, strong_fingerprint=strong_fingerprint)

    receptor_hint = _study_receptor_hint(root)
    manifest = StudyManifest(
        study_root=str(root),
        generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        simforge_version=_simforge_version(),
        schema_version=MANIFEST_SCHEMA_VERSION,
        discovery_settings={
            "strong_fingerprint": strong_fingerprint,
            "inspect_trajectories": inspect_trajectories,
            "receptor_hint": receptor_hint,
        },
        unassigned_files=disc.unassigned_files,
        warnings=list(disc.warnings),
    )

    records: list[SystemRecord] = []
    all_ambiguities: list[list[Ambiguity]] = []
    for cand in disc.candidates:
        rec, ambs = _candidate_to_record(
            cand, receptor_hint,
            inspect_trajectories=inspect_trajectories, gmx=gmx,
        )
        records.append(rec)
        all_ambiguities.append(ambs)

    # ── canonical ids (deterministic, collision-safe) ─────────────────────
    id_inputs: list[tuple[str, list[str]]] = []
    for rec in records:
        extra = {k: v for k, v in rec.condition_dimensions.items()
                 if k not in ("partner", "water_model")}
        cid = canonical_system_id(
            receptor=rec.receptor, partner=rec.partner,
            water_model=rec.water_model, replicate_id=rec.replicate_id,
            extra_dimensions=extra,
        )
        tokens = [
            *(fp.digest for fp in rec.source_fingerprints.values()),
            rec.topology_path or "", rec.structure_path or "",
            *rec.trajectory_paths,
        ]
        id_inputs.append((cid, tokens))
    final_ids = deduplicate_ids(id_inputs)
    for i, rec in enumerate(records):
        rec.system_id = final_ids[i]
        for amb in all_ambiguities[i]:
            amb.system_id = rec.system_id
            manifest.ambiguities.append(amb)

    # ── condition grouping (everything except replicate) ──────────────────
    conditions: dict[str, ConditionRecord] = {}
    dims_seen: set[str] = set()
    for rec in records:
        extra = {k: v for k, v in rec.condition_dimensions.items()
                 if k not in ("partner", "water_model")}
        cond_key = canonical_system_id(
            receptor=rec.receptor, partner=rec.partner,
            water_model=rec.water_model, replicate_id=None,
            extra_dimensions=extra,
        ).replace("__rep_unknown", "")
        cond = conditions.setdefault(cond_key, ConditionRecord(
            condition_id=cond_key,
            dimensions={
                "receptor": sanitize_token(rec.receptor),
                "partner": sanitize_token(rec.partner),
                "water_model": sanitize_token(rec.water_model),
                **{k: sanitize_token(v) for k, v in extra.items()},
            },
        ))
        cond.system_ids.append(rec.system_id)
        cond.replicate_ids.append(rec.replicate_id)
        rec.condition_id = cond_key
        dims_seen.update(rec.condition_dimensions.keys())

    manifest.systems = records
    manifest.conditions = sorted(conditions.values(), key=lambda c: c.condition_id)
    manifest.condition_dimensions = sorted({"receptor", "partner", "water_model"} | dims_seen)

    # ── study-level sanity warnings ──────────────────────────────────────
    rc = manifest.stats()["replicate_counts"]
    if rc and len(set(rc.values())) > 1:
        manifest.warnings.append(CampaignWarning(
            "unbalanced_replicates",
            f"conditions have unequal replicate counts: {rc}. Not necessarily invalid "
            f"— reported for transparency.",
            Severity.INFO,
        ))

    _enforce_assignment_invariant(manifest)
    return manifest


def assigned_paths(manifest: StudyManifest) -> set[str]:
    """Every file a system claims: production + workflow trajectories, and the
    selected/paired topology, structure and index."""
    out: set[str] = set()
    for s in manifest.systems:
        out.update(s.production_trajectory_paths)
        for a in s.trajectory_artifacts:
            out.add(a.path)
            if a.topology_path:
                out.add(a.topology_path)
            if a.structure_path:
                out.add(a.structure_path)
        for p in (s.topology_path, s.structure_path, s.index_path):
            if p:
                out.add(p)
    return out


def _enforce_assignment_invariant(manifest: StudyManifest) -> None:
    """Manifest contract: assigned ∩ unassigned == ∅.  Remove any leaked
    entries and record a warning if the invariant had been violated."""
    claimed = assigned_paths(manifest)
    leaked = [u for u in manifest.unassigned_files if u in claimed]
    if leaked:
        manifest.unassigned_files = [u for u in manifest.unassigned_files if u not in claimed]
        manifest.warnings.append(CampaignWarning(
            "assignment_invariant_repaired",
            f"{len(leaked)} file(s) appeared in both assigned and unassigned lists "
            f"and were removed from unassigned: {[Path(u).name for u in leaked]}",
            Severity.INFO,
        ))


# ═══════════════════════════════════════════════════════════════════════════════
# Serialisation
# ═══════════════════════════════════════════════════════════════════════════════

def write_manifest(manifest: StudyManifest, out_dir: str | Path) -> dict[str, str]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    d = manifest.to_dict()
    yaml_path = out_dir / "study_manifest.yaml"
    json_path = out_dir / "study_manifest.json"
    yaml_path.write_text(yaml.safe_dump(d, sort_keys=False, width=100))
    json_path.write_text(json.dumps(d, indent=2) + "\n")
    return {"yaml": str(yaml_path), "json": str(json_path)}


def load_manifest(path: str | Path) -> StudyManifest:
    p = Path(path)
    text = p.read_text()
    if p.suffix.lower() in (".yaml", ".yml"):
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    manifest = StudyManifest.from_dict(data)
    # honour manually-resolved ambiguities
    _apply_resolutions(manifest)
    return manifest


def _apply_resolutions(manifest: StudyManifest) -> None:
    for amb in manifest.ambiguities:
        if not amb.resolution:
            continue
        rec = manifest.system(amb.system_id)
        if not rec:
            continue
        if amb.kind == "component_classification":
            # resolution format: "receptor=A;peptide=B" (chain assignments)
            assigns = dict(
                part.split("=", 1) for part in amb.resolution.split(";") if "=" in part
            )
            for comp in rec.components:
                key = comp.component_type
                if key in assigns:
                    comp.chain_ids = [assigns[key].strip()]
                    comp.classification_state = ClassificationState.RESOLVED
                    comp.warnings.append("chain assignment set manually via manifest")
            rec.classification_state = ClassificationState.RESOLVED
            rec.user_overridden = True
