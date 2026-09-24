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
    study_root: Optional[Path] = None,
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

    # persistent scientific annotations (never components)
    ambiguities.extend(_resolve_annotations(rec, cand.sim_dir, study_root))

    # partner_type from resolved components
    peptide = det.by_type(ComponentType.PEPTIDE)
    ligand = det.by_type(ComponentType.LIGAND)
    if peptide and peptide.classification_state == ClassificationState.RESOLVED:
        rec.partner_type = "peptide"
    elif ligand and ligand.classification_state == ClassificationState.RESOLVED:
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
                options=_ambiguity_options(comp),
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

    from analysis.campaign.legacy import annotate_legacy
    rec.legacy_study = annotate_legacy(rec, cand.sim_dir)
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
            inspect_trajectories=inspect_trajectories, gmx=gmx, study_root=root,
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


def _ambiguity_options(comp) -> list[str]:
    if comp.component_type == ComponentType.LIGAND:
        for e in comp.evidence:
            if e.kind == "multiple_plausible_candidates" and isinstance(e.value, list):
                return [f"ligand=chain:{c['chain']},resname:{c['resname']},resid:{c['resid']}"
                        if c.get("chain") not in ("_", " ", None) else
                        f"ligand=resname:{c['resname']},resid:{c['resid']}" for c in e.value]
    return [f"chain {c}" for c in comp.chain_ids]


#: Non-polymer / nucleic components declarable explicitly in a manifest
#: resolution, e.g. "ligand=chain:A,resname:LIG,resid:301", "cofactor=resname:HEM",
#: "ligand=group:MyLig" (user index), "nucleic_acid=chain:N".
_EXPLICIT_TYPES = {
    "ligand": ComponentType.LIGAND,
    "cofactor": ComponentType.COFACTOR,
    "nucleic_acid": ComponentType.NUCLEIC_ACID,
}


def _resolve_annotations(rec: SystemRecord, sim_dir, study_root) -> list[Ambiguity]:
    """Resolve annotation declarations; conflicts / proposals become ambiguities."""
    from analysis.campaign.annotations import resolve_system_annotations
    from analysis.campaign.models import AnnotationState
    try:
        records, notes = resolve_system_annotations(rec, sim_dir=sim_dir, study_root=study_root)
    except Exception as exc:  # noqa: BLE001 — a broken annotation file must be visible
        rec.warnings.append(CampaignWarning(
            "annotations_invalid", f"annotations could not be loaded: {exc}", Severity.REVIEW,
            scope=str(sim_dir)))
        return []
    rec.annotations = records
    for n in notes:
        rec.warnings.append(CampaignWarning("annotations", n, Severity.INFO, scope=str(sim_dir)))
    ambs: list[Ambiguity] = []
    for r in records:
        if r.state == AnnotationState.REVIEW_REQUIRED and r.alternatives:
            ambs.append(Ambiguity(
                system_id="(pending)", kind="annotation_conflict",
                message=f"annotation {r.annotation_id!r}: " + "; ".join(r.reasons),
                options=[f"use:{r.annotation_id}:{i}" for i in range(len(r.alternatives))]))
        elif r.state == AnnotationState.PROPOSED:
            ambs.append(Ambiguity(
                system_id="(pending)", kind="annotation_proposal",
                message=f"proposed annotation {r.annotation_id!r} ({r.kind}, origin "
                        f"{r.origin}) is inactive until accepted",
                options=[f"accept:{r.annotation_id}", f"reject:{r.annotation_id}"]))
        if r.state in (AnnotationState.UNRESOLVED, AnnotationState.UNSUPPORTED) or (
                r.state == AnnotationState.REVIEW_REQUIRED and not r.alternatives):
            rec.warnings.append(CampaignWarning(
                "annotation_unresolved",
                f"annotation {r.annotation_id!r} is {r.state}: " + "; ".join(r.reasons),
                Severity.REVIEW, scope=str(sim_dir)))
    return ambs


def _apply_annotation_resolution(rec: SystemRecord, amb: Ambiguity) -> None:
    from analysis.campaign.annotations import GROUP_PREFIX
    from analysis.campaign.models import AnnotationCategory, AnnotationRecord, AnnotationState
    verb, _, rest = amb.resolution.strip().partition(":")
    ann_id, _, idx = rest.partition(":")
    for i, r in enumerate(rec.annotations):
        if r.annotation_id != ann_id:
            continue
        if amb.kind == "annotation_conflict" and verb == "use" and idx.isdigit() \
                and int(idx) < len(r.alternatives):
            chosen = AnnotationRecord.from_dict(r.alternatives[int(idx)])
            others = [a for j, a in enumerate(r.alternatives) if j != int(idx)]
            chosen.alternatives = others
            chosen.conflicts.append({"resolution": f"chosen via manifest ({amb.resolution})"})
            if chosen.state == AnnotationState.SUPERSEDED:
                chosen.state = AnnotationState.ACTIVE
            r = chosen
        elif amb.kind == "annotation_proposal" and verb == "accept" \
                and r.state == AnnotationState.PROPOSED:
            r.state = AnnotationState.ACTIVE
        elif amb.kind == "annotation_proposal" and verb == "reject":
            r.state = AnnotationState.REJECTED
        else:
            rec.warnings.append(CampaignWarning(
                "annotation_resolution_invalid",
                f"cannot apply resolution {amb.resolution!r} to annotation {ann_id!r}",
                Severity.REVIEW, scope=rec.system_id))
            return
        r.provenance["resolved_via"] = f"study_manifest: {amb.resolution}"
        r.group_name = (GROUP_PREFIX + r.annotation_id
                        if r.usable and r.category == AnnotationCategory.RESIDUE_SET else None)
        rec.annotations[i] = r
        rec.user_overridden = True
        return


def _apply_resolutions(manifest: StudyManifest) -> None:
    for amb in manifest.ambiguities:
        if not amb.resolution:
            continue
        rec = manifest.system(amb.system_id)
        if not rec:
            continue
        if amb.kind in ("annotation_conflict", "annotation_proposal"):
            _apply_annotation_resolution(rec, amb)
            continue
        if amb.kind == "component_classification":
            # resolution format: "receptor=A;peptide=B" (chain assignments), plus
            # explicit non-polymer declarations for the _EXPLICIT_TYPES keys
            assigns = dict(
                part.split("=", 1) for part in amb.resolution.split(";") if "=" in part
            )
            assigns = {k.strip(): v for k, v in assigns.items()}
            explicit = {k: v for k, v in assigns.items() if k in _EXPLICIT_TYPES}
            polymer = {k: v for k, v in assigns.items() if k not in _EXPLICIT_TYPES}
            for comp in rec.components:
                key = comp.component_type
                if key in polymer:
                    comp.chain_ids = [polymer[key].strip()]
                    comp.classification_state = ClassificationState.RESOLVED
                    comp.warnings.append("chain assignment set manually via manifest")
            for key, text in explicit.items():
                _apply_explicit_component(rec, _EXPLICIT_TYPES[key], key, text.strip())
            if polymer:
                rec.classification_state = ClassificationState.RESOLVED
            rec.user_overridden = True


def _apply_explicit_component(rec, ctype: str, key: str, text: str) -> None:
    """Explicit manifest intent → a RESOLVED component with an exact atom set.

    Automatic components of the same type are superseded (recorded, removed);
    automatic claims that disagree — a different resolved atom set, or another
    residue-name component covering the selected residues — are recorded as
    conflicts and the explicit declaration wins.
    """
    from analysis.campaign.models import ComponentEvidence, MolecularComponent
    from analysis.campaign.results import atom_set_hash
    from analysis.campaign.structure.entities import parse_selection, select_atoms
    from analysis.campaign.structure.spatial import _iter_atoms

    def warn(msg):
        rec.warnings.append(CampaignWarning("explicit_resolution_failed", msg, Severity.REVIEW,
                                            scope=rec.system_id))
    try:
        sel = parse_selection(text)
    except ValueError as exc:
        return warn(f"{key}={text}: {exc}")
    if not rec.structure_path or not Path(rec.structure_path).is_file():
        return warn(f"{key}={text}: no reference structure to resolve the selection against")
    try:
        ids = select_atoms(rec.structure_path, sel, rec.index_path)
    except ValueError as exc:
        return warn(f"{key}={text}: {exc}")
    if not ids:
        return warn(f"{key}={text}: selection matched no atoms in {Path(rec.structure_path).name}")

    idset = set(ids)
    atoms = [a for a in _iter_atoms(Path(rec.structure_path), first_model_only=True)
             if a.index in idset]
    resnames = sorted({a.resname for a in atoms})
    chains = sorted({a.chain for a in atoms if a.chain not in ("_", " ")})

    superseded, conflicts, kept = [], [], []
    for c in rec.components:
        if c.component_type == ctype:
            superseded.append({"component_type": c.component_type, "label": c.label,
                               "classification_state": c.classification_state,
                               "resnames": c.resnames, "chain_ids": c.chain_ids,
                               "n_atoms": len(c.atom_ids) if c.atom_ids else None})
            same = (set(c.atom_ids) == idset if c.atom_ids
                    else (set(c.chain_ids) == set(chains) if ctype == ComponentType.NUCLEIC_ACID
                          else set(c.resnames) == set(resnames)))
            if c.classification_state == ClassificationState.RESOLVED and not same:
                conflicts.append({"component_type": c.component_type, "label": c.label,
                                  "action": "superseded by explicit resolution"})
            continue
        overlap = sorted(set(c.resnames) & set(resnames)) if not c.chain_ids else []
        if overlap and c.component_type in (ComponentType.COFACTOR, ComponentType.WATER,
                                            ComponentType.IONS, ComponentType.MEMBRANE,
                                            ComponentType.LIGAND):
            c.resnames = [r for r in c.resnames if r not in overlap]
            conflicts.append({"component_type": c.component_type, "label": c.label,
                              "overlap_resnames": overlap,
                              "action": "residue names reassigned to the explicit component"
                                        + ("; component removed" if not c.resnames else "")})
            c.warnings.append(f"residue name(s) {overlap} reassigned by an explicit manifest "
                              f"resolution ({key})")
            if not c.resnames:
                continue
        kept.append(c)

    ev = [ComponentEvidence(
        "explicit_resolution", "declared in the study manifest", 1.0,
        value={"intent_source": "manifest", "selection": sel, "n_atoms": len(ids),
               "atoms_sha256": atom_set_hash(ids)})]
    if superseded:
        ev.append(ComponentEvidence("superseded_automatic_interpretation",
                                    "automatic interpretation(s) replaced by explicit intent",
                                    0.0, value=superseded))
    comp = MolecularComponent(
        component_type=ctype, label=f"{key.replace('_', ' ')} (explicit)",
        selection={"ligand": "Ligand", "cofactor": "Cofactor",
                   "nucleic_acid": "NucleicAcid"}[key],
        chain_ids=chains, resnames=resnames,
        residue_count=len({(a.chain, a.resnum, a.resname) for a in atoms}),
        atom_count=len(ids), atom_ids=sorted(ids), confidence=1.0,
        classification_state=ClassificationState.RESOLVED, evidence=ev,
    )
    if conflicts:
        comp.evidence.append(ComponentEvidence(
            "conflict_with_automatic",
            "automatic interpretation disagreed; explicit intent preserved", 0.0,
            value=conflicts))
        comp.warnings.append("explicit manifest resolution overrides automatic interpretation: "
                             + "; ".join(f"{c['component_type']} ({c['action']})"
                                         for c in conflicts))
    rec.components = kept + [comp]
