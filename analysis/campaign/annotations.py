"""Persistent scientific annotations → resolved atoms (Phase 6).

Declarations (``core.structural_annotation``) come from:

* ``simforge_annotations.yaml`` files from the system directory up to the
  study root (origin ``user_yaml`` unless declared otherwise);
* the build specification the run was compiled from (``metadata/run_info.json``
  → ``yaml_source`` → ``structural_annotation``; origin ``build_spec``);
* programmatic declarations (``explicit_api``).

Resolution happens once, against the system's reference structure, and yields
:class:`AnnotationRecord` objects with exact atom ids and an atom-set hash
(Phase 2 hashing).  Rules:

* residue numbers are matched only in the numbering scheme they were declared
  in — a mismatch is UNSUPPORTED, never silently re-interpreted;
* a selection that misses atoms, residues or chains is UNRESOLVED, never an
  empty or partial "valid" annotation;
* precedence: explicit (user_cli/explicit_api) > study-level (user_yaml,
  study_manifest, imported_ndx) > build_spec > derived; equal-rank conflicts are
  REVIEW_REQUIRED (recorded, not chosen); derived annotations stay PROPOSED
  until explicitly accepted.

Nothing here infers biology, scans trajectories, or queries databases.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Union

import yaml

from analysis.campaign.models import (
    AnnotationCategory, AnnotationRecord, AnnotationState, ClassificationState, ComponentType,
)
from core.structural_annotation import (
    ANNOTATION_FILENAME, AnnotationProvenance, AnnotationSet, AxisAnnotation, IndexGroupRef,
    ReferenceIntent, ResidueSelection, ResidueSetAnnotation, StructuralAnnotation,
    TransmembraneSegment, residues_in_range,
)

GROUP_PREFIX = "Ann_"

#: higher wins; equal-rank disagreements are conflicts
PRECEDENCE = {
    "user_cli": 5, "explicit_api": 5,
    "user_yaml": 4, "study_manifest": 4, "imported_ndx": 4,
    "build_spec": 3,
    "derived": 1,
}

Declaration = Union[ResidueSetAnnotation, AxisAnnotation, ReferenceIntent]


@dataclass
class _Decl:
    decl: Declaration
    origin: str
    base_dir: Optional[Path]           # for relative paths (ndx, reference files)
    source: Optional[str]              # file it came from


# ═══════════════════════════════════════════════════════════════════════════════
# Loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_annotation_file(path: str | Path) -> AnnotationSet:
    data = yaml.safe_load(Path(path).read_text()) or {}
    return AnnotationSet(**data)


def write_annotation_file(aset: AnnotationSet, path: str | Path) -> Path:
    p = Path(path)
    p.write_text(yaml.safe_dump(aset.model_dump(mode="json", exclude_defaults=False),
                                sort_keys=False, width=100))
    return p


def find_annotation_files(sim_dir: str | Path, study_root: Optional[str | Path]) -> list[Path]:
    """``simforge_annotations.yaml`` from the system directory up to the study root."""
    d = Path(sim_dir).resolve()
    root = Path(study_root).resolve() if study_root else d
    out: list[Path] = []
    while True:
        f = d / ANNOTATION_FILENAME
        if f.is_file():
            out.append(f)
        if d == root or d.parent == d or root not in d.parents:
            break
        d = d.parent
    return out


def _decls_from_set(aset: AnnotationSet, default_origin: str, base: Optional[Path],
                    source: Optional[str]) -> list[_Decl]:
    out: list[_Decl] = []
    sa = aset.structural_annotation
    for d in [*sa.residue_sets, *sa.axes, *aset.reference_intents]:
        # a declared provenance wins; otherwise the channel it arrived through decides
        origin = d.provenance.origin if "provenance" in d.model_fields_set else default_origin
        out.append(_Decl(d, origin, base, source))
    return out


def build_spec_declarations(sim_dir: str | Path, study_root: Optional[str | Path]
                            ) -> tuple[list[_Decl], list[str]]:
    """Explicit build-time annotations of the spec the run was compiled from.

    Build steps (orient_protein) apply these residue ranges to GRO residue
    numbers, so they are topology-numbered.  Nothing is inferred.
    """
    notes: list[str] = []
    d = Path(sim_dir).resolve()
    root = Path(study_root).resolve() if study_root else None
    info = None
    while True:
        cand = d / "metadata" / "run_info.json"
        if cand.is_file():
            info = cand
            break
        if d.parent == d or (root is not None and d == root):
            break
        d = d.parent
    if info is None:
        return [], notes
    try:
        src = json.loads(info.read_text()).get("yaml_source")
    except (OSError, ValueError):
        return [], [f"{info}: unreadable run_info.json"]
    if not src or not Path(src).is_file():
        return [], [f"build specification {src!r} referenced by {info} is not available"]
    raw = yaml.safe_load(Path(src).read_text()) or {}
    if "structural_annotation" not in raw:
        return [], notes
    sa = StructuralAnnotation(**raw["structural_annotation"])
    prov = AnnotationProvenance(origin="build_spec", loaded_from=str(Path(src).resolve()),
                                notes="build-time structural_annotation; ranges applied to GRO "
                                      "residue numbers by the builder (topology numbering)")
    decls: list[_Decl] = []

    def sel(rng: str) -> ResidueSelection:
        # the builder applies these ranges to protein Cα atoms only
        return ResidueSelection(numbering="topology", residues=rng, polymer_only=True)

    mt = sa.membrane_topology
    if mt is not None:
        for i, seg in enumerate(mt.transmembrane_segments, 1):
            rng = seg.residues if isinstance(seg, TransmembraneSegment) else seg
            label = seg.label if isinstance(seg, TransmembraneSegment) and seg.label else f"TM{i}"
            decls.append(_Decl(ResidueSetAnnotation(
                id=f"tm_{i}", kind="transmembrane_segment", selection=sel(rng), label=label,
                provenance=prov), "build_spec", Path(src).parent, src))
        for side, regions in (("extracellular", mt.extracellular_regions),
                              ("intracellular", mt.intracellular_regions)):
            for i, rng in enumerate(regions, 1):
                decls.append(_Decl(ResidueSetAnnotation(
                    id=f"{side}_{i}", kind="domain", selection=sel(rng), label=side,
                    provenance=prov), "build_spec", Path(src).parent, src))
    for i, dom in enumerate(sa.domains, 1):
        decls.append(_Decl(ResidueSetAnnotation(
            id=f"domain_{i}", kind="domain", selection=sel(dom.residues),
            label=dom.label or dom.role, provenance=prov), "build_spec", Path(src).parent, src))
    for d2 in [*sa.residue_sets, *sa.axes]:
        decls.append(_Decl(d2, "build_spec", Path(src).parent, src))
    return decls, notes


# ═══════════════════════════════════════════════════════════════════════════════
# Resolution
# ═══════════════════════════════════════════════════════════════════════════════

def reference_numbering(structure_path: Optional[str], declared: Optional[str]) -> Optional[str]:
    """.gro/.tpr are GROMACS-prepared → topology numbering; a .pdb needs a declaration."""
    if not structure_path:
        return None
    suffix = Path(structure_path).suffix.lower()
    if suffix in (".gro", ".tpr"):
        return "topology"
    return declared


def _atoms(structure_path: str):
    from analysis.campaign.structure.spatial import _iter_atoms
    return list(_iter_atoms(Path(structure_path), first_model_only=True))


def _resolve_selection(sel: ResidueSelection, atoms, ref_numbering: Optional[str]
                       ) -> tuple[str, list[int], list[str], list[str]]:
    """(state, atom_ids, residue identities, reasons)."""
    if ref_numbering is None:
        return (AnnotationState.UNSUPPORTED, [], [],
                ["the reference structure's residue numbering scheme is unknown (a .pdb "
                 "reference needs reference_numbering declared)"])
    if sel.numbering != ref_numbering:
        return (AnnotationState.UNSUPPORTED, [], [],
                [f"selection uses '{sel.numbering}' numbering but the reference structure is "
                 f"'{ref_numbering}'-numbered; no residue mapping is available — declare it in "
                 f"'{ref_numbering}' numbering"])
    chains = sorted({a.chain for a in atoms})
    if sel.chain is not None and sel.chain not in chains:
        return (AnnotationState.UNRESOLVED, [], [],
                [f"chain {sel.chain!r} is not present (available: {chains})"])
    wanted = residues_in_range(sel.residues) if sel.residues else None
    resn = set(sel.resnames)
    names = set(sel.atom_names)
    if sel.polymer_only:
        from analysis.campaign.structure.residue_classes import is_polymer_resname
        atoms = [a for a in atoms if is_polymer_resname(a.resname)]
    in_range = [a for a in atoms
                if (sel.chain is None or a.chain == sel.chain)
                and (wanted is None or a.resnum in wanted)]
    if wanted is not None:        # existence is checked before name filters narrow the set
        missing = sorted(wanted - {a.resnum for a in in_range})
        if missing:
            return (AnnotationState.UNRESOLVED, [], [],
                    [f"residue(s) {missing[:10]}{'…' if len(missing) > 10 else ''} not found"
                     + (f" in chain {sel.chain}" if sel.chain else "")])
    picked = [a for a in in_range if (not resn or a.resname in resn)
              and (not names or a.atomname in names)]
    if not picked:
        return AnnotationState.UNRESOLVED, [], [], ["selection matched no atoms"]
    residues = list(dict.fromkeys(f"{a.chain}:{a.resname}:{a.resnum}" for a in picked))
    reasons: list[str] = []
    state = AnnotationState.ACTIVE
    if sel.residues and sel.chain is None:
        spanned = sorted({a.chain for a in picked})
        if len(spanned) > 1:
            state = AnnotationState.REVIEW_REQUIRED
            reasons.append(f"residue numbers match in several chains {spanned}; qualify the "
                           f"selection with a chain")
    return state, [a.index for a in picked], residues, reasons


def _resolve_index_group(ref: IndexGroupRef, base: Optional[Path], atoms
                         ) -> tuple[str, list[int], list[str], list[str], dict]:
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.structure.index_groups import resolve_index_group
    path = Path(ref.path)
    if not path.is_absolute() and base is not None:
        path = base / path
    if not path.is_file():
        return AnnotationState.UNRESOLVED, [], [], [f"index file {ref.path} not found"], {}
    try:
        grp = resolve_index_group(path, ref.group)
    except ValueError as exc:
        return AnnotationState.UNRESOLVED, [], [], [str(exc)], {}
    info = {"index_file": str(path.resolve()), "index_digest": fingerprint_file(path).digest}
    if grp is None or not grp.atom_ids:
        return (AnnotationState.UNRESOLVED, [], [],
                [f"index group {ref.group!r} missing or empty in {path.name}"], info)
    ids = sorted(set(grp.atom_ids))
    by_index = {a.index: a for a in atoms}
    absent = [i for i in ids if i not in by_index]
    if absent:
        return (AnnotationState.UNRESOLVED, [], [],
                [f"group {ref.group!r} references atoms beyond the reference structure "
                 f"(e.g. {absent[:5]})"], info)
    residues = list(dict.fromkeys(
        f"{by_index[i].chain}:{by_index[i].resname}:{by_index[i].resnum}" for i in ids))
    return AnnotationState.ACTIVE, ids, residues, [], info


def _identity(rec: AnnotationRecord, canonical: dict, depends: list[str]) -> str:
    from analysis.campaign.results import definition_token
    return definition_token({"annotation": rec.annotation_id, "category": rec.category,
                             "kind": rec.kind, "definition": canonical,
                             "atoms_sha256": rec.atoms_sha256, "depends_on": depends})


def _canonical(decl: Declaration) -> dict:
    if isinstance(decl, ResidueSetAnnotation):
        if decl.selection is not None:
            return {"selection": decl.selection.canonical()}
        return {"index_group": decl.index_group.group}
    if isinstance(decl, AxisAnnotation):
        return {"axis": decl.definition.model_dump(mode="json")}
    return {"reference": {"purpose": decl.purpose, "source": decl.source, "file": decl.file}}


def _resolve_one(d: _Decl, *, atoms, ref_numbering, structure_path, rec_paths,
                 resolved: dict, components) -> AnnotationRecord:
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.results import atom_set_hash
    decl = d.decl
    prov = decl.provenance.model_dump(mode="json")
    prov["origin"] = d.origin
    if d.source:
        prov.setdefault("loaded_from", d.source)
    cat = (AnnotationCategory.RESIDUE_SET if isinstance(decl, ResidueSetAnnotation) else
           AnnotationCategory.AXIS if isinstance(decl, AxisAnnotation) else
           AnnotationCategory.REFERENCE)
    kind = decl.kind if cat != AnnotationCategory.REFERENCE else f"reference:{decl.purpose}"
    rec = AnnotationRecord(annotation_id=decl.id, category=cat, kind=kind, origin=d.origin,
                           state=AnnotationState.ACTIVE, definition=decl.model_dump(mode="json"),
                           provenance=prov, description=decl.description)
    rec.definition.pop("description", None)
    rec.definition.pop("provenance", None)
    depends: list[str] = []

    if cat == AnnotationCategory.RESIDUE_SET:
        ref = {"structure": str(Path(structure_path).resolve()) if structure_path else None,
               "structure_digest": fingerprint_file(Path(structure_path)).digest
               if structure_path and Path(structure_path).is_file() else None,
               "structure_numbering": ref_numbering}
        if decl.selection is not None:
            rec.numbering = decl.selection.numbering
            state, ids, residues, reasons = _resolve_selection(decl.selection, atoms, ref_numbering)
        else:
            state, ids, residues, reasons, info = _resolve_index_group(decl.index_group,
                                                                       d.base_dir, atoms)
            rec.numbering = "atom_index"
            ref.update(info)
        rec.reference = ref
        rec.state, rec.atom_ids, rec.residues, rec.reasons = state, ids, residues, reasons
        rec.n_atoms = len(ids)
        rec.atoms_sha256 = atom_set_hash(ids) if ids else None
    elif cat == AnnotationCategory.AXIS:
        dfn = decl.definition
        if dfn.type == "com_to_com":
            for ref_id in (dfn.from_annotation, dfn.to_annotation):
                other = resolved.get(ref_id)
                if other is None or not other.usable or other.category != AnnotationCategory.RESIDUE_SET:
                    rec.state = AnnotationState.UNRESOLVED
                    rec.reasons.append(f"axis endpoint annotation {ref_id!r} is not an active "
                                       f"residue-set annotation")
                else:
                    depends.append(other.identity)
            rec.depends_on = [dfn.from_annotation, dfn.to_annotation]
        elif dfn.type == "membrane_normal":
            mem = [c for c in components or [] if c.component_type == ComponentType.MEMBRANE]
            if components is not None and not (len(mem) == 1 and
                                               mem[0].classification_state == ClassificationState.RESOLVED):
                rec.state = AnnotationState.UNRESOLVED
                rec.reasons.append("membrane_normal requires a single RESOLVED membrane component")
        rec.reasons.append("axis definition only — not evaluated per frame in this phase")
    else:
        target = {"topology": rec_paths.get("topology"), "structure": rec_paths.get("structure"),
                  "file": None}[decl.source]
        if decl.source == "file":
            p = Path(decl.file)
            target = str(p if p.is_absolute() or d.base_dir is None else d.base_dir / p)
        if not target or not Path(target).is_file():
            rec.state = AnnotationState.UNRESOLVED
            rec.reasons.append(f"reference source {decl.source!r} is not available for this system")
        elif decl.source == "topology" and not target.lower().endswith(".tpr"):
            rec.state = AnnotationState.UNSUPPORTED
            rec.reasons.append("topology reference coordinates need a .tpr")
        else:
            rec.reference = {"file": str(Path(target).resolve()),
                             "digest": fingerprint_file(Path(target)).digest,
                             "source": decl.source, "purpose": decl.purpose}

    rec.identity = _identity(rec, _canonical(decl), depends)
    if rec.usable and cat == AnnotationCategory.RESIDUE_SET:
        rec.group_name = GROUP_PREFIX + rec.annotation_id
    return rec


def resolve_declarations(
    decls: Iterable[_Decl], *, structure_path: Optional[str], topology_path: Optional[str] = None,
    reference_numbering_declared: Optional[str] = None, accepted: Iterable[str] = (),
    components=None,
) -> list[AnnotationRecord]:
    """Resolve declarations with precedence, conflicts and acceptance."""
    decls = list(decls)
    accepted = set(accepted)
    atoms = _atoms(structure_path) if structure_path and Path(structure_path).is_file() else []
    ref_num = reference_numbering(structure_path, reference_numbering_declared)
    rec_paths = {"topology": topology_path, "structure": structure_path}
    by_id: dict[str, list[_Decl]] = {}
    for d in decls:
        by_id.setdefault(d.decl.id, []).append(d)

    order = sorted(by_id, key=lambda i: (
        0 if isinstance(by_id[i][0].decl, ResidueSetAnnotation) else
        1 if isinstance(by_id[i][0].decl, AxisAnnotation) else 2, i))
    resolved: dict[str, AnnotationRecord] = {}
    for ann_id in order:
        cands = sorted(by_id[ann_id], key=lambda d: -PRECEDENCE.get(d.origin, 0))
        top = PRECEDENCE.get(cands[0].origin, 0)
        top_c = [d for d in cands if PRECEDENCE.get(d.origin, 0) == top]
        kinds = {type(d.decl) for d in cands}
        distinct: list[_Decl] = []
        for d in top_c:
            if not any(_canonical(d.decl) == _canonical(x.decl) and type(d.decl) is type(x.decl)
                       for x in distinct):
                distinct.append(d)
        kw = dict(atoms=atoms, ref_numbering=ref_num, structure_path=structure_path,
                  rec_paths=rec_paths, resolved=resolved, components=components)
        if len(distinct) > 1 or len(kinds) > 1:
            alts = [_resolve_one(d, **kw) for d in cands]
            rec = AnnotationRecord(annotation_id=ann_id, category=alts[0].category,
                                   kind=alts[0].kind, origin=cands[0].origin,
                                   state=AnnotationState.REVIEW_REQUIRED,
                                   reasons=[f"{len(distinct)} different declarations of "
                                            f"{ann_id!r} at the same precedence — choose one"],
                                   alternatives=[a.to_dict() for a in alts])
            resolved[ann_id] = rec
            continue
        winner = distinct[0]
        rec = _resolve_one(winner, **kw)
        for loser in cands:
            if loser is winner or (_canonical(loser.decl) == _canonical(winner.decl)
                                   and loser.origin == winner.origin):
                continue
            alt = _resolve_one(loser, **kw)
            alt.state = AnnotationState.SUPERSEDED
            rec.alternatives.append(alt.to_dict())
            if _canonical(loser.decl) != _canonical(winner.decl):
                rec.conflicts.append({"origin": loser.origin,
                                      "resolution": f"{winner.origin} declaration takes precedence"})
        declared = winner.decl.state
        if rec.usable and winner.origin == "derived" and ann_id not in accepted:
            rec.state = AnnotationState.PROPOSED
            rec.reasons.append("derived annotation — inactive until explicitly accepted")
        elif rec.usable and declared in ("proposed", "review_required", "rejected") \
                and ann_id not in accepted:
            rec.state = {"proposed": AnnotationState.PROPOSED,
                         "review_required": AnnotationState.REVIEW_REQUIRED,
                         "rejected": AnnotationState.REJECTED}[declared]
        elif rec.usable and ann_id in accepted and winner.origin == "derived":
            rec.provenance["accepted"] = True
        if not rec.usable:
            rec.group_name = None
        resolved[ann_id] = rec
    return [resolved[i] for i in sorted(resolved)]


def resolve_system_annotations(
    rec, *, sim_dir: Optional[str | Path] = None, study_root: Optional[str | Path] = None,
    extra: Iterable[AnnotationSet] = (), extra_origin: str = "explicit_api",
) -> tuple[list[AnnotationRecord], list[str]]:
    """All declaration sources for one system → resolved records (+ notes)."""
    decls: list[_Decl] = []
    notes: list[str] = []
    accepted: set[str] = set()
    ref_decl: Optional[str] = None
    if sim_dir is not None:
        for f in find_annotation_files(sim_dir, study_root):
            aset = load_annotation_file(f)
            decls += _decls_from_set(aset, "user_yaml", f.parent, str(f.resolve()))
            accepted |= set(aset.accept)
            ref_decl = ref_decl or aset.reference_numbering
        bs, n = build_spec_declarations(sim_dir, study_root)
        decls += bs
        notes += n
    for aset in extra:
        decls += _decls_from_set(aset, extra_origin, None, None)
        accepted |= set(aset.accept)
        ref_decl = ref_decl or aset.reference_numbering
    if not decls:
        return [], notes
    records = resolve_declarations(
        decls, structure_path=rec.structure_path, topology_path=rec.topology_path,
        reference_numbering_declared=ref_decl, accepted=accepted, components=rec.components)
    return records, notes


def annotation_evidence(records: Iterable[AnnotationRecord], ann_id: str) -> Optional[dict]:
    """Identity-bearing evidence of one ACTIVE annotation (None otherwise)."""
    for r in records:
        if r.annotation_id == ann_id:
            return r.evidence() if r.usable else None
    return None


def evaluate_com_to_com_axis(axis: AnnotationRecord, com_series: dict, *, eps_nm: float = 1e-6):
    """Evaluate an ACTIVE ``com_to_com`` axis from per-frame centres.

    ``com_series`` maps annotation id -> (n_frames, 3) array (nm).  Returns
    ``(vectors, lengths, unit_vectors)``; frames whose endpoints coincide keep
    NaN unit vectors (undefined direction) — never silently normalised.
    Pure: computing the centres is the caller's (observable's) job.
    """
    import numpy as np
    if not axis.usable or axis.category != AnnotationCategory.AXIS:
        raise ValueError(f"axis {axis.annotation_id!r} is not an ACTIVE axis annotation")
    dfn = axis.definition.get("definition", {})
    if dfn.get("type") != "com_to_com":
        raise ValueError(f"axis {axis.annotation_id!r} is {dfn.get('type')!r}, not com_to_com")
    a, b = dfn["from_annotation"], dfn["to_annotation"]
    missing = [x for x in (a, b) if x not in com_series]
    if missing:
        raise ValueError(f"missing centre series for {missing}")
    pa, pb = np.asarray(com_series[a], float), np.asarray(com_series[b], float)
    if pa.shape != pb.shape:
        raise ValueError("endpoint series have different shapes")
    vec = pb - pa
    length = np.linalg.norm(vec, axis=1)
    unit = np.full_like(vec, np.nan)
    ok = length > eps_nm
    unit[ok] = vec[ok] / length[ok, None]
    return vec, length, unit
