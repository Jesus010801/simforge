"""Phase 6 — persistent scientific annotations: model, resolution, precedence,
activation, identity, semantic-index groups, observable evidence and policy."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from analysis.campaign.annotations import (
    annotation_evidence, build_spec_declarations, find_annotation_files,
    load_annotation_file, resolve_system_annotations, write_annotation_file,
)
from analysis.campaign.models import (
    AnnotationState as S, ClassificationState as CS, ComponentType as CT, MolecularComponent,
    SystemRecord,
)
from core.structural_annotation import (
    ANNOTATION_FILENAME, AnnotationProvenance, AnnotationSet, AxisAnnotation, AxisDefinition,
    IndexGroupRef, ReferenceIntent, ResidueSelection, ResidueSetAnnotation, StructuralAnnotation,
)
from tests.analysis.campaign.conftest import requires_gmx

AA = ["ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE"]
BACKBONE = ("N", "CA", "C")


def atoms_of_dimer():
    """Chains A and B, residues 1-30 each (N, CA, C), then a 10-atom LIG."""
    out = []
    for ci, chain in enumerate("AB"):
        for r in range(1, 31):
            for an in BACKBONE:
                out.append((chain, r, AA[(r - 1) % 10], an, 2.0 + 3.0 * ci, 0.1 * r, 1.0))
    out += [("L", 100, "LIG", f"C{k + 1}", 8.0, 1.0 + 0.1 * k, 1.0) for k in range(10)]
    return out


def write_gro(path: Path) -> Path:
    atoms = atoms_of_dimer()
    lines = ["dimer", str(len(atoms))]
    lines += [f"{rs:>5}{rn:<5}{an:>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}"
              for i, (_c, rs, rn, an, x, y, z) in enumerate(atoms, 1)]
    lines.append("  10.00000  10.00000  10.00000")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_pdb(path: Path) -> Path:
    lines = []
    for i, (c, rs, rn, an, x, y, z) in enumerate(atoms_of_dimer(), 1):
        rec = "HETATM" if rn == "LIG" else "ATOM"
        lines.append(f"{rec:<6}{i:>5} {an:<4} {rn:<3} {c:1}{rs:>4}    "
                     f"{10 * x:8.3f}{10 * y:8.3f}{10 * z:8.3f}  1.00  0.00           C")
    path.write_text("\n".join(lines + ["END"]) + "\n")
    return path


def ids(pred) -> list[int]:
    return [i for i, a in enumerate(atoms_of_dimer(), 1) if pred(a)]


def sel(numbering="topology", **kw) -> ResidueSelection:
    return ResidueSelection(numbering=numbering, **kw)


def rs(ann_id, selection=None, kind="binding_site", origin="explicit_api", **kw):
    return ResidueSetAnnotation(id=ann_id, kind=kind, selection=selection,
                                provenance=AnnotationProvenance(origin=origin), **kw)


def aset(*decls, **kw) -> AnnotationSet:
    res = [d for d in decls if isinstance(d, ResidueSetAnnotation)]
    axes = [d for d in decls if isinstance(d, AxisAnnotation)]
    refs = [d for d in decls if isinstance(d, ReferenceIntent)]
    return AnnotationSet(structural_annotation=StructuralAnnotation(residue_sets=res, axes=axes),
                         reference_intents=refs, **kw)


def record(tmp_path, structure="gro", components=None, topology=None):
    s = write_gro(tmp_path / "ref.gro") if structure == "gro" else write_pdb(tmp_path / "ref.pdb")
    return SystemRecord(system_id="s", condition_id="c", replicate_id="r", structure_path=str(s),
                        topology_path=topology, components=components or [])


def resolve(tmp_path, *decls, structure="gro", **kw):
    rec = record(tmp_path, structure, kw.pop("components", None), kw.pop("topology", None))
    records, _ = resolve_system_annotations(rec, extra=[aset(*decls, **kw)])
    return {r.annotation_id: r for r in records}


# ═══════════════════════════════════════════════════════════════════════════════
# Residue selections
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("selection,expected", [
    (dict(chain="A", residues="10"), lambda a: a[0] == "A" and a[1] == 10),
    (dict(chain="A", residues="10-12"), lambda a: a[0] == "A" and 10 <= a[1] <= 12),
    (dict(chain="B", residues="1-2,29-30"), lambda a: a[0] == "B" and a[1] in (1, 2, 29, 30)),
    (dict(chain="A", residues="1-30", resnames=["LYS"]), lambda a: a[0] == "A" and a[2] == "LYS"),
    (dict(chain="A", residues="5-6", atom_names=["CA"]),
     lambda a: a[0] == "A" and a[1] in (5, 6) and a[3] == "CA"),
    (dict(resnames=["LIG"]), lambda a: a[2] == "LIG"),
])
def test_selection_resolution(tmp_path, selection, expected):
    r = resolve(tmp_path, rs("site", sel(**selection)))["site"]
    assert r.state == S.ACTIVE and r.atom_ids == ids(expected)
    assert r.n_atoms == len(r.atom_ids) and r.atoms_sha256 and r.group_name == "Ann_site"


def test_chain_a_residue_differs_from_chain_b(tmp_path):
    out = resolve(tmp_path, rs("a100", sel(chain="A", residues="10")),
                  rs("b100", sel(chain="B", residues="10")))
    assert out["a100"].atoms_sha256 != out["b100"].atoms_sha256
    assert out["a100"].residues == ["A:ILE:10"] and out["b100"].residues == ["B:ILE:10"]


def test_unqualified_residues_across_chains_need_review(tmp_path):
    r = resolve(tmp_path, rs("site", sel(residues="10")))["site"]
    assert r.state == S.REVIEW_REQUIRED and "qualify" in r.reasons[0]
    assert r.group_name is None


def test_numbering_schemes_are_never_interchanged(tmp_path):
    # .gro reference is topology-numbered: an author_pdb selection is not re-interpreted
    r = resolve(tmp_path, rs("site", sel("author_pdb", chain="A", residues="10")))["site"]
    assert r.state == S.UNSUPPORTED and "numbering" in r.reasons[0]
    # a .pdb reference has unknown numbering unless declared
    r = resolve(tmp_path, rs("site", sel("author_pdb", chain="A", residues="10")),
                structure="pdb")["site"]
    assert r.state == S.UNSUPPORTED and "unknown" in r.reasons[0]
    ok = resolve(tmp_path, rs("site", sel("author_pdb", chain="A", residues="10")),
                 structure="pdb", reference_numbering="author_pdb")["site"]
    assert ok.state == S.ACTIVE and ok.numbering == "author_pdb"
    topo_on_author = resolve(tmp_path, rs("site", sel("topology", chain="A", residues="10")),
                             structure="pdb", reference_numbering="author_pdb")["site"]
    assert topo_on_author.state == S.UNSUPPORTED


def test_numbering_is_required():
    with pytest.raises(ValueError):
        ResidueSelection(chain="A", residues="10")


@pytest.mark.parametrize("selection,needle", [
    (dict(resnames=["ZZZ"]), "no atoms"),
    (dict(chain="A", residues="10-40"), "not found"),
    (dict(chain="Q", residues="1"), "not present"),
    (dict(chain="A", residues="10", atom_names=["XX"]), "no atoms"),
])
def test_bad_selections_are_unresolved_not_empty(tmp_path, selection, needle):
    r = resolve(tmp_path, rs("site", sel(**selection)))["site"]
    assert r.state == S.UNRESOLVED and r.atom_ids == [] and r.group_name is None
    assert any(needle in x for x in r.reasons)


# ═══════════════════════════════════════════════════════════════════════════════
# Index-group annotations
# ═══════════════════════════════════════════════════════════════════════════════

def test_index_group_annotation(tmp_path):
    ndx = tmp_path / "user.ndx"
    ndx.write_text("[ MySite ]\n4 5 6\n")
    decl = ResidueSetAnnotation(id="site", kind="binding_site",
                                index_group=IndexGroupRef(path=str(ndx), group="MySite"),
                                provenance=AnnotationProvenance(origin="imported_ndx"))
    r = resolve(tmp_path, decl)["site"]
    assert r.state == S.ACTIVE and r.atom_ids == [4, 5, 6] and r.residues == ["A:GLY:2"]
    first = r.identity
    ndx.write_text("[ MySite ]\n4 5 6 7\n")                      # same name, new atoms
    assert resolve(tmp_path, decl)["site"].identity != first
    moved = tmp_path / "sub"
    moved.mkdir()
    (moved / "user.ndx").write_text("[ MySite ]\n4 5 6\n")       # same atoms, other path
    decl2 = decl.model_copy(update={"index_group": IndexGroupRef(path=str(moved / "user.ndx"),
                                                                 group="MySite")})
    assert resolve(tmp_path, decl2)["site"].identity == first
    missing = decl.model_copy(update={"index_group": IndexGroupRef(path=str(ndx), group="Nope")})
    assert resolve(tmp_path, missing)["site"].state == S.UNRESOLVED
    gone = decl.model_copy(update={"index_group": IndexGroupRef(path="nowhere.ndx", group="X")})
    assert resolve(tmp_path, gone)["site"].state == S.UNRESOLVED


# ═══════════════════════════════════════════════════════════════════════════════
# Kinds, multiplicity, subunits, round trip
# ═══════════════════════════════════════════════════════════════════════════════

def test_multiple_annotations_per_kind_and_subunits(tmp_path):
    out = resolve(tmp_path,
                  rs("binding_site_A", sel(chain="A", residues="10-12")),
                  rs("binding_site_B", sel(chain="B", residues="10-12")),
                  rs("subunit_A", sel(chain="A"), kind="selected_subunit"),
                  rs("subunit_B", sel(chain="B"), kind="selected_subunit"))
    assert all(r.state == S.ACTIVE for r in out.values())
    assert out["subunit_A"].atoms_sha256 != out["subunit_B"].atoms_sha256
    assert out["subunit_A"].atom_ids == ids(lambda a: a[0] == "A")
    assert out["binding_site_A"].atom_ids != out["binding_site_B"].atom_ids


def test_nothing_is_inferred_without_declarations(tmp_path):
    rec = record(tmp_path, components=[MolecularComponent(CT.LIGAND, "ligand",
                                                          classification_state=CS.RESOLVED)])
    records, _ = resolve_system_annotations(rec, sim_dir=tmp_path, study_root=tmp_path)
    assert records == []                                           # no binding site invented


def test_annotation_file_round_trip(tmp_path):
    original = aset(
        rs("site", sel(chain="A", residues="10-12"), origin="user_yaml"),
        rs("cat", sel(chain="A", residues="11"), kind="catalytic_residues", origin="user_yaml"),
        rs("gate_upper", sel(chain="A", residues="1-2"), kind="gate", origin="user_yaml"),
        rs("lining", sel(chain="B", residues="5-8"), kind="pore_lining", origin="user_yaml"),
        rs("mine", sel(resnames=["LIG"]), kind="custom", origin="user_yaml",
           description="free text"),
        AxisAnnotation(id="pore_axis", kind="pore_axis", definition=AxisDefinition(
            type="com_to_com", from_annotation="gate_upper", to_annotation="lining")),
        ReferenceIntent(id="fit_ref", purpose="fit", source="structure"),
        accept=["later"])
    path = write_annotation_file(original, tmp_path / ANNOTATION_FILENAME)
    raw = yaml.safe_load(path.read_text())
    assert raw["schema_version"] == "simforge/annotations/v1"
    loaded = load_annotation_file(path)
    assert loaded == original
    assert [r.kind for r in loaded.structural_annotation.residue_sets] == [
        "binding_site", "catalytic_residues", "gate", "pore_lining", "custom"]
    assert loaded.structural_annotation.axes[0].definition.from_annotation == "gate_upper"


def test_bad_schema_rejected():
    with pytest.raises(ValueError):
        AnnotationSet(schema_version="simforge/annotations/v0")


# ═══════════════════════════════════════════════════════════════════════════════
# Axes
# ═══════════════════════════════════════════════════════════════════════════════

def test_axis_definitions(tmp_path):
    memb = [MolecularComponent(CT.MEMBRANE, "membrane", classification_state=CS.RESOLVED)]
    out = resolve(
        tmp_path,
        rs("gate_upper", sel(chain="A", residues="1-2"), kind="gate"),
        rs("gate_lower", sel(chain="A", residues="29-30"), kind="gate"),
        AxisAnnotation(id="z_axis", definition=AxisDefinition(type="explicit_vector",
                                                              vector=(0, 0, 1), frame="box")),
        AxisAnnotation(id="pore_axis", kind="pore_axis", definition=AxisDefinition(
            type="com_to_com", from_annotation="gate_upper", to_annotation="gate_lower")),
        AxisAnnotation(id="normal", definition=AxisDefinition(type="membrane_normal")),
        AxisAnnotation(id="broken", definition=AxisDefinition(
            type="com_to_com", from_annotation="gate_upper", to_annotation="missing")),
        components=memb)
    assert out["z_axis"].state == S.ACTIVE and out["z_axis"].group_name is None
    pa = out["pore_axis"]
    assert pa.state == S.ACTIVE and pa.depends_on == ["gate_upper", "gate_lower"]
    assert out["normal"].state == S.ACTIVE
    assert out["broken"].state == S.UNRESOLVED
    # the axis identity follows its endpoints' atoms
    moved = resolve(tmp_path,
                    rs("gate_upper", sel(chain="A", residues="1-3"), kind="gate"),
                    rs("gate_lower", sel(chain="A", residues="29-30"), kind="gate"),
                    AxisAnnotation(id="pore_axis", kind="pore_axis", definition=AxisDefinition(
                        type="com_to_com", from_annotation="gate_upper", to_annotation="gate_lower")))
    assert moved["pore_axis"].identity != pa.identity
    no_membrane = resolve(tmp_path, AxisAnnotation(
        id="normal", definition=AxisDefinition(type="membrane_normal")), components=[])
    assert no_membrane["normal"].state == S.UNRESOLVED
    with pytest.raises(ValueError):
        AxisDefinition(type="explicit_vector", vector=(0, 0, 1))          # frame required


# ═══════════════════════════════════════════════════════════════════════════════
# Precedence, conflicts, activation
# ═══════════════════════════════════════════════════════════════════════════════

def test_explicit_beats_derived(tmp_path):
    out = resolve(tmp_path,
                  rs("site", sel(chain="A", residues="20-22"), origin="derived"),
                  rs("site", sel(chain="A", residues="10-12"), origin="user_yaml"))
    r = out["site"]
    assert r.origin == "user_yaml" and r.state == S.ACTIVE
    assert r.atom_ids == ids(lambda a: a[0] == "A" and 10 <= a[1] <= 12)
    assert r.alternatives[0]["origin"] == "derived" and r.alternatives[0]["state"] == S.SUPERSEDED
    assert r.conflicts


def test_conflicting_explicit_annotations_need_review(tmp_path):
    r = resolve(tmp_path, rs("site", sel(chain="A", residues="10"), origin="user_yaml"),
                rs("site", sel(chain="B", residues="10"), origin="user_yaml"))["site"]
    assert r.state == S.REVIEW_REQUIRED and len(r.alternatives) == 2 and r.group_name is None


def test_derived_is_proposed_until_accepted(tmp_path):
    derived = rs("site", sel(chain="A", residues="10-12"), origin="derived")
    r = resolve(tmp_path, derived)["site"]
    assert r.state == S.PROPOSED and r.group_name is None and r.atom_ids   # resolved, inert
    assert resolve(tmp_path, derived, accept=["site"])["site"].state == S.ACTIVE


def _study_with_annotations(tmp_path, files: dict[str, AnnotationSet]):
    root = tmp_path / "study"
    d = root / "sys" / "rep1"
    d.mkdir(parents=True)
    write_gro(d / "reference.gro")
    (d / "production_run.xtc").write_bytes(b"x")
    (d / "topol.top").write_text("[ molecules ]\nProtein 2\n")
    for rel, a in files.items():
        write_annotation_file(a, root / rel / ANNOTATION_FILENAME)
    return root


def test_manifest_proposal_acceptance_and_conflict_resolution(tmp_path):
    from analysis.campaign.manifest import build_manifest, load_manifest, write_manifest
    from analysis.campaign.models import AnnotationRecord
    root = _study_with_annotations(tmp_path, {
        ".": aset(rs("proposed_site", sel(chain="A", residues="10"), origin="derived"),
                  rs("site", sel(chain="A", residues="10"), origin="user_yaml")),
        "sys/rep1": aset(rs("site", sel(chain="B", residues="10"), origin="user_yaml")),
    })
    m = build_manifest(root, inspect_trajectories=False)
    rec = m.systems[0]
    by = {a.annotation_id: a for a in rec.annotations}
    assert by["proposed_site"].state == S.PROPOSED
    assert by["site"].state == S.REVIEW_REQUIRED
    kinds = {a.kind: a for a in m.ambiguities if a.system_id == rec.system_id}
    assert kinds["annotation_proposal"].options == ["accept:proposed_site", "reject:proposed_site"]
    conflict = kinds["annotation_conflict"]
    kinds["annotation_proposal"].resolution = "accept:proposed_site"
    chosen = next(i for i, alt in enumerate(by["site"].alternatives)
                  if alt["residues"] == ["B:ILE:10"])
    conflict.resolution = f"use:site:{chosen}"
    loaded = load_manifest(write_manifest(m, tmp_path / "m")["yaml"]).system(rec.system_id)
    by2 = {a.annotation_id: a for a in loaded.annotations}
    assert by2["proposed_site"].state == S.ACTIVE and by2["proposed_site"].group_name == "Ann_proposed_site"
    assert by2["site"].state == S.ACTIVE and by2["site"].residues == ["B:ILE:10"]
    assert "study_manifest" in by2["site"].provenance["resolved_via"]
    assert loaded.user_overridden


def test_annotation_files_found_from_system_to_root(tmp_path):
    root = _study_with_annotations(tmp_path, {".": aset(), "sys": aset()})
    found = find_annotation_files(root / "sys" / "rep1", root)
    assert [f.parent.name for f in found] == ["sys", "study"]


# ═══════════════════════════════════════════════════════════════════════════════
# Identity
# ═══════════════════════════════════════════════════════════════════════════════

def test_identity_follows_scientific_definition_only(tmp_path):
    base = resolve(tmp_path, rs("site", sel(chain="A", residues="10-12")))["site"]
    def ident(**kw):
        d = tmp_path / str(len(list(tmp_path.iterdir())))
        d.mkdir()
        return resolve(d, rs("site", **kw))["site"]
    assert ident(selection=sel(chain="A", residues="10-13")).identity != base.identity
    assert ident(selection=sel(chain="B", residues="10-12")).identity != base.identity
    assert ident(selection=sel(chain="A", residues="10-12", atom_names=["CA"])).identity != base.identity
    same = ident(selection=sel(chain="A", residues="10-12"), description="new text")
    assert same.identity == base.identity and same.atoms_sha256 == base.atoms_sha256
    from analysis.campaign.results import atom_set_hash
    assert base.atoms_sha256 == atom_set_hash(base.atom_ids)          # Phase 2 hashing


# ═══════════════════════════════════════════════════════════════════════════════
# Downstream: observables and policy
# ═══════════════════════════════════════════════════════════════════════════════

def test_observable_definition_isolation(tmp_path):
    from analysis.campaign.observables.base import AnalysisContext, ObservableSpec
    from analysis.campaign.models import TrajectoryRequirements, TrajectoryView

    class UsesSite(ObservableSpec):
        id = "uses-site"
        required_annotations = ("binding_site",)
        def definition_evidence(self, ctx):
            ann = self.annotation_evidence(ctx)
            return None if ann is None else {"observable": self.id, "annotations": ann}

    class NoSite(ObservableSpec):
        id = "no-site"
        def definition_evidence(self, ctx):
            return {"observable": self.id}

    def ctx_for(residues):
        rec = record(tmp_path)
        rec.annotations, _ = resolve_system_annotations(
            rec, extra=[aset(rs("binding_site", sel(chain="A", residues=residues)))])
        return AnalysisContext(system=rec, semantic_index=None,
                               trajectory_view=TrajectoryView(kind="raw",
                                                              requirements=TrajectoryRequirements()),
                               topology_path="", structure_path=None, output_dir=tmp_path,
                               parameters={})
    a1, a2 = UsesSite().definition_token(ctx_for("10-12")), UsesSite().definition_token(ctx_for("10-13"))
    b1, b2 = NoSite().definition_token(ctx_for("10-12")), NoSite().definition_token(ctx_for("10-13"))
    assert a1 != a2 and b1 == b2
    rec = record(tmp_path)
    rec.annotations, _ = resolve_system_annotations(
        rec, extra=[aset(rs("binding_site", sel(chain="A", residues="10"), origin="derived"))])
    assert annotation_evidence(rec.annotations, "binding_site") is None     # proposed: unusable


def _policy_ctx(tmp_path, *decls, structure=True):
    from analysis.campaign.trajectory.policy import PolicyContext
    rec = record(tmp_path, topology=str(tmp_path / "md.tpr"))
    (tmp_path / "md.tpr").write_bytes(b"tpr")
    records, _ = resolve_system_annotations(rec, extra=[aset(*decls)])
    ndx = tmp_path / "semantic.ndx"
    ndx.write_text("[ System ]\n1\n" + "".join(f"[ {r.group_name} ]\n1\n" for r in records
                                               if r.group_name))
    return PolicyContext(components=[MolecularComponent(CT.RECEPTOR, "r",
                                                        classification_state=CS.RESOLVED)],
                         topology_path=rec.topology_path,
                         structure_path=rec.structure_path if structure else None,
                         index_path=str(ndx), annotations=records)


def test_policy_uses_annotation_state_not_group_existence(tmp_path):
    from analysis.campaign.models import DecisionClass as DC, ObservablePurpose as P
    from analysis.campaign.models import TrajectoryRequirements as TR
    from analysis.campaign.trajectory.policy import plan_preprocessing
    fit = TR(requires_whole_molecules=True, fit_selection="Ann_subunit_B")
    active = _policy_ctx(tmp_path, rs("subunit_B", sel(chain="B"), kind="selected_subunit"),
                         ReferenceIntent(id="fit_ref", purpose="fit", source="structure",
                                         provenance=AnnotationProvenance(origin="explicit_api")))
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=active).decisions[1]
    assert (d.classification, d.rule_id) == (DC.VALID, "fit.frame_alignment")
    assert d.system_context["component"]["annotation"]["id"] == "subunit_B"
    assert d.system_context["reference"]["intent"]["id"] == "fit_ref"

    (tmp_path / "p").mkdir()
    proposed = _policy_ctx(tmp_path / "p", rs("subunit_B", sel(chain="B"), kind="selected_subunit",
                                              origin="derived"),
                           ReferenceIntent(id="fit_ref", purpose="fit", source="structure",
                                           provenance=AnnotationProvenance(origin="explicit_api")))
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=proposed).decisions
    assert [x.classification for x in d if x.operation == "fit"] == [DC.UNSUPPORTED]  # no Ann_ group
    ctx = proposed
    ctx.index_path = active.index_path                                   # group exists, state proposed
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=ctx).decisions[1]
    assert (d.classification, d.rule_id) == (DC.VALID_WITH_INTENT, "fit.group_unresolved")


def test_reference_intent_must_match_what_executes(tmp_path):
    from analysis.campaign.models import DecisionClass as DC, ObservablePurpose as P
    from analysis.campaign.models import TrajectoryRequirements as TR
    from analysis.campaign.trajectory.policy import plan_preprocessing
    fit = TR(requires_whole_molecules=True, fit_selection="Ann_subunit_B")
    ambiguous = _policy_ctx(tmp_path, rs("subunit_B", sel(chain="B"), kind="selected_subunit"))
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=ambiguous).decisions[1]
    assert d.rule_id == "fit.reference_ambiguous"                         # implicit ambiguity
    (tmp_path / "t").mkdir()
    topo = _policy_ctx(tmp_path / "t", rs("subunit_B", sel(chain="B"), kind="selected_subunit"),
                       ReferenceIntent(id="fit_ref", purpose="fit", source="topology",
                                       provenance=AnnotationProvenance(origin="explicit_api")))
    # declared .tpr reference, but build_view would fit against the .gro -> not executable
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=topo).decisions[1]
    assert (d.classification, d.rule_id) == (DC.UNSUPPORTED, "fit.reference_intent_unexecutable")
    (tmp_path / "u").mkdir()
    consistent = _policy_ctx(tmp_path / "u", rs("subunit_B", sel(chain="B"), kind="selected_subunit"),
                             ReferenceIntent(id="fit_ref", purpose="fit", source="topology",
                                             provenance=AnnotationProvenance(origin="explicit_api")),
                             structure=False)
    # no structure file -> build_view fits against the .tpr, which is what was declared
    d = plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=consistent).decisions[1]
    assert (d.classification, d.rule_id) == (DC.VALID, "fit.frame_alignment")


# ═══════════════════════════════════════════════════════════════════════════════
# Build-time compatibility
# ═══════════════════════════════════════════════════════════════════════════════

def test_old_structural_annotation_yaml_still_loads():
    sa = StructuralAnnotation(**{"membrane_topology": {"extracellular_regions": ["1-10"],
                                                       "intracellular_regions": ["40-50"],
                                                       "transmembrane_segments": ["11-39"]}})
    assert sa.residue_sets == [] and sa.axes == [] and sa.is_partial_for_orient()


def test_build_spec_annotations_surface_in_campaign(tmp_path):
    spec = tmp_path / "spec.yaml"
    spec.write_text(yaml.safe_dump({"structural_annotation": {"membrane_topology": {
        "extracellular_regions": ["1-5"], "intracellular_regions": ["26-30"],
        "transmembrane_segments": ["6-25"]}}}))
    run = tmp_path / "run"
    (run / "metadata").mkdir(parents=True)
    (run / "metadata" / "run_info.json").write_text(json.dumps({"yaml_source": str(spec)}))
    step = run / "steps" / "11_production_md"
    step.mkdir(parents=True)
    decls, notes = build_spec_declarations(step, run)
    assert {d.decl.id for d in decls} == {"tm_1", "extracellular_1", "intracellular_1"}
    assert all(d.origin == "build_spec" and d.decl.selection.numbering == "topology" for d in decls)
    rec = record(step)
    records, _ = resolve_system_annotations(rec, sim_dir=step, study_root=run)
    tm = next(r for r in records if r.annotation_id == "tm_1")
    assert tm.kind == "transmembrane_segment" and tm.origin == "build_spec"
    # single residue range but a two-chain structure -> must be qualified, not guessed
    assert tm.state == S.REVIEW_REQUIRED


# ═══════════════════════════════════════════════════════════════════════════════
# Semantic index (gmx)
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
def test_active_annotations_become_ann_groups(tmp_path):
    from analysis.campaign.structure.component_detector import detect_components
    from analysis.campaign.structure.index_builder import build_semantic_index
    from analysis.campaign.structure.index_groups import parse_index_groups
    rec = record(tmp_path)
    rec.components = detect_components(structure_path=rec.structure_path).components
    rec.annotations, _ = resolve_system_annotations(rec, extra=[aset(
        rs("site_A", sel(chain="A", residues="10-12")),
        rs("proposed", sel(chain="B", residues="10-12"), origin="derived"),
        rs("broken", sel(chain="Q", residues="1")))])
    idx = build_semantic_index(structure_path=rec.structure_path, components=rec.components,
                               out_ndx=tmp_path / "semantic.ndx", annotations=rec.annotations)
    groups = {g.name: g.atom_ids for g in parse_index_groups(idx.path)}
    assert list(groups["Ann_site_A"]) == ids(lambda a: a[0] == "A" and 10 <= a[1] <= 12)
    assert "Ann_proposed" not in groups and "Ann_broken" not in groups
    assert list(groups)[-1] == "System"


# ═══════════════════════════════════════════════════════════════════════════════
# polymer_only and the real bundled membrane build spec
# ═══════════════════════════════════════════════════════════════════════════════

def test_polymer_only_excludes_lipid_residue_numbers(tmp_path):
    atoms = [("A", r, "ALA", "CA", 1.0, 0.1 * r, 1.0) for r in range(1, 11)]
    atoms += [("_", r, "POPC", "P", 5.0, 0.1 * r, 1.0) for r in range(11, 16)]
    lines = ["mix", str(len(atoms))]
    lines += [f"{rs:>5}{rn:<5}{an:>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}"
              for i, (_c, rs, rn, an, x, y, z) in enumerate(atoms, 1)]
    lines.append("  10.00000  10.00000  10.00000")
    gro = tmp_path / "mix.gro"
    gro.write_text("\n".join(lines) + "\n")
    rec = SystemRecord("s", "c", "r", structure_path=str(gro))
    def one(**kw):
        records, _ = resolve_system_annotations(rec, extra=[aset(rs("x", sel(**kw)))])
        return records[0]
    assert one(residues="12-14").state == S.ACTIVE                       # plain: lipids match
    r = one(residues="12-14", polymer_only=True)
    assert r.state == S.UNRESOLVED and "not found" in r.reasons[0]       # nothing silently kept
    assert one(residues="5-8", polymer_only=True).n_atoms == 4
    assert one(residues="5-8", polymer_only=True).identity != one(residues="5-8").identity


from tests.analysis.campaign.conftest import REAL_GRO, requires_real_traj  # noqa: E402


@requires_real_traj
def test_real_membrane_build_spec_annotations():
    run = REAL_GRO.parents[2]
    step = REAL_GRO.parent
    rec = SystemRecord("glp1r", "c", "r", structure_path=str(REAL_GRO))
    records, _ = resolve_system_annotations(rec, sim_dir=step, study_root=run)
    by = {r.annotation_id: r for r in records}
    for i in range(1, 6):
        r = by[f"tm_{i}"]
        assert r.state == S.ACTIVE and r.origin == "build_spec"
        assert {x.split(":")[0] for x in r.residues} == {"A"}            # protein only
    # the build spec declares residues the prepared protein does not contain
    assert by["tm_6"].state == S.UNRESOLVED and "[401, 402, 403, 404]" in by["tm_6"].reasons[0]
    assert by["intracellular_1"].state == S.UNRESOLVED
    assert by["extracellular_1"].state == S.ACTIVE and by["extracellular_1"].n_atoms == 1401
