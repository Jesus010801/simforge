"""Phase 13.5 — multi-molecule PBC semantics: topology partition, geometry
requirements, finite-assembly (clustered) views, policy, molecule-level
diagnostics, identity invalidation."""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from analysis.campaign.models import (
    AnalysisStatus, AnnotationRecord, ClassificationState, MolecularComponent, SemanticIndex,
    SemanticIndexGroup, SystemRecord, TrajectoryRequirements,
)
from analysis.campaign.structure.molecules import partition_from_dump
from tests.analysis.campaign._toy_assembly import build

requires_gmx = pytest.mark.skipif(shutil.which("gmx") is None, reason="needs gmx")
RES = ClassificationState.RESOLVED

DUMP = """   molblock (0):
      moltype              = 0 "Protein_chain_A"
      #molecules                     = 1
   molblock (1):
      moltype              = 1 "LIG"
      #molecules                     = 1
   molblock (2):
      moltype              = 2 "SOL"
      #molecules                     = 3
   moltype (0):
      name="Protein_chain_A"
      atoms:
         atom (5):
   moltype (1):
      name="LIG"
      atoms:
         atom (2):
   moltype (2):
      name="SOL"
      atoms:
         atom (3):
"""


# ═══════════════════════════════════════════════════════════════════════════════
# Topology partition (pure)
# ═══════════════════════════════════════════════════════════════════════════════

def test_molecule_partition_from_tpr_molblocks():
    p = partition_from_dump(DUMP.splitlines(), {"digest": "d"})
    assert (p.n_molecules, p.n_atoms) == (5, 16)
    assert p.starts == [1, 6, 8, 11, 14, 17] and p.moltypes[2:] == ["SOL"] * 3
    assert p.molecule_of(5) == 0 and p.molecule_of(6) == 1 and p.molecule_of(16) == 4
    one = p.select([1, 2, 3])
    assert one.n_molecules == 1 and not one.complete_molecules and not one.multi_molecule
    two = p.select([4, 5, 6, 7])                          # protein tail + ligand
    assert two.multi_molecule and two.complete_molecules is False
    assert p.closure([4, 6]) == [1, 2, 3, 4, 5, 6, 7]     # complete molecules touched
    with pytest.raises(ValueError):
        p.molecule_of(17)


# ═══════════════════════════════════════════════════════════════════════════════
# Molecule-level periodic image detector (pure geometry)
# ═══════════════════════════════════════════════════════════════════════════════

def test_lattice_decomposition_triclinic():
    from analysis.campaign.diagnostics.molecules import lattice_decomposition
    B = np.array([[5.0, 0, 0], [2.0, 5.0, 0], [1.0, 1.5, 5.0]])
    d = B[1] - B[2] + np.array([0.03, -0.02, 0.01])      # n = (0, 1, -1) + small motion
    n, lat, res = lattice_decomposition(d, B)
    assert n.tolist() == [0, 1, -1] and np.allclose(res, [0.03, -0.02, 0.01])


def test_diluted_aggregate_jump_is_found_per_molecule():
    """HMG case: 1 of 4 molecules changes image → aggregate centre moves 0.25 box
    (below the component detector's 0.30), the molecule detector sees it."""
    from analysis.campaign.diagnostics.molecules import image_changes, split_frames
    from analysis.campaign.diagnostics.motion import classify_displacements
    from analysis.campaign.diagnostics.params import DiagnosticParams
    L = 14.15
    B = np.diag([L, L, L])
    base = [np.array([10.5, 10.6, 4.5]), np.array([9.1, 11.3, 2.9]),
            np.array([5.5, 10.8, 4.4]), np.array([6.3, 9.6, 5.6])]
    frames = 6
    mols = [np.array([b + 0.01 * t for t in range(frames)]) for b in base]
    for t in (2, 3, 5):                                   # chain B changes image (−b)
        mols[1][t] = mols[1][t] - np.array([0, L, 0])
    boxes = [B] * frames
    aggregate = sum(mols) / 4
    assert classify_displacements(aggregate, boxes, DiagnosticParams()) == []   # H1 miss
    ev = image_changes(mols[1], boxes, 1.0)
    assert [e["frames"] for e in ev] == [[1, 2], [3, 4], [4, 5]]
    assert ev[0]["lattice_translation"] == [0, -1, 0]
    assert np.linalg.norm(ev[0]["minimum_image_displacement_nm"]) < 0.1
    assert split_frames(mols, boxes) == [2, 3, 5]
    assert image_changes(mols[0], boxes, 1.0) == []


# ═══════════════════════════════════════════════════════════════════════════════
# Policy
# ═══════════════════════════════════════════════════════════════════════════════

def _ctx(tmp_path, *components, tpr=True):
    from analysis.campaign.trajectory.policy import PolicyContext
    ndx = tmp_path / "p.ndx"
    ndx.write_text("[ Receptor ]\n1\n[ Ligand ]\n2\n[ Membrane ]\n3\n[ System ]\n1 2 3\n")
    top = tmp_path / ("md.tpr" if tpr else "md.gro")
    top.write_text("x")
    return PolicyContext(components=[MolecularComponent(t, t, classification_state=RES)
                                     for t in components],
                         topology_path=str(top), index_path=str(ndx))


def test_cluster_policy_finite_valid_extended_refused(tmp_path):
    from analysis.campaign.trajectory.policy import plan_preprocessing
    from analysis.campaign.models import ObservablePurpose as P
    ctx = _ctx(tmp_path, "receptor", "ligand", "membrane")
    ok = plan_preprocessing(TrajectoryRequirements(requires_whole_molecules=True,
                                                   cluster_groups=("Ligand", "Receptor")),
                            purpose=P.INTER_COMPONENT_GEOMETRY, context=ctx)
    d = {x.operation: x for x in ok.decisions}
    assert ok.executable and d["cluster_assembly"].classification == "valid"
    assert d["cluster_assembly"].intent_source == "auto"          # no intent needed: necessary
    bad = plan_preprocessing(TrajectoryRequirements(cluster_groups=("Membrane",)),
                             purpose=P.INTRAMOLECULAR_SHAPE, context=ctx)
    (c,) = [x for x in bad.decisions if x.operation == "cluster_assembly"]
    assert not bad.executable and c.classification == "refused"
    assert c.rule_id == "cluster.periodically_extended"
    diff = plan_preprocessing(TrajectoryRequirements(cluster_groups=("Receptor",)),
                              purpose=P.DIFFUSION, context=ctx)
    assert not diff.executable
    notpr = plan_preprocessing(TrajectoryRequirements(cluster_groups=("Receptor",)),
                               purpose=P.INTRAMOLECULAR_SHAPE,
                               context=_ctx(tmp_path, "receptor", tpr=False))
    assert [x.classification for x in notpr.decisions] == ["unsupported"]


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS on toy periodic systems
# ═══════════════════════════════════════════════════════════════════════════════

def _system(d: Path, membrane_group=None, scratch: Path = None) -> SystemRecord:
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.results import atom_set_hash
    rec = SystemRecord("toy", "c", "r", topology_path=str(d / "md.tpr"),
                       structure_path=str(d / "md.gro"),
                       production_trajectory_paths=[str(d / "md.xtc")])
    rec.source_fingerprints = {"trajectory": fingerprint_file(d / "md.xtc"),
                               "topology": fingerprint_file(d / "md.tpr")}
    rec.components = [MolecularComponent("receptor", "receptor", classification_state=RES),
                      MolecularComponent("ligand", "ligand", classification_state=RES)]
    ndx = d / "index.ndx"
    if membrane_group:                                   # a private copy: fixtures are shared
        own = scratch / "index_membrane.ndx"
        own.write_text(ndx.read_text() + f"[ Membrane ]\n{membrane_group}\n")
        ndx = own
        rec.components.append(MolecularComponent("membrane", "membrane", classification_state=RES))
    names = [l.strip("[] \n") for l in ndx.read_text().splitlines() if l.startswith("[")]
    rec.semantic_index = SemanticIndex(path=str(ndx), groups=[SemanticIndexGroup(n, 1) for n in names])
    for aid, atoms in (("span", [3, 4, 5, 6]), ("one", [1, 2, 3])):
        rec.annotations.append(AnnotationRecord(
            aid, "residue_set", "custom", "user_yaml", "active", atom_ids=atoms,
            n_atoms=len(atoms), atoms_sha256=atom_set_hash(atoms), identity=f"id-{aid}",
            group_name=f"Ann_{aid}"))
    return rec


def _run(rec, oid, params, tmp: Path, force_raw=False):
    from analysis.campaign.observables import registry
    from analysis.campaign.observables.base import AnalysisContext
    from analysis.campaign.orchestration.study_analyzer import (
        clear_view_cache, observable_requirements, plan_view,
    )
    registry.ensure_loaded()
    clear_view_cache()
    spec = registry.get(oid)
    req, why = observable_requirements(spec, rec, params)
    if force_raw:                                              # the pre-13.5 semantics
        req = spec.trajectory_requirements(params)
    assert why is None, why
    view, _ = plan_view(rec, req, tmp / "sys", purpose=spec.purpose)
    ctx = AnalysisContext(system=rec, semantic_index=rec.semantic_index, trajectory_view=view,
                          topology_path=rec.topology_path, structure_path=rec.structure_path,
                          output_dir=tmp / "out" / oid / ("raw" if force_raw else "new"),
                          parameters=params)
    return req, view, spec.execute(ctx)


def _values(res):
    from analysis.campaign.results import read_column
    return read_column(res.arrays[0].storage)


@pytest.fixture(scope="module", params=["rect", "tri"])
def toy(request, tmp_path_factory):
    if shutil.which("gmx") is None:
        pytest.skip("needs gmx")
    return build(tmp_path_factory.mktemp(request.param), request.param == "tri")


@requires_gmx
def test_multi_molecule_rg_is_image_invariant(toy, tmp_path):
    rec = _system(toy)
    params = {"rg": {"selection": "component:receptor"}}
    req, view, res = _run(rec, "rg", params, tmp_path)
    assert req.cluster_groups == ("Receptor",) and view.kind == "clustered" and view.safe
    decisions = {d["operation"]: d["classification"] for d in view.decisions}
    assert decisions["cluster_assembly"] == "valid"
    v = _values(res)
    assert res.status == AnalysisStatus.SUCCESS and len(v) == 3
    assert abs(v[1] - v[0]) < 1e-5 and abs(v[2] - v[0]) < 1e-5  # A == B (lattice image)
    assert res.definition_evidence["assembly"]["cluster_groups"] == ["Receptor"]
    t = res.arrays[0].axis("time")
    assert t.alignment.mode == "one_row_per_frame"               # clustered view time index
    _, _, old = _run(rec, "rg", params, tmp_path, force_raw=True)
    ov = _values(old)
    assert ov[1] > 5 * ov[0]                                     # the B1 artefact, reproduced


@requires_gmx
def test_pair_union_com_distance_is_image_invariant(toy, tmp_path):
    rec = _system(toy)
    params = {"com-distance": {"selection_a": "component:ligand",
                               "selection_b": "component:receptor"}}
    req, view, res = _run(rec, "com-distance", params, tmp_path)
    assert req.cluster_groups == ("Ligand", "Receptor")          # A ∪ B, not B alone
    v = _values(res)
    assert abs(v[1] - v[0]) < 2e-3 and abs(v[2] - v[0]) < 2e-3  # distance printed at 1e-3
    _, _, old = _run(rec, "com-distance", params, tmp_path, force_raw=True)
    assert abs(_values(old)[1] - v[0]) > 0.3


@requires_gmx
def test_pairwise_and_single_molecule_observables_are_unchanged(toy, tmp_path):
    rec = _system(toy)
    for oid, params in (
            ("min-distance", {"min-distance": {"selection_a": "component:ligand",
                                               "selection_b": "component:receptor"}}),
            ("sasa", {"sasa": {"selection": "component:receptor", "surface": "component:receptor"}}),
            ("rg", {"rg": {"selection": "annotation:one"}}),          # one molecule
            ("rg", {"rg": {"selection": "component:ligand"}})):
        req, view, res = _run(rec, oid, params, tmp_path)
        assert not req.cluster_groups and view.kind == "raw", oid
        assert "assembly" not in res.definition_evidence
        _, _, old = _run(rec, oid, params, tmp_path, force_raw=True)
        assert _values(res) == _values(old), oid                  # identical numerics
        assert res.definition_token == old.definition_token      # identical identity


@requires_gmx
def test_multi_molecule_rmsd_is_assembled_before_fit(toy, tmp_path):
    rec = _system(toy)
    req, view, res = _run(rec, "rmsd-receptor", {}, tmp_path)
    assert req.cluster_groups == ("Receptor_Backbone",) and view.kind == "clustered"
    assert res.status == AnalysisStatus.SUCCESS, res.message
    assert max(_values(res)) < 1e-3                              # A, B, A: same structure
    _, _, old = _run(rec, "rmsd-receptor", {}, tmp_path, force_raw=True)
    assert max(_values(old)) > 0.5                               # whole-per-molecule only


@requires_gmx
def test_annotation_spanning_molecules_assembles_complete_molecules(toy, tmp_path):
    from analysis.campaign.trajectory.preprocessor import _assembly_index
    rec = _system(toy)
    req, view, res = _run(rec, "rg", {"rg": {"selection": "annotation:span"}}, tmp_path)
    assert req.cluster_groups == ("Ann_span",) and res.status == AnalysisStatus.SUCCESS
    v = _values(res)
    assert abs(v[1] - v[0]) < 1e-5
    ndx, err = _assembly_index(("Ann_span",), rec.topology_path, rec.semantic_index,
                               tmp_path / "asm")
    body = ndx.read_text()
    assert "[ Assembly ]\n1 2 3 4 5 6 7 8\n" in body            # closure: both whole subunits


@requires_gmx
def test_membrane_global_geometry_is_refused(toy, tmp_path):
    from analysis.campaign.observables import registry
    from analysis.campaign.orchestration.study_analyzer import observable_requirements
    registry.ensure_loaded()
    rec = _system(toy, membrane_group="1 2 3 4", scratch=tmp_path)  # pretend subunit 1 is membrane
    rec.components[0].classification_state = RES
    _, why = observable_requirements(registry.get("rg"), rec, {"rg": {"selection": "component:membrane"}})
    assert why and "periodically extended" in why
    _, why = observable_requirements(
        registry.get("com-distance"), rec,
        {"com-distance": {"selection_a": "component:membrane", "selection_b": "component:ligand"}})
    assert why and "membrane" in why
    req, why = observable_requirements(                   # pairwise geometry stays available
        registry.get("min-distance"), rec,
        {"min-distance": {"selection_a": "component:membrane", "selection_b": "component:ligand"}})
    assert why is None and not req.cluster_groups


@requires_gmx
def test_cluster_view_identity_and_old_result_is_incompatible(toy, tmp_path):
    from analysis.campaign.compatibility import Candidate, ProvenanceTier, evaluate
    rec = _system(toy)
    params = {"rg": {"selection": "component:receptor"}}
    _, view, new = _run(rec, "rg", params, tmp_path)
    _, raw_view, old = _run(rec, "rg", params, tmp_path, force_raw=True)
    assert view.cache_key != raw_view.cache_key
    from analysis.campaign.results import atom_set_hash
    assert "cluster=Receptor" in view.requirements.cache_token()
    rep = evaluate(Candidate(source="old", tier=ProvenanceTier.NATIVE, result=old),
                   requested_evidence=new.definition_evidence, requested_view_ref=view.cache_key,
                   requested_schema=[new.arrays[0]], requested_timeline=[0.0, 10.0, 20.0])
    assert rep.status == "incompatible"
    failed = {c.dimension for c in rep.checks if c.outcome == "mismatch"}
    assert {"definition", "view"} <= failed and "assembly" in rep.reason


@requires_gmx
def test_molecule_detector_on_toy_and_not_applicable_for_membrane(toy, tmp_path):
    from analysis.campaign.diagnostics.context import DiagnosticContext
    from analysis.campaign.diagnostics.registry import run_diagnostics
    ctx = DiagnosticContext(trajectory_path=str(toy / "md.xtc"), structure_path=str(toy / "md.tpr"),
                            index_path=str(toy / "index.ndx"),
                            groups={"receptor": "Receptor", "ligand": "Ligand", "system": "System"},
                            cache_dir=str(tmp_path / "cache"))
    rep = run_diagnostics(ctx, only=["molecule_periodic_image_change", "component_periodic_jump"])
    runs = {r.detector_id: r for r in rep.detectors}
    assert runs["molecule_periodic_image_change"].status == "ran"
    (d,) = [x for x in rep.diagnostics if x.code == "molecule_periodic_image_change"]
    assert d.scope == "component:receptor" and d.measured["n_frames_split"] == 1
    assert d.interpretation == "periodic_image_change_of_intact_molecules"
    ev = d.measured["events"][0]
    assert ev["molecule"] == 2 and ev["frames"] == [0, 1] and "box_nm" in ev
    assert sum(ev["lattice_translation"]) != 0
    ctx2 = DiagnosticContext(trajectory_path=str(toy / "md.xtc"), structure_path=str(toy / "md.tpr"),
                             index_path=str(toy / "index.ndx"), groups={"membrane": "Receptor"},
                             cache_dir=str(tmp_path / "cache2"))
    rep2 = run_diagnostics(ctx2, only=["molecule_periodic_image_change"])
    (r,) = rep2.detectors
    assert r.status == "not_applicable" and "periodically extended" in r.reason
