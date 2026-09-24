"""Phase 4.6 — nucleic-acid detection fix, conservative ligand resolution, and
explicit manifest resolution of non-polymer components.

Structural evidence only (offline); no binding-site or functional inference.
All fixtures are synthetic — exact-behaviour tests, not real-world validation.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from analysis.campaign.gmx import run_gmx
from analysis.campaign.manifest import load_manifest, write_manifest
from analysis.campaign.models import (
    Ambiguity, CampaignRunResult, ClassificationState, ComponentType, DetectorRunStatus,
    StudyManifest, SystemRecord,
)
from analysis.campaign.orchestration.study_analyzer import _ensure_semantic_index
from analysis.campaign.structure.component_detector import (
    LIGAND_MIN_HEAVY_ATOMS, detect_components,
)
from analysis.campaign.structure.index_builder import build_semantic_index
from analysis.campaign.structure.index_groups import parse_index_groups
from tests.analysis.campaign.conftest import requires_gmx

AA = ["ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE"]


# ═══════════════════════════════════════════════════════════════════════════════
# Fixture builder: records are (record, atom, resname, chain, resseq, x, y, z) in nm
# ═══════════════════════════════════════════════════════════════════════════════

def protein(chain="A", n=30, x0=8.9):
    return [("ATOM", "CA", AA[i % 10], chain, i + 1, x0 + 0.2 * (i % 5), 4.0 + 0.1 * i, 5.0)
            for i in range(n)]


def small_molecule(resname, chain="A", resseq=301, heavy=12, hydrogens=3, x=9.7):
    atoms = [("HETATM", f"C{k + 1}", resname, chain, resseq, x, 5.0 + 0.05 * k, 5.2)
             for k in range(heavy)]
    atoms += [("HETATM", f"H{k + 1}", resname, chain, resseq, x, 5.0 + 0.05 * k, 5.3)
              for k in range(hydrogens)]
    return atoms


def nucleic(chain, names, resseq0=1, x0=2.0):
    out = []
    for i, rn in enumerate(names):
        for k, an in enumerate(("P", "C4'", "N1")):
            out.append(("ATOM", an, rn, chain, resseq0 + i, x0 + 0.1 * i, 2.0, 2.0 + 0.1 * k))
    return out


def solvent():
    return [("HETATM", "OW", "SOL", "W", 500, 3.0, 3.0, 3.0),
            ("HETATM", "NA", "NA", "I", 600, 3.5, 3.0, 3.0),
            ("HETATM", "CL", "CL", "I", 601, 3.6, 3.0, 3.0)]


def write_pdb(path: Path, atoms) -> Path:
    lines = ["CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1"]
    for s, (rec, an, rn, ch, rs, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{rec:<6}{s:>5} {an:<4} {rn:<3} {ch:1}{rs:>4}    "
                     f"{10 * x:8.3f}{10 * y:8.3f}{10 * z:8.3f}  1.00  0.00           C")
    path.write_text("\n".join(lines + ["END"]) + "\n")
    return path


def atom_ids_of(atoms, pred) -> list[int]:
    return [i for i, a in enumerate(atoms, 1) if pred(a)]


def by_type(det, ctype):
    return [c for c in det.components if c.component_type == ctype]


# ═══════════════════════════════════════════════════════════════════════════════
# Nucleic acids
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("names", [["DA", "DT", "DG", "DC"] * 3, ["A", "U", "G", "C"] * 3],
                         ids=["dna", "rna"])
def test_pure_nucleic_acid_is_nucleic_acid(tmp_path, names):
    atoms = nucleic("N", names)
    det = detect_components(structure_path=write_pdb(tmp_path / "na.pdb", atoms))
    (na,) = by_type(det, ComponentType.NUCLEIC_ACID)
    assert na.classification_state == ClassificationState.RESOLVED
    assert na.chain_ids == ["N"] and na.residue_count == len(names)
    assert not by_type(det, ComponentType.RECEPTOR)
    assert any(e.kind == "nucleotide_polymer" for e in na.evidence)


def test_double_stranded_dna_is_one_nucleic_component(tmp_path):
    atoms = nucleic("M", ["DA", "DT", "DG"] * 3) + nucleic("N", ["DC", "DA", "DT"] * 3, x0=4.0)
    det = detect_components(structure_path=write_pdb(tmp_path / "ds.pdb", atoms))
    (na,) = by_type(det, ComponentType.NUCLEIC_ACID)
    assert na.chain_ids == ["M", "N"] and na.classification_state == ClassificationState.RESOLVED
    assert not by_type(det, ComponentType.RECEPTOR)


@pytest.mark.parametrize("names", [["DA", "DT", "DG", "DC"] * 3, ["A", "U", "G", "C"] * 3],
                         ids=["dna", "rna"])
def test_protein_plus_nucleic_acid_keeps_both_roles(tmp_path, names):
    atoms = protein() + nucleic("N", names)
    det = detect_components(structure_path=write_pdb(tmp_path / "pn.pdb", atoms))
    (rec,) = by_type(det, ComponentType.RECEPTOR)
    (na,) = by_type(det, ComponentType.NUCLEIC_ACID)
    assert rec.chain_ids == ["A"] and rec.classification_state == ClassificationState.RESOLVED
    assert na.chain_ids == ["N"]
    assert not by_type(det, ComponentType.PEPTIDE)            # DNA is not scored as a peptide
    (cx,) = by_type(det, ComponentType.COMPLEX)
    assert cx.chain_ids == ["A", "N"]                         # Complex = union of polymer chains
    assert not [c for c in det.components if c.label == "no partner"]


def test_single_nucleotide_is_not_a_polymer_nor_ligand(tmp_path):
    atoms = protein() + [("HETATM", "P", "DA", "X", 900, 2.0, 2.0, 2.0),
                         ("HETATM", "N1", "DA", "X", 900, 2.1, 2.0, 2.0)]
    det = detect_components(structure_path=write_pdb(tmp_path / "n1.pdb", atoms))
    assert not by_type(det, ComponentType.NUCLEIC_ACID)
    assert not by_type(det, ComponentType.LIGAND)


def test_nucleotide_cofactor_stays_cofactor(tmp_path):
    atoms = protein() + small_molecule("ATP", heavy=31)
    det = detect_components(structure_path=write_pdb(tmp_path / "atp.pdb", atoms))
    assert by_type(det, ComponentType.COFACTOR)[0].resnames == ["ATP"]
    assert not by_type(det, ComponentType.NUCLEIC_ACID) and not by_type(det, ComponentType.LIGAND)


# ═══════════════════════════════════════════════════════════════════════════════
# Ligand — conservative automatic resolution
# ═══════════════════════════════════════════════════════════════════════════════

def test_single_unambiguous_ligand_resolved(tmp_path):
    atoms = protein() + small_molecule("LIG") + solvent()
    det = detect_components(structure_path=write_pdb(tmp_path / "pl.pdb", atoms))
    (lig,) = by_type(det, ComponentType.LIGAND)
    assert lig.classification_state == ClassificationState.RESOLVED
    assert lig.resnames == ["LIG"] and lig.chain_ids == ["A"]
    assert lig.atom_ids == atom_ids_of(atoms, lambda a: a[2] == "LIG")
    kinds = {e.kind for e in lig.evidence}
    assert {"non_polymer_candidate", "excluded_known_classes",
            "sole_plausible_candidate"} <= kinds
    excl = next(e for e in lig.evidence if e.kind == "excluded_known_classes").value
    assert excl["water"] == 1 and excl["ion"] == 2
    cand = next(e for e in lig.evidence if e.kind == "non_polymer_candidate").value
    assert cand["heavy_atoms"] == 12 and cand["n_atoms"] == 15


def test_unknown_resname_ligand_unl(tmp_path):
    atoms = protein() + small_molecule("UNL", chain="B", resseq=1)
    det = detect_components(structure_path=write_pdb(tmp_path / "unl.pdb", atoms))
    (lig,) = by_type(det, ComponentType.LIGAND)
    assert lig.resnames == ["UNL"] and lig.classification_state == ClassificationState.RESOLVED


def test_known_classes_never_become_ligand(tmp_path):
    atoms = (protein() + solvent() + small_molecule("HEM", heavy=43)
             + [("HETATM", "P", "POPC", "L", 700, 1.0, 1.0, 1.0)]
             + nucleic("N", ["DA", "DT", "DG"] * 2))
    det = detect_components(structure_path=write_pdb(tmp_path / "known.pdb", atoms))
    assert not by_type(det, ComponentType.LIGAND)
    assert by_type(det, ComponentType.COFACTOR)


def test_small_additive_not_resolved_and_reported(tmp_path):
    atoms = protein() + small_molecule("GOL", heavy=LIGAND_MIN_HEAVY_ATOMS - 1)
    det = detect_components(structure_path=write_pdb(tmp_path / "gol.pdb", atoms))
    assert not by_type(det, ComponentType.LIGAND)
    assert any("GOL" in w and "below" in w for w in det.warnings)


def test_small_additive_does_not_block_a_real_ligand(tmp_path):
    atoms = protein() + small_molecule("LIG") + small_molecule("GOL", resseq=302, heavy=6)
    det = detect_components(structure_path=write_pdb(tmp_path / "lg.pdb", atoms))
    (lig,) = by_type(det, ComponentType.LIGAND)
    assert lig.classification_state == ClassificationState.RESOLVED and lig.resnames == ["LIG"]
    small = next(e for e in lig.evidence if e.kind == "small_molecules_not_considered").value
    assert small[0]["resname"] == "GOL"


@pytest.mark.parametrize("second", [("LI2", 302), ("LIG", 302)], ids=["distinct", "two_copies"])
def test_multiple_plausible_candidates_ambiguous_not_merged(tmp_path, second):
    atoms = protein() + small_molecule("LIG") + small_molecule(second[0], resseq=second[1])
    det = detect_components(structure_path=write_pdb(tmp_path / "two.pdb", atoms))
    (lig,) = by_type(det, ComponentType.LIGAND)
    assert lig.classification_state == ClassificationState.AMBIGUOUS
    assert lig.atom_ids is None
    cands = next(e for e in lig.evidence if e.kind == "multiple_plausible_candidates").value
    assert len(cands) == 2


def test_internal_modified_residue_not_a_candidate(tmp_path):
    prot = protein()
    atoms = prot[:10] + [("ATOM", "CA", "MLY", "A", 11, 9.0, 5.0, 5.0)] + \
        [(r, an, rn, ch, rs + 1, x, y, z) for (r, an, rn, ch, rs, x, y, z) in prot[10:]]
    det = detect_components(structure_path=write_pdb(tmp_path / "mly.pdb", atoms))
    assert not by_type(det, ComponentType.LIGAND)


def test_no_polymer_no_ligand(tmp_path):
    det = detect_components(structure_path=write_pdb(tmp_path / "lonely.pdb",
                                                     small_molecule("LIG") + solvent()))
    assert not by_type(det, ComponentType.LIGAND)
    assert any("no polymer" in w for w in det.warnings)


def test_existing_protein_peptide_outputs_unchanged(tmp_path):
    from tests.analysis.campaign.conftest import write_pdb as write_ca_pdb
    det = detect_components(structure_path=write_ca_pdb(tmp_path / "rp.pdb", {"A": 120, "B": 25}))
    summary = [(c.component_type, c.classification_state, c.chain_ids) for c in det.components]
    assert summary == [
        ("receptor", "ambiguous", ["A"]), ("peptide", "ambiguous", ["B"]),
        ("complex", "resolved", ["A", "B"])]


# ═══════════════════════════════════════════════════════════════════════════════
# Explicit manifest resolution
# ═══════════════════════════════════════════════════════════════════════════════

def _manifest(tmp_path, atoms, resolution, *, index_path=None):
    pdb = write_pdb(tmp_path / "sys.pdb", atoms)
    det = detect_components(structure_path=pdb)
    rec = SystemRecord(system_id="s1", condition_id="c", replicate_id="r",
                       structure_path=str(pdb), index_path=index_path,
                       components=det.components)
    m = StudyManifest(study_root=str(tmp_path), generated_utc="", simforge_version="t",
                      systems=[rec],
                      ambiguities=[Ambiguity(system_id="s1", kind="component_classification",
                                             message="declared", resolution=resolution)])
    paths = write_manifest(m, tmp_path / "m")
    return load_manifest(paths["yaml"]).system("s1"), atoms


def test_explicit_resolution_resolves_ambiguous_ligand(tmp_path):
    atoms = protein() + small_molecule("LIG") + small_molecule("LI2", resseq=302)
    rec, atoms = _manifest(tmp_path, atoms, "ligand=resname:LI2")
    (lig,) = [c for c in rec.components if c.component_type == ComponentType.LIGAND]
    assert lig.classification_state == ClassificationState.RESOLVED
    assert lig.atom_ids == atom_ids_of(atoms, lambda a: a[2] == "LI2")
    ev = next(e for e in lig.evidence if e.kind == "explicit_resolution").value
    assert ev["intent_source"] == "manifest" and ev["selection"] == {"resname": "LI2"}
    assert ev["n_atoms"] == 15 and ev["atoms_sha256"]
    sup = next(e for e in lig.evidence if e.kind == "superseded_automatic_interpretation").value
    assert sup[0]["classification_state"] == ClassificationState.AMBIGUOUS
    assert rec.user_overridden


def test_explicit_resolution_wins_conflict_with_auto(tmp_path):
    atoms = protein() + small_molecule("LIG") + small_molecule("HEM", resseq=400, heavy=43)
    rec, atoms = _manifest(tmp_path, atoms, "ligand=chain:A,resname:HEM,resid:400")
    ligs = [c for c in rec.components if c.component_type == ComponentType.LIGAND]
    assert len(ligs) == 1                                      # auto LIG superseded, not kept
    lig = ligs[0]
    assert lig.atom_ids == atom_ids_of(atoms, lambda a: a[2] == "HEM")
    conflict = next(e for e in lig.evidence if e.kind == "conflict_with_automatic").value
    assert any(c["component_type"] == "ligand" for c in conflict)
    assert any(c["component_type"] == "cofactor" for c in conflict)
    assert not [c for c in rec.components if c.component_type == ComponentType.COFACTOR]
    assert any("explicit" in w for w in lig.warnings)


def test_explicit_resolution_by_index_group(tmp_path):
    atoms = protein() + small_molecule("GOL", heavy=6)       # too small for auto resolution
    ids = atom_ids_of(atoms, lambda a: a[2] == "GOL")
    ndx = tmp_path / "user.ndx"
    ndx.write_text("[ MyLig ]\n" + " ".join(map(str, ids)) + "\n")
    rec, _ = _manifest(tmp_path, atoms, "ligand=group:MyLig", index_path=str(ndx))
    (lig,) = [c for c in rec.components if c.component_type == ComponentType.LIGAND]
    assert lig.atom_ids == ids
    assert next(e for e in lig.evidence if e.kind == "explicit_resolution").value["selection"] \
        == {"group": "MyLig"}


def test_explicit_nucleic_acid_and_bad_selection(tmp_path):
    atoms = protein() + nucleic("N", ["DA", "DT", "DG"] * 2)
    rec, _ = _manifest(tmp_path, atoms, "nucleic_acid=chain:N;cofactor=resname:ZZZ")
    (na,) = [c for c in rec.components if c.component_type == ComponentType.NUCLEIC_ACID]
    assert na.classification_state == ClassificationState.RESOLVED
    assert any(e.kind == "explicit_resolution" for e in na.evidence)
    assert not [c for c in rec.components if c.component_type == ComponentType.COFACTOR]
    assert any("ZZZ" in w.message for w in rec.warnings)       # matched no atoms: reported


def test_legacy_chain_resolution_unchanged(tmp_path):
    from tests.analysis.campaign.conftest import write_pdb as write_ca_pdb
    pdb = write_ca_pdb(tmp_path / "rp.pdb", {"A": 120, "B": 25})
    det = detect_components(structure_path=pdb)
    rec = SystemRecord(system_id="s1", condition_id="c", replicate_id="r",
                       structure_path=str(pdb), components=det.components)
    m = StudyManifest(study_root=str(tmp_path), generated_utc="", simforge_version="t",
                      systems=[rec], ambiguities=[Ambiguity(
                          system_id="s1", kind="component_classification", message="x",
                          resolution="receptor=B;peptide=A")])
    loaded = load_manifest(write_manifest(m, tmp_path / "m")["yaml"]).system("s1")
    r = next(c for c in loaded.components if c.component_type == ComponentType.RECEPTOR)
    assert r.chain_ids == ["B"] and r.classification_state == ClassificationState.RESOLVED


# ═══════════════════════════════════════════════════════════════════════════════
# Semantic index + diagnostics through the study-style path (gmx)
# ═══════════════════════════════════════════════════════════════════════════════

X_WRAP = [9.70, 9.80, 0.10, 0.15, 0.20]


def _trajectory(tmp_path, atoms_for_frame) -> Path:
    gl = []
    for f in range(len(X_WRAP)):
        a = atoms_for_frame(f)
        gl += [f"fixture t= {20.0 * f:.5f} step= {f}", str(len(a))]
        for s, (rec, an, rn, ch, rs, x, y, z) in enumerate(a, 1):
            gl.append(f"{rs % 100000:>5}{rn:<5}{an:>5}{s:>5}{x:8.3f}{y:8.3f}{z:8.3f}")
        gl.append("  10.00000  10.00000  10.00000")
    gro = tmp_path / "frames.gro"
    gro.write_text("\n".join(gl) + "\n")
    xtc = tmp_path / "md.xtc"
    res = run_gmx(["trjconv", "-f", str(gro), "-s", str(gro), "-o", str(xtc)], stdin="0\n")
    assert res.ok, res.stderr[-300:]
    gro.unlink()
    return xtc


def _study(tmp_path, atoms_for_frame, components=None):
    atoms = atoms_for_frame(0)
    pdb = write_pdb(tmp_path / "reference.pdb", atoms)
    xtc = _trajectory(tmp_path, atoms_for_frame)
    comps = components if components is not None else detect_components(structure_path=pdb).components
    rec = SystemRecord(system_id="s", condition_id="c", replicate_id="r",
                       structure_path=str(pdb), production_trajectory_paths=[str(xtc)],
                       components=comps)
    _ensure_semantic_index(rec, tmp_path / "out", "gmx",
                           CampaignRunResult(study_root="x", output_dir="y"))
    return rec, atoms


@requires_gmx
def test_auto_ligand_study_path_index_and_diagnostics(tmp_path):
    from analysis.campaign.diagnostics import diagnose_system
    rec, atoms = _study(tmp_path, lambda f: protein() + small_molecule("LIG", x=X_WRAP[f]))
    groups = {g.name: g.atom_ids for g in parse_index_groups(rec.semantic_index.path)}
    assert list(groups["Ligand"]) == atom_ids_of(atoms, lambda a: a[2] == "LIG")
    rep = diagnose_system(rec)
    run = next(r for r in rep.detectors if r.detector_id == "partner_receptor_separation")
    assert run.status == DetectorRunStatus.RAN
    sep = [d for d in rep.diagnostics if d.code == "ligand_receptor_periodic_separation"]
    assert sep and sep[0].provenance["partner_group_evidence"]["n_atoms"] == 15


@requires_gmx
def test_explicit_ligand_study_path_detector_identity(tmp_path):
    from analysis.campaign.diagnostics import diagnose_system

    def frame(f):
        return (protein() + small_molecule("LIG", x=X_WRAP[f])
                + small_molecule("LI2", resseq=302, x=3.0))
    pdb = write_pdb(tmp_path / "sys.pdb", frame(0))
    det = detect_components(structure_path=pdb)
    rec0 = SystemRecord(system_id="s1", condition_id="c", replicate_id="r",
                        structure_path=str(pdb), components=det.components)
    m = StudyManifest(study_root=str(tmp_path), generated_utc="", simforge_version="t",
                      systems=[rec0], ambiguities=[Ambiguity(
                          system_id="s1", kind="component_classification", message="x",
                          resolution="ligand=resname:LIG")])
    comps = load_manifest(write_manifest(m, tmp_path / "m")["yaml"]).system("s1").components
    rec, atoms = _study(tmp_path, frame, components=comps)
    rep = diagnose_system(rec)
    run = next(r for r in rep.detectors if r.detector_id == "partner_receptor_separation")
    assert run.status == DetectorRunStatus.RAN
    lig_ev = rep.groups["ligand"]["evidence"]
    from analysis.campaign.results import atom_set_hash
    assert lig_ev["atoms_sha256"] == atom_set_hash(atom_ids_of(atoms, lambda a: a[2] == "LIG"))


@requires_gmx
def test_nucleic_acid_study_path_index(tmp_path):
    atoms = protein() + nucleic("N", ["DA", "DT", "DG"] * 2)
    pdb = write_pdb(tmp_path / "pn.pdb", atoms)
    det = detect_components(structure_path=pdb)
    idx = build_semantic_index(structure_path=pdb, components=det.components,
                               out_ndx=tmp_path / "i.ndx")
    groups = {g.name: g.atom_ids for g in parse_index_groups(idx.path)}
    assert list(groups["NucleicAcid"]) == atom_ids_of(atoms, lambda a: a[3] == "N")
    assert list(groups["Receptor"]) == atom_ids_of(atoms, lambda a: a[3] == "A")


# ═══════════════════════════════════════════════════════════════════════════════
# Architectural invariants (documented, not redesigned)
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
def test_group_presence_is_not_resolution(tmp_path):
    """An AMBIGUOUS receptor still gets a ``Receptor`` group (pre-existing
    behaviour).  Policy must read ``classification_state``, never infer
    resolution from the presence of an index group."""
    from tests.analysis.campaign.conftest import write_pdb as write_ca_pdb
    pdb = write_ca_pdb(tmp_path / "rp.pdb", {"A": 120, "B": 25})
    det = detect_components(structure_path=pdb)
    rec = by_type(det, ComponentType.RECEPTOR)[0]
    assert rec.classification_state == ClassificationState.AMBIGUOUS
    idx = build_semantic_index(structure_path=pdb, components=det.components,
                               out_ndx=tmp_path / "i.ndx")
    assert "Receptor" in idx.group_names()


@requires_gmx
def test_complex_membership_is_chain_union_not_a_binding_relation(tmp_path):
    """``Complex`` unions polymer *chains*: a cofactor sharing chain A falls
    inside it, one on its own chain does not.  Complex membership is not
    evidence of a receptor–partner relation."""
    atoms = protein() + small_molecule("HEM", chain="A", resseq=400, heavy=10) \
        + small_molecule("HEM", chain="H", resseq=401, heavy=10)
    pdb = write_pdb(tmp_path / "cx.pdb", atoms)
    det = detect_components(structure_path=pdb)
    idx = build_semantic_index(structure_path=pdb, components=det.components,
                               out_ndx=tmp_path / "i.ndx")
    cx = set(next(g for g in parse_index_groups(idx.path) if g.name == "Complex").atom_ids)
    hem_a = set(atom_ids_of(atoms, lambda a: a[2] == "HEM" and a[3] == "A"))
    hem_h = set(atom_ids_of(atoms, lambda a: a[2] == "HEM" and a[3] == "H"))
    assert hem_a <= cx and not (hem_h & cx)
