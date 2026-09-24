"""Phase 4.5 — resolved Ligand / Cofactor / NucleicAcid components become
semantic index groups, so downstream diagnostics/preprocessing can use them.

Before the fix, a RESOLVED cofactor (the detector emits these from its
residue-class table) produced no group and ``partner_receptor_separation`` was
NOT_APPLICABLE through the study path.  Goldens were captured pre-fix.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from analysis.campaign.diagnostics import (
    DiagnosticContext, detector_identity, diagnose_system, groups_from_semantic_index,
)
from analysis.campaign.diagnostics import registry as dreg
from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    CampaignRunResult, ClassificationState, ComponentType, DetectorRunStatus,
    MolecularComponent, SemanticIndex, SemanticIndexGroup, SystemRecord,
)
from analysis.campaign.orchestration.study_analyzer import _ensure_semantic_index
from analysis.campaign.structure.component_detector import detect_components
from analysis.campaign.structure.index_builder import build_semantic_index
from analysis.campaign.structure.index_groups import parse_index_groups
from analysis.campaign.trajectory import requirements as reqs
from tests.analysis.campaign.conftest import requires_gmx

AA = ["ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE"]
X_WRAP = [9.70, 9.80, 0.10, 0.15, 0.20]


def _ordered(ids) -> str:
    return hashlib.sha256(",".join(map(str, ids)).encode()).hexdigest()[:24]


def _groups(ndx) -> dict[str, tuple]:
    return {g.name: g.atom_ids for g in parse_index_groups(ndx)}


def _atoms(frame: int, het_resname: str, extra_chain: list | None = None):
    out = [("A", i + 1, AA[i % 10], "CA", 8.9 + 0.8 * ((i % 5) / 4), 4.0 + 0.1 * i, 5.0, "ATOM")
           for i in range(30)]
    for k, n in enumerate(("C1", "C2", "C3", "C4")):
        out.append(("B", 31, het_resname, n, X_WRAP[frame], 5.0 + 0.1 * k, 5.0, "HETATM"))
    out.append(("W", 32, "SOL", "OW", 3.0, 3.0, 3.0, "HETATM"))
    for a in extra_chain or []:
        out.append(a)
    return out


def _write_system(d: Path, het_resname: str, *, extra_chain=None, trajectory=True):
    d.mkdir(parents=True, exist_ok=True)
    pdb = d / "reference.pdb"
    lines = ["CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1"]
    for s, (ch, rs, rn, an, x, y, z, rec) in enumerate(_atoms(0, het_resname, extra_chain), 1):
        lines.append(f"{rec:<6}{s:>5} {an:<4} {rn:<3} {ch:1}{rs:>4}    "
                     f"{10 * x:8.3f}{10 * y:8.3f}{10 * z:8.3f}  1.00  0.00           C")
    pdb.write_text("\n".join(lines + ["END"]) + "\n")
    if not trajectory:
        return pdb, None
    gl = []
    for f in range(len(X_WRAP)):
        a = _atoms(f, het_resname, extra_chain)
        gl += [f"fixture t= {20.0 * f:.5f} step= {f}", str(len(a))]
        for s, (ch, rs, rn, an, x, y, z, rec) in enumerate(a, 1):
            gl.append(f"{rs:>5}{rn:<5}{an:>5}{s:>5}{x:8.3f}{y:8.3f}{z:8.3f}")
        gl.append("  10.00000  10.00000  10.00000")
    gro = d / "frames.gro"
    gro.write_text("\n".join(gl) + "\n")
    xtc = d / "md.xtc"
    res = run_gmx(["trjconv", "-f", str(gro), "-s", str(gro), "-o", str(xtc)], stdin="0\n")
    assert res.ok, res.stderr[-300:]
    gro.unlink()
    return pdb, xtc


def _record(pdb, xtc, components) -> SystemRecord:
    return SystemRecord(system_id="s", condition_id="c", replicate_id="r",
                        structure_path=str(pdb),
                        production_trajectory_paths=[str(xtc)] if xtc else [],
                        components=components)


def _study_index(rec: SystemRecord, out: Path) -> SemanticIndex:
    """The exact index path `simforge study analyze` takes."""
    _ensure_semantic_index(rec, out, "gmx", CampaignRunResult(study_root="x", output_dir="y"))
    return rec.semantic_index


def _ligand(state=ClassificationState.RESOLVED, resnames=("LIG",), label="ligand"):
    return MolecularComponent(component_type=ComponentType.LIGAND, label=label,
                              resnames=list(resnames), classification_state=state)


# Captured from the pre-fix builder on the receptor + HEM fixture.
_GOLDEN_COFACTOR_FIXTURE = {
    "Water": "9f14025af0065b30e47e23eb",
    "Receptor": "a8da6dc1099b8b38805d26f0", "Receptor_Backbone": "a8da6dc1099b8b38805d26f0",
    "Receptor_CA": "a8da6dc1099b8b38805d26f0", "Complex": "a8da6dc1099b8b38805d26f0",
    "Complex_Backbone": "a8da6dc1099b8b38805d26f0", "Complex_CA": "a8da6dc1099b8b38805d26f0",
    "System": "0620dc9876a92272bf8358b2",
}


# ═══════════════════════════════════════════════════════════════════════════════
# Cofactor — the reproducible study-path gap
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
def test_resolved_cofactor_becomes_group_and_diagnostics_run(tmp_path):
    pdb, xtc = _write_system(tmp_path / "sys", "HEM")
    det = detect_components(structure_path=pdb)
    cof = [c for c in det.components if c.component_type == ComponentType.COFACTOR]
    assert cof and cof[0].classification_state == ClassificationState.RESOLVED
    rec = _record(pdb, xtc, det.components)
    idx = _study_index(rec, tmp_path / "out")
    groups = _groups(idx.path)

    assert groups["Cofactor"] == (31, 32, 33, 34)                # exactly the HEM atoms
    assert {k: _ordered(v) for k, v in groups.items() if k != "Cofactor"} == \
        _GOLDEN_COFACTOR_FIXTURE                                  # existing groups unchanged
    assert list(groups)[-2:] == ["Cofactor", "System"]           # appended before System

    rep = diagnose_system(rec)
    run = next(r for r in rep.detectors if r.detector_id == "partner_receptor_separation")
    assert run.status == DetectorRunStatus.RAN
    codes = {d.code for d in rep.diagnostics}
    assert "cofactor_receptor_periodic_separation" in codes
    assert "component_periodic_wrap" in codes


# ═══════════════════════════════════════════════════════════════════════════════
# Ligand
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
def test_resolved_ligand_group_exact_membership_and_detector_runs(tmp_path):
    pdb, xtc = _write_system(tmp_path / "sys", "LIG")
    det = detect_components(structure_path=pdb)
    # the detector never emits LIGAND (see the documented gap below); a resolved
    # ligand component is supplied as a SystemRecord would carry it
    rec = _record(pdb, xtc, det.components + [_ligand()])
    idx = _study_index(rec, tmp_path / "out")
    groups = _groups(idx.path)
    assert groups["Ligand"] == (31, 32, 33, 34)
    assert len(set(groups["Ligand"])) == len(groups["Ligand"])

    rep = diagnose_system(rec)
    run = next(r for r in rep.detectors if r.detector_id == "partner_receptor_separation")
    assert run.status == DetectorRunStatus.RAN
    (sep,) = [d for d in rep.diagnostics if d.code == "ligand_receptor_periodic_separation"]
    assert sep.measured["events"][0]["frames"] == [1, 2]
    assert sep.provenance["partner_group"] == "Ligand"
    assert sep.provenance["partner_group_evidence"]["n_atoms"] == 4


@requires_gmx
def test_no_ligand_no_group(tmp_path):
    pdb, _ = _write_system(tmp_path / "sys", "LIG", trajectory=False)
    det = detect_components(structure_path=pdb)
    idx = build_semantic_index(structure_path=pdb, components=det.components,
                               out_ndx=tmp_path / "i.ndx")
    assert "Ligand" not in _groups(idx.path)


@requires_gmx
@pytest.mark.parametrize("state", [ClassificationState.AMBIGUOUS, ClassificationState.UNKNOWN,
                                   ClassificationState.REVIEW_REQUIRED])
def test_unresolved_ligand_not_guessed(tmp_path, state):
    pdb, _ = _write_system(tmp_path / "sys", "LIG", trajectory=False)
    det = detect_components(structure_path=pdb)
    idx = build_semantic_index(structure_path=pdb, components=det.components + [_ligand(state)],
                               out_ndx=tmp_path / "i.ndx")
    assert "Ligand" not in _groups(idx.path)
    assert any("no 'Ligand' group" in w for w in idx.warnings)


@requires_gmx
def test_multiple_resolved_ligands_not_merged(tmp_path):
    pdb, _ = _write_system(tmp_path / "sys", "LIG", trajectory=False)
    det = detect_components(structure_path=pdb)
    idx = build_semantic_index(
        structure_path=pdb,
        components=det.components + [_ligand(label="ligand 1"), _ligand(label="ligand 2")],
        out_ndx=tmp_path / "i.ndx")
    assert not any(n.startswith("Ligand") for n in _groups(idx.path))
    assert any("not merged" in w for w in idx.warnings)


# ═══════════════════════════════════════════════════════════════════════════════
# Nucleic acid
# ═══════════════════════════════════════════════════════════════════════════════

_DNA = [("N", 40 + i, rn, an, 2.0 + 0.1 * i, 2.0, 2.0 + 0.1 * k, "ATOM")
        for i, rn in enumerate(("DA", "DT", "DG")) for k, an in enumerate(("P", "C4'", "N1"))]


@requires_gmx
def test_resolved_nucleic_acid_group(tmp_path):
    pdb, _ = _write_system(tmp_path / "sys", "LIG", extra_chain=_DNA, trajectory=False)
    det = detect_components(structure_path=pdb)
    # since Phase 4.6 the detector resolves the DNA chain itself
    (na,) = [c for c in det.components if c.component_type == ComponentType.NUCLEIC_ACID]
    assert na.chain_ids == ["N"] and na.classification_state == ClassificationState.RESOLVED
    idx = build_semantic_index(structure_path=pdb, components=det.components,
                               out_ndx=tmp_path / "i.ndx")
    groups = _groups(idx.path)
    assert groups["NucleicAcid"] == tuple(range(36, 45))       # the 9 DNA atoms, in order
    assert "NucleicAcid_Backbone" not in groups                  # no protein conventions


# ═══════════════════════════════════════════════════════════════════════════════
# Identity (Phase 3 views, Phase 4 detectors)
# ═══════════════════════════════════════════════════════════════════════════════

def _view_key(tmp_path, ndx_text, req):
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.trajectory.preprocessor import (
        _identity_evidence, _plan_steps, view_identity,
    )
    ndx = tmp_path / "semantic.ndx"
    ndx.write_text(ndx_text)
    traj, tpr = tmp_path / "md.xtc", tmp_path / "md.tpr"
    traj.write_bytes(b"t")
    tpr.write_bytes(b"p")
    names = [l.strip("[] \n") for l in ndx_text.splitlines() if l.startswith("[")]
    sidx = SemanticIndex(path=str(ndx), groups=[SemanticIndexGroup(n, 0) for n in names])
    steps, err = _plan_steps(req, str(traj), str(tpr), None, sidx, tmp_path / "o")
    assert err is None
    ev, groups = _identity_evidence(req, req.view_kind(), steps, [fingerprint_file(traj)],
                                    {str(tpr): fingerprint_file(tpr)}, str(ndx), "v")
    return view_identity(ev), groups


_BASE = "[ Receptor ]\n1 2 3\n[ Receptor_Backbone ]\n1 2\n[ System ]\n1 2 3 4 5\n"


def test_unused_ligand_group_keeps_protein_view_identity(tmp_path):
    for req in (reqs.whole_only("w"), reqs.fit_to("Receptor_Backbone", "f")):
        a, _ = _view_key(tmp_path, _BASE, req)
        b, _ = _view_key(tmp_path, _BASE + "[ Ligand ]\n4\n", req)
        assert a == b


def test_view_using_ligand_hashes_its_atoms(tmp_path):
    from analysis.campaign.results import atom_set_hash
    req = reqs.centered_on("Ligand", "c")
    k1, groups = _view_key(tmp_path, _BASE + "[ Ligand ]\n4\n", req)
    lig = [g for g in groups if g["name"] == "Ligand"]
    assert lig and lig[0]["atoms_sha256"] == atom_set_hash([4])
    k2, _ = _view_key(tmp_path, _BASE + "[ Ligand ]\n4 5\n", req)
    assert k1 != k2


def test_detector_identity_includes_ligand_atoms(tmp_path):
    ndx = tmp_path / "semantic.ndx"
    traj = tmp_path / "md.xtc"
    traj.write_bytes(b"t")
    spec = dreg.get("partner_receptor_separation")

    def ident(ligand_atoms):
        ndx.write_text(_BASE + f"[ Ligand ]\n{ligand_atoms}\n")
        ctx = DiagnosticContext(trajectory_path=str(traj), index_path=str(ndx),
                                groups=groups_from_semantic_index(str(ndx)))
        assert ctx.groups["ligand"] == "Ligand"
        return detector_identity(spec, ctx)

    assert ident("4") == ident("4") != ident("4 5")


# ═══════════════════════════════════════════════════════════════════════════════
# Former upstream detector gaps — fixed in Phase 4.6, locked as regressions
# ═══════════════════════════════════════════════════════════════════════════════

def test_detector_resolves_small_molecule_ligand(tmp_path):
    """Was a strict xfail (the detector never emitted LIGAND)."""
    d = tmp_path / "lig"
    d.mkdir()
    lines = [f"ATOM  {i + 1:>5} CA   {AA[i % 10]:<3} A{i + 1:>4}    {10.0:8.3f}{10.0:8.3f}{i:8.3f}"
             f"  1.00  0.00           C" for i in range(30)]
    lines += [f"HETATM{31 + k:>5} C{k + 1:<3} LIG B   1    {20.0:8.3f}{10.0 + k:8.3f}{0.0:8.3f}"
              f"  1.00  0.00           C" for k in range(10)]
    pdb = d / "pl.pdb"
    pdb.write_text("\n".join(lines + ["END"]) + "\n")
    det = detect_components(structure_path=pdb)
    (lig,) = [c for c in det.components if c.component_type == ComponentType.LIGAND]
    assert lig.classification_state == ClassificationState.RESOLVED
    assert lig.atom_ids == list(range(31, 41))


def test_tiny_unknown_molecule_stays_unresolved(tmp_path):
    """The 4-heavy-atom fixture molecule is below the conservative ligand threshold."""
    pdb, _ = _write_system(tmp_path / "sys", "LIG", trajectory=False)
    det = detect_components(structure_path=pdb)
    assert not any(c.component_type == ComponentType.LIGAND for c in det.components)


def test_detector_resolves_single_nucleic_acid_chain(tmp_path):
    """Was a strict xfail (dominant_class() 'nucleotide' vs 'nucleic_acid' mismatch)."""
    d = tmp_path / "dna"
    d.mkdir()
    pdb = d / "dna.pdb"
    lines = [f"ATOM  {i + 1:>5} P    {rn:<3} N{i + 1:>4}    {10.0:8.3f}{10.0:8.3f}{i:8.3f}"
             f"  1.00  0.00           P" for i, rn in enumerate(["DA", "DT", "DG", "DC"] * 5)]
    pdb.write_text("\n".join(lines + ["END"]) + "\n")
    det = detect_components(structure_path=pdb)
    assert any(c.component_type == ComponentType.NUCLEIC_ACID for c in det.components)
    assert not any(c.component_type == ComponentType.RECEPTOR for c in det.components)
