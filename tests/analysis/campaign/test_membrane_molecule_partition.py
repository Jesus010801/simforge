"""Phase 13.6 — membrane molecule-partition audit.

A protein-membrane topology can represent lipids either as their own
moleculetypes (``DPP  N``) or — legacy SimForge builds that ran pdb2gmx on the
mixed protein+bilayer coordinates — merged into the protein moleculetype with
no bonded term between protein and lipid atoms.  The partition reader must
report exactly what the .tpr says, and the finite-assembly refusal for the
membrane must come from the *component semantics* (the Membrane group), never
from moltype names or molecule counts.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from analysis.campaign.models import SemanticIndex, SemanticIndexGroup, SystemRecord
from analysis.campaign.structure.molecules import molecule_partition, partition_from_dump
from analysis.campaign.trajectory.assembly import decide_assembly

requires_gmx = pytest.mark.skipif(shutil.which("gmx") is None, reason="needs gmx")

# protein (1) + two lipid species with repeated moltypes + water + two ion types
DUMP_SEPARATE = """   molblock (0):
      moltype              = 0 "Protein"
      #molecules                     = 1
   molblock (1):
      moltype              = 1 "DPP"
      #molecules                     = 3
   molblock (2):
      moltype              = 2 "POP"
      #molecules                     = 2
   molblock (3):
      moltype              = 3 "SOL"
      #molecules                     = 4
   molblock (4):
      moltype              = 4 "NA"
      #molecules                     = 2
   molblock (5):
      moltype              = 5 "CL"
      #molecules                     = 2
   moltype (0):
      name="Protein"
      atoms:
         atom (4):
   moltype (1):
      name="DPP"
      atoms:
         atom (2):
   moltype (2):
      name="POP"
      atoms:
         atom (3):
   moltype (3):
      name="SOL"
      atoms:
         atom (3):
   moltype (4):
      name="NA"
      atoms:
         atom (1):
   moltype (5):
      name="CL"
      atoms:
         atom (1):
"""

# the same atoms, lipids merged into the protein moleculetype (GLP-1R legacy build)
DUMP_MERGED = """   molblock (0):
      moltype              = 0 "Protein"
      #molecules                     = 1
   molblock (1):
      moltype              = 1 "SOL"
      #molecules                     = 4
   molblock (2):
      moltype              = 2 "NA"
      #molecules                     = 2
   molblock (3):
      moltype              = 3 "CL"
      #molecules                     = 2
   moltype (0):
      name="Protein"
      atoms:
         atom (16):
   moltype (1):
      name="SOL"
      atoms:
         atom (3):
   moltype (2):
      name="NA"
      atoms:
         atom (1):
   moltype (3):
      name="CL"
      atoms:
         atom (1):
"""

RECEPTOR = list(range(1, 5))
MEMBRANE = list(range(5, 17))
WATER = list(range(17, 29))
IONS = list(range(29, 33))


def test_partition_separate_lipids_repeated_moltypes():
    p = partition_from_dump(DUMP_SEPARATE.splitlines())
    assert (p.n_molecules, p.n_atoms) == (14, 32)
    assert p.moltypes == ["Protein"] + ["DPP"] * 3 + ["POP"] * 2 + ["SOL"] * 4 + ["NA"] * 2 + ["CL"] * 2
    assert p.starts[:7] == [1, 5, 7, 9, 11, 14, 17]
    rec, mem, wat, ion = (p.select(x) for x in (RECEPTOR, MEMBRANE, WATER, IONS))
    assert (rec.n_molecules, rec.complete_molecules) == (1, True)
    assert (mem.n_molecules, mem.complete_molecules, mem.moltypes) == (5, True, ["DPP", "POP"])
    assert mem.atoms_per_molecule == [2, 2, 2, 3, 3]
    assert (wat.n_molecules, ion.n_molecules, ion.moltypes) == (4, 4, ["CL", "NA"])
    assert not set(p.closure(RECEPTOR)) & set(MEMBRANE)


def test_partition_merged_lipids_is_reported_as_one_molecule():
    """The reader does not 'repair' a merged topology: it reports it."""
    p = partition_from_dump(DUMP_MERGED.splitlines())
    assert (p.n_molecules, p.n_atoms) == (9, 32)
    rec, mem = p.select(RECEPTOR), p.select(MEMBRANE)
    assert rec.molecules == mem.molecules == [0]
    assert not rec.complete_molecules and not mem.complete_molecules
    assert p.closure(RECEPTOR) == RECEPTOR + MEMBRANE          # molecule 0 = Receptor ∪ Membrane
    assert p.select(WATER).n_molecules == 4 and p.select(IONS).n_molecules == 4


# ═══════════════════════════════════════════════════════════════════════════════
# Real grompp: both representations of the same coordinates
# ═══════════════════════════════════════════════════════════════════════════════

_HEAD = """[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
CX 12.011 0.0 A 0.34 0.36
"""


def _moltype(name, residues, bonds):
    """residues: [(resnr, resname, n_atoms)]; bonds: [(i, j)] within the moltype."""
    rows, k = [], 0
    for r, rn, n in residues:
        for _ in range(n):
            k += 1
            rows.append(f"{k} CX {r} {rn} A{k} {k} 0 12.011")
    out = [f"[ moleculetype ]\n{name} 3\n[ atoms ]", *rows]
    if bonds:
        out += ["[ bonds ]", *(f"{i} {j} 1 0.15 1000" for i, j in bonds)]
    return "\n".join(out) + "\n"


def _chain(start, n):
    return [(start + i, start + i + 1) for i in range(n - 1)]


def _topologies() -> dict[str, str]:
    prot = [(1, "PRO", 4)]
    dpp = [(2 + i, "DPP", 2) for i in range(3)]
    pop = [(5 + i, "POP", 3) for i in range(2)]
    tail = (_moltype("SOL", [(1, "SOL", 3)], _chain(1, 3)) + _moltype("NA", [(1, "NA", 1)], [])
            + _moltype("CL", [(1, "CL", 1)], []))
    separate = (_HEAD + _moltype("Protein", prot, _chain(1, 4))
                + _moltype("DPP", [(1, "DPP", 2)], _chain(1, 2))
                + _moltype("POP", [(1, "POP", 3)], _chain(1, 3)) + tail
                + "[ system ]\nsep\n[ molecules ]\nProtein 1\nDPP 3\nPOP 2\nSOL 4\nNA 2\nCL 2\n")
    # merged: lipid fragments inside the protein moltype, no protein-lipid bond
    bonds = _chain(1, 4) + _chain(5, 2) + _chain(7, 2) + _chain(9, 2) + _chain(11, 3) + _chain(14, 3)
    merged = (_HEAD + _moltype("Protein", prot + dpp + pop, bonds) + tail
              + "[ system ]\nmerged\n[ molecules ]\nProtein 1\nSOL 4\nNA 2\nCL 2\n")
    return {"separate": separate, "merged": merged}


def _gro(path: Path):
    names = ["PRO"] * 4 + ["DPP"] * 6 + ["POP"] * 6 + ["SOL"] * 12 + ["NA"] * 2 + ["CL"] * 2
    lines = ["toy", f"{len(names):5d}"]
    for i, rn in enumerate(names, 1):
        x, y, z = 0.3 * (i % 8) + 0.5, 0.3 * ((i // 8) % 8) + 0.5, 0.3 * (i // 64) + 0.5
        lines.append(f"{i:5d}{rn:<5s}{'A' + str(i):>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append("   5.00000   5.00000   5.00000")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def tprs(tmp_path_factory):
    if shutil.which("gmx") is None:
        pytest.skip("needs gmx")
    d = tmp_path_factory.mktemp("membrane_partition")
    _gro(d / "conf.gro")
    (d / "md.mdp").write_text("integrator=md\nnsteps=0\ncutoff-scheme=Verlet\nrvdw=1.0\n"
                              "rcoulomb=1.0\ncoulombtype=cut-off\n")
    (d / "index.ndx").write_text(
        "".join(f"[ {n} ]\n{' '.join(map(str, a))}\n" for n, a in
                (("Receptor", RECEPTOR), ("Membrane", MEMBRANE), ("Water", WATER),
                 ("Ions", IONS), ("Complex", RECEPTOR + MEMBRANE))))
    (d / "unlabelled.ndx").write_text(                      # same atoms, no Membrane semantics
        f"[ Receptor ]\n{' '.join(map(str, RECEPTOR))}\n"
        f"[ Lipids ]\n{' '.join(map(str, MEMBRANE))}\n")
    out = {}
    for name, top in _topologies().items():
        (d / f"{name}.top").write_text(top)
        r = subprocess.run(["gmx", "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", f"{name}.top",
                            "-o", f"{name}.tpr", "-maxwarn", "10"], cwd=d, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-800:]
        out[name] = d / f"{name}.tpr"
    out["dir"] = d
    return out


def _rec(tpr: Path, ndx: Path) -> tuple[SystemRecord, SemanticIndex]:
    rec = SystemRecord("toy", "c", "r", topology_path=str(tpr))
    names = [l.strip("[] \n") for l in ndx.read_text().splitlines() if l.startswith("[")]
    return rec, SemanticIndex(path=str(ndx), groups=[SemanticIndexGroup(n, 1) for n in names])


@requires_gmx
def test_real_tpr_partition_matches_topology(tprs):
    sep = molecule_partition(str(tprs["separate"]))
    assert sep.moltypes == partition_from_dump(DUMP_SEPARATE.splitlines()).moltypes
    assert sep.starts == partition_from_dump(DUMP_SEPARATE.splitlines()).starts
    mer = molecule_partition(str(tprs["merged"]))
    assert mer.starts == partition_from_dump(DUMP_MERGED.splitlines()).starts
    assert mer.moltypes == ["Protein"] + ["SOL"] * 4 + ["NA"] * 2 + ["CL"] * 2


@requires_gmx
@pytest.mark.parametrize("rep", ["separate", "merged"])
def test_receptor_alone_needs_no_assembly_in_either_representation(tprs, rep):
    """A merged topology must not drag the membrane into a receptor cluster."""
    rec, si = _rec(tprs[rep], tprs["dir"] / "index.ndx")
    d = decide_assembly(rec, si, ("Receptor",))
    assert not d.blocked and not d.needed
    assert d.evidence["molecules"]["n_molecules"] == 1
    assert d.evidence["molecules"]["complete_molecules"] is (rep == "separate")


@requires_gmx
@pytest.mark.parametrize("rep", ["separate", "merged"])
@pytest.mark.parametrize("groups", [("Membrane",), ("Complex",), ("Receptor", "Membrane")])
def test_membrane_refusal_comes_from_component_semantics(tprs, rep, groups):
    rec, si = _rec(tprs[rep], tprs["dir"] / "index.ndx")
    d = decide_assembly(rec, si, groups)
    assert not d.needed and "periodically extended membrane component" in d.blocked
    assert d.geometry == "periodic_extended"


@requires_gmx
def test_without_membrane_semantics_the_lipids_are_just_molecules(tprs):
    """Counterfactual: the same lipid atoms without a Membrane label are a finite
    multi-molecule set (separate) or part of one molecule (merged) — the refusal
    above is therefore semantic, not topological."""
    rec, si = _rec(tprs["separate"], tprs["dir"] / "unlabelled.ndx")
    d = decide_assembly(rec, si, ("Lipids",))
    assert d.needed and not d.blocked and d.evidence["molecules"]["n_molecules"] == 5
    rec, si = _rec(tprs["merged"], tprs["dir"] / "unlabelled.ndx")
    d = decide_assembly(rec, si, ("Lipids",))
    assert not d.needed and not d.blocked and d.evidence["molecules"]["n_molecules"] == 1
