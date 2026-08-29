"""Shared fixtures for the ``analysis/campaign`` test suite.

All structure/topology files here are tiny, hand-crafted and synthetic.  The
only real MD data used anywhere in the suite is the trajectory under
``simforge_runs/protein-membrane/...`` which a handful of gmx-guarded
integration tests symlink into a temporary tree.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]

REAL_RUN_DIR = (
    REPO_ROOT / "simforge_runs" / "protein-membrane" / "runs"
    / "2026-06-25_14-38-57" / "steps" / "11_production_md"
)
REAL_XTC = REAL_RUN_DIR / "md.xtc"
REAL_TPR = REAL_RUN_DIR / "md.tpr"
REAL_GRO = REAL_RUN_DIR / "md.gro"

HAVE_GMX = shutil.which("gmx") is not None
HAVE_REAL_TRAJ = REAL_XTC.is_file() and REAL_TPR.is_file() and REAL_GRO.is_file()

requires_gmx = pytest.mark.skipif(not HAVE_GMX, reason="needs gmx")
requires_real_traj = pytest.mark.skipif(
    not HAVE_REAL_TRAJ, reason="needs the bundled protein-membrane trajectory"
)

_AA_CYCLE = [
    "ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE",
    "PHE", "TYR", "TRP", "ARG", "ASN", "GLN", "HIS", "PRO", "MET", "CYS",
]


# ── low-level writers ────────────────────────────────────────────────────────

def pdb_atom(serial: int, name: str, resname: str, chain: str, resseq: int,
             x: float, y: float, z: float) -> str:
    """One ATOM record with columns the campaign parsers rely on."""
    return (
        f"ATOM  {serial:>5} {name:<4} {resname:<3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
    )


def gro_atom(resnum: int, resname: str, atomname: str, atomnum: int,
             x: float, y: float, z: float) -> str:
    return (
        f"{resnum % 100000:>5}{resname:<5}{atomname:>5}{atomnum % 100000:>5}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def write_pdb(path: Path, chains: dict[str, int], *, start_z: float = 1.0,
              dz: float = 0.5) -> Path:
    """Write a CA-only PDB. ``chains`` maps chain id -> residue count."""
    lines: list[str] = ["REMARK  synthetic test structure"]
    serial = 1
    z = start_z
    for cid, nres in chains.items():
        for i in range(1, nres + 1):
            rn = _AA_CYCLE[(i - 1) % len(_AA_CYCLE)]
            lines.append(pdb_atom(serial, "CA", rn, cid, i, 10.0, 10.0, z))
            serial += 1
            z += dz
        lines.append(f"TER   {serial:>5}      {rn} {cid}{nres:>4}")
        serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_gro(path: Path, *, n_protein: int = 50, n_sol: int = 2000,
              n_na: int = 5, n_cl: int = 5, box: float = 8.0) -> Path:
    atoms: list[str] = []
    n = 0
    z = 1.0
    for i in range(1, n_protein + 1):
        rn = _AA_CYCLE[(i - 1) % len(_AA_CYCLE)]
        n += 1
        atoms.append(gro_atom(i, rn, "CA", n, 4.0, 4.0, z))
        z += 0.05
    resnum = n_protein
    for i in range(n_sol):
        resnum += 1
        n += 1
        atoms.append(gro_atom(resnum, "SOL", "OW", n, 1.0 + (i % 5) * 0.3, 1.0, 1.0))
    for i in range(n_na):
        resnum += 1
        n += 1
        atoms.append(gro_atom(resnum, "NA", "NA", n, 2.0, 2.0, 2.0))
    for i in range(n_cl):
        resnum += 1
        n += 1
        atoms.append(gro_atom(resnum, "CL", "CL", n, 3.0, 3.0, 3.0))
    text = "synthetic solvated protein\n"
    text += f"{n}\n"
    text += "\n".join(atoms) + "\n"
    text += f"{box:10.5f}{box:10.5f}{box:10.5f}\n"
    path.write_text(text)
    return path


def write_top(path: Path, molecules: list[tuple[str, int]] | None = None,
              *, with_moleculetype: bool = True) -> Path:
    molecules = molecules or [
        ("Protein_chain_A", 1), ("Peptide", 1), ("SOL", 2000), ("NA", 10),
    ]
    lines: list[str] = ['#include "amber99sb.ff/forcefield.itp"', ""]
    if with_moleculetype:
        lines += [
            "[ moleculetype ]",
            "; name   nrexcl",
            "Protein_chain_A   3",
            "",
            "[ atoms ]",
            "  1  N  1  ALA  N  1  -0.4  14.01",
            "",
        ]
    lines += ["[ system ]", "Synthetic", "", "[ molecules ]",
              "; Compound   #mols"]
    for name, count in molecules:
        lines.append(f"{name}   {count}")
    path.write_text("\n".join(lines) + "\n")
    return path


# ── structure fixtures ──────────────────────────────────────────────────────

@pytest.fixture
def minimal_pdb_receptor_peptide(tmp_path) -> Path:
    return write_pdb(tmp_path / "receptor_peptide.pdb", {"A": 120, "B": 25})


@pytest.fixture
def minimal_pdb_single_chain(tmp_path) -> Path:
    return write_pdb(tmp_path / "single_chain.pdb", {"A": 150})


@pytest.fixture
def minimal_pdb_homodimer(tmp_path) -> Path:
    return write_pdb(tmp_path / "homodimer.pdb", {"A": 100, "B": 98})


@pytest.fixture
def minimal_pdb_three_chains(tmp_path) -> Path:
    """A long chain + two short chains (used with a monkeypatched embedding)."""
    return write_pdb(tmp_path / "three_chains.pdb", {"A": 140, "B": 22, "C": 24})


@pytest.fixture
def minimal_pdb_unequal(tmp_path) -> Path:
    return write_pdb(tmp_path / "unequal.pdb", {"A": 200, "B": 60})


@pytest.fixture
def minimal_gro_solvated_protein(tmp_path) -> Path:
    return write_gro(tmp_path / "solvated.gro")


@pytest.fixture
def minimal_top(tmp_path) -> Path:
    return write_top(tmp_path / "topol.top")


# ── real-trajectory fixtures ────────────────────────────────────────────────

@pytest.fixture
def real_membrane_gro(tmp_path) -> Path:
    if not HAVE_REAL_TRAJ:
        pytest.skip("needs the bundled protein-membrane trajectory")
    link = tmp_path / "membrane_ref.gro"
    link.symlink_to(REAL_GRO)
    return link


# ── study-tree factory ──────────────────────────────────────────────────────

@pytest.fixture
def synthetic_study_tree(tmp_path):
    """Build a directory tree of synthetic MD "systems".

    Each spec is ``(partner, water, replicate, traj_name, top_name, struct_name)``.
    Trajectory files are tiny stubs — discovery only looks at paths/extensions.
    """
    default_specs = [
        ("GLP1", "TIP3P", "rep1", "final.xtc", "topol.top", "conf.pdb"),
        ("GLP1", "TIP3P", "rep2", "production.xtc", "system.top", "md.pdb"),
        ("GLP1", "SPCE", "rep1", "md_100ns.xtc", "prod.top", "confout.pdb"),
        ("Semaglutide", "TIP3P", "rep1", "run3.trr", "topol.top", "conf.pdb"),
    ]

    def _make(specs=None, *, root_name="GLP1R_study", structure="single",
              traj_bytes=b"FAKEXTCDATA", make_structure=True, make_topology=True):
        specs = default_specs if specs is None else specs
        root = tmp_path / root_name
        for partner, water, rep, traj, top, struct in specs:
            d = root / partner / water / rep
            d.mkdir(parents=True, exist_ok=True)
            (d / traj).write_bytes(traj_bytes)
            if make_topology and top:
                write_top(d / top)
            if make_structure and struct:
                if structure == "single":
                    write_pdb(d / struct, {"A": 120})
                elif structure == "receptor_peptide":
                    write_pdb(d / struct, {"A": 120, "B": 25})
                elif structure == "gro":
                    write_gro(d / struct)
        return root

    return _make


@pytest.fixture
def study_tree(synthetic_study_tree):
    return synthetic_study_tree()


@pytest.fixture(autouse=True)
def _reset_observable_registry():
    """Keep the global observable registry pristine between tests.

    NOTE: we snapshot/restore ``registry._REGISTRY`` directly rather than call
    ``registry._reset_for_tests()`` — the latter is broken (it clears the dict
    and then relies on a one-time import side-effect to repopulate it, so the
    registry ends up permanently empty).  See the suite's bug report.
    """
    from analysis.campaign.observables import registry
    snapshot = dict(registry._REGISTRY)
    yield
    registry._REGISTRY.clear()
    registry._REGISTRY.update(snapshot)


@pytest.fixture(autouse=True)
def _clear_view_cache():
    from analysis.campaign.orchestration import study_analyzer
    study_analyzer.clear_view_cache()
    yield
    study_analyzer.clear_view_cache()


def file_stat(path: Path) -> tuple[int, int]:
    st = Path(path).stat()
    return st.st_size, st.st_mtime_ns
