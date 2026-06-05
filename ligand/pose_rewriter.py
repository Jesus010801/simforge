"""
LigandPoseRewriter — transfers pose coordinates from a PDB into the
normalized LigParGen GRO topology.

Inputs
------
pose_pdb        : PDB with the ligand in its docked/bound pose (Å)
reference_gro   : Normalized L01.gro — canonical atom order, names, residue
reference_itp   : Normalized L01.itp — used for element inference via mass

Output
------
{output_dir}/{residue_name}_pose.gro  (coordinates from pose, topology from GRO)

Does not assemble systems, edit topol.top, or run grompp.
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np

from core.ligand_workflow_models import LigandPoseRewriteResult
from utils.geometry_alignment import align
from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro
from utils.itp_parser import parse_itp
from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

logger = logging.getLogger(__name__)

# Å → nm
_ANG_TO_NM = 0.1

# Reference masses (Da) and tolerance for element lookup
_MASS_TABLE: list[tuple[float, str]] = [
    (1.008,    "H"),
    (12.011,   "C"),
    (14.007,   "N"),
    (15.999,   "O"),
    (18.998,   "F"),
    (30.974,   "P"),
    (32.060,   "S"),
    (35.453,   "Cl"),
    (79.904,   "Br"),
    (126.904,  "I"),
]
_MASS_TOL = 0.5  # Da


def _element_from_mass(mass: float) -> str:
    for ref, symbol in _MASS_TABLE:
        if abs(mass - ref) < _MASS_TOL:
            return symbol
    return "X"


def _match_atoms(
    ref_elements: list[str],
    ref_coords: np.ndarray,   # (N, 3) nm
    pose_elements: list[str],
    pose_coords: np.ndarray,  # (N, 3) nm
) -> list[int]:
    """Greedy nearest-neighbour matching per element (centered coordinates).

    Returns mapping[i] = pose index that corresponds to ref index i.
    Raises ValueError on composition mismatch or impossible assignment.
    """
    from collections import Counter
    if Counter(ref_elements) != Counter(pose_elements):
        ref_cnt = dict(Counter(ref_elements))
        pose_cnt = dict(Counter(pose_elements))
        raise ValueError(
            f"Element composition mismatch — "
            f"reference: {ref_cnt}, pose: {pose_cnt}"
        )

    ref_c = ref_coords.mean(axis=0)
    pose_c = pose_coords.mean(axis=0)
    ref_cen = ref_coords - ref_c
    pose_cen = pose_coords - pose_c

    # Pool available pose indices per element
    pose_pool: dict[str, list[int]] = defaultdict(list)
    for i, el in enumerate(pose_elements):
        pose_pool[el].append(i)

    mapping: list[int] = [-1] * len(ref_elements)
    used: set[int] = set()

    for ref_i, el in enumerate(ref_elements):
        candidates = [j for j in pose_pool[el] if j not in used]
        if not candidates:
            raise ValueError(
                f"No available {el!r} atom in pose for reference atom {ref_i}"
            )
        dists = [
            float(np.linalg.norm(ref_cen[ref_i] - pose_cen[j]))
            for j in candidates
        ]
        best = candidates[int(np.argmin(dists))]
        mapping[ref_i] = best
        used.add(best)

    return mapping


def _heavy_atom_rmsd(
    gro_a: GroFile,
    gro_b: GroFile,
    elements: list[str],
) -> float:
    """Kabsch-aligned RMSD on heavy atoms (nm) between two GroFiles."""
    heavy = [i for i, el in enumerate(elements) if el != "H"]
    coords_a = np.array([[gro_a.atoms[i].x, gro_a.atoms[i].y, gro_a.atoms[i].z] for i in heavy])
    coords_b = np.array([[gro_b.atoms[i].x, gro_b.atoms[i].y, gro_b.atoms[i].z] for i in heavy])
    _, _, rmsd = align(coords_a, coords_b)
    return rmsd


class LigandPoseRewriter:
    """Transfers pose coordinates from a PDB into the normalized GRO topology."""

    def rewrite(
        self,
        pose_pdb: Path | str,
        reference_gro: Path | str,
        reference_itp: Path | str,
        output_dir: Path | str,
    ) -> LigandPoseRewriteResult:
        pose_pdb = Path(pose_pdb)
        reference_gro = Path(reference_gro)
        reference_itp = Path(reference_itp)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # ── Parse inputs ───────────────────────────────────────────────────────
        if not pose_pdb.exists():
            return _fail(f"Pose PDB not found: {pose_pdb}")
        if not reference_gro.exists():
            return _fail(f"Reference GRO not found: {reference_gro}")
        if not reference_itp.exists():
            return _fail(f"Reference ITP not found: {reference_itp}")

        try:
            ref_gro = parse_gro(reference_gro)
        except Exception as exc:
            return _fail(f"Cannot parse reference GRO: {exc}")

        try:
            ref_itp = parse_itp(reference_itp)
        except Exception as exc:
            return _fail(f"Cannot parse reference ITP: {exc}")

        pdb_raw = parse_pdb_ligand(pose_pdb)
        if pdb_raw.get("error"):
            return _fail(f"Cannot parse pose PDB: {pdb_raw['error']}")

        # ── Consistency checks ─────────────────────────────────────────────────
        n_ref = ref_gro.atom_count
        n_pose = pdb_raw["n_atoms"]
        if n_ref != n_pose:
            return _fail(
                f"Atom count mismatch: reference GRO has {n_ref} atoms, "
                f"pose PDB has {n_pose} atoms"
            )

        if len(ref_itp.atoms) != n_ref:
            return _fail(
                f"ITP atom count ({len(ref_itp.atoms)}) does not match "
                f"GRO atom count ({n_ref})"
            )

        # ── Element inference ──────────────────────────────────────────────────
        ref_elements = [_element_from_mass(a.mass) for a in ref_itp.atoms]
        pose_elements = [a["element"] for a in pdb_raw["atoms"]]

        # ── Coordinate arrays (both in nm) ─────────────────────────────────────
        ref_coords = np.array(
            [[a.x, a.y, a.z] for a in ref_gro.atoms], dtype=float
        )
        pose_coords = np.array(
            [[a["x"] * _ANG_TO_NM, a["y"] * _ANG_TO_NM, a["z"] * _ANG_TO_NM]
             for a in pdb_raw["atoms"]],
            dtype=float,
        )

        # ── Atom matching ──────────────────────────────────────────────────────
        try:
            mapping = _match_atoms(ref_elements, ref_coords, pose_elements, pose_coords)
        except ValueError as exc:
            return _fail(str(exc))

        # ── Build output GRO ───────────────────────────────────────────────────
        pose_coords_mapped = pose_coords[mapping]
        residue_name = ref_gro.atoms[0].residue_name

        out_atoms: list[GroAtom] = []
        for i, ref_atom in enumerate(ref_gro.atoms):
            x, y, z = pose_coords_mapped[i]
            out_atoms.append(
                GroAtom(
                    residue_number=ref_atom.residue_number,
                    residue_name=ref_atom.residue_name,
                    atom_name=ref_atom.atom_name,
                    atom_number=ref_atom.atom_number,
                    x=float(x),
                    y=float(y),
                    z=float(z),
                )
            )

        out_gro = GroFile(
            title=f"{residue_name} pose",
            atoms=out_atoms,
            box=ref_gro.box,
        )

        output_path = output_dir / f"{residue_name}_pose.gro"
        try:
            write_gro(out_gro, output_path)
        except Exception as exc:
            return _fail(f"Cannot write output GRO: {exc}")

        # ── Heavy-atom RMSD ────────────────────────────────────────────────────
        try:
            rmsd = _heavy_atom_rmsd(out_gro, ref_gro, ref_elements)
        except Exception as exc:
            logger.warning("RMSD computation failed: %s", exc)
            rmsd = None

        logger.info(
            "Pose rewrite: %d atoms → %s (heavy-atom RMSD %.4f nm)",
            n_ref, output_path.name, rmsd or 0.0,
        )

        return LigandPoseRewriteResult(
            success=True,
            output_path=output_path,
            ligand_residue_name=residue_name,
            atoms_written=n_ref,
            rmsd_from_reference=rmsd,
        )


def _fail(message: str) -> LigandPoseRewriteResult:
    logger.error("LigandPoseRewriter: %s", message)
    return LigandPoseRewriteResult(success=False, error=message)
