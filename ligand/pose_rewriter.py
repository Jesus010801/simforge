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
from utils.geometry_alignment import align, kabsch_rotation
from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro
from utils.itp_parser import parse_itp
from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

logger = logging.getLogger(__name__)

# Å → nm
_ANG_TO_NM = 0.1

# Minimum number of local points (center + neighbours) needed for a
# well-conditioned per-atom Kabsch rotation during hydrogen reconstruction.
_MIN_LOCAL_FRAME_POINTS = 3

# Sanity limit (nm) on centroid-aligned heavy-atom deviation for the
# geometric nearest-neighbour correspondence fallback. Above this, the
# match is treated as ambiguous/unreliable rather than silently accepted.
_MAX_GEOMETRIC_MATCH_DEVIATION = 0.5

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

    def rewrite_heavy_atom_transfer(
        self,
        reference_gro: Path | str,
        reference_itp: Path | str,
        target_pose_pdb: Path | str,
        output_dir: Path | str,
    ) -> LigandPoseRewriteResult:
        """
        Transfer receptor-specific heavy-atom coordinates onto the normalized
        LigParGen topology, then reconstruct hydrogen positions (mode
        "heavy_atom_transfer").

        Unlike ``rewrite()``, the target pose does not need to already
        contain hydrogens or match the reference atom count — only the
        heavy atoms need to correspond. This is the correct mode when the
        same ligand adopts different torsions/conformations across
        receptors, where a single global (Kabsch) rotation cannot align
        every atom simultaneously.

        Correspondence strategy (first that succeeds is used):
          1. Positional — heavy atoms compared in file order, accepted only
             when every position's element matches exactly. Robust to
             arbitrary rotation/translation between reference and target,
             which is the expected case when the same ligand file was
             docked independently into different receptors.
          2. Geometric fallback — centroid-centered nearest-neighbour
             matching per element (as used by ``rewrite()``), accepted only
             when every matched pair falls within a sanity deviation
             threshold after centroid alignment.

        Hydrogens are repositioned using each hydrogen's bonded heavy atom's
        local frame: a Kabsch rotation fit over that heavy atom's 1- and
        2-bond heavy neighbours (falling back to the whole-molecule
        rotation when too few neighbours exist for a well-conditioned fit).
        Heavy atoms are never moved from the exact target-pose coordinates
        and are never included in any global optimization.

        Fails explicitly (returns success=False, does not guess) when:
          - heavy-atom counts or element composition disagree,
          - neither correspondence strategy is confident,
          - the reference ITP has no ``[ bonds ]`` section, or a hydrogen's
            bonded heavy atom cannot be determined from it.
        """
        reference_gro = Path(reference_gro)
        reference_itp = Path(reference_itp)
        target_pose_pdb = Path(target_pose_pdb)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not reference_gro.exists():
            return _fail(f"Reference GRO not found: {reference_gro}")
        if not reference_itp.exists():
            return _fail(f"Reference ITP not found: {reference_itp}")
        if not target_pose_pdb.exists():
            return _fail(f"Target pose PDB not found: {target_pose_pdb}")

        try:
            ref_gro = parse_gro(reference_gro)
        except Exception as exc:
            return _fail(f"Cannot parse reference GRO: {exc}")
        try:
            ref_itp = parse_itp(reference_itp)
        except Exception as exc:
            return _fail(f"Cannot parse reference ITP: {exc}")

        n_ref = ref_gro.atom_count
        if len(ref_itp.atoms) != n_ref:
            return _fail(
                f"ITP atom count ({len(ref_itp.atoms)}) does not match "
                f"GRO atom count ({n_ref})"
            )
        if not ref_itp.bonds:
            return _fail(
                "Reference ITP has no [ bonds ] section — cannot determine "
                "hydrogen connectivity for pose reconstruction."
            )

        pdb_raw = parse_pdb_ligand(target_pose_pdb)
        if pdb_raw.get("error"):
            return _fail(f"Cannot parse target pose PDB: {pdb_raw['error']}")

        elements = [_element_from_mass(a.mass) for a in ref_itp.atoms]
        heavy_idx = [i for i, e in enumerate(elements) if e != "H"]
        h_idx = [i for i, e in enumerate(elements) if e == "H"]

        target_heavy = [a for a in pdb_raw["atoms"] if a["element"].upper() != "H"]
        if len(heavy_idx) != len(target_heavy):
            return _fail(
                f"Heavy-atom count mismatch: reference has {len(heavy_idx)} heavy "
                f"atoms, target pose has {len(target_heavy)} heavy atoms."
            )

        ref_coords = np.array([[a.x, a.y, a.z] for a in ref_gro.atoms], dtype=float)
        target_heavy_coords = np.array(
            [[a["x"] * _ANG_TO_NM, a["y"] * _ANG_TO_NM, a["z"] * _ANG_TO_NM] for a in target_heavy],
            dtype=float,
        )
        target_heavy_elements = [a["element"] for a in target_heavy]
        ref_heavy_elements = [elements[i] for i in heavy_idx]

        warnings: list[str] = []
        positional_ok = len(ref_heavy_elements) == len(target_heavy_elements) and all(
            r.upper() == t.upper() for r, t in zip(ref_heavy_elements, target_heavy_elements)
        )

        if positional_ok:
            heavy_target_for_ref = {
                heavy_idx[k]: target_heavy_coords[k] for k in range(len(heavy_idx))
            }
            mapping_method = "atom_order"
        else:
            try:
                mapping = _match_atoms(
                    ref_heavy_elements, ref_coords[heavy_idx],
                    target_heavy_elements, target_heavy_coords,
                )
            except ValueError as exc:
                return _fail(f"Heavy-atom correspondence could not be determined: {exc}")

            matched = target_heavy_coords[mapping]
            centered_ref = ref_coords[heavy_idx] - ref_coords[heavy_idx].mean(axis=0)
            centered_new = matched - matched.mean(axis=0)
            pair_dists = np.linalg.norm(centered_ref - centered_new, axis=1)
            max_dev = float(pair_dists.max()) if len(pair_dists) else 0.0
            if max_dev > _MAX_GEOMETRIC_MATCH_DEVIATION:
                return _fail(
                    "Ambiguous heavy-atom mapping: geometric nearest-neighbour "
                    f"matching produced a deviation of {max_dev:.3f} nm (limit "
                    f"{_MAX_GEOMETRIC_MATCH_DEVIATION} nm) after centroid alignment. "
                    "Reference and target pose do not confidently share atom "
                    "ordering or orientation — provide a chemically explicit "
                    "reference or verify the poses manually."
                )
            heavy_target_for_ref = {heavy_idx[k]: matched[k] for k in range(len(heavy_idx))}
            mapping_method = "geometric_nearest_neighbor"
            warnings.append(
                "Heavy-atom order differed between reference topology and target "
                "pose; used geometric nearest-neighbour matching instead "
                "(lower confidence than atom-order correspondence)."
            )

        new_coords = ref_coords.copy()
        for i, coord in heavy_target_for_ref.items():
            new_coords[i] = coord

        # ── Hydrogen reconstruction via local rigid-frame transfer ───────────
        heavy_neighbors: dict[int, list[int]] = defaultdict(list)
        h_parent: dict[int, int] = {}
        for a1, a2 in ref_itp.bonds:
            i, j = a1 - 1, a2 - 1
            if not (0 <= i < n_ref and 0 <= j < n_ref):
                continue
            el_i, el_j = elements[i], elements[j]
            if el_i != "H" and el_j != "H":
                heavy_neighbors[i].append(j)
                heavy_neighbors[j].append(i)
            elif el_i == "H" and el_j != "H":
                h_parent.setdefault(i, j)
            elif el_j == "H" and el_i != "H":
                h_parent.setdefault(j, i)

        missing_parent = [i for i in h_idx if i not in h_parent]
        if missing_parent:
            names = [ref_itp.atoms[i].atom_name for i in missing_parent]
            return _fail(
                "Cannot reconstruct hydrogen positions: no bonded heavy atom "
                f"found in the reference ITP [ bonds ] section for hydrogen(s) {names}."
            )

        try:
            _, global_R, _ = align(ref_coords[heavy_idx], new_coords[heavy_idx])
        except Exception as exc:
            return _fail(
                f"Could not compute whole-molecule alignment for hydrogen placement: {exc}"
            )

        for h in h_idx:
            parent = h_parent[h]
            R = _local_rotation(parent, heavy_neighbors, ref_coords, new_coords, global_R)
            offset_ref = ref_coords[h] - ref_coords[parent]
            new_coords[h] = new_coords[parent] + offset_ref @ R

        # ── Build output GRO ──────────────────────────────────────────────
        residue_name = ref_gro.atoms[0].residue_name
        out_atoms: list[GroAtom] = []
        for i, ref_atom in enumerate(ref_gro.atoms):
            x, y, z = new_coords[i]
            out_atoms.append(
                GroAtom(
                    residue_number=ref_atom.residue_number,
                    residue_name=ref_atom.residue_name,
                    atom_name=ref_atom.atom_name,
                    atom_number=ref_atom.atom_number,
                    x=float(x), y=float(y), z=float(z),
                )
            )

        out_gro = GroFile(
            title=f"{residue_name} pose (heavy-atom transfer)",
            atoms=out_atoms,
            box=ref_gro.box,
        )
        output_path = output_dir / f"{residue_name}_pose_parameterized.gro"
        try:
            write_gro(out_gro, output_path)
        except Exception as exc:
            return _fail(f"Cannot write output GRO: {exc}")

        try:
            rmsd = _heavy_atom_rmsd(out_gro, ref_gro, elements)
        except Exception as exc:
            logger.warning("RMSD computation failed: %s", exc)
            rmsd = None

        logger.info(
            "Heavy-atom pose transfer: %d heavy + %d H -> %s (mapping=%s)",
            len(heavy_idx), len(h_idx), output_path.name, mapping_method,
        )

        return LigandPoseRewriteResult(
            success=True,
            output_path=output_path,
            ligand_residue_name=residue_name,
            atoms_written=n_ref,
            rmsd_from_reference=rmsd,
            hydrogens_reconstructed=len(h_idx) > 0,
            heavy_atom_mapping_method=mapping_method,
            warnings=warnings,
        )


def _local_rotation(
    center_idx: int,
    heavy_neighbors: dict[int, list[int]],
    ref_coords: np.ndarray,
    new_coords: np.ndarray,
    fallback_R: np.ndarray,
) -> np.ndarray:
    """
    Kabsch rotation derived from a heavy atom's local bonded neighbourhood
    (1- and 2-bond heavy neighbours). Falls back to the whole-molecule
    rotation when too few neighbours are available for a well-conditioned
    local frame (e.g. a terminal heavy atom bonded only to hydrogens).
    """
    pts = {center_idx}
    one_bond = set(heavy_neighbors.get(center_idx, []))
    pts |= one_bond
    two_bond: set[int] = set()
    for n in one_bond:
        two_bond |= set(heavy_neighbors.get(n, []))
    pts |= two_bond

    if len(pts) < _MIN_LOCAL_FRAME_POINTS:
        return fallback_R

    idx = sorted(pts)
    ref_local = ref_coords[idx] - ref_coords[center_idx]
    new_local = new_coords[idx] - new_coords[center_idx]
    try:
        return kabsch_rotation(ref_local, new_local)
    except Exception:
        return fallback_R


def _fail(message: str) -> LigandPoseRewriteResult:
    logger.error("LigandPoseRewriter: %s", message)
    return LigandPoseRewriteResult(success=False, error=message)
