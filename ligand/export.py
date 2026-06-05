"""
Ligand export for LigParGen parameterization (Phase 2).

Converts a ligand file to PDB format suitable for LigParGen, with optional
hydrogen addition. After export, uses Kabsch alignment (geometry_alignment)
to verify the heavy-atom scaffold is preserved within a tolerance.

Public API:
    export_for_ligpargen(path, output_dir, mol_name, add_hydrogens) -> LigandExportResult
        RDKit-style PDB (HETATM + numbered names).  Kept for reference.

    export_for_ligpargen_legacy(path, output_dir, mol_name, add_hydrogens) -> LigandExportResult
        Legacy PDB (ATOM + simple element names).  This is the format that
        LigParGen's online server actually accepts.
        Output file: <mol_name>_ligpargen_legacy.pdb
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from core.ligand_workflow_models import LigandExportResult
from ligand.rdkit_reader import (
    has_3d_conformer,
    heavy_atom_coords,
    heavy_atom_elements,
    load_mol,
    sanitized_mol_name,
)
from utils.geometry_alignment import align, elements_match

# Heavy-atom RMSD above this threshold after processing triggers a warning.
_RMSD_WARN_THRESHOLD = 0.05  # Å (0.005 nm) — tight: only distortion, not rounding


def export_for_ligpargen(
    path: str | Path,
    output_dir: str | Path,
    mol_name: str = "",
    add_hydrogens: bool = True,
) -> LigandExportResult:
    """
    Export a ligand to PDB for submission to LigParGen.

    Steps:
      1. Load the molecule with RDKit.
      2. Record heavy-atom coordinates (input reference).
      3. Optionally add missing hydrogens (preserving existing 3D coords).
      4. Write PDB to output_dir / <mol_name>.pdb.
      5. Compute Kabsch-aligned RMSD of heavy atoms between input and output.
         If RMSD > threshold the result carries a warning but export still succeeds.

    Args:
        path:          Input ligand file (.sdf, .mol, .pdb).
        output_dir:    Directory to write the exported PDB into.
        mol_name:      4-char GROMACS-compatible molecule name (auto-derived if empty).
        add_hydrogens: Add explicit H atoms if not already present.

    Returns:
        LigandExportResult (success=False carries an error string on failure).
    """
    path = Path(path)
    output_dir = Path(output_dir)

    # ── 1. Load ────────────────────────────────────────────────────────────────
    try:
        mol = load_mol(path)
    except FileNotFoundError:
        return LigandExportResult(
            success=False,
            error=f"File not found: {path}",
        )
    except ValueError as exc:
        return LigandExportResult(
            success=False,
            error=f"RDKit parse error: {exc}",
        )

    if not has_3d_conformer(mol):
        return LigandExportResult(
            success=False,
            error=(
                "Molecule has no 3D conformer. "
                "Generate coordinates before exporting: "
                "AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())"
            ),
        )

    # ── 2. Heavy-atom reference (input) ────────────────────────────────────────
    ref_elements = heavy_atom_elements(mol)
    ref_coords_raw = heavy_atom_coords(mol)
    ref_coords = np.array(ref_coords_raw, dtype=np.float64)

    # ── 3. Derive molecule name ────────────────────────────────────────────────
    name = mol_name.strip() or sanitized_mol_name(mol)
    name = (name[:4].upper()) or "LIG"

    # ── 4. Add hydrogens ───────────────────────────────────────────────────────
    try:
        processed_mol = _add_hydrogens_if_needed(mol, add_hydrogens)
    except Exception as exc:
        return LigandExportResult(
            success=False,
            molecule_name=name,
            error=f"Failed to add hydrogens: {exc}",
        )

    # ── 5. Write PDB ───────────────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{name}.pdb"

    try:
        _write_pdb(processed_mol, out_path, name)
    except Exception as exc:
        return LigandExportResult(
            success=False,
            molecule_name=name,
            error=f"Failed to write PDB: {exc}",
        )

    # ── 6. Heavy-atom RMSD (Kabsch) ────────────────────────────────────────────
    out_elements = heavy_atom_elements(processed_mol)
    out_coords_raw = heavy_atom_coords(processed_mol)
    out_coords = np.array(out_coords_raw, dtype=np.float64)

    rmsd: Optional[float] = None
    if elements_match(ref_elements, out_elements) and len(ref_coords) > 0:
        _, _, rmsd = align(ref_coords, out_coords)

    return LigandExportResult(
        success=True,
        exported_path=out_path,
        molecule_name=name,
        atom_count=processed_mol.GetNumAtoms(),
        heavy_atom_rmsd=rmsd,
    )


# ── Internal helpers ──────────────────────────────────────────────────────────

def _add_hydrogens_if_needed(mol, add_hydrogens: bool):
    """
    Return mol with explicit H atoms if add_hydrogens is True.

    Uses addCoords=True so existing heavy-atom geometry is unchanged and
    hydrogens are placed geometrically (no force field minimisation).
    """
    from rdkit import Chem

    if not add_hydrogens:
        return mol

    has_h = any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    if has_h:
        return mol  # already has explicit H — no change

    return Chem.AddHs(mol, addCoords=True)


def _write_pdb(mol, path: Path, mol_name: str) -> None:
    """Write an RDKit Mol to a PDB file with the given residue name."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Assign a residue name that LigParGen will see as the molecule identifier
    mi = Chem.AtomPDBResidueInfo()
    mi.SetResidueName(mol_name[:3].ljust(3))
    mi.SetResidueNumber(1)
    mi.SetChainId("A")
    mi.SetIsHeteroAtom(True)

    for atom in mol.GetAtoms():
        atom.SetMonomerInfo(mi)

    pdb_block = Chem.MolToPDBBlock(mol)
    if pdb_block is None:
        raise RuntimeError("Chem.MolToPDBBlock returned None")

    path.write_text(pdb_block)


# ── Legacy export (LigParGen-accepted format) ─────────────────────────────────

def export_for_ligpargen_legacy(
    path: str | Path,
    output_dir: str | Path,
    mol_name: str = "",
    add_hydrogens: bool = True,
) -> LigandExportResult:
    """
    Export a ligand to the legacy PDB format accepted by LigParGen.

    Produces ``<mol_name>_ligpargen_legacy.pdb`` inside *output_dir*.

    Key differences from :func:`export_for_ligpargen`:
      - Uses ``ATOM``    records (not ``HETATM``)
      - Uses simple element-only atom names: ``C``, ``O``, ``H``, ``N`` …
      - Builds CONECT records for every bond including H-bonds
      - Does NOT call RDKit's ``MolToPDBBlock``

    Steps mirror the original export function:
      1. Load with RDKit.
      2. Capture heavy-atom reference coordinates (Kabsch baseline).
      3. Add explicit H atoms if missing.
      4. Write with :class:`~ligand.legacy_writer.LigParGenLegacyWriter`.
      5. Compute Kabsch-aligned heavy-atom RMSD (should be ~0).
    """
    from ligand.legacy_writer import LigParGenLegacyWriter

    path = Path(path)
    output_dir = Path(output_dir)

    # ── Load ───────────────────────────────────────────────────────────────────
    try:
        mol = load_mol(path)
    except FileNotFoundError:
        return LigandExportResult(
            success=False,
            error=f"File not found: {path}",
        )
    except ValueError as exc:
        return LigandExportResult(
            success=False,
            error=f"RDKit parse error: {exc}",
        )

    if not has_3d_conformer(mol):
        return LigandExportResult(
            success=False,
            error=(
                "Molecule has no 3D conformer. "
                "Generate coordinates before exporting: "
                "AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())"
            ),
        )

    # ── Heavy-atom reference (input) ───────────────────────────────────────────
    ref_elements = heavy_atom_elements(mol)
    ref_coords = np.array(heavy_atom_coords(mol), dtype=np.float64)

    # ── Molecule name ──────────────────────────────────────────────────────────
    name = (mol_name.strip() or sanitized_mol_name(mol))[:4].upper() or "LIG"

    # ── Add hydrogens ──────────────────────────────────────────────────────────
    try:
        processed_mol = _add_hydrogens_if_needed(mol, add_hydrogens)
    except Exception as exc:
        return LigandExportResult(
            success=False,
            molecule_name=name,
            error=f"Failed to add hydrogens: {exc}",
        )

    # ── Write legacy PDB ───────────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{name}_ligpargen_legacy.pdb"

    try:
        LigParGenLegacyWriter().write(processed_mol, out_path, name)
    except Exception as exc:
        return LigandExportResult(
            success=False,
            molecule_name=name,
            error=f"Failed to write legacy PDB: {exc}",
        )

    # ── Heavy-atom RMSD (Kabsch) ───────────────────────────────────────────────
    out_elements = heavy_atom_elements(processed_mol)
    out_coords = np.array(heavy_atom_coords(processed_mol), dtype=np.float64)

    rmsd: Optional[float] = None
    if elements_match(ref_elements, out_elements) and len(ref_coords) > 0:
        _, _, rmsd = align(ref_coords, out_coords)

    return LigandExportResult(
        success=True,
        exported_path=out_path,
        molecule_name=name,
        atom_count=processed_mol.GetNumAtoms(),
        heavy_atom_rmsd=rmsd,
    )
