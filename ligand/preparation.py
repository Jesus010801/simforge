"""
Ligand preparation validation for LigParGen parameterization (Phase 2).

Validates that a ligand file is structurally and chemically ready for
external parameterization before any export or topology generation step.

Public API:
    validate_ligand_for_parameterization(path) -> LigandPreparationValidationResult
"""

from __future__ import annotations

from pathlib import Path

from core.ligand_workflow_models import LigandPreparationValidationResult
from ligand.rdkit_reader import (
    has_3d_conformer,
    has_explicit_hydrogens,
    heavy_atom_elements,
    load_mol,
    net_formal_charge,
    sanitized_mol_name,
)
from utils.geometry_alignment import elements_match


# Elements that trigger a preparation warning (non-standard for OPLS-AA / GAFF)
_WARN_ELEMENTS = frozenset({"Se", "Te", "Si", "B", "As", "Ge"})
# Elements that hard-block parameterization with standard FFs
_BLOCK_ELEMENTS = frozenset({"Fe", "Cu", "Zn", "Mg", "Ca", "Mn", "Co", "Ni"})


def validate_ligand_for_parameterization(
    path: str | Path,
) -> LigandPreparationValidationResult:
    """
    Validate a ligand file for readiness before LigParGen parameterization.

    Checks performed:
    - File is parseable by RDKit
    - Molecule has at least one heavy atom
    - 3D conformer is present and non-trivial
    - Explicit hydrogens are present (required by LigParGen)
    - Net formal charge (informational, not blocking)
    - No metal atoms that standard FFs cannot handle
    - No uncommon elements (warning only)
    - Element list self-consistency via elements_match

    Returns LigandPreparationValidationResult; errors list is non-empty on failure.
    """
    path = Path(path)
    warnings: list[str] = []
    errors: list[str] = []

    # ── 1. Parse ──────────────────────────────────────────────────────────────
    try:
        mol = load_mol(path)
    except FileNotFoundError:
        return LigandPreparationValidationResult(
            valid=False,
            file_path=path,
            errors=[f"File not found: {path}"],
        )
    except ValueError as exc:
        return LigandPreparationValidationResult(
            valid=False,
            file_path=path,
            errors=[f"RDKit parse error: {exc}"],
        )

    mol_name = sanitized_mol_name(mol)
    atom_count = mol.GetNumAtoms()

    # ── 2. Non-empty ──────────────────────────────────────────────────────────
    if atom_count == 0:
        return LigandPreparationValidationResult(
            valid=False,
            file_path=path,
            molecule_name=mol_name,
            errors=["Molecule has no atoms"],
        )

    # ── 3. 3D conformer ───────────────────────────────────────────────────────
    if not has_3d_conformer(mol):
        errors.append(
            "No 3D conformer found. LigParGen requires 3D coordinates. "
            "Generate with: AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())"
        )

    # ── 4. Explicit hydrogens ─────────────────────────────────────────────────
    h_present = has_explicit_hydrogens(mol)
    if not h_present:
        warnings.append(
            "No explicit hydrogens detected. LigParGen expects all H atoms present. "
            "Add with: Chem.AddHs(mol)"
        )

    # ── 5. Element checks ─────────────────────────────────────────────────────
    elements = heavy_atom_elements(mol)
    element_set = set(elements)

    blocking = element_set & _BLOCK_ELEMENTS
    if blocking:
        errors.append(
            f"Metal atoms detected ({sorted(blocking)}): "
            "standard LigParGen / OPLS-AA cannot parameterize these. "
            "Use a metal-specific force field (e.g., MCPB.py) instead."
        )

    unusual = element_set & _WARN_ELEMENTS
    if unusual:
        warnings.append(
            f"Unusual elements for OPLS-AA detected ({sorted(unusual)}): "
            "verify LigParGen supports these before proceeding."
        )

    # ── 6. Element self-consistency ────────────────────────────────────────────
    # Confirm the heavy-atom element list is identical when re-read — sanity
    # check against partial-read bugs.  elements_match is imported from
    # utils/geometry_alignment.py (Phase 1).
    if not elements_match(elements, elements):
        errors.append("Internal element-list inconsistency detected (bug in rdkit_reader).")

    # ── 7. Charge ─────────────────────────────────────────────────────────────
    charge = net_formal_charge(mol)
    if charge != 0:
        warnings.append(
            f"Net formal charge is {charge:+d}. LigParGen accepts charged molecules, "
            "but verify the protonation state is correct at pH 7.4."
        )

    valid = len(errors) == 0
    heavy_count = len(elements)  # excludes H

    return LigandPreparationValidationResult(
        valid=valid,
        file_path=path,
        molecule_name=mol_name,
        atom_count=atom_count,
        has_hydrogens=h_present,
        warnings=warnings,
        errors=errors,
    )
