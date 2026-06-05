"""
Pydantic models for the protein-ligand LigParGen workflow (Phase 1).

These models represent the data contracts between workflow stages.
No builder, compiler, or decision engine integration yet.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


class LigandExportResult(BaseModel):
    """Outcome of exporting a ligand structure for external parameterization."""

    success: bool
    exported_path: Optional[Path] = None
    molecule_name: str = ""
    atom_count: int = 0
    error: Optional[str] = None
    # RMSD of heavy atoms between input and processed output (via Kabsch alignment)
    heavy_atom_rmsd: Optional[float] = None


class LigandPreparationValidationResult(BaseModel):
    """Validation of a ligand structure file before parameterization."""

    valid: bool
    file_path: Optional[Path] = None
    molecule_name: str = ""
    atom_count: int = 0
    has_hydrogens: bool = False
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class LigParGenImportValidationResult(BaseModel):
    """Validation of topology files returned by LigParGen (not yet implemented)."""

    valid: bool
    itp_path: Optional[Path] = None
    gro_path: Optional[Path] = None
    molecule_name: str = ""
    atom_count: int = 0
    total_charge: float = 0.0
    charge_integer: bool = False
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class LigandPoseRewriteResult(BaseModel):
    """Result of rewriting ligand coordinates into the system .gro file."""

    success: bool
    output_path: Optional[Path] = None
    ligand_residue_name: str = ""
    atoms_written: int = 0
    rmsd_from_reference: Optional[float] = None
    error: Optional[str] = None


class ProteinLigandAssemblyReport(BaseModel):
    """Summary of assembling protein + ligand into a single simulation system."""

    success: bool
    output_gro: Optional[Path] = None
    output_top: Optional[Path] = None
    protein_atom_count: int = 0
    ligand_atom_count: int = 0
    total_atom_count: int = 0
    ligand_residue_name: str = ""
    ligand_itp_included: bool = False
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class GromppPreflightReport(BaseModel):
    """Pre-flight check before running grompp on a protein-ligand system."""

    ready: bool
    mdp_path: Optional[Path] = None
    gro_path: Optional[Path] = None
    top_path: Optional[Path] = None
    ligand_itp_present: bool = False
    ligand_charge_integer: bool = False
    ligand_total_charge: Optional[float] = None
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)
