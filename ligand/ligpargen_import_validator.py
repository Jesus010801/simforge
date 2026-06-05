"""
Validator for LigParGen topology output files (.gro + .itp).

Validates before any system assembly:
  1. File existence and parseability
  2. Required ITP sections: [ moleculetype ], [ atoms ]
  3. GRO / ITP atom count consistency
  4. Total charge near-integrality
  5. Generic source name detection
  6. Normalization round-trip consistency

Does not assemble systems, edit topol.top, run grompp, or touch the compiler.
"""

from __future__ import annotations

import logging
from pathlib import Path

from core.ligand_workflow_models import LigParGenImportValidationResult
from ligand.normalization import (
    GENERIC_LIGPARGEN_NAMES,
    LigandIdentity,
    normalize_ligpargen_outputs,
)
from utils.gro_parser import parse_gro
from utils.itp_parser import parse_itp

logger = logging.getLogger(__name__)


class LigParGenImportValidator:
    """Validates LigParGen .gro + .itp output before simulation assembly.

    Runs normalization as a consistency check: files that cannot be
    normalized to the given identity are considered invalid.
    """

    def validate(
        self,
        ligand_gro: Path | str,
        ligand_itp: Path | str,
        identity: LigandIdentity,
        work_dir: Path | str,
    ) -> LigParGenImportValidationResult:
        """Validate LigParGen output and normalize to identity.

        Normalized files are written under work_dir/normalized/ when
        all earlier checks pass.  On success, result.gro_path and
        result.itp_path point to the normalized files.
        """
        ligand_gro = Path(ligand_gro)
        ligand_itp = Path(ligand_itp)
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        errors: list[str] = []
        warnings: list[str] = []

        # ── 1. File existence ──────────────────────────────────────────────────
        if not ligand_gro.exists():
            errors.append(f"GRO file not found: {ligand_gro}")
        if not ligand_itp.exists():
            errors.append(f"ITP file not found: {ligand_itp}")
        if errors:
            return _early_fail(ligand_gro, ligand_itp, errors, warnings)

        # ── 2. Parse ───────────────────────────────────────────────────────────
        try:
            gro = parse_gro(ligand_gro)
        except Exception as exc:
            return _early_fail(ligand_gro, ligand_itp, [f"Cannot parse GRO: {exc}"], warnings)

        try:
            itp = parse_itp(ligand_itp)
        except Exception as exc:
            return _early_fail(ligand_gro, ligand_itp, [f"Cannot parse ITP: {exc}"], warnings)

        # ── 3. Required ITP sections ───────────────────────────────────────────
        if itp.moleculetype is None:
            errors.append("ITP is missing required [ moleculetype ] section")
        if not itp.atoms:
            errors.append("ITP is missing required [ atoms ] section or contains no atoms")
        if errors:
            return _early_fail(ligand_gro, ligand_itp, errors, warnings)

        source_moleculetype = itp.moleculetype.name

        # ── 4. Generic name detection ──────────────────────────────────────────
        if source_moleculetype in GENERIC_LIGPARGEN_NAMES:
            warnings.append(
                f"Generic LigParGen moleculetype name {source_moleculetype!r} detected; "
                f"normalizing to {identity.moleculetype!r}"
            )

        # ── 5. GRO / ITP atom count consistency ───────────────────────────────
        gro_count = gro.atom_count
        itp_count = len(itp.atoms)
        if gro_count != itp_count:
            errors.append(
                f"Atom count mismatch: GRO has {gro_count} atoms, "
                f"ITP [ atoms ] has {itp_count} atoms"
            )

        # ── 6. Total charge ────────────────────────────────────────────────────
        total_charge = itp.total_charge
        charge_integer = itp.charge_is_integer()
        if not charge_integer:
            warnings.append(
                f"Non-integer total charge {total_charge:.4f} — "
                "verify parametrization or check for truncation errors"
            )

        if errors:
            return LigParGenImportValidationResult(
                valid=False,
                gro_path=ligand_gro,
                itp_path=ligand_itp,
                molecule_name=source_moleculetype,
                atom_count=gro_count,
                total_charge=total_charge,
                charge_integer=charge_integer,
                errors=errors,
                warnings=warnings,
            )

        # ── 7. Normalize and round-trip consistency check ─────────────────────
        norm_dir = work_dir / "normalized"
        try:
            norm = normalize_ligpargen_outputs(ligand_gro, ligand_itp, identity, norm_dir)
        except Exception as exc:
            errors.append(f"Normalization failed: {exc}")
            return LigParGenImportValidationResult(
                valid=False,
                gro_path=ligand_gro,
                itp_path=ligand_itp,
                molecule_name=source_moleculetype,
                atom_count=gro_count,
                total_charge=total_charge,
                charge_integer=charge_integer,
                errors=errors,
                warnings=warnings,
            )

        norm_gro = parse_gro(norm.normalized_gro)
        norm_itp = parse_itp(norm.normalized_itp)

        if norm_gro.atom_count != gro_count:
            errors.append(
                f"Normalization altered GRO atom count: "
                f"{gro_count} → {norm_gro.atom_count}"
            )
        if len(norm_itp.atoms) != itp_count:
            errors.append(
                f"Normalization altered ITP atom count: "
                f"{itp_count} → {len(norm_itp.atoms)}"
            )

        unexpected_gro = {a.residue_name for a in norm_gro.atoms} - {identity.residue_name}
        if unexpected_gro:
            errors.append(
                f"Normalized GRO contains unexpected residue names: {unexpected_gro}"
            )

        unexpected_itp = {a.residue_name for a in norm_itp.atoms} - {identity.residue_name}
        if unexpected_itp:
            errors.append(
                f"Normalized ITP [ atoms ] contains unexpected residue names: {unexpected_itp}"
            )

        if (
            norm_itp.moleculetype is not None
            and norm_itp.moleculetype.name != identity.moleculetype
        ):
            errors.append(
                f"Normalized ITP [ moleculetype ] is {norm_itp.moleculetype.name!r}, "
                f"expected {identity.moleculetype!r}"
            )

        return LigParGenImportValidationResult(
            valid=len(errors) == 0,
            gro_path=norm.normalized_gro,
            itp_path=norm.normalized_itp,
            molecule_name=source_moleculetype,
            atom_count=gro_count,
            total_charge=total_charge,
            charge_integer=charge_integer,
            errors=errors,
            warnings=warnings,
        )


def _early_fail(
    gro: Path,
    itp: Path,
    errors: list[str],
    warnings: list[str],
) -> LigParGenImportValidationResult:
    return LigParGenImportValidationResult(
        valid=False,
        gro_path=gro,
        itp_path=itp,
        errors=errors,
        warnings=warnings,
    )
