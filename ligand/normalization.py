"""
LigParGen identity normalization layer.

LigParGen often returns generic molecule names (H, UNK, LIG, MOL, UNL).
This module normalizes those to deterministic SimForge internal identifiers
such as L01, L02, L03, while preserving all atomic data and bonded terms.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field, field_validator

from utils.gro_parser import parse_gro, write_gro
from utils.itp_parser import parse_itp

logger = logging.getLogger(__name__)

GENERIC_LIGPARGEN_NAMES: frozenset[str] = frozenset({"H", "UNK", "LIG", "MOL", "UNL"})

# GROMACS residue/molecule names: 1-5 uppercase alphanumeric characters
_GROMACS_NAME_RE = re.compile(r"^[A-Z0-9]{1,5}$")


class LigandIdentity(BaseModel):
    component_id: str
    display_name: str
    source_filename: str
    internal_id: str
    residue_name: str
    moleculetype: str
    source_moleculetype: Optional[str] = None

    @field_validator("moleculetype", "residue_name")
    @classmethod
    def _gromacs_safe(cls, v: str) -> str:
        if not _GROMACS_NAME_RE.match(v):
            raise ValueError(
                f"Name {v!r} must be 1-5 uppercase alphanumeric characters (GROMACS-safe)"
            )
        return v


class LigandNormalizationResult(BaseModel):
    identity: LigandIdentity
    normalized_gro: Path
    normalized_itp: Path
    identity_json: Path
    normalization_report: Path
    source_was_generic: bool
    warnings: list[str] = Field(default_factory=list)

    model_config = {"arbitrary_types_allowed": True}


def normalize_ligpargen_outputs(
    ligand_gro: Path | str,
    ligand_itp: Path | str,
    identity: LigandIdentity,
    output_dir: Path | str,
) -> LigandNormalizationResult:
    """Normalize LigParGen output files to SimForge internal naming conventions.

    Rewrites residue/molecule names to identity.residue_name / identity.moleculetype
    while preserving all atomic data, coordinates, bonded terms, and comments.
    """
    ligand_gro = Path(ligand_gro)
    ligand_itp = Path(ligand_itp)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []

    gro = parse_gro(ligand_gro)
    itp = parse_itp(ligand_itp)

    # Detect source moleculetype from parsed ITP
    source_moleculetype = itp.moleculetype.name if itp.moleculetype else None
    identity = identity.model_copy(update={"source_moleculetype": source_moleculetype})

    source_was_generic = source_moleculetype in GENERIC_LIGPARGEN_NAMES
    if source_was_generic:
        msg = (
            f"Source moleculetype {source_moleculetype!r} is a generic LigParGen name; "
            f"normalizing to {identity.moleculetype!r}"
        )
        warnings.append(msg)
        logger.warning(msg)

    # Rewrite GRO: update residue names, preserve everything else
    for atom in gro.atoms:
        atom.residue_name = identity.residue_name
    out_gro = output_dir / f"{identity.internal_id}.gro"
    write_gro(gro, out_gro)

    # Rewrite ITP: update moleculetype name and residue names in [ atoms ]
    out_itp = output_dir / f"{identity.internal_id}.itp"
    _rewrite_itp(
        ligand_itp,
        out_itp,
        new_moleculetype_name=identity.moleculetype,
        new_residue_name=identity.residue_name,
    )

    # Write identity.json
    identity_json = output_dir / "identity.json"
    identity_json.write_text(identity.model_dump_json(indent=2))

    # Write normalization_report.json
    report = {
        "source_gro": str(ligand_gro),
        "source_itp": str(ligand_itp),
        "source_moleculetype": source_moleculetype,
        "normalized_moleculetype": identity.moleculetype,
        "normalized_residue_name": identity.residue_name,
        "internal_id": identity.internal_id,
        "atom_count": gro.atom_count,
        "total_charge": itp.total_charge,
        "source_was_generic": source_was_generic,
        "warnings": warnings,
    }
    normalization_report = output_dir / "normalization_report.json"
    normalization_report.write_text(json.dumps(report, indent=2))

    return LigandNormalizationResult(
        identity=identity,
        normalized_gro=out_gro,
        normalized_itp=out_itp,
        identity_json=identity_json,
        normalization_report=normalization_report,
        source_was_generic=source_was_generic,
        warnings=warnings,
    )


def _replace_nth_field(line: str, n: int, new_value: str) -> str:
    """Replace the nth (0-indexed) whitespace-separated token in line, preserving spacing."""
    # re.split with a capturing group keeps the separators in the result list
    tokens = re.split(r"(\s+)", line)
    # tokens alternates: [possible_empty_prefix, sep, tok, sep, tok, ...]
    field_indices = [i for i, t in enumerate(tokens) if t and not t[0].isspace()]
    if n < len(field_indices):
        tokens[field_indices[n]] = new_value
    return "".join(tokens)


def _rewrite_itp(
    source: Path,
    dest: Path,
    new_moleculetype_name: str,
    new_residue_name: str,
) -> None:
    """Rewrite ITP replacing only the moleculetype name and residue names in [ atoms ]."""
    raw_lines = source.read_text().splitlines()
    out: list[str] = []
    current_section: Optional[str] = None
    mol_data_done = False

    for line in raw_lines:
        data_part = line.split(";")[0]
        stripped = data_part.strip()

        if stripped.startswith("["):
            section_name = stripped.strip("[]").strip().lower()
            current_section = section_name
            mol_data_done = False
            out.append(line)
            continue

        if current_section == "moleculetype" and not mol_data_done:
            if stripped:
                parts = stripped.split()
                if len(parts) >= 2:
                    # First data line: "NAME   nrexcl" — replace just the name field
                    comment_idx = line.find(";")
                    if comment_idx >= 0:
                        data_section = line[:comment_idx]
                        comment_section = line[comment_idx:]
                    else:
                        data_section = line
                        comment_section = ""
                    new_data = _replace_nth_field(data_section, 0, new_moleculetype_name)
                    out.append(new_data + comment_section)
                    mol_data_done = True
                    continue

        if current_section == "atoms" and stripped:
            parts = stripped.split()
            # Atom lines have at least 7 fields: nr type resnr resname atomname cgnr charge
            if len(parts) >= 7:
                try:
                    int(parts[0])
                    comment_idx = line.find(";")
                    if comment_idx >= 0:
                        data_section = line[:comment_idx]
                        comment_section = line[comment_idx:]
                    else:
                        data_section = line
                        comment_section = ""
                    # residue name is the 4th field (index 3)
                    new_data = _replace_nth_field(data_section, 3, new_residue_name)
                    out.append(new_data + comment_section)
                    continue
                except ValueError:
                    pass

        out.append(line)

    dest.write_text("\n".join(out) + "\n")
