"""
Parser for GROMACS .itp topology include files.

Extracts [ moleculetype ] and [ atoms ] sections.
Handles ; comments and blank lines correctly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class MoleculetypeRecord:
    name: str
    nrexcl: int


@dataclass
class AtomRecord:
    number: int          # atom index (1-based in GROMACS)
    atom_type: str
    residue_number: int
    residue_name: str
    atom_name: str
    charge_group: int
    charge: float
    mass: Optional[float] = None


@dataclass
class ItpFile:
    sections: dict[str, list[str]] = field(default_factory=dict)
    moleculetype: Optional[MoleculetypeRecord] = None
    atoms: list[AtomRecord] = field(default_factory=list)

    @property
    def total_charge(self) -> float:
        return sum(a.charge for a in self.atoms)

    def charge_is_integer(self, tolerance: float = 1e-3) -> bool:
        q = self.total_charge
        return abs(q - round(q)) < tolerance


def _strip_comment(line: str) -> str:
    """Remove inline ; comments and strip whitespace."""
    return line.split(";")[0].strip()


def parse_itp(path: str | Path) -> ItpFile:
    """Parse a GROMACS .itp file into an ItpFile record."""
    raw_lines = Path(path).read_text().splitlines()

    sections: dict[str, list[str]] = {}
    current_section: Optional[str] = None

    for line in raw_lines:
        stripped = _strip_comment(line)
        if not stripped:
            continue

        if stripped.startswith("["):
            section_name = stripped.strip("[]").strip().lower()
            current_section = section_name
            sections.setdefault(current_section, [])
        elif current_section is not None:
            sections[current_section].append(stripped)

    result = ItpFile(sections=sections)

    if "moleculetype" in sections:
        result.moleculetype = _parse_moleculetype(sections["moleculetype"])

    if "atoms" in sections:
        result.atoms = _parse_atoms(sections["atoms"])

    return result


def _parse_moleculetype(lines: list[str]) -> Optional[MoleculetypeRecord]:
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                return MoleculetypeRecord(name=parts[0], nrexcl=int(parts[1]))
            except ValueError:
                continue
    return None


def _parse_atoms(lines: list[str]) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            number = int(parts[0])
            atom_type = parts[1]
            residue_number = int(parts[2])
            residue_name = parts[3]
            atom_name = parts[4]
            charge_group = int(parts[5])
            charge = float(parts[6])
            mass = float(parts[7]) if len(parts) >= 8 else None
        except (ValueError, IndexError):
            continue

        atoms.append(AtomRecord(
            number=number,
            atom_type=atom_type,
            residue_number=residue_number,
            residue_name=residue_name,
            atom_name=atom_name,
            charge_group=charge_group,
            charge=charge,
            mass=mass,
        ))
    return atoms
