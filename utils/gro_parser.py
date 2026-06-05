"""
Parser and writer for GROMACS .gro coordinate files.

Format reference (columns are fixed-width, 5+5+5+8+8+8[+8+8+8]):
  Line 1 : title (arbitrary text)
  Line 2 : atom count
  Lines 3..(N+2): atom records
  Last line: box vectors (3 or 9 floats)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class GroAtom:
    residue_number: int
    residue_name: str
    atom_name: str
    atom_number: int
    x: float  # nm
    y: float  # nm
    z: float  # nm
    vx: Optional[float] = None
    vy: Optional[float] = None
    vz: Optional[float] = None


@dataclass
class GroFile:
    title: str
    atoms: list[GroAtom] = field(default_factory=list)
    box: list[float] = field(default_factory=list)

    @property
    def atom_count(self) -> int:
        return len(self.atoms)


def parse_gro(path: str | Path) -> GroFile:
    """Parse a GROMACS .gro file into a GroFile record."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Too few lines in .gro file: {path}")

    title = lines[0]

    try:
        declared_count = int(lines[1].strip())
    except ValueError:
        raise ValueError(f"Line 2 must be atom count, got: {lines[1]!r}")

    atom_lines = lines[2 : 2 + declared_count]
    if len(atom_lines) != declared_count:
        raise ValueError(
            f"Declared {declared_count} atoms but found {len(atom_lines)} atom lines"
        )

    atoms: list[GroAtom] = []
    for lineno, line in enumerate(atom_lines, start=3):
        atom = _parse_atom_line(line, lineno)
        atoms.append(atom)

    box_line = lines[2 + declared_count]
    box = [float(v) for v in box_line.split()]

    return GroFile(title=title, atoms=atoms, box=box)


def _parse_atom_line(line: str, lineno: int) -> GroAtom:
    """Parse one atom record from a .gro file line."""
    # Fixed-width format: resnum(5) resname(5) atomname(5) atomnum(5) x(8) y(8) z(8) [vx(8) vy(8) vz(8)]
    if len(line) < 44:
        raise ValueError(f"Line {lineno} too short for a .gro atom record: {line!r}")

    try:
        residue_number = int(line[0:5])
        residue_name = line[5:10].strip()
        atom_name = line[10:15].strip()
        atom_number = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Cannot parse atom line {lineno}: {line!r}") from exc

    vx = vy = vz = None
    if len(line) >= 68:
        try:
            vx = float(line[44:52])
            vy = float(line[52:60])
            vz = float(line[60:68])
        except ValueError:
            pass

    return GroAtom(
        residue_number=residue_number,
        residue_name=residue_name,
        atom_name=atom_name,
        atom_number=atom_number,
        x=x, y=y, z=z,
        vx=vx, vy=vy, vz=vz,
    )


def write_gro(gro: GroFile, path: str | Path) -> None:
    """Write a GroFile back to disk in standard .gro format."""
    lines: list[str] = []
    lines.append(gro.title)
    lines.append(f"{gro.atom_count:5d}")

    for atom in gro.atoms:
        line = (
            f"{atom.residue_number % 100000:5d}"
            f"{atom.residue_name:<5s}"
            f"{atom.atom_name:>5s}"
            f"{atom.atom_number % 100000:5d}"
            f"{atom.x:8.3f}"
            f"{atom.y:8.3f}"
            f"{atom.z:8.3f}"
        )
        if atom.vx is not None and atom.vy is not None and atom.vz is not None:
            line += f"{atom.vx:8.4f}{atom.vy:8.4f}{atom.vz:8.4f}"
        lines.append(line)

    lines.append("   ".join(f"{v:8.5f}" for v in gro.box))
    Path(path).write_text("\n".join(lines) + "\n")
