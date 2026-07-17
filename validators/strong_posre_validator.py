"""
validators/strong_posre_validator.py
─────────────────────────────────────
Pre-grompp validator for [ position_restraints ] atom index bounds.

Walks the topology include graph keeping track of the current active
moleculetype and its atom count.  For every [ position_restraints ]
section it checks that every atom index satisfies:

    1 <= index <= n_atoms_in_active_moleculetype

On violation it raises a descriptive error message prefixed with
TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Optional


def _count_atoms_in_itp(path: Path) -> int:
    """Count atoms listed in the first [ atoms ] section of a topology file."""
    in_atoms = False
    count = 0
    for line in path.read_text().splitlines():
        ls = line.strip()
        if re.match(r'^\[\s*atoms\s*\]', ls, re.IGNORECASE):
            in_atoms = True
            continue
        if in_atoms:
            if ls.startswith('['):
                break
            if ls and not ls.startswith(';'):
                parts = ls.split()
                if parts and parts[0].isdigit():
                    count += 1
    return count


def _active_moleculetype_at(path: Path, visited: Optional[set] = None) -> list[tuple[str, int, Path]]:
    """
    Walk include graph from `path`.  Return list of (mol_name, n_atoms, itp_path)
    for each moleculetype discovered.
    """
    if visited is None:
        visited = set()
    key = str(path.resolve())
    if key in visited:
        return []
    visited.add(key)

    result: list[tuple[str, int, Path]] = []
    try:
        text = path.read_text()
    except OSError:
        return result

    current_mol: Optional[str] = None
    in_atoms = False
    atom_count = 0

    for line in text.splitlines():
        ls = line.strip()

        # Follow #include directives
        if ls.startswith('#include'):
            pts = ls.split('"')
            if len(pts) >= 2:
                inc_name = pts[1]
                cand = (path.parent / inc_name).resolve()
                if cand.exists():
                    result += _active_moleculetype_at(cand, visited)
            continue

        sec = re.match(r'^\[\s*(\w+)\s*\]', ls, re.IGNORECASE)
        if sec:
            sec_name = sec.group(1).lower()
            if sec_name == 'moleculetype':
                if current_mol is not None and atom_count > 0:
                    result.append((current_mol, atom_count, path))
                current_mol = None
                atom_count = 0
                in_atoms = False
            elif sec_name == 'atoms':
                in_atoms = True
            else:
                in_atoms = False
            continue

        # Inside [ moleculetype ] — grab the name (first non-comment data line)
        if current_mol is None and not in_atoms:
            if ls and not ls.startswith(';') and not ls.startswith('['):
                parts = ls.split()
                if parts:
                    current_mol = parts[0]
            continue

        # Inside [ atoms ] — count rows
        if in_atoms and ls and not ls.startswith(';'):
            parts = ls.split()
            if parts and parts[0].isdigit():
                atom_count += 1

    if current_mol is not None and atom_count > 0:
        result.append((current_mol, atom_count, path))

    return result


def validate_posre_bounds(topol_path: Path) -> list[str]:
    """
    Walk *topol_path* and every included file.  For each [ position_restraints ]
    section, verify that all atom indices are within the local moleculetype
    atom count.

    Returns a list of error strings (empty = OK).
    Each error is prefixed with TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR.
    """
    errors: list[str] = []
    _walk_posre(topol_path, None, 0, errors, visited=set())
    return errors


def _walk_posre(
    path: Path,
    active_mol: Optional[str],
    active_natoms: int,
    errors: list[str],
    visited: set,
) -> tuple[Optional[str], int]:
    """Recursive walk; returns (active_mol, active_natoms) after processing path."""
    key = str(path.resolve())
    if key in visited:
        return active_mol, active_natoms
    visited.add(key)

    try:
        text = path.read_text()
    except OSError:
        return active_mol, active_natoms

    in_posre = False
    in_atoms = False
    in_mol   = False
    local_atom_count = 0

    mol_name   = active_mol
    n_atoms    = active_natoms

    for line in text.splitlines():
        ls = line.strip()

        # #include: recurse, then pick up updated state
        if ls.startswith('#include'):
            pts = ls.split('"')
            if len(pts) >= 2:
                inc_name = pts[1]
                cand = (path.parent / inc_name).resolve()
                if cand.exists():
                    mol_name, n_atoms = _walk_posre(
                        cand, mol_name, n_atoms, errors, visited
                    )
            in_posre = False  # section context resets at include boundaries
            in_atoms = False
            in_mol   = False
            continue

        # Skip preprocessor conditionals
        if ls.startswith('#'):
            continue

        sec = re.match(r'^\[\s*(\w+)\s*\]', ls, re.IGNORECASE)
        if sec:
            sec_name = sec.group(1).lower()
            if sec_name == 'moleculetype':
                mol_name = None
                n_atoms  = 0
                local_atom_count = 0
                in_mol  = True
                in_atoms = False
                in_posre = False
            elif sec_name == 'atoms':
                in_atoms = True
                in_mol   = False
                in_posre = False
            elif sec_name == 'position_restraints':
                in_posre = True
                in_atoms = False
                in_mol   = False
                if mol_name and n_atoms == 0:
                    n_atoms = local_atom_count
            else:
                if sec_name not in ('bonds', 'angles', 'dihedrals', 'pairs',
                                    'exclusions', 'constraints', 'settles',
                                    'virtual_sites2', 'virtual_sites3',
                                    'virtual_sitesn', 'cmap'):
                    in_posre = False
                in_atoms = False
                in_mol   = False
            continue

        if ls.startswith(';') or not ls:
            continue

        if in_mol and mol_name is None:
            parts = ls.split()
            if parts:
                mol_name = parts[0]
            continue

        if in_atoms:
            parts = ls.split()
            if parts and parts[0].isdigit():
                local_atom_count += 1
                n_atoms = local_atom_count
            continue

        if in_posre:
            parts = ls.split()
            if not parts or not parts[0].isdigit():
                continue
            try:
                idx = int(parts[0])
            except ValueError:
                continue
            if n_atoms > 0 and idx > n_atoms:
                errors.append(
                    f"TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR: "
                    f"{path.name} contains atom index {idx}, but "
                    f"{mol_name or 'active moleculetype'} has only {n_atoms} atoms. "
                    f"Position restraint indices must be local to the moleculetype "
                    f"(1-{n_atoms}), not global system indices."
                )

    return mol_name, n_atoms
