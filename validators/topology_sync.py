"""
validators/topology_sync.py

Synchronizes a topology's [ molecules ] lipid count with the actual lipid
residue count present in a coordinate (.gro) file.

Problem this fixes: inflategro mutates work.gro in place, potentially
removing protein-overlapping lipid residues (hard-overlap or per-iteration
deflation cleanup). If the topology consumed by the next `gmx grompp` call
still lists the pre-removal lipid count, grompp fails with a hard atom-count
mismatch — the coordinate file has fewer atoms than the topology expects.

Real-run evidence (membrane_embedding, GLP-1R):
    work.gro:            30421 atoms, 477 DPP residues
    bootstrap_topol.top:  478 DPP  ->  implies 30471 atoms
    grompp: "number of coordinates in coordinate file does not match topology"

This module counts lipid residues directly from the coordinate file (never
trusts a previously-recorded count or inflategro's log text) and rewrites
only the matching [ molecules ] line, leaving every other molecule entry
(Protein, SOL, ions, ...) untouched. If the lipid entry is absent from the
topology, no line is invented — the sync is reported as unsynchronized.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

_SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\]')
_MOL_LINE_RE = re.compile(r'^(\S+)\s+(\d+)\s*(;.*)?$')


def _parse_gro_atoms_and_lipids(gro_path: Path, lipid_resname: str) -> tuple[int, int, int]:
    """Return (total_atoms, lipid_atom_count, lipid_residue_count)."""
    return _parse_gro_atoms_and_molecules(gro_path, lipid_resname)


def _parse_gro_atoms_and_molecules(gro_path: Path, resname: str) -> tuple[int, int, int]:
    """Return (total_atoms, matching_atom_count, matching_residue_count)
    for any resname (SOL, DPP, ...). Never trusts a prior recorded count —
    always recounts directly from the coordinate file."""
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        return 0, 0, 0
    try:
        n = int(lines[1].strip())
    except ValueError:
        return 0, 0, 0

    atom_count = 0
    resids: set[int] = set()
    for line in lines[2: 2 + n]:
        if len(line) < 15:
            continue
        line_resname = line[5:10].strip()
        if line_resname != resname:
            continue
        atom_count += 1
        try:
            resids.add(int(line[0:5]))
        except ValueError:
            continue
    return n, atom_count, len(resids)


def _find_molecule_count(topol_text: str, mol_name: str) -> Optional[int]:
    """Find mol_name's count, restricted to the [ molecules ] section only —
    a bare '<name> <count>' pattern can otherwise collide with an unrelated
    [ moleculetype ] header line (e.g. 'DPP   3' for nrexcl) if that block
    happens to be inlined in the same file."""
    in_molecules = False
    for raw_line in topol_text.splitlines():
        line = raw_line.strip()
        sec_m = _SECTION_RE.match(line)
        if sec_m:
            in_molecules = sec_m.group(1).lower() == "molecules"
            continue
        if not in_molecules or not line or line.startswith(";"):
            continue
        m = _MOL_LINE_RE.match(line)
        if m and m.group(1) == mol_name:
            return int(m.group(2))
    return None


def _update_molecule_count(topol_text: str, mol_name: str, new_count: int) -> tuple[str, bool]:
    """Rewrite only mol_name's line inside [ molecules ]. Returns (new_text, updated)."""
    out_lines: list[str] = []
    in_molecules = False
    updated = False
    for raw_line in topol_text.splitlines():
        stripped = raw_line.strip()
        sec_m = _SECTION_RE.match(stripped)
        if sec_m:
            in_molecules = sec_m.group(1).lower() == "molecules"
            out_lines.append(raw_line)
            continue
        if in_molecules and not updated and stripped and not stripped.startswith(";"):
            m = _MOL_LINE_RE.match(stripped)
            if m and m.group(1) == mol_name:
                comment = m.group(3)
                new_line = f"{mol_name:<16}{new_count}"
                if comment:
                    new_line += f"   {comment}"
                out_lines.append(new_line)
                updated = True
                continue
        out_lines.append(raw_line)
    return "\n".join(out_lines), updated


def sync_topology_molecule_count(
    gro_path:      Path | str,
    topol_path:    Path | str,
    resname:       str,
    report_path:   Path | str | None = None,
    atoms_per_molecule: Optional[int] = None,
    molecule_count_override: Optional[int] = None,
) -> dict:
    """
    Synchronize `resname`'s molecule count in topol_path's [ molecules ]
    section with the actual count in gro_path.

    Generalizes the lipid-only sync this module originally shipped with
    (sync_topology_lipid_count, now a thin wrapper around this function) to
    any resname — in particular SOL, replacing the three duplicated
    unscoped `^SOL\\s+\\d+` regex updaters that used to live in
    validators/membrane_protein_spatial_classifier.py, validators/pore_hydration.py,
    and the generated run_clean_water.py template.

    - Counts `resname` molecules directly from the coordinate file — never
      from a prior expected count, a removed-water delta, or log text.
    - If `atoms_per_molecule` is given (e.g. 3 for SOL), the molecule count
      is derived as atom_count // atoms_per_molecule instead of counting
      distinct residue numbers — this avoids undercounting on a .gro file
      whose 5-digit residue-number field has wrapped past 99999, which a
      distinct-resid count cannot detect.
    - Rewrites only the matching line inside [ molecules ]; every other
      molecule entry (Protein, other solvents, ions, ...) is left untouched.
    - If the entry is absent from the topology, no entry is invented; the
      report is written with synchronized=False.

    Returns the report dict (also written to report_path, default
    "topology_sync_report.json" next to topol_path). Report keys:
        coordinate_file, topology_file, resname,
        old_count, new_count, removed_count,
        old_expected_atoms, new_expected_atoms, coordinate_atoms,
        synchronized, warnings
    """
    gro_path   = Path(gro_path)
    topol_path = Path(topol_path)
    report_path = (
        Path(report_path) if report_path is not None
        else topol_path.parent / "topology_sync_report.json"
    )

    report: dict = {
        "coordinate_file":    str(gro_path),
        "topology_file":      str(topol_path),
        "resname":            resname,
        "old_count":          None,
        "new_count":          None,
        "removed_count":      None,
        "old_expected_atoms": None,
        "new_expected_atoms": None,
        "coordinate_atoms":   None,
        "synchronized":       False,
        "warnings":           [],
    }

    if not gro_path.exists():
        report["warnings"].append(f"coordinate file not found: {gro_path}")
        report_path.write_text(json.dumps(report, indent=2))
        return report

    if not topol_path.exists():
        report["warnings"].append(f"topology file not found: {topol_path}")
        report_path.write_text(json.dumps(report, indent=2))
        return report

    coordinate_atoms, matching_atom_count, new_count_by_resid = _parse_gro_atoms_and_molecules(
        gro_path, resname
    )
    if molecule_count_override is not None:
        if molecule_count_override < 0:
            raise ValueError("molecule_count_override must be non-negative")
        new_count = int(molecule_count_override)
    elif atoms_per_molecule:
        new_count, remainder = divmod(matching_atom_count, atoms_per_molecule)
        if remainder != 0:
            report["warnings"].append(
                f"{resname} atom count ({matching_atom_count}) not evenly divisible by "
                f"atoms_per_molecule ({atoms_per_molecule}); falling back to distinct-resid count"
            )
            new_count = new_count_by_resid
    else:
        new_count = new_count_by_resid

    report["coordinate_atoms"] = coordinate_atoms
    report["new_count"]        = new_count

    topol_text = topol_path.read_text()
    old_count  = _find_molecule_count(topol_text, resname)
    report["old_count"] = old_count

    if old_count is None:
        report["warnings"].append(
            f"'{resname}' entry not found in {topol_path} [ molecules ] — "
            "not inventing a new entry; topology left unchanged"
        )
        report_path.write_text(json.dumps(report, indent=2))
        return report

    atoms_per_molecule_est: Optional[int] = atoms_per_molecule
    if not atoms_per_molecule_est and new_count > 0:
        atoms_per_molecule_est, remainder = divmod(matching_atom_count, new_count)
        if remainder != 0:
            report["warnings"].append(
                f"{resname} atom count ({matching_atom_count}) not evenly divisible by "
                f"molecule count ({new_count}); atoms-per-molecule estimate may be inaccurate"
            )

    if atoms_per_molecule_est:
        base_atoms = coordinate_atoms - new_count * atoms_per_molecule_est
        report["new_expected_atoms"] = base_atoms + new_count * atoms_per_molecule_est
        report["old_expected_atoms"] = base_atoms + old_count * atoms_per_molecule_est
    else:
        report["new_expected_atoms"] = coordinate_atoms
        report["warnings"].append(
            f"no {resname} atoms found in coordinate file — cannot derive atoms-per-molecule; "
            "old_expected_atoms left unknown"
        )

    report["removed_count"] = old_count - new_count

    if old_count == new_count:
        report["synchronized"] = True
    else:
        new_text, updated = _update_molecule_count(topol_text, resname, new_count)
        if updated:
            topol_path.write_text(new_text)
            report["synchronized"] = True
        else:
            report["warnings"].append(
                f"failed to rewrite '{resname}' line in {topol_path}"
            )
            report["synchronized"] = False

    report_path.write_text(json.dumps(report, indent=2))
    return report


def sync_topology_lipid_count(
    gro_path:      Path | str,
    topol_path:    Path | str,
    lipid_resname: str,
    report_path:   Path | str | None = None,
) -> dict:
    """Lipid-specific wrapper around sync_topology_molecule_count(), kept for
    backward compatibility. See that function's docstring for behavior.
    Report keys (both the returned dict and the JSON written to disk) are
    translated back to this function's original names (old_dpp_count,
    new_dpp_count, removed_dpp_count, lipid_resname)."""
    resolved_report_path = (
        Path(report_path) if report_path is not None
        else Path(topol_path).parent / "topology_sync_report.json"
    )
    generic = sync_topology_molecule_count(gro_path, topol_path, lipid_resname, resolved_report_path)
    legacy = {
        "coordinate_file":    generic["coordinate_file"],
        "topology_file":      generic["topology_file"],
        "lipid_resname":      generic["resname"],
        "old_dpp_count":      generic["old_count"],
        "new_dpp_count":      generic["new_count"],
        "removed_dpp_count":  generic["removed_count"],
        "old_expected_atoms": generic["old_expected_atoms"],
        "new_expected_atoms": generic["new_expected_atoms"],
        "coordinate_atoms":   generic["coordinate_atoms"],
        "synchronized":       generic["synchronized"],
        "warnings":           generic["warnings"],
    }
    resolved_report_path.write_text(json.dumps(legacy, indent=2))
    return legacy
