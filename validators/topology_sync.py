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
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        return 0, 0, 0
    try:
        n = int(lines[1].strip())
    except ValueError:
        return 0, 0, 0

    lipid_atom_count = 0
    lipid_resids: set[int] = set()
    for line in lines[2: 2 + n]:
        if len(line) < 15:
            continue
        resname = line[5:10].strip()
        if resname != lipid_resname:
            continue
        lipid_atom_count += 1
        try:
            lipid_resids.add(int(line[0:5]))
        except ValueError:
            continue
    return n, lipid_atom_count, len(lipid_resids)


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


def sync_topology_lipid_count(
    gro_path:      Path | str,
    topol_path:    Path | str,
    lipid_resname: str,
    report_path:   Path | str | None = None,
) -> dict:
    """
    Synchronize the lipid molecule count in topol_path's [ molecules ]
    section with the actual lipid residue count in gro_path.

    - Counts DPP (or whichever lipid_resname) molecules directly from the
      coordinate file — never from a prior expected count or log text.
    - Rewrites only the matching lipid line; Protein/SOL/ion counts are left
      untouched.
    - If the lipid entry is absent from the topology, no entry is invented;
      the report is written with synchronized=False.

    Returns the report dict (also written to report_path, default
    "topology_sync_report.json" next to topol_path). Report keys:
        coordinate_file, topology_file, lipid_resname,
        old_dpp_count, new_dpp_count, removed_dpp_count,
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
        "lipid_resname":      lipid_resname,
        "old_dpp_count":      None,
        "new_dpp_count":      None,
        "removed_dpp_count":  None,
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

    coordinate_atoms, lipid_atom_count, new_dpp_count = _parse_gro_atoms_and_lipids(
        gro_path, lipid_resname
    )
    report["coordinate_atoms"] = coordinate_atoms
    report["new_dpp_count"]    = new_dpp_count

    topol_text    = topol_path.read_text()
    old_dpp_count = _find_molecule_count(topol_text, lipid_resname)
    report["old_dpp_count"] = old_dpp_count

    if old_dpp_count is None:
        report["warnings"].append(
            f"'{lipid_resname}' entry not found in {topol_path} [ molecules ] — "
            "not inventing a new entry; topology left unchanged"
        )
        report_path.write_text(json.dumps(report, indent=2))
        return report

    atoms_per_lipid: Optional[int] = None
    if new_dpp_count > 0:
        atoms_per_lipid, remainder = divmod(lipid_atom_count, new_dpp_count)
        if remainder != 0:
            report["warnings"].append(
                f"lipid atom count ({lipid_atom_count}) not evenly divisible by "
                f"lipid residue count ({new_dpp_count}); atoms-per-lipid estimate "
                "may be inaccurate"
            )

    if atoms_per_lipid:
        base_atoms = coordinate_atoms - new_dpp_count * atoms_per_lipid
        report["new_expected_atoms"] = base_atoms + new_dpp_count * atoms_per_lipid
        report["old_expected_atoms"] = base_atoms + old_dpp_count * atoms_per_lipid
    else:
        report["new_expected_atoms"] = coordinate_atoms
        report["warnings"].append(
            "no lipid atoms found in coordinate file — cannot derive atoms-per-lipid; "
            "old_expected_atoms left unknown"
        )

    report["removed_dpp_count"] = old_dpp_count - new_dpp_count

    if old_dpp_count == new_dpp_count:
        report["synchronized"] = True
    else:
        new_text, updated = _update_molecule_count(topol_text, lipid_resname, new_dpp_count)
        if updated:
            topol_path.write_text(new_text)
            report["synchronized"] = True
        else:
            report["warnings"].append(
                f"failed to rewrite '{lipid_resname}' line in {topol_path}"
            )
            report["synchronized"] = False

    report_path.write_text(json.dumps(report, indent=2))
    return report
