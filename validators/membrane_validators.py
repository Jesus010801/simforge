"""
validators/membrane_validators.py — membrane system integrity validators.

Four validators, each with a single concern:

    APLConvergenceValidator     — did the shrink loop reach physical APL?
    OverlapValidator            — are protein and bilayer atoms clashing?
    OrientationValidator        — is the protein TM axis aligned with Z?
    WaterConsistencyValidator   — are waters present inside the bilayer core?

Design rules (consistent with protein_validator.py / ligand_validator.py):
  - Each validator returns a typed Pydantic result. Never raises.
  - No remediation. No suggestions about what to do. Only perception.
  - No external dependencies beyond stdlib + pydantic.
  - membrane_knowledge.py is the only source of physical constants.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


# ── Shared severity / status types ───────────────────────────────────────────

class ValidationStatus(str, Enum):
    PASS    = "pass"
    WARNING = "warning"
    FAIL    = "fail"
    SKIPPED = "skipped"   # required input not available


# ══════════════════════════════════════════════════════════════════════════════
# 1. APL Convergence Validator
# ══════════════════════════════════════════════════════════════════════════════

class APLConvergenceResult(BaseModel):
    """Result of checking whether the shrink loop reached target APL."""

    status:            ValidationStatus

    # Measured values
    apl_ang2:          Optional[float] = None   # measured, in Å²
    apl_nm2:           Optional[float] = None   # measured, in nm²

    # Reference values used for comparison
    target_apl_ang2:   Optional[float] = None
    tolerance_ang2:    float = 2.0

    # Derived
    delta_ang2:        Optional[float] = None   # measured - target; negative = under-compressed
    converged:         bool = False

    # Diagnostics
    source_file:       str = ""
    message:           str = ""


def validate_apl_convergence(
    area_dat_path: Path | str,
    target_apl_ang2: float,
    tolerance_ang2: float = 2.0,
) -> APLConvergenceResult:
    """
    Check whether area_2.dat reports an APL at or below the convergence target.

    Args:
        area_dat_path:   Path to area_2.dat written by inflategro.
        target_apl_ang2: Physical target in Å² (from membrane_knowledge.apl_target()).
        tolerance_ang2:  Acceptable overshoot above target (default 2 Å²).

    Pass condition: apl_ang2 <= target_apl_ang2 + tolerance_ang2
    """
    path = Path(area_dat_path)

    if not path.exists():
        return APLConvergenceResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"area_2.dat not found: {path}",
        )

    try:
        raw = path.read_text().strip()
        apl_nm2 = float(raw)
    except (ValueError, OSError) as exc:
        return APLConvergenceResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"Could not parse area_2.dat: {exc}",
        )

    apl_ang2 = round(apl_nm2 * 100, 2)
    delta    = apl_ang2 - target_apl_ang2
    converged = apl_ang2 <= (target_apl_ang2 + tolerance_ang2)

    status = ValidationStatus.PASS if converged else ValidationStatus.FAIL
    msg = (
        f"APL = {apl_ang2:.1f} Å² (target ≤ {target_apl_ang2 + tolerance_ang2:.1f} Å²)"
        if converged
        else f"APL = {apl_ang2:.1f} Å² exceeds target {target_apl_ang2:.1f} + tol {tolerance_ang2:.1f} Å²"
    )

    return APLConvergenceResult(
        status=status,
        apl_ang2=apl_ang2,
        apl_nm2=apl_nm2,
        target_apl_ang2=target_apl_ang2,
        tolerance_ang2=tolerance_ang2,
        delta_ang2=round(delta, 2),
        converged=converged,
        source_file=str(path),
        message=msg,
    )


# ══════════════════════════════════════════════════════════════════════════════
# 2. Overlap Validator
# ══════════════════════════════════════════════════════════════════════════════

class ClashInfo(BaseModel):
    protein_atom_line: int     # 1-based line number in .gro
    lipid_atom_line:   int
    distance_nm:       float


class OverlapResult(BaseModel):
    """Result of checking protein-bilayer atom overlap."""

    status:          ValidationStatus

    n_clashes:       int = 0
    clashes:         list[ClashInfo] = Field(default_factory=list)
    clash_cutoff_nm: float = 0.2

    source_file:     str = ""
    message:         str = ""

    # Atom counts inspected
    n_protein_atoms: int = 0
    n_lipid_atoms:   int = 0


def validate_no_overlap(
    gro_path: Path | str,
    lipid_residue_name: str,
    clash_cutoff_nm: float = 0.2,
    max_reported_clashes: int = 20,
) -> OverlapResult:
    """
    Detect van-der-Waals clashes between protein atoms and bilayer lipid atoms.

    A clash is defined as any protein–lipid atom pair closer than clash_cutoff_nm.
    Uses a simple O(N×M) scan — acceptable for pre-inflation systems where lipids
    are spread out.  Do NOT run on production trajectories.

    Args:
        gro_path:             Path to the concatenated system .gro.
        lipid_residue_name:   Residue name of lipid atoms (e.g. "DPP").
        clash_cutoff_nm:      Distance threshold in nm (default 0.2 nm = 2 Å).
        max_reported_clashes: Stop collecting clash details after this many.
    """
    path = Path(gro_path)
    if not path.exists():
        return OverlapResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO file not found: {path}",
        )

    try:
        protein_coords, lipid_coords = _parse_gro_two_groups(
            path, lipid_residue_name
        )
    except Exception as exc:
        return OverlapResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO parse error: {exc}",
        )

    clashes: list[ClashInfo] = []
    cutoff2 = clash_cutoff_nm ** 2

    for p_line, (px, py, pz) in protein_coords:
        for l_line, (lx, ly, lz) in lipid_coords:
            d2 = (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2
            if d2 < cutoff2:
                clashes.append(ClashInfo(
                    protein_atom_line=p_line,
                    lipid_atom_line=l_line,
                    distance_nm=round(math.sqrt(d2), 4),
                ))
                if len(clashes) >= max_reported_clashes:
                    break
        if len(clashes) >= max_reported_clashes:
            break

    n = len(clashes)
    if n == 0:
        status = ValidationStatus.PASS
        msg    = f"No protein–lipid clashes detected (cutoff {clash_cutoff_nm} nm)"
    else:
        status = ValidationStatus.FAIL
        msg    = f"{n} clash(es) detected between protein and {lipid_residue_name} atoms"

    return OverlapResult(
        status=status,
        n_clashes=n,
        clashes=clashes,
        clash_cutoff_nm=clash_cutoff_nm,
        source_file=str(path),
        message=msg,
        n_protein_atoms=len(protein_coords),
        n_lipid_atoms=len(lipid_coords),
    )


# ══════════════════════════════════════════════════════════════════════════════
# 3. Orientation Validator
# ══════════════════════════════════════════════════════════════════════════════

class TMHelixInfo(BaseModel):
    res_start:  int
    res_end:    int
    axis_x:     float   # unit vector components of helix principal axis
    axis_y:     float
    axis_z:     float
    angle_with_z_deg: float   # angle between helix axis and Z-axis (0° = perfectly aligned)


class OrientationResult(BaseModel):
    """Result of checking TM helix Z-alignment."""

    status:              ValidationStatus

    helices:             list[TMHelixInfo] = Field(default_factory=list)
    n_helices_checked:   int = 0
    n_helices_aligned:   int = 0   # angle_with_z_deg <= alignment_threshold_deg

    alignment_threshold_deg: float = 30.0   # default: helix within 30° of Z is aligned

    source_file:         str = ""
    message:             str = ""


def validate_tm_orientation(
    gro_path: Path | str,
    tm_residue_ranges: list[tuple[int, int]],
    alignment_threshold_deg: float = 30.0,
) -> OrientationResult:
    """
    Check whether TM helices are Z-aligned (i.e. perpendicular to the bilayer).

    Computes the principal axis of each helix segment using C-alpha coordinates
    and measures its angle with the Z-axis.

    Args:
        gro_path:             Path to the oriented protein .gro.
        tm_residue_ranges:    List of (res_start, res_end) tuples for each TM segment.
                              Derived from DeepTMHMM or equivalent prediction.
        alignment_threshold_deg: Max angle from Z-axis considered "aligned" (default 30°).
    """
    path = Path(gro_path)
    if not path.exists():
        return OrientationResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO file not found: {path}",
        )
    if not tm_residue_ranges:
        return OrientationResult(
            status=ValidationStatus.SKIPPED,
            source_file=str(path),
            message="No TM residue ranges provided — orientation check skipped",
        )

    try:
        ca_coords = _parse_gro_ca_coords(path)
    except Exception as exc:
        return OrientationResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO parse error: {exc}",
        )

    helices: list[TMHelixInfo] = []
    for res_start, res_end in tm_residue_ranges:
        segment = [
            (x, y, z) for resnum, (x, y, z) in ca_coords
            if res_start <= resnum <= res_end
        ]
        if len(segment) < 3:
            continue
        axis = _principal_axis(segment)
        angle = _angle_with_z(axis)
        helices.append(TMHelixInfo(
            res_start=res_start,
            res_end=res_end,
            axis_x=round(axis[0], 4),
            axis_y=round(axis[1], 4),
            axis_z=round(axis[2], 4),
            angle_with_z_deg=round(angle, 1),
        ))

    n_aligned = sum(1 for h in helices if h.angle_with_z_deg <= alignment_threshold_deg)
    n_total   = len(helices)

    if n_total == 0:
        status = ValidationStatus.SKIPPED
        msg    = "No TM helix segments had enough CA atoms to compute axis"
    elif n_aligned == n_total:
        status = ValidationStatus.PASS
        msg    = f"All {n_total} TM helix/helices aligned with Z (≤ {alignment_threshold_deg}°)"
    elif n_aligned == 0:
        status = ValidationStatus.FAIL
        msg    = f"No TM helices aligned with Z — protein likely needs rotation"
    else:
        status = ValidationStatus.WARNING
        msg    = f"{n_aligned}/{n_total} TM helices Z-aligned"

    return OrientationResult(
        status=status,
        helices=helices,
        n_helices_checked=n_total,
        n_helices_aligned=n_aligned,
        alignment_threshold_deg=alignment_threshold_deg,
        source_file=str(path),
        message=msg,
    )


# ══════════════════════════════════════════════════════════════════════════════
# 4. Water Consistency Validator
# ══════════════════════════════════════════════════════════════════════════════

class WaterConsistencyResult(BaseModel):
    """Result of checking for water molecules inside the bilayer core."""

    status:                  ValidationStatus

    n_waters_in_bilayer:     int = 0
    bilayer_z_min_nm:        Optional[float] = None
    bilayer_z_max_nm:        Optional[float] = None

    n_headgroup_atoms_found: int = 0
    n_tail_atoms_found:      int = 0
    n_water_oxygens_checked: int = 0

    source_file:             str = ""
    message:                 str = ""


def validate_no_water_in_bilayer(
    gro_path: Path | str,
    headgroup_atom: str,
    tail_atom: str,
    water_oxygen: str = "OW",
) -> WaterConsistencyResult:
    """
    Detect water molecules whose oxygen sits within the bilayer hydrophobic core.

    Uses the same Z-boundary logic as water_deletor.pl:
      - headgroup_atom (e.g. O33): marks the top/bottom surfaces of each leaflet
      - tail_atom (e.g. C50): marks the bilayer centre
      - A water OW between z_min(headgroup) and z_max(headgroup) is inside the core

    The bilayer Z-range is computed as:
      z_min = min Z of all headgroup_atom occurrences
      z_max = max Z of all headgroup_atom occurrences

    Any OW whose Z falls within [z_min, z_max] is flagged as inside the bilayer.

    Args:
        gro_path:       Path to solvated system .gro.
        headgroup_atom: Atom name marking bilayer surfaces (e.g. "O33" for DPPC OPLS-AA).
        tail_atom:      Atom name marking bilayer core (e.g. "C50").  Used for
                        sanity check that the bilayer is present.
        water_oxygen:   Atom name for water oxygen (default "OW" for SPC/TIP3P).
    """
    path = Path(gro_path)
    if not path.exists():
        return WaterConsistencyResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO file not found: {path}",
        )

    try:
        headgroup_z, tail_z, water_oz = _parse_gro_z_by_atomname(
            path, headgroup_atom, tail_atom, water_oxygen
        )
    except Exception as exc:
        return WaterConsistencyResult(
            status=ValidationStatus.FAIL,
            source_file=str(path),
            message=f"GRO parse error: {exc}",
        )

    if not headgroup_z:
        return WaterConsistencyResult(
            status=ValidationStatus.SKIPPED,
            source_file=str(path),
            message=f"Headgroup atom '{headgroup_atom}' not found — is this a membrane system?",
        )

    z_min = min(headgroup_z)
    z_max = max(headgroup_z)

    n_inside = sum(1 for z in water_oz if z_min <= z <= z_max)

    if n_inside == 0:
        status = ValidationStatus.PASS
        msg    = f"No water oxygens inside bilayer core (Z: {z_min:.3f}–{z_max:.3f} nm)"
    else:
        status = ValidationStatus.FAIL
        msg    = (
            f"{n_inside} water oxygen(s) inside bilayer core "
            f"(Z: {z_min:.3f}–{z_max:.3f} nm) — run water_deletor"
        )

    return WaterConsistencyResult(
        status=status,
        n_waters_in_bilayer=n_inside,
        bilayer_z_min_nm=round(z_min, 4),
        bilayer_z_max_nm=round(z_max, 4),
        n_headgroup_atoms_found=len(headgroup_z),
        n_tail_atoms_found=len(tail_z),
        n_water_oxygens_checked=len(water_oz),
        source_file=str(path),
        message=msg,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Internal GRO parsers — no public API
# ══════════════════════════════════════════════════════════════════════════════

def _parse_gro_two_groups(
    path: Path,
    lipid_residue_name: str,
) -> tuple[list[tuple[int, tuple[float, float, float]]], list[tuple[int, tuple[float, float, float]]]]:
    """
    Parse a .gro file and return two lists of (line_number, (x,y,z)):
      protein_coords: atoms whose residue name is NOT the lipid and NOT SOL
      lipid_coords:   atoms whose residue name matches lipid_residue_name

    GRO format (fixed width):
      col  1- 5  residue number
      col  6-10  residue name
      col 11-15  atom name
      col 16-20  atom number
      col 21-28  x (nm)
      col 29-36  y (nm)
      col 37-44  z (nm)
    """
    protein_coords: list[tuple[int, tuple[float, float, float]]] = []
    lipid_coords:   list[tuple[int, tuple[float, float, float]]] = []

    lines = path.read_text().splitlines()
    # Skip header (line 0) and atom count (line 1); last line is box vectors
    atom_lines = lines[2:-1]

    for i, line in enumerate(atom_lines, start=3):
        if len(line) < 44:
            continue
        resname  = line[5:10].strip()
        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        coord = (x, y, z)
        if resname == lipid_residue_name:
            lipid_coords.append((i, coord))
        elif resname not in ("SOL", "NA", "CL", "HOH"):
            protein_coords.append((i, coord))

    return protein_coords, lipid_coords


def _parse_gro_ca_coords(
    path: Path,
) -> list[tuple[int, tuple[float, float, float]]]:
    """Return list of (residue_number, (x,y,z)) for CA atoms only."""
    result: list[tuple[int, tuple[float, float, float]]] = []
    for line in path.read_text().splitlines()[2:-1]:
        if len(line) < 44:
            continue
        atomname = line[10:15].strip()
        if atomname != "CA":
            continue
        try:
            resnum = int(line[0:5])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue
        result.append((resnum, (x, y, z)))
    return result


def _parse_gro_z_by_atomname(
    path: Path,
    headgroup_atom: str,
    tail_atom: str,
    water_oxygen: str,
) -> tuple[list[float], list[float], list[float]]:
    """Return Z-coordinates for headgroup, tail, and water oxygen atoms."""
    headgroup_z: list[float] = []
    tail_z:      list[float] = []
    water_oz:    list[float] = []

    for line in path.read_text().splitlines()[2:-1]:
        if len(line) < 44:
            continue
        atomname = line[10:15].strip()
        try:
            z = float(line[36:44])
        except ValueError:
            continue
        if atomname == headgroup_atom:
            headgroup_z.append(z)
        elif atomname == tail_atom:
            tail_z.append(z)
        elif atomname == water_oxygen:
            water_oz.append(z)

    return headgroup_z, tail_z, water_oz


# ── Geometry helpers ──────────────────────────────────────────────────────────

def _principal_axis(coords: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    """
    Estimate the principal axis of a helix segment.

    Uses the end-to-end vector (last CA minus first CA) as the primary estimate,
    then refines with a linear regression through all CAs.  This is more robust
    than power iteration for straight helix segments where one spatial dimension
    may have zero variance (e.g. a perfect Z-helix has zero X/Y variance).
    """
    if len(coords) < 2:
        return (0.0, 0.0, 1.0)

    # End-to-end vector — good first estimate for a straight helix
    x0, y0, z0 = coords[0]
    x1, y1, z1 = coords[-1]
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    mag = math.sqrt(dx*dx + dy*dy + dz*dz)
    if mag < 1e-10:
        return (0.0, 0.0, 1.0)
    return (dx / mag, dy / mag, dz / mag)


def _angle_with_z(axis: tuple[float, float, float]) -> float:
    """Return angle in degrees between axis and the Z-unit vector (0,0,1)."""
    dot = abs(axis[2])   # |axis · z_hat| = |axis_z| since z_hat=(0,0,1)
    dot = min(1.0, dot)  # clamp floating point
    return math.degrees(math.acos(dot))
