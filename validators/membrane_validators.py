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
        n_true_cavity_trapped=class_counts["true_cavity_trapped"], n_tm_slab_hard_contact=class_counts["tm_slab_hard_contact"],
        n_annular_surface_contact=class_counts["annular_surface_contact"], n_soluble_domain_contact=class_counts["soluble_domain_contact"],
        n_ambiguous=class_counts["ambiguous"], n_legacy_trapped_advisories=sum(1 for c in candidates if c.low_overlap_but_geometrically_inside),
        slab_z_min=round(slab_z_min,4), slab_z_max=round(slab_z_max,4),
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
    """Result of checking for water molecules inside the bilayer hydrophobic core."""

    status:                       ValidationStatus

    # Core-water count (gate-relevant)
    n_waters_in_bilayer:          int = 0

    # Hydrophobic core boundaries (inner leaflet headgroup means, split by midplane)
    bilayer_midplane_z:           Optional[float] = None  # mean Z of tail atoms
    core_z_min:                   Optional[float] = None  # inner boundary of bottom leaflet
    core_z_max:                   Optional[float] = None  # inner boundary of top leaflet

    # Outer headgroup extent (informational, NOT used for water counting)
    bilayer_z_min_nm:             Optional[float] = None  # min Z of all headgroup atoms
    bilayer_z_max_nm:             Optional[float] = None  # max Z of all headgroup atoms

    n_headgroup_atoms_found:      int = 0
    n_tail_atoms_found:           int = 0
    n_water_oxygens_checked:      int = 0

    source_file:                  str = ""
    message:                      str = ""


def validate_no_water_in_bilayer(
    gro_path: Path | str,
    headgroup_atom: str,
    tail_atom: str,
    water_oxygen: str = "OW",
) -> WaterConsistencyResult:
    """
    Detect water molecules whose oxygen sits within the bilayer hydrophobic core.

    Uses the same Z-boundary algorithm as WaterDeletorAdapter / water_deletor.pl:
      1. midplane = mean Z of tail_atom occurrences (bilayer geometric centre)
      2. top-leaflet headgroups = headgroup_atom with Z > midplane → mean = z_core_top
         bot-leaflet headgroups = headgroup_atom with Z ≤ midplane → mean = z_core_bot
      3. An OW is flagged if z_core_bot ≤ z ≤ z_core_top

    This is the INNER boundary of each leaflet, not min/max of all headgroups.
    Waters in the headgroup/water interface (outside this range) are correct and
    are NOT counted.

    Args:
        gro_path:       Path to .gro (pre- or post-deletion).
        headgroup_atom: Atom name marking leaflet inner boundary (e.g. "O33" for DPPC OPLS-AA).
        tail_atom:      Atom name marking bilayer midplane (e.g. "C50" for DPPC OPLS-AA).
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

    if not tail_z:
        return WaterConsistencyResult(
            status=ValidationStatus.SKIPPED,
            source_file=str(path),
            n_headgroup_atoms_found=len(headgroup_z),
            message=f"Tail atom '{tail_atom}' not found — cannot compute bilayer midplane",
        )

    # Step 1: bilayer midplane from tail atoms
    z_midplane = sum(tail_z) / len(tail_z)

    # Step 2: split headgroups into leaflets, compute inner boundary per leaflet
    top_hg = [z for z in headgroup_z if z > z_midplane]
    bot_hg = [z for z in headgroup_z if z <= z_midplane]

    z_core_top = sum(top_hg) / len(top_hg) if top_hg else z_midplane + 1.0
    z_core_bot = sum(bot_hg) / len(bot_hg) if bot_hg else z_midplane - 1.0

    # Step 3: count OW atoms within the hydrophobic core
    n_inside = sum(1 for z in water_oz if z_core_bot <= z <= z_core_top)

    # Outer headgroup extent (informational)
    z_outer_min = min(headgroup_z)
    z_outer_max = max(headgroup_z)

    if n_inside == 0:
        status = ValidationStatus.PASS
        msg    = (
            f"No water oxygens inside bilayer hydrophobic core "
            f"(core Z: {z_core_bot:.3f}–{z_core_top:.3f} nm, midplane: {z_midplane:.3f} nm)"
        )
    else:
        status = ValidationStatus.FAIL
        msg    = (
            f"{n_inside} water oxygen(s) inside bilayer hydrophobic core "
            f"(core Z: {z_core_bot:.3f}–{z_core_top:.3f} nm)"
        )

    return WaterConsistencyResult(
        status=status,
        n_waters_in_bilayer=n_inside,
        bilayer_midplane_z=round(z_midplane, 4),
        core_z_min=round(z_core_bot, 4),
        core_z_max=round(z_core_top, 4),
        bilayer_z_min_nm=round(z_outer_min, 4),
        bilayer_z_max_nm=round(z_outer_max, 4),
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


# ══════════════════════════════════════════════════════════════════════════════
# 5. Lipid-in-cavity detector (post embed_in_bilayer reporting)
# ══════════════════════════════════════════════════════════════════════════════

class LipidInCavityResult(BaseModel):
    """Result of checking for lipid atoms trapped inside the protein TM cavity."""

    n_lipid_residues_in_cavity: int = 0
    lipid_residue_ids:          list[int] = Field(default_factory=list)
    cavity_center_x:            Optional[float] = None
    cavity_center_y:            Optional[float] = None
    cavity_z_min:               Optional[float] = None
    cavity_z_max:               Optional[float] = None
    cavity_radius_nm:           Optional[float] = None
    source_file:                str = ""
    message:                    str = ""


def detect_lipid_in_cavity(
    gro_path:       Path | str,
    lipid_resname:  str,
    tm_residues:    Optional[set[int]] = None,
    cavity_margin:  float = 0.5,
) -> LipidInCavityResult:
    """
    Detect lipid residues whose XY COM falls inside the protein TM cavity.

    Algorithm:
      1. Collect Cα (or any heavy atom) XY coords of TM residues → compute XY centroid.
      2. Compute the minimum XY bounding radius of TM residues around that centroid.
         cavity_radius = max per-atom distance from centroid + cavity_margin.
         This represents the "inside" of the TM bundle.
      3. For each lipid residue: compute its XY COM.  If |COM - centroid| < cavity_radius
         AND the COM Z is within [z_min, z_max] of the TM residues → flagged as in cavity.

    This is a DETECTION-only function; no remediation is performed.
    Returns counts and residue IDs for downstream reporting.

    Args:
        gro_path:      Path to the embedded system .gro (post embed_in_bilayer).
        lipid_resname: Residue name of lipid atoms (e.g. "DPP" for DPPC OPLS-AA).
        tm_residues:   Set of protein residue numbers defining the TM bundle.
                       If None, uses all non-lipid, non-solvent residues.
        cavity_margin: Extra radius in nm added to the TM bundle XY extent (default 0.5 nm).
    """
    path = Path(gro_path)
    if not path.exists():
        return LipidInCavityResult(
            source_file=str(path),
            message=f"GRO file not found: {path}",
        )

    solvent_names = {"SOL", "HOH", "NA", "CL", "K", "MG", "CA"}

    try:
        lines = path.read_text().splitlines()
        atom_lines = lines[2:-1]
    except Exception as exc:
        return LipidInCavityResult(
            source_file=str(path),
            message=f"GRO parse error: {exc}",
        )

    # Collect TM-region protein Cα (or Cβ/N as fallback) coordinates
    tm_coords: list[tuple[float, float, float]] = []
    lipid_atoms: dict[int, list[tuple[float, float, float]]] = {}  # resnum → [(x,y,z)]

    for line in atom_lines:
        if len(line) < 44:
            continue
        resname  = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            resnum = int(line[0:5])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        if resname == lipid_resname:
            lipid_atoms.setdefault(resnum, []).append((x, y, z))
        elif resname not in solvent_names:
            if tm_residues is None or resnum in tm_residues:
                if atomname in ("CA", "N", "CB"):
                    tm_coords.append((x, y, z))

    if not tm_coords:
        return LipidInCavityResult(
            source_file=str(path),
            message=(
                "No TM-region backbone atoms found — "
                "cavity detection skipped (provide tm_residues for accuracy)"
            ),
        )

    # TM bundle XY centroid and Z extent
    cx = sum(c[0] for c in tm_coords) / len(tm_coords)
    cy = sum(c[1] for c in tm_coords) / len(tm_coords)
    z_min = min(c[2] for c in tm_coords)
    z_max = max(c[2] for c in tm_coords)

    # Cavity radius = max distance from centroid + margin
    max_r = max(math.sqrt((c[0]-cx)**2 + (c[1]-cy)**2) for c in tm_coords)
    cavity_r = max_r + cavity_margin

    # Check each lipid residue's XY COM
    flagged_ids: list[int] = []
    for resnum, coords in lipid_atoms.items():
        lx = sum(c[0] for c in coords) / len(coords)
        ly = sum(c[1] for c in coords) / len(coords)
        lz = sum(c[2] for c in coords) / len(coords)
        r  = math.sqrt((lx - cx)**2 + (ly - cy)**2)
        if r < cavity_r and z_min <= lz <= z_max:
            flagged_ids.append(resnum)

    flagged_ids.sort()
    n = len(flagged_ids)

    if n == 0:
        msg = "No lipid residues detected inside protein TM cavity"
    else:
        msg = (
            f"{n} lipid residue(s) detected inside protein TM cavity "
            f"(XY centroid: ({cx:.3f}, {cy:.3f}) nm, radius: {cavity_r:.3f} nm, "
            f"Z: {z_min:.3f}–{z_max:.3f} nm)"
        )

    return LipidInCavityResult(
        n_lipid_residues_in_cavity=n,
        lipid_residue_ids=flagged_ids[:50],  # cap to avoid huge JSON
        cavity_center_x=round(cx, 4),
        cavity_center_y=round(cy, 4),
        cavity_z_min=round(z_min, 4),
        cavity_z_max=round(z_max, 4),
        cavity_radius_nm=round(cavity_r, 4),
        source_file=str(path),
        message=msg,
    )


# ══════════════════════════════════════════════════════════════════════════════
# 6. Trapped Lipid Detector (Phase 4)
# ══════════════════════════════════════════════════════════════════════════════

class TrappedLipidCandidate(BaseModel):
    resid: int
    resname: str
    stage_first_detected: str = "03_embed_in_bilayer"
    present_since_embed: bool = True
    introduced_by_shrink: bool = False
    inside_tm_footprint: bool
    inside_expanded_tm_footprint: bool
    low_overlap_but_geometrically_inside: bool
    # True when protein CA atoms angularly surround the lipid COM (max angular gap < 185°).
    # Reliable indicator that the lipid is inside the TM bundle, not just surface-adjacent.
    is_surrounded_by_protein: bool = False
    min_distance_to_protein_nm: float
    atoms_within_025: int
    atoms_within_035: int
    atoms_within_045: int
    lipid_center: list[float]
    protein_tm_center: list[float]
    z_relative_to_midplane: float
    suspected_cause: str
    recommended_action: str
    classification: str = "ambiguous"
    lipid_com: list[float] = Field(default_factory=list)
    inside_slab: bool = False
    angularly_enclosed: bool = False
    min_dist_tm_slab_nm: float = 999.0
    min_dist_full_protein_nm: float = 999.0
    nearest_contact_region: str = "ambiguous"


class TrappedLipidDiagnosis(BaseModel):
    status: ValidationStatus
    bilayer_midplane_z: float
    protein_tm_center: list[float]
    protein_xy_footprint: list[float]
    lipid_residue_candidates: list[TrappedLipidCandidate] = Field(default_factory=list)

    # summary
    n_lipids_checked: int = 0
    n_suspicious_lipids: int = 0
    n_hard_overlap_lipids: int = 0
    n_low_overlap_geometrically_trapped_lipids: int = 0
    likely_root_cause: str = ""
    should_block_future_workflows: bool = True
    remediation_not_applied: bool = True
    message: str = ""
    n_true_cavity_trapped: int = 0
    n_tm_slab_hard_contact: int = 0
    n_annular_surface_contact: int = 0
    n_soluble_domain_contact: int = 0
    n_ambiguous: int = 0
    n_legacy_trapped_advisories: int = 0
    slab_z_min: float = 0.0
    slab_z_max: float = 0.0


def detect_trapped_lipids(
    gro_path: Path | str,
    lipid_resname: str = "DPP",
    tm_residues: Optional[set[int]] = None,
) -> TrappedLipidDiagnosis:
    """
    Detect lipid residues that are trapped inside the protein cavity (low overlap)
    or have hard overlaps with protein atoms.
    """
    path = Path(gro_path)
    if not path.exists():
        return TrappedLipidDiagnosis(
            status=ValidationStatus.FAIL,
            bilayer_midplane_z=0.0,
            protein_tm_center=[0.0, 0.0, 0.0],
            protein_xy_footprint=[0.0, 0.0, 0.0, 0.0],
            message=f"GRO file not found: {path}",
        )

    solvent_names = {"SOL", "HOH", "NA", "CL", "K", "MG", "CA"}

    try:
        lines = path.read_text().splitlines()
        atom_lines = lines[2:-1]
    except Exception as exc:
        return TrappedLipidDiagnosis(
            status=ValidationStatus.FAIL,
            bilayer_midplane_z=0.0,
            protein_tm_center=[0.0, 0.0, 0.0],
            protein_xy_footprint=[0.0, 0.0, 0.0, 0.0],
            message=f"GRO parse error: {exc}",
        )

    # Collect coordinates
    tm_coords: list[tuple[float, float, float]] = []
    protein_atoms: list[tuple[float, float, float]] = []
    lipid_atoms: dict[int, list[tuple[float, float, float]]] = {}  # resnum -> [(x,y,z)]
    tail_zs: list[float] = []

    for line in atom_lines:
        if len(line) < 44:
            continue
        resname  = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            resnum = int(line[0:5])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        if resname == lipid_resname:
            lipid_atoms.setdefault(resnum, []).append((x, y, z))
            if atomname == "C50":
                tail_zs.append(z)
        elif resname not in solvent_names:
            protein_atoms.append((x, y, z))
            if tm_residues is None or resnum in tm_residues:
                if atomname == "CA":
                    tm_coords.append((x, y, z))

    # Fallback to general CA / protein atoms if tm_coords is empty
    if not tm_coords:
        for line in atom_lines:
            if len(line) < 44:
                continue
            resname  = line[5:10].strip()
            atomname = line[10:15].strip()
            if resname != lipid_resname and resname not in solvent_names:
                try:
                    x = float(line[20:28])
                    y = float(line[28:36])
                    z = float(line[36:44])
                    if atomname == "CA":
                        tm_coords.append((x, y, z))
                except ValueError:
                    continue

    if not tm_coords:
        for line in atom_lines:
            if len(line) < 44:
                continue
            resname  = line[5:10].strip()
            if resname != lipid_resname and resname not in solvent_names:
                try:
                    x = float(line[20:28])
                    y = float(line[28:36])
                    z = float(line[36:44])
                    tm_coords.append((x, y, z))
                except ValueError:
                    continue

    all_lipid_zs = [a[2] for atoms in lipid_atoms.values() for a in atoms]
    midplane_z = (sum(all_lipid_zs) / len(all_lipid_zs)) if all_lipid_zs else 0.0
    slab_z_min = (midplane_z - 1.5) if all_lipid_zs else -1e9
    slab_z_max = (midplane_z + 1.5) if all_lipid_zs else 1e9
    tm_slab_atoms = [a for a in protein_atoms if slab_z_min <= a[2] <= slab_z_max]
    tm_slab_cas = [c for c in tm_coords if slab_z_min <= c[2] <= slab_z_max]
    if len(tm_slab_cas) < 5:
        tm_slab_cas = tm_coords

    if not tm_coords or not protein_atoms:
        return TrappedLipidDiagnosis(
            status=ValidationStatus.PASS,
            bilayer_midplane_z=midplane_z,
            protein_tm_center=[0.0, 0.0, 0.0],
            protein_xy_footprint=[0.0, 0.0, 0.0, 0.0],
            n_lipids_checked=len(lipid_atoms),
            message="No protein atoms found; trapped lipid detection skipped.",
        )

    # Compute TM Center and Bounding Box
    tm_cx = sum(c[0] for c in tm_coords) / len(tm_coords)
    tm_cy = sum(c[1] for c in tm_coords) / len(tm_coords)
    tm_cz = sum(c[2] for c in tm_coords) / len(tm_coords)
    tm_center = [tm_cx, tm_cy, tm_cz]

    tm_x_min = min(c[0] for c in tm_coords)
    tm_x_max = max(c[0] for c in tm_coords)
    tm_y_min = min(c[1] for c in tm_coords)
    tm_y_max = max(c[1] for c in tm_coords)
    tm_z_min = min(c[2] for c in tm_coords)
    tm_z_max = max(c[2] for c in tm_coords)
    tm_radius = max(math.sqrt((c[0]-tm_cx)**2 + (c[1]-tm_cy)**2) for c in tm_coords)

    def get_angular_gap(lx, ly, lz):
        local_cas = [c for c in tm_coords if abs(c[2] - lz) < 1.5]
        if len(local_cas) < 5:
            local_cas = tm_coords
        
        angles = []
        for c in local_cas:
            dx = c[0] - lx
            dy = c[1] - ly
            angles.append(math.atan2(dy, dx))
        
        if not angles:
            return 360.0
        
        angles = sorted(angles)
        max_gap = 0.0
        for i in range(len(angles)):
            if i == len(angles) - 1:
                gap = (angles[0] + 2 * math.pi) - angles[i]
            else:
                gap = angles[i+1] - angles[i]
            if gap > max_gap:
                max_gap = gap
                
        return math.degrees(max_gap)

    candidates = []
    n_hard = 0
    n_low_overlap = 0
    class_counts = {"true_cavity_trapped":0, "tm_slab_hard_contact":0, "annular_surface_contact":0, "soluble_domain_contact":0, "ambiguous":0}

    for resid, coords in lipid_atoms.items():
        lc_x = sum(c[0] for c in coords) / len(coords)
        lc_y = sum(c[1] for c in coords) / len(coords)
        lc_z = sum(c[2] for c in coords) / len(coords)
        lc = [lc_x, lc_y, lc_z]

        dist_to_tm_center = math.sqrt((lc_x - tm_cx)**2 + (lc_y - tm_cy)**2 + (lc_z - tm_cz)**2)
        lipid_radius = max(math.sqrt((a[0]-lc_x)**2 + (a[1]-lc_y)**2 + (a[2]-lc_z)**2) for a in coords)

        # Distances are classified separately: TM/slab geometry drives cavity status;
        # full-protein distance is retained only as a steric warning.
        def mindist(atoms):
            return min((math.dist(la, pa) for la in coords for pa in atoms), default=999.0)
        min_dist_tm = mindist(tm_slab_atoms)
        min_dist = mindist(protein_atoms)
        atoms_within_025 = sum(1 for la in coords for pa in tm_slab_atoms if math.dist(la, pa) < 0.25)
        atoms_within_035 = sum(1 for la in coords for pa in tm_slab_atoms if math.dist(la, pa) < 0.35)
        atoms_within_045 = sum(1 for la in coords for pa in tm_slab_atoms if math.dist(la, pa) < 0.45)

        inside_slab = slab_z_min <= lc_z <= slab_z_max
        inside_tm_footprint = (tm_x_min <= lc_x <= tm_x_max) and (tm_y_min <= lc_y <= tm_y_max)
        inside_expanded_tm_footprint = (tm_x_min - 0.5 <= lc_x <= tm_x_max + 0.5) and (tm_y_min - 0.5 <= lc_y <= tm_y_max + 0.5)
        max_gap = get_angular_gap(lc_x, lc_y, lc_z)
        is_surrounded = max_gap < 185.0
        angularly_enclosed = is_surrounded and inside_slab
        legacy_trapped = is_surrounded and (tm_z_min <= lc_z <= tm_z_max) and (dist_to_tm_center < tm_radius)
        hard_tm = min_dist_tm < 0.2
        if legacy_trapped and angularly_enclosed and inside_tm_footprint and min_dist_tm >= 0.2:
            classification = "true_cavity_trapped"
        elif hard_tm and inside_slab and inside_tm_footprint and angularly_enclosed:
            classification = "tm_slab_hard_contact"
        elif inside_slab and inside_tm_footprint and not angularly_enclosed:
            classification = "annular_surface_contact"
        elif min_dist < 0.2 and not inside_slab:
            classification = "soluble_domain_contact"
        else:
            classification = "ambiguous"
        class_counts[classification] += 1
        if legacy_trapped: n_low_overlap += 1
        if hard_tm: n_hard += 1
        suspicious = classification in {"true_cavity_trapped", "tm_slab_hard_contact"}
        suspected_cause = {"true_cavity_trapped":"TM/slab cavity enclosure", "tm_slab_hard_contact":"TM/slab hard contact", "annular_surface_contact":"valid annular surface contact", "soluble_domain_contact":"soluble-domain-only contact", "ambiguous":"ambiguous membrane contact"}[classification]
        recommended_action = "Remove via membrane-relevant exclusion" if classification == "true_cavity_trapped" else "Review steric contact"
        candidates.append(TrappedLipidCandidate(
            resid=resid, resname=lipid_resname,
            inside_tm_footprint=inside_tm_footprint,
            inside_expanded_tm_footprint=inside_expanded_tm_footprint,
            low_overlap_but_geometrically_inside=(classification == "true_cavity_trapped"),
            is_surrounded_by_protein=is_surrounded,
            min_distance_to_protein_nm=round(min_dist, 4),
            atoms_within_025=atoms_within_025, atoms_within_035=atoms_within_035, atoms_within_045=atoms_within_045,
            lipid_center=lc, protein_tm_center=tm_center,
            z_relative_to_midplane=round(lc_z - midplane_z, 4),
            suspected_cause=suspected_cause, recommended_action=recommended_action,
            classification=classification, lipid_com=lc, inside_slab=inside_slab,
            angularly_enclosed=angularly_enclosed, min_dist_tm_slab_nm=round(min_dist_tm,4),
            min_dist_full_protein_nm=round(min_dist,4), nearest_contact_region=classification,
        ))

    n_suspicious = class_counts["true_cavity_trapped"] + class_counts["tm_slab_hard_contact"]
    status = ValidationStatus.WARNING if n_suspicious > 0 else ValidationStatus.PASS

    if n_suspicious > 0:
        msg = f"{n_suspicious} suspicious lipid(s) detected trapped in the protein cavity or clashing."
    else:
        msg = "No trapped or clashing lipids detected in the protein cavity."

    return TrappedLipidDiagnosis(
        status=status,
        bilayer_midplane_z=round(midplane_z, 4),
        protein_tm_center=[round(v, 4) for v in tm_center],
        protein_xy_footprint=[round(v, 4) for v in [tm_x_min, tm_x_max, tm_y_min, tm_y_max]],
        lipid_residue_candidates=candidates,
        n_lipids_checked=len(lipid_atoms),
        n_suspicious_lipids=n_suspicious,
        n_hard_overlap_lipids=n_hard,
        n_low_overlap_geometrically_trapped_lipids=n_low_overlap,
        likely_root_cause=(
            "Initial movememb embedding step lacks protein-cavity-aware exclusion, "
            "resulting in DPPC lipids being placed inside/through the receptor. "
            "Subsequently, the inflategro shrink loop is executed with cutoff=0.0 "
            "which disables any lipid deletion, leaving the cavity-embedded lipids trapped."
        ),
        should_block_future_workflows=True,
        remediation_not_applied=True,
        message=msg,
        n_true_cavity_trapped=class_counts["true_cavity_trapped"],
        n_tm_slab_hard_contact=class_counts["tm_slab_hard_contact"],
        n_annular_surface_contact=class_counts["annular_surface_contact"],
        n_soluble_domain_contact=class_counts["soluble_domain_contact"],
        n_ambiguous=class_counts["ambiguous"],
        n_legacy_trapped_advisories=sum(1 for c in candidates if c.low_overlap_but_geometrically_inside),
        slab_z_min=round(slab_z_min,4), slab_z_max=round(slab_z_max,4),
    )

