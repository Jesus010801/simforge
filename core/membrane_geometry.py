# core/membrane_geometry.py
"""
Membrane orientation geometry for TM protein alignment.

Given a GRO file produced by gmx editconf -princ (principal axes aligned)
and semantic annotations (which residues are EC vs IC), this module
computes the gmx editconf -rotate angles needed to place the EC face at +Z.

Rotation table (after -princ, axes aligned to box):
    TM dominant axis   EC direction   editconf -rotate args
    X                  +X             0  270  0
    X                  -X             0   90  0
    Y                  +Y           -90    0  0
    Y                  -Y            90    0  0
    Z                  +Z             0    0  0   (no rotation needed)
    Z                  -Z           180    0  0
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple


# ── Public types ──────────────────────────────────────────────────────────────

class OrientRotation(NamedTuple):
    rx:          float   # degrees around X
    ry:          float   # degrees around Y
    rz:          float   # degrees around Z
    description: str


# ── Residue range parser ──────────────────────────────────────────────────────

def parse_residue_range(s: str) -> list[int]:
    """
    Parse a residue range string into a sorted list of residue numbers.

    Accepts:
        "1-50"           → [1, 2, ..., 50]
        "1-20,45-60"     → [1..20] + [45..60]
        "5,10,15"        → [5, 10, 15]
        "1-50,75,90-100" → mixed forms
    """
    if not s or not s.strip():
        return []

    residues: list[int] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            residues.extend(range(int(lo.strip()), int(hi.strip()) + 1))
        else:
            residues.append(int(part))

    return sorted(set(residues))


# ── GRO reader ────────────────────────────────────────────────────────────────

def _read_gro_ca_positions(
    gro_path: Path,
    residue_ids: set[int],
) -> list[tuple[float, float, float]]:
    """
    Return Cα positions (nm) for the given residue numbers from a GRO file.

    GRO atom line layout (fixed-width):
        cols  0-4   : residue number (5 chars)
        cols  5-9   : residue name (5 chars)
        cols 10-14  : atom name    (5 chars)
        cols 15-19  : atom number  (5 chars)
        cols 20-27  : x (nm)
        cols 28-35  : y (nm)
        cols 36-43  : z (nm)
    """
    positions: list[tuple[float, float, float]] = []
    lines = Path(gro_path).read_text().splitlines()

    if len(lines) < 3:
        return positions

    try:
        n_atoms = int(lines[1].strip())
    except ValueError:
        return positions

    for line in lines[2:2 + n_atoms]:
        if len(line) < 44:
            continue
        try:
            resnum   = int(line[0:5])
            atomname = line[10:15].strip()
        except ValueError:
            continue

        if resnum not in residue_ids or atomname != "CA":
            continue

        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        positions.append((x, y, z))

    return positions


def _mean_pos(
    positions: list[tuple[float, float, float]],
) -> tuple[float, float, float]:
    n = len(positions)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (
        sum(p[0] for p in positions) / n,
        sum(p[1] for p in positions) / n,
        sum(p[2] for p in positions) / n,
    )


# ── Core rotation computation ─────────────────────────────────────────────────

def compute_orient_rotation(
    gro_path:  Path,
    ec_resids: list[int],
    ic_resids: list[int],
) -> OrientRotation:
    """
    Compute editconf rotation angles to orient a TM protein.

    Reads Cα positions for EC and IC residue groups from gro_path
    (assumed to be the output of editconf -princ), computes their
    centres of mass, finds the dominant TM axis, and returns the
    rotation that places the EC face at +Z.

    Raises ValueError if Cα atoms cannot be found for either group.
    """
    ec_pos = _read_gro_ca_positions(gro_path, set(ec_resids))
    ic_pos = _read_gro_ca_positions(gro_path, set(ic_resids))

    if not ec_pos:
        raise ValueError(
            f"No Cα atoms found for EC residues {ec_resids[:5]}{'...' if len(ec_resids) > 5 else ''} "
            f"in {gro_path}"
        )
    if not ic_pos:
        raise ValueError(
            f"No Cα atoms found for IC residues {ic_resids[:5]}{'...' if len(ic_resids) > 5 else ''} "
            f"in {gro_path}"
        )

    ec_com = _mean_pos(ec_pos)
    ic_com = _mean_pos(ic_pos)

    dx = ec_com[0] - ic_com[0]
    dy = ec_com[1] - ic_com[1]
    dz = ec_com[2] - ic_com[2]

    # Dominant axis: the coordinate with the largest absolute displacement
    axis = max(range(3), key=lambda i: abs((dx, dy, dz)[i]))
    delta = (dx, dy, dz)[axis]

    if axis == 0:
        if delta > 0:
            return OrientRotation(0.0, 270.0, 0.0, "TM along +X → rotate 0 270 0")
        else:
            return OrientRotation(0.0, 90.0, 0.0, "TM along -X → rotate 0 90 0")
    elif axis == 1:
        if delta > 0:
            return OrientRotation(-90.0, 0.0, 0.0, "TM along +Y → rotate -90 0 0")
        else:
            return OrientRotation(90.0, 0.0, 0.0, "TM along -Y → rotate 90 0 0")
    else:
        if delta > 0:
            return OrientRotation(0.0, 0.0, 0.0, "TM along +Z, EC at +Z → no rotation needed")
        else:
            return OrientRotation(180.0, 0.0, 0.0, "TM along -Z → rotate 180 0 0")


# ═══════════════════════════════════════════════════════════════════════════════
# Phase 4 — Geometric validation of an oriented GRO
# ═══════════════════════════════════════════════════════════════════════════════

class GeometricValidationResult(NamedTuple):
    level:   str   # "ok", "warning", "error"
    code:    str   # machine-readable tag
    message: str   # human-readable description


def validate_oriented_gro(
    gro_path:             Path,
    ec_resids:            list[int],
    ic_resids:            list[int],
    tm_resids:            list[int] | None = None,
    extracellular_side:   str = "+z",
) -> list[GeometricValidationResult]:
    """
    Validate that an oriented GRO file is geometrically consistent with the
    structural annotation.

    Checks:
        1. EC COM is on the expected Z side (extracellular_side).
        2. IC COM is on the opposite side.
        3. EC and IC COMs are on opposite sides of each other.
        4. If TM residues are provided, their COM is between EC and IC Z positions.

    Returns a list of GeometricValidationResult. An empty list means no issues.
    """
    results: list[GeometricValidationResult] = []

    ec_pos = _read_gro_ca_positions(gro_path, set(ec_resids))
    ic_pos = _read_gro_ca_positions(gro_path, set(ic_resids))

    if not ec_pos:
        results.append(GeometricValidationResult(
            "warning", "ec_no_ca",
            f"No Cα atoms found for EC residues in {gro_path.name}. "
            "Validation skipped for EC check."
        ))
    if not ic_pos:
        results.append(GeometricValidationResult(
            "warning", "ic_no_ca",
            f"No Cα atoms found for IC residues in {gro_path.name}. "
            "Validation skipped for IC check."
        ))

    if not ec_pos or not ic_pos:
        return results

    ec_com = _mean_pos(ec_pos)
    ic_com = _mean_pos(ic_pos)
    ec_z   = ec_com[2]
    ic_z   = ic_com[2]

    # Check 1 & 2: EC and IC are on opposite sides
    ec_above_ic = ec_z > ic_z
    if not ec_above_ic and not (ec_z < ic_z):
        results.append(GeometricValidationResult(
            "warning", "ec_ic_same_z",
            f"EC COM z={ec_z:.3f} nm and IC COM z={ic_z:.3f} nm are identical — "
            "cannot determine orientation."
        ))
    else:
        # Check: does EC land on the declared side?
        want_ec_above = (extracellular_side == "+z")
        if ec_above_ic != want_ec_above:
            expected = "+Z" if want_ec_above else "-Z"
            actual   = "+Z" if ec_above_ic else "-Z"
            results.append(GeometricValidationResult(
                "error", "orientation_inverted",
                f"Orientation appears inverted: EC COM at {actual} "
                f"(expected {expected} per structural_annotation.orientation.extracellular_side='{extracellular_side}'). "
                f"EC z={ec_z:.3f} nm, IC z={ic_z:.3f} nm. "
                "Check residue numbers or flip extracellular_side in the YAML."
            ))
        else:
            results.append(GeometricValidationResult(
                "ok", "orientation_correct",
                f"EC COM at {'+'  if ec_above_ic else '-'}Z as expected "
                f"(EC z={ec_z:.3f} nm, IC z={ic_z:.3f} nm, Δz={abs(ec_z - ic_z):.3f} nm)."
            ))

    # Check 3: TM COM is between EC and IC z-coordinates
    if tm_resids:
        tm_pos = _read_gro_ca_positions(gro_path, set(tm_resids))
        if tm_pos:
            tm_com = _mean_pos(tm_pos)
            tm_z   = tm_com[2]
            z_min  = min(ec_z, ic_z)
            z_max  = max(ec_z, ic_z)
            margin = (z_max - z_min) * 0.15  # 15% tolerance at boundaries

            if z_min - margin <= tm_z <= z_max + margin:
                results.append(GeometricValidationResult(
                    "ok", "tm_in_membrane_zone",
                    f"TM COM z={tm_z:.3f} nm is within the EC–IC span "
                    f"[{z_min:.3f}, {z_max:.3f}] nm."
                ))
            else:
                results.append(GeometricValidationResult(
                    "warning", "tm_outside_membrane_zone",
                    f"TM COM z={tm_z:.3f} nm is outside the EC–IC span "
                    f"[{z_min:.3f}, {z_max:.3f}] nm. "
                    "Check TM segment residue numbers."
                ))

    return results
