"""
adapters/movememb_adapter.py — Python reimplementation of MoveMemb.f

MoveMemb.f is an interactive Fortran program that translates the entire
bilayer along Z so its midplane coincides with the protein Z-centre.
The original is interactive (reads stdin), making it unsuitable for
scripted workflows.

This adapter reimplements the algorithm directly in Python:
  1. compute_z_shift(): reads protein + bilayer GRO, returns the nm shift
     needed to align bilayer midplane with protein Z-centre.
  2. run(): applies the shift, combines protein + bilayer into one GRO.

No gfortran or MoveMemb.f required.

Guaranteed metadata keys on success:
    z_shift_nm        float  shift applied to all bilayer atoms (+ = upward)
    protein_z_min     float  protein Z minimum (nm)
    protein_z_max     float  protein Z maximum (nm)
    bilayer_z_min     float  original bilayer Z minimum (nm)
    bilayer_z_max     float  original bilayer Z maximum (nm)
    overlap_before_nm float  Z overlap before shift (negative means gap)
    atoms_protein     int    protein atom count
    atoms_bilayer     int    bilayer atom count
    atoms_total       int    combined atom count
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

from adapters.base import (
    AdapterResult,
    AvailabilityResult,
    ExternalToolAdapter,
    PreconditionViolation,
)


# ── GRO I/O helpers ───────────────────────────────────────────────────────────

def _parse_gro(path: Path) -> tuple[str, list[str], list[float]]:
    """Return (title, atom_lines, box_xyz) from a .gro file."""
    text = path.read_text().splitlines()
    title = text[0]
    n = int(text[1].strip())
    atom_lines = text[2: 2 + n]
    box_parts  = text[2 + n].split()
    box        = [float(v) for v in box_parts[:3]]
    return title, atom_lines, box


def _atom_z(line: str) -> float:
    return float(line[36:44])


def _shift_atom_z(line: str, dz: float) -> str:
    """Return the atom line with Z coordinate shifted by dz nm."""
    z_new = _atom_z(line) + dz
    return line[:36] + f"{z_new:8.3f}" + line[44:]


def _renumber_atoms(atom_lines: list[str], start: int = 1) -> list[str]:
    """Renumber the atom-number field (chars 15-19) starting from `start`."""
    out = []
    for i, line in enumerate(atom_lines, start=start):
        out.append(line[:15] + f"{i % 100000:5d}" + line[20:])
    return out


def _z_extents(atom_lines: list[str]) -> tuple[float, float]:
    zs = [_atom_z(l) for l in atom_lines]
    return min(zs), max(zs)


# ── Adapter ───────────────────────────────────────────────────────────────────

class MoveMembAdapter(ExternalToolAdapter):
    """
    Aligns a lipid bilayer with a membrane protein along the Z axis.

    Computes the Z shift needed so the bilayer midplane coincides with
    the protein Z-centre, then writes a combined .gro with both systems.

    Python reimplementation — no gfortran or MoveMemb.f required.

    Args:
        fortran_source: Ignored (kept for API compatibility).
        prefer_python:  Always True in the current implementation.
    """

    tool_name = "movememb"

    def __init__(
        self,
        fortran_source: Optional[Path | str] = None,
        prefer_python:  bool = True,
    ) -> None:
        self._fortran_source = Path(fortran_source).resolve() if fortran_source else None
        self._prefer_python  = prefer_python

    # ── Availability ──────────────────────────────────────────────────────────

    def check_availability(self) -> AvailabilityResult:
        return AvailabilityResult(
            available=True,
            tool_name=self.tool_name,
            binary_path="python3 (built-in reimplementation)",
        )

    # ── Preconditions ─────────────────────────────────────────────────────────

    def validate_preconditions(self, **kwargs) -> list[PreconditionViolation]:
        violations: list[PreconditionViolation] = []
        for field in ("protein_gro", "bilayer_gro"):
            p = kwargs.get(field)
            if p is None or not Path(p).exists():
                violations.append(PreconditionViolation(field, f"file not found: {p}"))
        gro_out = kwargs.get("gro_out")
        if gro_out is not None and not Path(gro_out).parent.exists():
            violations.append(PreconditionViolation("gro_out", f"parent directory does not exist: {Path(gro_out).parent}"))
        return violations

    # ── compute_z_shift ───────────────────────────────────────────────────────

    def compute_z_shift(
        self,
        protein_gro:  Path | str,
        bilayer_gro:  Path | str,
        clearance_nm: float = 0.0,
    ) -> float:
        """
        Compute the Z shift needed to centre the bilayer at the protein Z-midpoint.

        The bilayer midplane is moved to match the protein Z-centre so that
        the TM region (assumed to span the full protein height in v1) sits
        in the hydrophobic core of the bilayer.

        Args:
            protein_gro:  Oriented protein .gro.
            bilayer_gro:  Pre-built bilayer .gro.
            clearance_nm: Extra gap added above the bilayer top (nm).  Use 0
                          for TM proteins — the protein should sit inside the
                          bilayer, not above it.

        Returns:
            z_shift_nm (positive = bilayer moves upward in Z).
        """
        _, prot_lines, _ = _parse_gro(Path(protein_gro))
        _, bil_lines,  _ = _parse_gro(Path(bilayer_gro))

        prot_z_min, prot_z_max = _z_extents(prot_lines)
        bil_z_min,  bil_z_max  = _z_extents(bil_lines)

        prot_z_centre  = (prot_z_min + prot_z_max) / 2.0
        bil_z_midplane = (bil_z_min  + bil_z_max)  / 2.0

        return prot_z_centre - bil_z_midplane + clearance_nm

    # ── run ───────────────────────────────────────────────────────────────────

    def run(  # type: ignore[override]
        self,
        *,
        protein_gro:  Path | str,
        bilayer_gro:  Path | str,
        gro_out:      Path | str,
        z_shift_nm:   Optional[float] = None,
        clearance_nm: float = 0.0,
    ) -> AdapterResult:
        """
        Shift the bilayer in Z and write a combined protein+bilayer .gro.

        The output GRO contains all protein atoms first, followed by all
        bilayer atoms with Z shifted by z_shift_nm.  Box dimensions are
        taken from the protein GRO (which should already be sized to match
        the bilayer XY footprint).

        Args:
            protein_gro:  Oriented + box-sized protein .gro.
            bilayer_gro:  Pre-built equilibrated bilayer .gro.
            gro_out:      Output combined system .gro.
            z_shift_nm:   Explicit shift (nm).  None = auto-compute via
                          compute_z_shift().
            clearance_nm: Passed to compute_z_shift() when z_shift_nm is None.
        """
        started_at = datetime.now()
        self.assert_available()
        self.assert_preconditions(
            protein_gro=protein_gro, bilayer_gro=bilayer_gro, gro_out=gro_out,
        )

        protein_gro = Path(protein_gro).resolve()
        bilayer_gro = Path(bilayer_gro).resolve()
        gro_out     = Path(gro_out).resolve()

        try:
            prot_title, prot_lines, prot_box = _parse_gro(protein_gro)
            _,          bil_lines,  _         = _parse_gro(bilayer_gro)

            prot_z_min, prot_z_max = _z_extents(prot_lines)
            bil_z_min,  bil_z_max  = _z_extents(bil_lines)

            bil_z_mid   = (bil_z_min + bil_z_max) / 2.0
            prot_z_mid  = (prot_z_min + prot_z_max) / 2.0
            overlap_nm  = min(prot_z_max, bil_z_max) - max(prot_z_min, bil_z_min)

            if z_shift_nm is None:
                z_shift_nm = prot_z_mid - bil_z_mid + clearance_nm

            # Apply shift to bilayer atoms
            shifted_bil = [_shift_atom_z(l, z_shift_nm) for l in bil_lines]

            # Combine: protein first, then bilayer
            combined = _renumber_atoms(prot_lines + shifted_bil, start=1)
            n_total  = len(combined)

            title = f"Protein + bilayer (dz={z_shift_nm:+.3f} nm)"
            box_str = "  ".join(f"{v:.5f}" for v in prot_box)
            lines   = [title, f"{n_total}", *combined, box_str]
            gro_out.write_text("\n".join(lines) + "\n")

        except Exception as exc:
            return self._make_result(
                tool_name=self.tool_name,
                adapter_type=type(self).__name__,
                success=False,
                started_at=started_at,
                error_message=str(exc),
                stderr=str(exc),
            )

        return self._make_result(
            tool_name=self.tool_name,
            adapter_type=type(self).__name__,
            success=True,
            exit_code=0,
            stdout=(
                f"Z-shift applied: {z_shift_nm:+.3f} nm\n"
                f"Protein Z: [{prot_z_min:.3f}, {prot_z_max:.3f}] nm  "
                f"centre={prot_z_mid:.3f}\n"
                f"Bilayer Z (original): [{bil_z_min:.3f}, {bil_z_max:.3f}] nm  "
                f"midplane={bil_z_mid:.3f}\n"
                f"Combined atoms: {n_total}"
            ),
            started_at=started_at,
            outputs={"gro_out": str(gro_out)},
            metadata={
                "z_shift_nm":        z_shift_nm,
                "protein_z_min":     prot_z_min,
                "protein_z_max":     prot_z_max,
                "bilayer_z_min":     bil_z_min,
                "bilayer_z_max":     bil_z_max,
                "overlap_before_nm": overlap_nm,
                "atoms_protein":     len(prot_lines),
                "atoms_bilayer":     len(bil_lines),
                "atoms_total":       n_total,
            },
        )
