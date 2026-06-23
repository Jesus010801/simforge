"""
adapters/movememb_adapter.py — Python reimplementation of MoveMemb.f

MoveMemb.f is an interactive Fortran program that translates the entire
bilayer along Z so its midplane coincides with the protein Z-centre.
The original is interactive (reads stdin), making it unsuitable for
scripted workflows.

This adapter reimplements the algorithm directly in Python:
  1. compute_z_shift(): reads protein + bilayer GRO, returns the nm shift
     needed to align bilayer midplane with the protein alignment Z-centre.
  2. run(): applies the shift, combines protein + bilayer into one GRO.

TM-aware alignment (Phase 2):
  When tm_residues is provided, the bilayer midplane is aligned to the
  Z-centre of CA atoms in those residues rather than the full protein
  Z-extent centre.  This is critical for proteins with large EC/IC
  domains — using the full extent would mis-position the hydrophobic core.

No gfortran or MoveMemb.f required.

Guaranteed metadata keys on success:
    z_shift_nm          float       shift applied to all bilayer atoms (+ = upward)
    protein_z_min       float       protein Z minimum (nm)
    protein_z_max       float       protein Z maximum (nm)
    protein_center_z    float       (z_min + z_max) / 2 (nm)
    tm_center_z         float|None  mean Z of TM CA atoms (nm); None if no TM annotation
    bilayer_midplane_z  float       original bilayer midplane Z (nm)
    bilayer_z_min       float       original bilayer Z minimum (nm)
    bilayer_z_max       float       original bilayer Z maximum (nm)
    alignment_method    str         "tm_center" | "protein_center_fallback"
    tm_annotation_used  bool        True when TM CA atoms were found and used
    overlap_before_nm   float       Z overlap before shift (negative means gap)
    atoms_protein       int         protein atom count
    atoms_bilayer       int         bilayer atom count
    atoms_total         int         combined atom count
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


def _tm_ca_center_z(
    atom_lines:     list[str],
    tm_residue_set: set[int],
) -> float | None:
    """
    Return the mean Z of CA atoms belonging to TM residues.

    GRO fixed-width columns:
        [0:5]   residue number
        [5:10]  residue name
        [10:15] atom name
        [36:44] Z coordinate (nm)

    Returns None if no matching CA atoms are found (empty TM set, or
    the GRO lacks backbone atoms — e.g. coarse-grained models).
    """
    zs: list[float] = []
    for line in atom_lines:
        try:
            resnum = int(line[0:5].strip())
        except (ValueError, IndexError):
            continue
        if resnum not in tm_residue_set:
            continue
        if line[10:15].strip() == "CA":
            zs.append(_atom_z(line))
    return sum(zs) / len(zs) if zs else None


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
        protein_gro:    Path | str,
        bilayer_gro:    Path | str,
        clearance_nm:   float = 0.0,
        tm_residues:    Optional[set[int]] = None,
    ) -> float:
        """
        Compute the Z shift needed to align the bilayer midplane with the protein
        TM-region centre (or full protein Z-centre as a fallback).

        When tm_residues is provided the shift is computed from the mean Z of
        CA atoms in those residues.  When absent, the full protein Z-extent
        centre is used (original behaviour, but potentially inaccurate for
        proteins with large EC/IC domains).

        Args:
            protein_gro:  Oriented protein .gro.
            bilayer_gro:  Pre-built bilayer .gro.
            clearance_nm: Extra gap added to the shift (nm); normally 0 for TM.
            tm_residues:  Set of residue numbers annotated as TM; if None the
                          full protein Z-centre is used as fallback.

        Returns:
            z_shift_nm (positive = bilayer moves upward in Z).
        """
        _, prot_lines, _ = _parse_gro(Path(protein_gro))
        _, bil_lines,  _ = _parse_gro(Path(bilayer_gro))

        prot_z_min, prot_z_max = _z_extents(prot_lines)
        bil_z_min,  bil_z_max  = _z_extents(bil_lines)
        bil_z_midplane = (bil_z_min + bil_z_max) / 2.0

        if tm_residues is not None:
            tm_z = _tm_ca_center_z(prot_lines, tm_residues)
        else:
            tm_z = None

        alignment_z = tm_z if tm_z is not None else (prot_z_min + prot_z_max) / 2.0
        return alignment_z - bil_z_midplane + clearance_nm

    # ── run ───────────────────────────────────────────────────────────────────

    def run(  # type: ignore[override]
        self,
        *,
        protein_gro:  Path | str,
        bilayer_gro:  Path | str,
        gro_out:      Path | str,
        z_shift_nm:   Optional[float] = None,
        clearance_nm: float = 0.0,
        tm_residues:  Optional[set[int]] = None,
    ) -> AdapterResult:
        """
        Shift the bilayer in Z and write a combined protein+bilayer .gro.

        The output GRO contains all protein atoms first, followed by all
        bilayer atoms with Z shifted by z_shift_nm.  Box dimensions are
        taken from the protein GRO (which should already be sized to match
        the bilayer XY footprint).

        When tm_residues is provided the bilayer is aligned to the mean Z of
        CA atoms in those residues, not to the full protein Z-extent centre.
        This is the correct behaviour for proteins with large EC/IC domains.

        Args:
            protein_gro:  Oriented + box-sized protein .gro.
            bilayer_gro:  Pre-built equilibrated bilayer .gro.
            gro_out:      Output combined system .gro.
            z_shift_nm:   Explicit shift (nm).  None = auto-compute.
            clearance_nm: Extra gap added to the shift (nm); normally 0.
            tm_residues:  Set of residue numbers annotated as TM.  When None
                          the full protein Z-centre is used as fallback.
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

            # ── TM-aware Z-centre computation ─────────────────────────────────
            tm_center_z: float | None = None
            if tm_residues is not None:
                tm_center_z = _tm_ca_center_z(prot_lines, tm_residues)

            if tm_center_z is not None:
                alignment_z        = tm_center_z
                alignment_method   = "tm_center"
                tm_annotation_used = True
            else:
                alignment_z        = prot_z_mid
                alignment_method   = "protein_center_fallback"
                tm_annotation_used = False

            if z_shift_nm is None:
                z_shift_nm = alignment_z - bil_z_mid + clearance_nm

            # Apply shift to bilayer atoms
            shifted_bil = [_shift_atom_z(l, z_shift_nm) for l in bil_lines]

            # Combine: protein first, then bilayer
            combined = _renumber_atoms(prot_lines + shifted_bil, start=1)
            n_total  = len(combined)

            title = f"Protein + bilayer (dz={z_shift_nm:+.3f} nm, method={alignment_method})"
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

        stdout_lines = [
            f"Z-shift applied:    {z_shift_nm:+.3f} nm",
            f"Alignment method:   {alignment_method}",
        ]
        if tm_annotation_used and tm_center_z is not None:
            stdout_lines.append(f"TM center Z:        {tm_center_z:.3f} nm")
        else:
            stdout_lines.append(
                "WARNING: TM annotation absent — bilayer aligned to full protein Z centre"
            )
        stdout_lines += [
            f"Protein Z centre:   {prot_z_mid:.3f} nm  "
            f"[{prot_z_min:.3f}, {prot_z_max:.3f}]",
            f"Bilayer midplane:   {bil_z_mid:.3f} nm  "
            f"[{bil_z_min:.3f}, {bil_z_max:.3f}]",
            f"Combined atoms:     {n_total}",
        ]

        return self._make_result(
            tool_name=self.tool_name,
            adapter_type=type(self).__name__,
            success=True,
            exit_code=0,
            stdout="\n".join(stdout_lines),
            started_at=started_at,
            outputs={"gro_out": str(gro_out)},
            metadata={
                "z_shift_nm":          z_shift_nm,
                "protein_z_min":       prot_z_min,
                "protein_z_max":       prot_z_max,
                "protein_center_z":    prot_z_mid,
                "tm_center_z":         tm_center_z,
                "bilayer_midplane_z":  bil_z_mid,
                "bilayer_z_min":       bil_z_min,
                "bilayer_z_max":       bil_z_max,
                "alignment_method":    alignment_method,
                "tm_annotation_used":  tm_annotation_used,
                "overlap_before_nm":   overlap_nm,
                "atoms_protein":       len(prot_lines),
                "atoms_bilayer":       len(bil_lines),
                "atoms_total":         n_total,
            },
        )
