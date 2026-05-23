"""
core/geometry_advisor.py — Compile-time geometry analysis and scientific advisories.

Reads PDB coordinates, computes bounding-box dimensions and aspect ratio,
estimates system size (atom count, solvent overhead), and emits structured
advisories that the CLI can display before building a workspace.

Advisory levels:
  INFO    — neutral observations (system size estimate)
  WARNING — potentially expensive or problematic geometry
  SUGGEST — actionable suggestions for the user

All functions are pure (no I/O side effects beyond reading the PDB).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


AdvisoryLevel = Literal["INFO", "WARNING", "SUGGEST"]

# Water density for TIP3P/SPC at ~300 K: ~33 molecules/nm³
_WATER_DENSITY_NM3 = 33.0
# Approximate average mass of a residue (amino acid) for rough atom count
_AVG_ATOMS_PER_RESIDUE = 8.0  # Cα + heavy atoms, no hydrogens
# Aspect ratio thresholds
_ELONGATED_THRESHOLD      = 1.8
_HIGHLY_ELONGATED_THRESHOLD = 2.8


@dataclass
class Advisory:
    level:   AdvisoryLevel
    message: str
    detail:  str = ""


@dataclass
class GeometryReport:
    """
    Output of GeometryAdvisor.analyze().

    All lengths in nanometres; atom counts are estimates (not GROMACS-exact).
    """
    # Bounding box of the raw PDB coordinates (no padding)
    dim_x: float = 0.0
    dim_y: float = 0.0
    dim_z: float = 0.0

    # Effective simulation box after adding padding on all sides
    box_x: float = 0.0
    box_y: float = 0.0
    box_z: float = 0.0

    # Geometry classification
    aspect_ratio: float = 1.0   # max_dim / min_dim of the raw protein
    geometry:     str   = "compact"
    # "compact"          → aspect < 1.8
    # "elongated"        → 1.8 ≤ aspect < 2.8
    # "highly_elongated" → aspect ≥ 2.8

    # Size estimates
    n_protein_atoms:   int = 0
    n_water_estimated: int = 0
    n_total_estimated: int = 0

    # Advisories generated
    advisories: list[Advisory] = field(default_factory=list)

    def has_warnings(self) -> bool:
        return any(a.level in ("WARNING", "SUGGEST") for a in self.advisories)


class GeometryAdvisor:
    """
    Analyzes PDB geometry and generates scientific advisories.

    Usage:
        report = GeometryAdvisor().analyze(pdb_path, padding_nm=1.0)
    """

    def analyze(
        self,
        pdb_path: str | Path,
        padding_nm: float = 1.0,
        box_type:   str   = "cubic",
    ) -> GeometryReport:
        """
        Parse *pdb_path*, compute geometry, return GeometryReport with advisories.

        padding_nm : solvent shell thickness (–d flag in gmx editconf), in nm.
        box_type   : "cubic" | "dodecahedron" | "triclinic" — affects volume estimate.
        """
        pdb_path = Path(pdb_path)
        report   = GeometryReport()

        coords = self._parse_coords(pdb_path)
        if not coords:
            report.advisories.append(Advisory(
                level="WARNING",
                message="Could not parse atom coordinates from PDB.",
                detail=f"File: {pdb_path}",
            ))
            return report

        # ── Bounding box ─────────────────────────────────────────────────────
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]

        # PDB coordinates are in Å; convert to nm
        dim_x = (max(xs) - min(xs)) / 10.0
        dim_y = (max(ys) - min(ys)) / 10.0
        dim_z = (max(zs) - min(zs)) / 10.0

        report.dim_x = round(dim_x, 2)
        report.dim_y = round(dim_y, 2)
        report.dim_z = round(dim_z, 2)

        # ── Simulation box (padding on each side) ────────────────────────────
        report.box_x = round(dim_x + 2 * padding_nm, 2)
        report.box_y = round(dim_y + 2 * padding_nm, 2)
        report.box_z = round(dim_z + 2 * padding_nm, 2)

        # ── Aspect ratio / geometry ───────────────────────────────────────────
        dims = sorted([dim_x, dim_y, dim_z])
        min_d, max_d = dims[0], dims[2]
        aspect = max_d / min_d if min_d > 0 else 1.0

        report.aspect_ratio = round(aspect, 2)
        if aspect >= _HIGHLY_ELONGATED_THRESHOLD:
            report.geometry = "highly_elongated"
        elif aspect >= _ELONGATED_THRESHOLD:
            report.geometry = "elongated"
        else:
            report.geometry = "compact"

        # ── Atom count estimates ──────────────────────────────────────────────
        report.n_protein_atoms = len(coords)

        box_vol_nm3 = self._box_volume(
            report.box_x, report.box_y, report.box_z, box_type
        )
        # Subtract protein volume (rough: ~0.5 Å³/Da × 110 Da/residue × residues)
        protein_vol_nm3 = (len(coords) / _AVG_ATOMS_PER_RESIDUE) * 0.012
        solvent_vol_nm3 = max(box_vol_nm3 - protein_vol_nm3, 0.0)

        report.n_water_estimated = int(solvent_vol_nm3 * _WATER_DENSITY_NM3)
        # Ions: roughly 0.15 M NaCl in ~1 ion per ~370 water molecules
        n_ions = max(report.n_water_estimated // 370, 0)
        report.n_total_estimated = (
            report.n_protein_atoms
            + report.n_water_estimated * 3   # TIP3P: 3 atoms per water
            + n_ions * 2                      # Na⁺ + Cl⁻ pairs
        )

        # ── Generate advisories ───────────────────────────────────────────────
        self._advise_geometry(report, padding_nm, box_type)
        self._advise_size(report)

        return report

    # ── Advice generators ─────────────────────────────────────────────────────

    def _advise_geometry(
        self,
        report:     GeometryReport,
        padding_nm: float,
        box_type:   str,
    ) -> None:
        dims_str = (
            f"{report.dim_x:.1f} × {report.dim_y:.1f} × {report.dim_z:.1f} nm"
        )

        if report.geometry == "highly_elongated":
            report.advisories.append(Advisory(
                level="WARNING",
                message=(
                    f"Highly elongated structure detected  "
                    f"(aspect ratio {report.aspect_ratio:.1f}:1,  {dims_str})."
                ),
                detail=(
                    "Non-cubic boxes waste a large fraction of the volume as solvent. "
                    "Simulating highly elongated systems is expensive."
                ),
            ))
            report.advisories.append(Advisory(
                level="SUGGEST",
                message="Align the longest axis with Z and use a smaller padding.",
                detail=(
                    "In your YAML, set:  assembly → principal_axis_alignment: yes\n"
                    f"and reduce:         assembly → padding: {max(padding_nm - 0.2, 0.6):.1f}  "
                    f"(currently {padding_nm:.1f} nm)"
                ),
            ))
            if box_type in ("cubic", "rectangular"):
                report.advisories.append(Advisory(
                    level="SUGGEST",
                    message="Consider a dodecahedral box to reduce solvent volume.",
                    detail=(
                        "In your YAML, set:  assembly → box_type: dodecahedron\n"
                        "A dodecahedral box uses ~29% less volume than a cubic box "
                        "for the same minimum image distance."
                    ),
                ))

        elif report.geometry == "elongated":
            report.advisories.append(Advisory(
                level="INFO",
                message=(
                    f"Elongated structure  "
                    f"(aspect ratio {report.aspect_ratio:.1f}:1,  {dims_str})."
                ),
                detail=(
                    "This is common for rod-like proteins or domains. "
                    "Consider principal-axis alignment if solvent overhead is high."
                ),
            ))

        else:
            report.advisories.append(Advisory(
                level="INFO",
                message=f"Compact structure  ({dims_str}).",
            ))

    def _advise_size(self, report: GeometryReport) -> None:
        n = report.n_total_estimated
        report.advisories.append(Advisory(
            level="INFO",
            message=f"Projected system size: ~{n:,} atoms.",
            detail=(
                f"  protein atoms : {report.n_protein_atoms:,}\n"
                f"  water (TIP3P) : ~{report.n_water_estimated:,} molecules\n"
                f"  box           : "
                f"{report.box_x:.1f} × {report.box_y:.1f} × {report.box_z:.1f} nm"
            ),
        ))

        if n > 300_000:
            report.advisories.append(Advisory(
                level="WARNING",
                message=f"Very large system (~{n // 1000}k atoms) — production will be slow.",
                detail=(
                    "Consider:\n"
                    "  · Reducing padding (assembly → padding)\n"
                    "  · Using a dodecahedral box (assembly → box_type: dodecahedron)\n"
                    "  · Principal-axis alignment to minimize bounding box volume"
                ),
            ))
        elif n > 150_000:
            report.advisories.append(Advisory(
                level="INFO",
                message=f"Moderate system size (~{n // 1000}k atoms). GPU recommended.",
            ))

    # ── PDB parser ────────────────────────────────────────────────────────────

    @staticmethod
    def _parse_coords(pdb_path: Path) -> list[tuple[float, float, float]]:
        """
        Extract (x, y, z) in Å from ATOM/HETATM records.
        Skips hydrogen atoms (element H or name starting with H).
        Returns empty list on any parse error.
        """
        coords: list[tuple[float, float, float]] = []
        try:
            for line in pdb_path.read_text(errors="replace").splitlines():
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                # PDB column format: cols 77-78 = element symbol
                element = line[76:78].strip() if len(line) >= 78 else ""
                atom_name = line[12:16].strip()
                if element.upper() == "H" or atom_name.startswith("H"):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
        except Exception:
            return []
        return coords

    @staticmethod
    def _box_volume(x: float, y: float, z: float, box_type: str) -> float:
        """Approximate box volume in nm³ for the given box type."""
        if box_type == "dodecahedron":
            # Rhombic dodecahedron: V ≈ 0.707 × a³ where a = min(x,y,z)
            a = min(x, y, z)
            return 0.707 * a ** 3
        # cubic / rectangular / triclinic: V = x × y × z
        return x * y * z
