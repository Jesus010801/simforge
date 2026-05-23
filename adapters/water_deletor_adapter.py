"""
adapters/water_deletor_adapter.py — Python reimplementation of water_deletor.pl

water_deletor.pl (Lemkul, 2017) removes water molecules whose oxygen
falls within the bilayer hydrophobic core, defined by:
  - reference atom (e.g. O33 = headgroup phosphate oxygen in DPPC OPLS-AA)
  - middle atom   (e.g. C50 = terminal tail carbon)

The script splits the bilayer into top/bottom leaflets using z(ref) vs
z(middle), computes average z-boundaries per leaflet, and deletes any
SOL molecule whose OW falls within those boundaries.

This adapter reimplements the algorithm directly in Python; no Perl
interpreter or external script is required.

Guaranteed metadata keys on success:
    waters_removed   int    number of water molecules deleted
    atoms_out        int    total atoms in output .gro
    atoms_in         int    total atoms in input .gro
    z_top_nm         float  top leaflet headgroup boundary (nm)
    z_bot_nm         float  bottom leaflet headgroup boundary (nm)
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


# ── GRO parsing helpers ───────────────────────────────────────────────────────

def _parse_gro(path: Path) -> tuple[str, list[str], str]:
    """Return (title, atom_lines, box_line) from a .gro file."""
    text = path.read_text().splitlines()
    title = text[0]
    n = int(text[1].strip())
    return title, text[2: 2 + n], text[2 + n]


def _atom_fields(line: str) -> dict:
    """Parse one GRO atom line into a dict.  Keeps original line for output."""
    return {
        "resnum":   int(line[0:5]),
        "resname":  line[5:10].strip(),
        "atomname": line[10:15].strip(),
        "z":        float(line[36:44]),
        "raw":      line,
    }


# ── Core deletion algorithm ───────────────────────────────────────────────────

def _delete_bilayer_waters(
    atom_lines:      list[str],
    ref_atom:        str,
    middle_atom:     str,
    water_resname:   str = "SOL",
    water_oxygen:    str = "OW",
) -> tuple[list[str], int, float, float]:
    """
    Remove water molecules whose OW falls inside the bilayer hydrophobic core.

    Algorithm (mirrors water_deletor.pl):
      1. Collect Z of ref_atom from lipid residues → headgroup boundaries.
      2. Collect Z of middle_atom from lipid residues → bilayer midplane.
      3. Split headgroups into top/bottom leaflets by midplane.
      4. z_top = mean Z of top-leaflet headgroups.
         z_bot = mean Z of bottom-leaflet headgroups.
      5. Delete any SOL residue whose OW Z ∈ [z_bot, z_top].

    Returns (filtered_atom_lines, n_waters_removed, z_top, z_bot).
    """
    atoms = [_atom_fields(l) for l in atom_lines]

    # ── Step 1-2: bilayer geometry ────────────────────────────────────────────
    ref_z = [a["z"] for a in atoms if a["atomname"] == ref_atom and a["resname"] != water_resname]
    mid_z = [a["z"] for a in atoms if a["atomname"] == middle_atom and a["resname"] != water_resname]

    if not ref_z:
        raise ValueError(f"No atoms named '{ref_atom}' found in non-solvent residues")
    if not mid_z:
        raise ValueError(f"No atoms named '{middle_atom}' found in non-solvent residues")

    z_midplane = sum(mid_z) / len(mid_z)

    top_z = [z for z in ref_z if z > z_midplane]
    bot_z = [z for z in ref_z if z <= z_midplane]

    z_top = sum(top_z) / len(top_z) if top_z else z_midplane + 1.0
    z_bot = sum(bot_z) / len(bot_z) if bot_z else z_midplane - 1.0

    # ── Step 3-5: mark SOL residues for deletion ──────────────────────────────
    delete_indices: set[int] = set()
    n_mols_removed = 0
    i = 0

    while i < len(atoms):
        a = atoms[i]
        if a["resname"] != water_resname:
            i += 1
            continue

        # Walk to end of this SOL residue (same resnum + resname)
        j = i + 1
        while (
            j < len(atoms)
            and atoms[j]["resnum"] == a["resnum"]
            and atoms[j]["resname"] == a["resname"]
        ):
            j += 1

        # Check OW position
        ow = next((atoms[k] for k in range(i, j) if atoms[k]["atomname"] == water_oxygen), None)
        if ow is not None and z_bot <= ow["z"] <= z_top:
            delete_indices.update(range(i, j))
            n_mols_removed += 1

        i = j

    # ── Rebuild atom lines with renumbered atom indices ───────────────────────
    out_lines: list[str] = []
    atom_num = 0
    for k, a in enumerate(atoms):
        if k in delete_indices:
            continue
        atom_num += 1
        line = a["raw"]
        out_lines.append(line[:15] + f"{atom_num % 100000:5d}" + line[20:])

    return out_lines, n_mols_removed, z_top, z_bot


# ── Adapter ───────────────────────────────────────────────────────────────────

class WaterDeletorAdapter(ExternalToolAdapter):
    """
    Removes water molecules embedded in a lipid bilayer.

    Pure Python reimplementation of water_deletor.pl — no Perl interpreter
    or external script required.

    Args:
        script_path: Ignored (kept for API compatibility).
        timeout_s:   Unused (pure Python, no subprocess).
    """

    tool_name    = "water_deletor"
    _WATER_OX    = "OW"
    _WATER_RESN  = "SOL"

    def __init__(
        self,
        script_path: Optional[Path | str] = None,
        timeout_s:   int = 120,
    ) -> None:
        self._script    = Path(script_path).resolve() if script_path else None
        self._timeout_s = timeout_s

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
        gro_in      = kwargs.get("gro_in")
        gro_out     = kwargs.get("gro_out")
        ref_atom    = kwargs.get("ref_atom")
        middle_atom = kwargs.get("middle_atom")
        nwater      = kwargs.get("nwater", 3)

        if gro_in is None or not Path(gro_in).exists():
            violations.append(PreconditionViolation("gro_in", f"file not found: {gro_in}"))
        if gro_out is not None and not Path(gro_out).parent.exists():
            violations.append(PreconditionViolation("gro_out", f"parent directory does not exist: {Path(gro_out).parent}"))
        if not ref_atom:
            violations.append(PreconditionViolation("ref_atom", "required (e.g. 'O33' for DPPC OPLS-AA)"))
        if not middle_atom:
            violations.append(PreconditionViolation("middle_atom", "required (e.g. 'C50' for DPPC OPLS-AA)"))
        if int(nwater) < 1:
            violations.append(PreconditionViolation("nwater", "must be ≥ 1"))
        return violations

    # ── Execution ─────────────────────────────────────────────────────────────

    def run(  # type: ignore[override]
        self,
        *,
        gro_in:      Path | str,
        gro_out:     Path | str,
        ref_atom:    str,
        middle_atom: str,
        nwater:      int  = 3,
        verbose:     bool = False,
    ) -> AdapterResult:
        """
        Remove bilayer-embedded water molecules.

        Args:
            gro_in:      Input .gro (solvated system).
            gro_out:     Output .gro (bilayer-interior waters removed).
            ref_atom:    Headgroup reference atom name (e.g. "O33").
            middle_atom: Tail atom marking bilayer centre (e.g. "C50").
            nwater:      Atoms per water molecule (3 for TIP3P/SPC; unused in
                         Python impl — residue boundaries are detected automatically).
            verbose:     If True, include per-residue Z values in stdout.
        """
        started_at = datetime.now()
        self.assert_available()
        self.assert_preconditions(
            gro_in=gro_in, gro_out=gro_out,
            ref_atom=ref_atom, middle_atom=middle_atom, nwater=nwater,
        )

        gro_in  = Path(gro_in).resolve()
        gro_out = Path(gro_out).resolve()

        try:
            title, atom_lines, box_line = _parse_gro(gro_in)
            atoms_in = len(atom_lines)

            kept, n_removed, z_top, z_bot = _delete_bilayer_waters(
                atom_lines, ref_atom, middle_atom,
                water_resname=self._WATER_RESN,
                water_oxygen=self._WATER_OX,
            )
            atoms_out = len(kept)

            lines = [title, f"{atoms_out}", *kept, box_line]
            gro_out.write_text("\n".join(lines) + "\n")

            stdout = (
                f"Bilayer boundaries: z_bot={z_bot:.3f} nm  z_top={z_top:.3f} nm\n"
                f"Removed {n_removed} water molecules  "
                f"({atoms_in} → {atoms_out} atoms)"
            )

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
            stdout=stdout,
            started_at=started_at,
            outputs={"gro_out": str(gro_out)},
            metadata={
                "waters_removed": n_removed,
                "atoms_in":       atoms_in,
                "atoms_out":      atoms_out,
                "z_top_nm":       z_top,
                "z_bot_nm":       z_bot,
            },
        )
