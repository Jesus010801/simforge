"""
LigParGen legacy PDB writer (Phase 2 follow-up).

Background
----------
RDKit's MolToPDBBlock() produces HETATM records with numbered atom names
(C1, C2, H1, H2 …) which LigParGen's online server rejects.
The server accepts the older PDB convention:

  - ATOM    records  (not HETATM)
  - Simple element-symbol atom names: " C  ", " O  ", " H  ", "CL  "
  - Full CONECT block for every bond (including H–X bonds)
  - Coordinates preserved exactly from the RDKit conformer

The standard fixed-width PDB column layout (1-indexed):
  1–6   record name   ("ATOM  ")
  7–11  atom serial   (right-justified integer)
  12    space
  13–16 atom name     (4 chars; see _pdb_atom_name)
  17    alt loc       (space)
  18–20 residue name  (3 chars)
  21    space
  22    chain ID      (space = no chain)
  23–26 residue seq   (right-justified integer)
  27–30 spaces        (iCode + 3 blanks)
  31–38 x coord       (%8.3f)
  39–46 y coord       (%8.3f)
  47–54 z coord       (%8.3f)
  55–60 occupancy     (%6.2f)
  61–66 B-factor      (%6.2f)

Public API
----------
    LigParGenLegacyWriter().write(mol, path, mol_name)
    write_legacy_pdb(mol, path, mol_name)           ← convenience alias
"""

from __future__ import annotations

from pathlib import Path


class LigParGenLegacyWriter:
    """
    Write an RDKit Mol to the legacy PDB format accepted by LigParGen.

    The output file is named  ``<mol_name>_ligpargen_legacy.pdb``  when
    called via :func:`ligand.export.export_for_ligpargen_legacy`, or the
    caller supplies the full path directly.
    """

    def write(self, mol, path: str | Path, mol_name: str = "LIG") -> None:
        """
        Write *mol* to *path*.

        Parameters
        ----------
        mol      : RDKit Mol with at least one 3D conformer and explicit Hs.
        path     : Destination .pdb file (parents must exist or be created by
                   the caller).
        mol_name : Residue name written into the PDB (≤3 chars, auto-truncated
                   and uppercased).  Appears in both ATOM records and as the
                   stem of any default filename.
        """
        path = Path(path)
        resname = mol_name[:3].upper().ljust(3)

        conf = mol.GetConformer()
        lines: list[str] = []

        # ── ATOM records ──────────────────────────────────────────────────────
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            serial = idx + 1
            name_field = _pdb_atom_name(atom.GetSymbol())
            pos = conf.GetAtomPosition(idx)

            # Exact column layout verified against PDB standard (see module docstring)
            line = (
                "ATOM  "
                f"{serial:5d}"
                " "
                f"{name_field}"   # 4 chars: " C  ", " O  ", "CL  ", …
                " "               # col 17 — alt location indicator (blank)
                f"{resname}"      # cols 18-20
                " "               # col 21 — blank chain ID separator
                " "               # col 22 — chain ID (blank = no chain)
                f"{1:4d}"         # cols 23-26 — residue sequence number
                "    "            # cols 27-30 — iCode + 3 blanks
                f"{pos.x:8.3f}"
                f"{pos.y:8.3f}"
                f"{pos.z:8.3f}"
                "  1.00  0.00"   # occupancy + B-factor
            )
            lines.append(line)

        # ── CONECT records ────────────────────────────────────────────────────
        # One CONECT entry per atom listing all bonded neighbours.
        # Up to 4 neighbours per line; wrap to additional lines beyond that.
        for atom in mol.GetAtoms():
            serial = atom.GetIdx() + 1
            neighbours = [
                bond.GetOtherAtomIdx(atom.GetIdx()) + 1
                for bond in atom.GetBonds()
            ]
            if not neighbours:
                continue
            for chunk_start in range(0, len(neighbours), 4):
                chunk = neighbours[chunk_start : chunk_start + 4]
                lines.append(
                    f"CONECT{serial:5d}" + "".join(f"{n:5d}" for n in chunk)
                )

        lines.append("END")
        path.write_text("\n".join(lines) + "\n")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _pdb_atom_name(element: str) -> str:
    """
    Format element symbol as the 4-char PDB atom-name field (cols 13–16).

    PDB convention:
      1-char element (C, N, O, S, H …) → " X  "  (leading space, then element, 2 trailing)
      2-char element (Cl, Br, Fe …)    → "XX  "  (element at col 13, 2 trailing spaces)

    This produces simple names with no trailing digits, which is exactly what
    LigParGen requires.
    """
    el = element.strip().upper()
    if len(el) == 1:
        return f" {el}  "
    elif len(el) == 2:
        return f"{el}  "
    # Fallback: truncate to 4 chars (unusual elements)
    return f"{el[:4]:<4s}"


# Convenience alias so callers do not need to instantiate the class
def write_legacy_pdb(mol, path: str | Path, mol_name: str = "LIG") -> None:
    """Thin wrapper — calls ``LigParGenLegacyWriter().write(mol, path, mol_name)``."""
    LigParGenLegacyWriter().write(mol, path, mol_name)
