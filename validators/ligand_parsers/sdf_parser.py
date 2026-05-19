# validators/ligand_parsers/sdf_parser.py
"""
Parser SDF mínimo (implementación actual, línea por línea).

Interfaz pública:
    parse_sdf(path: Path) -> dict

Retorna:
    {
        "mol_name": str,
        "n_atoms":  int,
        "n_bonds":  int,
        "atoms":    list[dict],   # index, element, x, y, z, charge
        "bonds":    list[dict],   # atom1, atom2, bond_type
        "error":    str | None,
    }

Reemplazar esta implementación por RDKit = solo tocar este archivo.
"""

from pathlib import Path


def parse_sdf(path: Path) -> dict:
    lines = path.read_text().splitlines()

    if len(lines) < 4:
        return {"error": "Archivo SDF demasiado corto o vacío"}

    # ─── Header ──────────────────────────────────────────────────────────────
    mol_name = lines[0].strip()

    # ─── Counts line (línea 4, índice 3) ─────────────────────────────────────
    counts_line = lines[3]
    try:
        n_atoms = int(counts_line[0:3].strip())
        n_bonds = int(counts_line[3:6].strip())
    except (ValueError, IndexError):
        return {"error": f"Counts line inválida: '{counts_line}'"}

    # ─── Atom block ──────────────────────────────────────────────────────────
    atoms = []
    for i in range(n_atoms):
        line_idx = 4 + i
        if line_idx >= len(lines):
            return {"error": f"Atom block incompleto en línea {line_idx}"}
        line = lines[line_idx]
        try:
            x       = float(line[0:10].strip())
            y       = float(line[10:20].strip())
            z       = float(line[20:30].strip())
            element = line[31:34].strip()

            charge_code = 0
            try:
                charge_code = int(line[36:39].strip()) if len(line) > 36 else 0
            except ValueError:
                charge_code = 0

            # SDF charge codes → carga formal real
            charge_map = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
            formal_charge = charge_map.get(charge_code, 0)

            atoms.append({
                "index":   i + 1,
                "element": element,
                "x": x, "y": y, "z": z,
                "charge":  formal_charge,
            })
        except (ValueError, IndexError) as e:
            return {"error": f"Error parseando átomo {i+1}: {e}"}

    # ─── Bond block ──────────────────────────────────────────────────────────
    bonds = []
    for i in range(n_bonds):
        line_idx = 4 + n_atoms + i
        if line_idx >= len(lines):
            return {"error": f"Bond block incompleto en línea {line_idx}"}
        line = lines[line_idx]
        try:
            atom1     = int(line[0:3].strip())
            atom2     = int(line[3:6].strip())
            bond_type = int(line[6:9].strip())
            bonds.append({
                "atom1":     atom1,
                "atom2":     atom2,
                "bond_type": bond_type,
            })
        except (ValueError, IndexError) as e:
            return {"error": f"Error parseando enlace {i+1}: {e}"}

    return {
        "mol_name": mol_name,
        "n_atoms":  n_atoms,
        "n_bonds":  n_bonds,
        "atoms":    atoms,
        "bonds":    bonds,
        "error":    None,
    }
