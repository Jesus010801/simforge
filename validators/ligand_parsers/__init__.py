# validators/ligand_parsers/__init__.py
"""
Parsers modulares para archivos de ligandos.

Importar directamente el parser que necesites:
    from validators.ligand_parsers.sdf_parser import parse_sdf
    from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

Contrato compartido de salida:
    {
        "mol_name": str,
        "n_atoms":  int,
        "n_bonds":  int,
        "atoms":    list[dict],   # index, element, x, y, z, charge
        "bonds":    list[dict],   # atom1, atom2, bond_type (1/2/3/4)
        "error":    str | None,
    }
"""

from validators.ligand_parsers.sdf_parser import parse_sdf
from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

__all__ = ["parse_sdf", "parse_pdb_ligand"]
