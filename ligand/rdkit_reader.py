"""
Shared RDKit molecule loading for the ligand package.

Supports .sdf / .mol (V2000 and V3000) and .pdb inputs.
RDKit is imported lazily so this module can be imported in environments
where RDKit is absent without raising ImportError at collection time —
functions that need RDKit will raise a clear error when called.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def _rdkit_chem():
    try:
        from rdkit import Chem
        return Chem
    except ImportError:
        raise ImportError(
            "RDKit is required for ligand processing. "
            "Install with: conda install -c conda-forge rdkit"
        )


def load_mol(path: str | Path):
    """
    Load a molecule from an SDF/MOL or PDB file using RDKit.

    Returns an RDKit Mol object (with explicit Hs preserved from the file),
    or raises ValueError if the file cannot be parsed.
    """
    Chem = _rdkit_chem()
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Ligand file not found: {path}")

    ext = path.suffix.lower()

    if ext in (".sdf", ".mol"):
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        mol = next(iter(supplier), None)
    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
    else:
        raise ValueError(
            f"Unsupported format '{ext}'. Supported: .sdf, .mol, .pdb"
        )

    if mol is None:
        raise ValueError(f"RDKit could not parse '{path.name}'")

    return mol


def heavy_atom_coords(mol) -> list[tuple[float, float, float]]:
    """
    Return (x, y, z) tuples for every heavy atom (non-hydrogen) in mol.

    Requires mol to have an embedded conformer.
    """
    conf = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append((pos.x, pos.y, pos.z))
    return coords


def heavy_atom_elements(mol) -> list[str]:
    """Return element symbols for every heavy atom in mol."""
    return [
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() != 1
    ]


def net_formal_charge(mol) -> int:
    """Sum of formal charges across all atoms."""
    Chem = _rdkit_chem()
    return Chem.GetFormalCharge(mol)


def has_explicit_hydrogens(mol) -> bool:
    """Return True if the molecule has at least one explicit hydrogen atom."""
    return any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())


def has_3d_conformer(mol) -> bool:
    """
    Return True if mol has a genuine 3D conformer.

    RDKit's Compute2DCoords() sets all z=0 while x/y vary, so a conformer
    with max|z| < 1e-3 Å across all atoms is classified as 2D-only.
    """
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    if n == 0:
        return False
    max_z = max(abs(conf.GetAtomPosition(i).z) for i in range(n))
    return max_z > 1e-3


def sanitized_mol_name(mol, fallback: str = "LIG") -> str:
    """
    Extract a molecule name from the '_Name' property, or return fallback.

    Truncates to 4 characters (GROMACS residue name limit).
    """
    name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    name = name or fallback
    # Keep only alphanumeric and strip to 4 chars for GROMACS compatibility
    cleaned = "".join(c for c in name if c.isalnum()).upper()
    return cleaned[:4] or fallback[:4].upper()
