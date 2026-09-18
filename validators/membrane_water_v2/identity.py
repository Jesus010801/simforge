"""File identities are atom-record identities, never wrapped GRO numbers."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from utils.gro_parser import parse_gro, GroFile
from .models import SystemGeometry, WaterMolecule
from .cell import PeriodicCell

PROTEIN = 1
LIPID = 2
OTHER = 3
ION = 4
WATERS = frozenset("SOL HOH WAT TIP3 TIP4 TIP5 TIP3P TIP4P SPC SPCE".split())
AMINO = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP HSD HSE HSP CYX ASH GLH LYN ACE NME".split())
LIPIDS = frozenset("DPP DPPC POP POPC POPE POPG POPS CHOL PALM OLEO PALC SM CER".split())
IONS = frozenset("NA CL SOD CLA K POT CA MG ZN".split())
RADII = {"C": .170, "N": .155, "O": .152, "P": .180, "S": .180, "F": .147,
         "CL": .175, "BR": .185, "I": .198, "NA": .227, "K": .275, "CA": .231, "MG": .173, "ZN": .139}


@dataclass
class MolecularInput:
    geometry: SystemGeometry
    gro: GroFile
    waters: list[WaterMolecule]
    source: Path
    lines: list[str]


def load_system(path, lipid_resnames=None, protein_resnames=None):
    path = Path(path)
    gro = parse_gro(path)
    lines = path.read_text().splitlines(keepends=True)
    lipids = set(lipid_resnames or LIPIDS)
    proteins = set(protein_resnames or AMINO)
    waters, groups, current = [], [], []
    last_key = None
    for i, a in enumerate(gro.atoms):
        key = (a.residue_number, a.residue_name)
        oxygen = a.residue_name in WATERS and a.atom_name.upper().startswith("O")
        if current and (key != last_key or oxygen):
            groups.append(current)
            current = []
        current.append(i)
        last_key = key
    if current:
        groups.append(current)
    coords, radii, species, names, ids, known = [], [], [], [], [], []
    warnings = []
    obstacle_group = 0
    for group in groups:
        a = gro.atoms[group[0]]
        if a.residue_name in WATERS:
            ox = [i for i in group if gro.atoms[i].atom_name.upper().startswith("O")]
            nh = sum(gro.atoms[i].atom_name.upper().lstrip("0123456789").startswith("H") for i in group)
            valid = len(ox) == 1 and nh == 2 and len(group) in (3, 4, 5)
            if not valid:
                warnings.append(f"Malformed water instance {len(waters)} preserved unless input is repaired")
            waters.append(WaterMolecule(len(waters), tuple(group), ox[0] if ox else group[0], a.residue_number, a.residue_name, valid))
            continue
        role = LIPID if a.residue_name in lipids else PROTEIN if a.residue_name in proteins else ION if a.residue_name in IONS else OTHER
        for i in group:
            atom = gro.atoms[i]
            name = atom.atom_name.upper().lstrip("0123456789")
            if name.startswith("H") or name.startswith(("MW", "VS")):
                continue
            element = {"SOD": "NA", "CLA": "CL", "POT": "K"}.get(a.residue_name, a.residue_name) if role == ION else name[:2] if name[:2] in {"CL", "BR"} else name[:1]
            coords.append((atom.x, atom.y, atom.z))
            radii.append(RADII.get(element, .170))
            known.append(element in RADII)
            species.append(role)
            names.append(atom.atom_name)
            ids.append(obstacle_group)
        obstacle_group += 1
    # Water presence is deliberately excluded from geometry warnings and IDs.
    geometry = SystemGeometry(PeriodicCell.from_gro(gro.box), np.asarray(coords, float).reshape(-1, 3),
                              np.asarray(radii), np.asarray(species, np.uint8), tuple(names),
                              np.asarray(ids), np.asarray(known, bool), [])
    if not all(known):
        geometry.warnings.append("Unknown atomic radii: affected geometry is uncertain")
    return MolecularInput(geometry, gro, waters, path, lines)
