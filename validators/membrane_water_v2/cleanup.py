"""Coordinate-preserving GRO filtering by unique water-instance identity."""
from pathlib import Path
from .identity import load_system
from .policy import decide
from utils.gro_parser import parse_gro


def clean_gro(input_path, output_path, mappings, policy="conservative"):
    system = load_system(input_path)
    decisions = [decide(m, policy) for m in mappings]
    remove = {d.molecule_id for d in decisions if d.action == "REMOVE"}
    remove_atoms = {i for w in system.waters if w.uid in remove for i in w.atom_indices}
    lines = system.lines
    atom_lines = lines[2:2+len(system.gro.atoms)]
    retained = [line for i, line in enumerate(atom_lines) if i not in remove_atoms]
    target = Path(output_path); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(lines[0] + f"{len(retained)}\n" + "".join(retained) + "".join(lines[2+len(atom_lines):]))
    result = parse_gro(target)
    if len(result.atoms) != len(system.gro.atoms) - len(remove_atoms):
        raise RuntimeError("GRO atom-count validation failed")
    return decisions

