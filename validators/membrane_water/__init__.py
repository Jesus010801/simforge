"""Public membrane-water adapters; v2 owns all classification and cleanup semantics."""
from pathlib import Path
from validators.membrane_water_v2 import classify_membrane_water, run_cleanup
from validators.membrane_water_v2.identity import load_system
from validators.membrane_water_v2.models import ExclusionModel
from validators.membrane_water_v2.clearance import ClearanceIndex

CLEANUP_MODE_CONSERVATIVE = "conservative"
CLEANUP_MODE_AGGRESSIVE = "aggressive"

def determine_removal_set(report, cleanup_mode=None, remove_membrane_core_water=True, remove_isolated_internal_water=False):
    mode = cleanup_mode or report.get("cleanup_mode", CLEANUP_MODE_CONSERVATIVE)
    remove = set()
    for item in report.get("water_mappings", []):
        cls = item["classification"]; uid = item["molecule_uid"]
        if item.get("hard_clash"):
            remove.add(uid)
        elif remove_membrane_core_water and cls == "membrane_defect_candidate" and not item.get("uncertainty"):
            remove.add(uid)
        elif mode == CLEANUP_MODE_AGGRESSIVE and cls == "steric_clash" and not item.get("uncertainty"):
            remove.add(uid)
        elif remove_isolated_internal_water and cls in {"internal_cavity", "buried_pocket"} and not item.get("uncertainty"):
            remove.add(uid)
    return remove

def remove_waters_and_write(gro_in, gro_out, remove_resids, lipid_resnames=None):
    system = load_system(gro_in, lipid_resnames); remove = set(remove_resids)
    remove_atoms = {i for water in system.waters if water.uid in remove for i in water.atom_indices}
    atoms = system.lines[2:2 + len(system.gro.atoms)]
    kept = [line for i, line in enumerate(atoms) if i not in remove_atoms]
    target = Path(gro_out); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(system.lines[0] + f"{len(kept)}\n" + "".join(kept) + "".join(system.lines[2 + len(atoms):]))
    check = load_system(target, lipid_resnames)
    if len(check.gro.atoms) != len(system.gro.atoms) - len(remove_atoms):
        raise RuntimeError("GRO atom-count validation failed")
    return {"n_water_molecules_removed": len(remove & {w.uid for w in system.waters}),
            "n_water_atoms_removed": len(remove_atoms), "n_atoms_before": len(system.gro.atoms),
            "n_atoms_after": len(check.gro.atoms), "output_gro_path": str(target)}

def validate_cleanup_artifacts(gro_in, gro_out, pre_report=None, post_report=None, removed_resids=(), lipid_resnames=None):
    source = load_system(gro_in, lipid_resnames); output = load_system(gro_out, lipid_resnames)
    remove = set(removed_resids); errors = []; warnings = []
    atom_lines = source.lines[2:2 + len(source.gro.atoms)]
    remove_atoms = {i for w in source.waters if w.uid in remove for i in w.atom_indices}
    expected = [line for i, line in enumerate(atom_lines) if i not in remove_atoms]
    output_lines = output.lines[2:2 + len(output.gro.atoms)]
    if output_lines != expected: errors.append("retained atom records changed or water-removal set was not applied atomically")
    if any(not w.valid for w in output.waters): errors.append("output contains incomplete or malformed water molecules")
    if output.waters:
        import numpy as np
        xyz = np.array([[output.gro.atoms[w.oxygen_index].x, output.gro.atoms[w.oxygen_index].y, output.gro.atoms[w.oxygen_index].z] for w in output.waters])
        _, margin, _ = ClearanceIndex(source.geometry).query(xyz, ExclusionModel())
        if np.any(margin < 0): errors.append(f"{int(np.sum(margin < 0))} retained water oxygen(s) overlap certified excluded volume")
    return {"passed": not errors, "errors": errors, "warnings": warnings}

def build_human_report(report, removed_resids):
    return "Membrane-water region-atlas v2 cleanup\n" + f"Regions: {len(report.get("regions", []))}\nRemoved water instance IDs: {len(removed_resids)}\n"

def write_human_report(report, removed_resids, path):
    Path(path).write_text(build_human_report(report, removed_resids))
