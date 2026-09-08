"""Safe, non executing command plans for legacy campaign protocols.

The command plans deliberately contain argv and stdin separately.  This keeps
them suitable for JSON export and prevents a renderer from accidentally
executing shell text.
"""
from __future__ import annotations

import shlex
from pathlib import Path


_ANALYSES = (
    ("protein_rmsd", "gmx rms", ("protein", "protein"), "rmsd_protein.xvg"),
    ("ligand_rmsd", "gmx rms", ("protein", "ligand"), "rmsd_ligand.xvg"),
    ("protein_rmsf", "gmx rmsf", ("protein",), "rmsf_protein.xvg"),
    ("sasa", "gmx sasa", ("protein",), "sasa_protein.xvg"),
    ("hbonds", "gmx hbond", ("protein", "ligand"), "hbonds_protein_ligand.xvg"),
    ("active_site_mindist", "gmx mindist", ("ligand", "active_site"), "mindist_lig_active.xvg"),
    ("active_site_contacts", "gmx mindist", ("ligand", "active_site"), "contacts_lig_active.xvg"),
    ("catalytic_com_distance", "gmx distance", ("ligand", "catalytic_site"), "dist_lig_catalytic.xvg"),
    ("interaction_energy", "gmx energy", ("protein", "ligand"), "interaction_energy.xvg"),
    ("mmpbsa", "g_mmpbsa run", ("protein", "ligand"), "binding_energy.xvg"),
)


def _file(directory: Path, patterns: tuple[str, ...]) -> str | None:
    for pattern in patterns:
        candidates = sorted(p for p in directory.glob(pattern) if p.is_file())
        if candidates:
            return str(candidates[0])
    return None


def _shell(argv: list[str], stdin: str) -> str:
    """Render for humans; callers should use ``argv``/``stdin`` to execute."""
    return f"printf %s {shlex.quote(stdin)} | " + " ".join(shlex.quote(x) for x in argv)


def build_xanthone_commands(rec, legacy: dict, *, gmx: str = "gmx") -> dict:
    """Build a deterministic dry-run report; this function never invokes GROMACS."""
    directory = Path(legacy["directory"])
    selections = legacy.get("selections", {})
    existing = legacy.get("existing_results", {})
    tpr = _file(directory, ("md.tpr", "*.tpr"))
    trajectory = _file(directory, ("md.xtc", "*.xtc", "md.trr", "*.trr"))
    index = legacy.get("index")
    protein_id = selections.get("protein", {}).get("group_id")
    preprocessing = [
        {"step": "center", "argv": [gmx, "trjconv", "-f", trajectory or "md.xtc", "-s", tpr or "md.tpr", "-pbc", "res", "-ur", "compact", "-center", "-o", str(directory / "mdcenter.xtc")], "stdin": str(protein_id) + "\n" if protein_id is not None else "", "selection_role": "protein", "selection_confidence": 0.7, "status": "PLANNED" if protein_id is not None else "REVIEW_REQUIRED"},
        {"step": "fit", "argv": [gmx, "trjconv", "-f", str(directory / "mdcenter.xtc"), "-s", tpr or "md.tpr", "-fit", "rot+trans", "-o", str(directory / "mdfit.xtc")], "stdin": str(protein_id) + "\n" if protein_id is not None else "", "selection_role": "protein", "selection_confidence": 0.7, "status": "PLANNED" if protein_id is not None else "REVIEW_REQUIRED"},
    ]
    rows = []
    for analysis, tool, roles, output in _ANALYSES:
        missing = [r for r in roles if not selections.get(r, {}).get("group_id") is not None]
        if legacy.get("profile", {}).get("value") != "xanthone_short":
            status, reason = "REVIEW_REQUIRED", "system is not safely classified as xanthone_short"
        elif legacy.get("status") == "review_required":
            status, reason = "REVIEW_REQUIRED", "legacy system requires review before execution"
        elif analysis == "interaction_energy":
            status, reason = "REVIEW_REQUIRED", "energy-term indices and compatible EDR terms require manual confirmation"
        elif analysis in existing:
            status, reason = "SKIP_EXISTING", "validated existing result"
        elif analysis == "mmpbsa":
            status, reason = "REVIEW_REQUIRED", "historical g_mmpbsa parameters are incomplete; no guessing"
        elif missing or not tpr or not trajectory or not index:
            status = "REVIEW_REQUIRED"
            reason = "missing " + ", ".join(missing or [x for x, v in (("md.tpr", tpr), ("trajectory", trajectory), ("index.ndx", index)) if not v])
        else:
            status, reason = "PLANNED", "dry-run command only"
        argv = [gmx, tool.split()[1], "-s", tpr or "md.tpr", "-f", trajectory or "md.xtc", "-n", index or "index.ndx"]
        output_flag = "-o"
        if analysis == "active_site_mindist": output_flag = "-od"
        elif analysis == "active_site_contacts":
            output_flag = "-on"
        elif analysis == "catalytic_com_distance": output_flag = "-ox"
        elif analysis == "hbonds": output_flag = "-num"
        argv += [output_flag, str(directory / output)]
        if analysis == "mmpbsa":
            argv = [gmx, "run", "-f", str(directory / "mdfit.xtc"), "-s", tpr or "md.tpr", "-i", str(directory / "mmpbsa.mdp"), "-n", index or "index.ndx", "-unit1", 'group "Protein"', "-unit2", 'group "LIG"']
        if analysis == "protein_rmsf": argv += ["-res"]
        if analysis == "active_site_contacts": argv += ["-d", "0.6"]
        if analysis == "sasa": argv += ["-tu", "ns"]
        stdin = "\n".join(str(selections.get(r, {}).get("group_id", "")) for r in roles) + "\n"
        rows.append({"analysis": analysis, "status": status, "roles": list(roles),
                     "group_ids": {r: selections.get(r, {}).get("group_id") for r in roles},
                     "argv": argv, "stdin": stdin, "command": _shell(argv, stdin), "reason": reason})
    return {"system": legacy.get("system", directory.name), "profile": "xanthone_short",
            "directory": str(directory), "preprocessing": preprocessing,
            "analysis_window": legacy.get("analysis_window", {"status": "REVIEW_REQUIRED", "reason": "historical thesis interval not established for this system"}),
            "commands": rows,
            "execution": "disabled; dry-run only"}

