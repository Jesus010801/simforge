"""
validators/membrane_water/constants.py

Single source of truth for the membrane-water topology engine: element van
der Waals radii, solvent/lipid/ion resname sets, and the region/category
enums shared by every module in this package and by the two legacy adapters
(validators/membrane_protein_spatial_classifier.py, validators/pore_hydration.py).

No other module in this package (or its adapters) should redefine these.
"""
from __future__ import annotations

import re

# ── Solvent / ion resnames ───────────────────────────────────────────────────

SOLVENT_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
})
ION_RESNAMES: frozenset[str] = frozenset({
    "NA", "CL", "SOD", "MG", "K", "CA", "ZN", "POT", "CLA",
})
NON_PROTEIN_SOLVENT_LIKE: frozenset[str] = SOLVENT_RESNAMES | ION_RESNAMES

# ── Lipid resnames (fallback set; core.membrane_knowledge is consulted first
#    for the canonical single headgroup/tail atom name of a registered
#    lipid+forcefield pair — see species.py) ─────────────────────────────────

DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset({
    "DPP", "DPPC", "POP", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER",
})

# Broader headgroup/tail atom-name heuristic, used as a fallback for lipids
# not registered in core.membrane_knowledge.lipid_atom_names(). Matches the
# set the legacy classifier used, kept here as the one copy.
HEADGROUP_ATOM_NAMES: frozenset[str] = frozenset({"P", "P1", "N", "NZ", "O33", "O11", "O13"})
TAIL_ATOM_NAMES: frozenset[str] = frozenset({
    "C50", "C49", "C48", "C47", "C46", "C45", "C34", "C35", "C218",
})

# ── Element van der Waals radii (Bondi 1964 + standard extensions), nm ──────

VDW_RADII_NM: dict[str, float] = {
    "H": 0.120, "C": 0.170, "N": 0.155, "O": 0.152, "P": 0.180, "S": 0.180,
    "F": 0.147, "CL": 0.175, "BR": 0.185, "I": 0.198,
    "NA": 0.227, "MG": 0.173, "K": 0.275, "CA": 0.231, "ZN": 0.139,
}
DEFAULT_VDW_RADIUS_NM = 0.170  # carbon-like fallback for unrecognized elements

WATER_OXYGEN_VDW_RADIUS_NM = VDW_RADII_NM["O"]
DEFAULT_SOLVENT_PROBE_RADIUS_NM = 0.140  # standard water-probe SASA radius

# Hard steric-clash test scales the raw VDW radius sum down (standard
# structural-biology "bump check" convention, e.g. ~0.7-0.8x is typical in
# PyMOL/Rosetta-style clash checks) rather than using it unscaled. A real
# equilibrated MD snapshot has hydrogen-bonded/electrostatically-stabilized
# water sitting closer to polar/charged heavy atoms than a pure repulsive
# hard-sphere model would allow; an unscaled sum over-flags legitimate,
# stable hydration-shell water as "clashing".
CLASH_VDW_SCALE = 0.8

_TWO_LETTER_ELEMENTS = frozenset({"CL", "BR", "NA", "MG", "ZN"})
_LEADING_DIGIT_RE = re.compile(r"^\d+")


def infer_element(atom_name: str, resname: str) -> str:
    """Best-effort element inference from a GRO atom/residue name.

    Ions are single-atom residues whose atom name usually equals the element
    (e.g. resname "NA", atomname "NA"); everything else follows the common
    PDB/GRO convention of an optional leading digit followed by 1-2 letters.
    """
    name = _LEADING_DIGIT_RE.sub("", atom_name.strip().upper())
    resname_u = resname.strip().upper()
    if resname_u in ION_RESNAMES:
        if resname_u in VDW_RADII_NM:
            return resname_u
        if resname_u == "SOD":
            return "NA"
        if resname_u == "CLA":
            return "CL"
        if resname_u == "POT":
            return "K"
    if len(name) >= 2 and name[:2] in _TWO_LETTER_ELEMENTS:
        return name[:2]
    if name:
        return name[0]
    return "C"


def vdw_radius_nm(atom_name: str, resname: str) -> float:
    element = infer_element(atom_name, resname)
    return VDW_RADII_NM.get(element, DEFAULT_VDW_RADIUS_NM)


# ── Region / water-category enums (the ONE definition; adapters map onto
#    these, they never redefine them) ────────────────────────────────────────

REGION_POOL = "bulk"
REGION_PORE = "pore"
REGION_VESTIBULE = "vestibule"
REGION_ISOLATED_CAVITY = "isolated_cavity"
REGION_MEMBRANE_DEFECT = "membrane_defect"

WATER_CATEGORIES: tuple[str, ...] = (
    "bulk",
    "pore",
    "vestibule",
    "isolated_cavity",
    "membrane_defect",
    "steric_clash",
    "ambiguous",
)

CLEANUP_MODE_CONSERVATIVE = "conservative"
CLEANUP_MODE_AGGRESSIVE = "aggressive"

# Categories removed under each cleanup mode (before per-region confidence
# gating, which can still force preservation of an ambiguous water).
CONSERVATIVE_REMOVE_CATEGORIES: frozenset[str] = frozenset({"membrane_defect", "steric_clash"})
AGGRESSIVE_PRESERVE_CATEGORIES: frozenset[str] = frozenset({"pore", "vestibule", "isolated_cavity"})

# ── Tunable defaults ─────────────────────────────────────────────────────────

DEFAULT_GRID_SPACING_NM = 0.20
DEFAULT_PROTEIN_WALL_THRESHOLD = 0.55   # pore vs. membrane_defect (two-sided)
DEFAULT_VESTIBULE_WALL_THRESHOLD = 0.40  # vestibule vs. membrane_defect (one-sided)
DEFAULT_MIN_TM_CONFIDENCE = 0.5
DEFAULT_MIN_REGION_CONFIDENCE = 0.5
