"""
membrane_knowledge.py — single source of truth for lipid physical constants.

No builder, validator, or adapter should hardcode these values.
All data derived from experimental literature and validated protocols.

Supported lipids: DPPC
Supported force fields: opls-aa (oplsaa_membrane.ff)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


# ── APL targets (Å²/lipid) ────────────────────────────────────────────────────
# Area Per Lipid at experimental/simulation convergence.
# Values are lipid + force field + temperature specific.
# Source: Kandt et al. 2007 (InflateGRO paper); Berger DPPC reference simulations.

_APL_TARGETS: dict[tuple[str, str, float], float] = {
    ("DPPC", "opls-aa", 298.0): 62.0,
    ("DPPC", "opls-aa", 309.65): 64.0,   # estimated from thermal expansion ~0.08 Å²/K above 298K
    ("DPPC", "opls-aa", 323.0): 68.0,    # above main phase transition; liquid-disordered phase
}

# Tolerance used during shrink-loop convergence check (Å²)
APL_CONVERGENCE_TOLERANCE: float = 2.0

# Deflation factor per shrink iteration (< 1.0 = compress)
SHRINK_DEFLATION_FACTOR: float = 0.95

# Initial inflation factor for InflateGRO (> 1.0 = expand lipid XY)
_INFLATION_FACTORS: dict[str, float] = {
    "single_pass_tm":  4.0,   # 1 TM helix, e.g. single-pass receptor fragment
    "multi_pass_tm":   5.0,   # GPCRs, channels (>4 TM helices)
    "beta_barrel":     4.5,   # beta-barrel membrane proteins
    "default":         4.0,
}

# Max shrink iterations before aborting with a non-convergence error
SHRINK_MAX_ITERATIONS: int = 200


# ── Lipid residue names per force field ───────────────────────────────────────
# GROMACS .gro / topology residue name; tool arguments use this name (e.g. InflateGRO -name)

_LIPID_RESIDUE_NAMES: dict[tuple[str, str], str] = {
    ("DPPC", "opls-aa"): "DPP",
    ("POPC", "opls-aa"): "POP",
    ("POPE", "opls-aa"): "POE",
    ("DPPC", "charmm36"): "DPPC",
    ("POPC", "charmm36"): "POPC",
}


# ── Atom name mappings per force field ────────────────────────────────────────
# Used by water_deletor and bilayer Z-range detection.
# ref_atom: headgroup phosphate oxygen (marks top/bottom of bilayer)
# middle_atom: tail carbon (marks bilayer core center)
# These are OPLS-AA DPPC atom names from oplsaa_membrane.ff

@dataclass(frozen=True)
class LipidAtomNames:
    headgroup_ref: str    # used as -ref in water_deletor; splits top/bottom leaflets
    tail_middle: str      # used as -middle in water_deletor; marks hydrophobic core


_LIPID_ATOM_NAMES: dict[tuple[str, str], LipidAtomNames] = {
    ("DPPC", "opls-aa"): LipidAtomNames(headgroup_ref="O33", tail_middle="C50"),
    # charmm36 DPPC uses different naming convention
    ("DPPC", "charmm36"): LipidAtomNames(headgroup_ref="O13", tail_middle="C218"),
}


# ── Standard bilayer geometries ───────────────────────────────────────────────
# Pre-built equilibrated bilayer files available in the project.
# box_xy in nm. n_lipids = total lipids (both leaflets).

@dataclass(frozen=True)
class BilayerGeometry:
    filename: str
    lipid: str
    n_lipids: int
    box_x_nm: float
    box_y_nm: float
    box_z_nm: float    # bilayer thickness including water layers


BILAYER_LIBRARY: dict[str, BilayerGeometry] = {
    "dppc512": BilayerGeometry(
        filename="dppc512_whole.gro",
        lipid="DPPC",
        n_lipids=512,
        box_x_nm=12.8368,
        box_y_nm=12.8870,
        box_z_nm=12.0,
    ),
    "dppc128": BilayerGeometry(
        filename="dppc128.gro",
        lipid="DPPC",
        n_lipids=128,
        box_x_nm=6.44,
        box_y_nm=6.44,
        box_z_nm=12.0,
    ),
}


# ── Equilibration protocol defaults ───────────────────────────────────────────
# MDP parameter defaults for membrane systems.
# These differ from standard protein-in-water defaults.

@dataclass
class MembraneEquilibrationDefaults:
    # NVT
    nvt_nsteps: int           # steps at dt=0.002 → ps
    nvt_dt: float
    nvt_tcoupl: str
    nvt_tc_grps: str          # "system" not "Protein Non-Protein"
    nvt_tau_t: float
    nvt_ref_t: float
    nvt_gen_vel: bool
    # NPT
    npt_nsteps: int
    npt_dt: float
    npt_pcoupl: str
    npt_pcoupltype: str       # must be semiisotropic
    npt_tau_p: float
    npt_ref_p_xy: float
    npt_ref_p_z: float
    npt_compressibility: float
    # Production
    prod_nsteps: int
    prod_dt: float            # 0.001 for OPLS-AA lipids
    prod_tcoupl: str
    prod_pcoupl: str
    prod_pcoupltype: str
    prod_tau_p: float
    prod_ref_p_xy: float
    prod_ref_p_z: float
    prod_constraints: str     # h-bonds (not all-bonds) in production


MEMBRANE_EQUILIBRATION_DEFAULTS = MembraneEquilibrationDefaults(
    # NVT — 100 ps
    nvt_nsteps=25000,
    nvt_dt=0.002,
    nvt_tcoupl="V-rescale",
    nvt_tc_grps="system",
    nvt_tau_t=0.1,
    nvt_ref_t=298.0,          # overridden by state.environment.temperature_K
    nvt_gen_vel=True,
    # NPT — 1 ns
    npt_nsteps=150000,
    npt_dt=0.002,
    npt_pcoupl="Berendsen",
    npt_pcoupltype="semiisotropic",
    npt_tau_p=5.0,
    npt_ref_p_xy=0.5,
    npt_ref_p_z=0.5,
    npt_compressibility=4.5e-5,
    # Production — 500 ns default
    prod_nsteps=50000000,
    prod_dt=0.001,            # OPLS-AA lipid stability: must not exceed 0.001
    prod_tcoupl="Nose-Hoover",
    prod_pcoupl="Parrinello-Rahman",
    prod_pcoupltype="semiisotropic",
    prod_tau_p=2.0,
    prod_ref_p_xy=1.0,
    prod_ref_p_z=1.0,
    prod_constraints="h-bonds",
)

# Shrink-loop minimization uses stricter emtol than standard EM
SHRINK_LOOP_EMTOL: float = 1000.0

# Force constants for strong position restraints (kJ/mol/nm²)
STRONG_POSRES_FC: tuple[int, int, int] = (100000, 100000, 100000)


# ── Public API ────────────────────────────────────────────────────────────────

def apl_target(lipid: str, forcefield: str, temperature_K: float) -> float:
    """Return the APL convergence target (Å²) for a lipid/FF/temperature combo.

    Rounds temperature to nearest reference point if exact match not found.
    Raises KeyError if lipid+FF combination is unknown.
    """
    key = (lipid.upper(), _normalise_ff(forcefield), float(temperature_K))
    if key in _APL_TARGETS:
        return _APL_TARGETS[key]

    # Find closest temperature for same lipid+FF
    candidates = [
        (abs(t - temperature_K), apl)
        for (lip, ff, t), apl in _APL_TARGETS.items()
        if lip == lipid.upper() and ff == _normalise_ff(forcefield)
    ]
    if not candidates:
        raise KeyError(f"Unknown lipid/FF combination: {lipid!r} + {forcefield!r}")
    _, best_apl = min(candidates, key=lambda x: x[0])
    return best_apl


def lipid_residue_name(lipid: str, forcefield: str) -> str:
    """Return the GROMACS residue name for a lipid in a given force field."""
    key = (lipid.upper(), _normalise_ff(forcefield))
    if key not in _LIPID_RESIDUE_NAMES:
        raise KeyError(f"Unknown lipid residue name for {lipid!r} + {forcefield!r}")
    return _LIPID_RESIDUE_NAMES[key]


def lipid_atom_names(lipid: str, forcefield: str) -> LipidAtomNames:
    """Return headgroup and tail atom names for water_deletor and Z-range detection."""
    key = (lipid.upper(), _normalise_ff(forcefield))
    if key not in _LIPID_ATOM_NAMES:
        raise KeyError(f"Unknown atom names for {lipid!r} + {forcefield!r}")
    return _LIPID_ATOM_NAMES[key]


def inflation_factor(protein_topology: str = "default") -> float:
    """Return the InflateGRO initial scale factor for a protein topology class."""
    return _INFLATION_FACTORS.get(protein_topology, _INFLATION_FACTORS["default"])


def bilayer_for_box(target_x_nm: float, target_y_nm: float, lipid: str = "DPPC") -> Optional[BilayerGeometry]:
    """Return the best pre-built bilayer geometry matching the target box XY dimensions.

    Returns None if no bilayer is close enough (within 0.5 nm on each axis).
    """
    for geom in BILAYER_LIBRARY.values():
        if (
            geom.lipid == lipid.upper()
            and abs(geom.box_x_nm - target_x_nm) < 0.5
            and abs(geom.box_y_nm - target_y_nm) < 0.5
        ):
            return geom
    return None


def _normalise_ff(forcefield: str) -> str:
    """Normalise force field name to canonical form used as dict key."""
    return forcefield.lower().replace("_", "-")
