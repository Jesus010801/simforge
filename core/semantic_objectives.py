"""
core/semantic_objectives.py — Semantic objective system.

Provides:
  - CANONICAL_OBJECTIVES  — authoritative set of recognized objectives
  - OBJECTIVE_ALIASES     — informal strings → canonical list
  - SIMULATION_PRESETS    — named profiles that expand to objectives + hints
  - normalize_objective() — main normalization function
  - suggest_objectives()  — fuzzy suggestions for unknowns
"""
from __future__ import annotations

import difflib


# ── Canonical objectives ───────────────────────────────────────────────────────

CANONICAL_OBJECTIVES: dict[str, str] = {
    "stability":               "Protein structural stability (RMSD, RMSF, Rg)",
    "binding":                 "Ligand binding and free energy estimation",
    "competitive_binding":     "Competitive ligand vs substrate at active site",
    "membrane_insertion":      "Membrane insertion / translocation mechanism",
    "membrane_perturbation":   "Bilayer perturbation, lipid order, thickness",
    "conformational_sampling": "Conformational landscape exploration",
    "aggregation":             "Protein aggregation propensity",
    "permeability":            "Membrane permeability / ion channel dynamics",
    "allosteric_effect":       "Allosteric communication and signal propagation",
    "active_site_dynamics":    "Active site residue dynamics and flexibility",
    "active_site_stability":   "Active site structural conservation during MD",
}

# ── Alias → canonical list ─────────────────────────────────────────────────────

_ALIASES: dict[str, list[str]] = {
    # stability
    "protein_stability":            ["stability"],
    "structural_stability":         ["stability"],
    "thermal_stability":            ["stability"],
    "thermostability":              ["stability"],
    "folding_stability":            ["stability"],
    "protein_folding":              ["conformational_sampling", "stability"],
    "unfolding":                    ["conformational_sampling", "stability"],
    "denaturation":                 ["conformational_sampling", "stability"],

    # binding
    "binding_affinity":             ["binding"],
    "free_energy":                  ["binding"],
    "drug_binding":                 ["binding"],
    "ligand_binding":               ["binding"],
    "docking":                      ["binding"],
    "mmgbsa":                       ["binding"],
    "mmpbsa":                       ["binding"],
    "hit_optimization":             ["binding"],

    # competitive
    "competitive_inhibition":       ["competitive_binding"],
    "inhibitor_binding":            ["competitive_binding", "binding"],

    # membrane
    "membrane_protein_dynamics":    ["membrane_perturbation", "stability"],
    "membrane_dynamics":            ["membrane_perturbation", "conformational_sampling"],
    "protein_membrane_interaction": ["membrane_perturbation", "membrane_insertion"],
    "lipid_protein_interaction":    ["membrane_perturbation"],
    "membrane_fluidity":            ["membrane_perturbation"],
    "lipid_dynamics":               ["membrane_perturbation"],
    "bilayer_dynamics":             ["membrane_perturbation"],
    "lipid_bilayer":                ["membrane_perturbation"],

    # permeability / channels
    "pore_formation":               ["membrane_perturbation", "permeability"],
    "ion_channel":                  ["membrane_perturbation", "permeability"],
    "transporter":                  ["membrane_perturbation", "permeability"],
    "channel_gating":               ["membrane_perturbation", "permeability"],

    # conformational
    "protein_motion":               ["conformational_sampling"],
    "conformational_flexibility":   ["conformational_sampling"],
    "conformational_change":        ["conformational_sampling"],
    "loop_dynamics":                ["conformational_sampling"],
    "domain_motion":                ["conformational_sampling"],
    "idr":                          ["conformational_sampling"],
    "intrinsically_disordered":     ["conformational_sampling"],
    "enhanced_sampling":            ["conformational_sampling"],

    # allosteric
    "allostery":                    ["allosteric_effect"],
    "allosteric_modulation":        ["allosteric_effect"],
    "signal_propagation":           ["allosteric_effect"],
    "network_analysis":             ["allosteric_effect"],

    # active site
    "enzyme_dynamics":              ["active_site_dynamics", "active_site_stability"],
    "catalytic_site":               ["active_site_dynamics"],
    "active_site":                  ["active_site_dynamics"],
    "catalysis":                    ["active_site_dynamics", "active_site_stability"],
}


# ── Scientific presets ─────────────────────────────────────────────────────────

SIMULATION_PRESETS: dict[str, dict] = {
    "membrane_protein": {
        "description": "Integral membrane protein — bilayer + semiisotropic coupling",
        "objectives":  ["membrane_perturbation", "stability"],
        "hints": {
            "membrane_required":       True,
            "semiisotropic_coupling":  True,
            "conservative_timestep":   True,
            "lipid_aware_restraints":  True,
            "membrane_equilibration":  True,
        },
    },
    "soluble_protein": {
        "description": "Soluble globular protein stability",
        "objectives":  ["stability"],
        "hints": {},
    },
    "protein_ligand_binding": {
        "description": "Protein-ligand binding free energy / pose stability",
        "objectives":  ["binding", "stability"],
        "hints": {
            "extended_equilibration": True,
        },
    },
    "idr_sampling": {
        "description": "Intrinsically disordered protein / region",
        "objectives":  ["conformational_sampling"],
        "hints": {
            "enhanced_sampling":  True,
            "long_production":    True,
        },
    },
    "competitive_inhibition": {
        "description": "Competitive inhibitor displacing substrate",
        "objectives":  ["competitive_binding", "binding"],
        "hints": {},
    },
    "membrane_lipid_dynamics": {
        "description": "Pure bilayer / lipid-only dynamics",
        "objectives":  ["membrane_perturbation"],
        "hints": {
            "membrane_required":      True,
            "semiisotropic_coupling": True,
        },
    },
    "allosteric_modulation": {
        "description": "Allosteric communication and signal propagation",
        "objectives":  ["allosteric_effect", "conformational_sampling"],
        "hints": {
            "extended_production": True,
        },
    },
    "ion_channel": {
        "description": "Ion channel / transporter permeability",
        "objectives":  ["permeability", "membrane_perturbation"],
        "hints": {
            "membrane_required":      True,
            "semiisotropic_coupling": True,
        },
    },
}


# ── Normalization ──────────────────────────────────────────────────────────────

def _canonical_key(text: str) -> str:
    """Normalize input text to a lookup key."""
    return text.lower().strip().replace(" ", "_").replace("-", "_")


def normalize_objective(text: str) -> tuple[list[str], str | None]:
    """
    Map an informal objective string to a list of canonical objectives.

    Returns:
        (canonicals, note)
        canonicals  — list of canonical objective names (may be empty if unknown)
        note        — human-readable normalization note, or None if already canonical
    """
    key = _canonical_key(text)

    # Already canonical
    if key in CANONICAL_OBJECTIVES:
        return [key], None

    # Exact alias match
    if key in _ALIASES:
        return _ALIASES[key], f"'{text}' → {_ALIASES[key]}"

    # Fuzzy match against all known keys
    all_keys = list(CANONICAL_OBJECTIVES.keys()) + list(_ALIASES.keys())
    matches = difflib.get_close_matches(key, all_keys, n=3, cutoff=0.65)

    if matches:
        top = matches[0]
        if top in CANONICAL_OBJECTIVES:
            return [top], f"fuzzy '{text}' → '{top}'"
        if top in _ALIASES:
            return _ALIASES[top], f"fuzzy '{text}' → {_ALIASES[top]}"

    return [], None


def suggest_objectives(text: str) -> list[str]:
    """
    Return the top 3 canonical objective names that are closest to *text*.
    Used for error messages / autocomplete suggestions.
    """
    key = _canonical_key(text)
    all_keys = list(CANONICAL_OBJECTIVES.keys()) + list(_ALIASES.keys())
    raw = difflib.get_close_matches(key, all_keys, n=5, cutoff=0.35)

    # Resolve aliases → canonical, deduplicate
    seen: set[str] = set()
    suggestions: list[str] = []
    for m in raw:
        targets = CANONICAL_OBJECTIVES.keys() if m in CANONICAL_OBJECTIVES else _ALIASES.get(m, [])
        for t in ([m] if m in CANONICAL_OBJECTIVES else targets):
            if t not in seen:
                seen.add(t)
                suggestions.append(t)
            if len(suggestions) == 3:
                break
        if len(suggestions) == 3:
            break

    return suggestions
