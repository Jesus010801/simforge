# core/bilayer_geometry.py
"""
Per-lipid geometric parameters for box sizing and bilayer embedding.

Used by MatchBoxBuilder at compile time; constants are baked into
box_match_helper.py so the standalone script needs no SimForge imports.

Sources:
  DPPC thickness: Kandt et al. 2007; Pabst et al. 2010 (SAXS)
  POPC thickness: Pan et al. 2012
  POPE thickness: Rappolt et al. 2003
  DMPC thickness: Nagle & Tristram-Nagle 2000
  APL values: membrane_knowledge.DPPC@opls-aa + analogues
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class LipidBoxParams:
    bilayer_thickness_nm: float     # headgroup-to-headgroup bilayer thickness
    apl_A2: float                   # area per lipid (Å²) at reference temperature
    lateral_padding_nm: float       # XY padding added around protein footprint
    solvent_z_padding_nm: float     # water layer on each Z face (above EC, below IC)
    protein_coverage_warn: float    # warn if protein XY / box XY exceeds this fraction


LIPID_BOX_PARAMS: dict[str, LipidBoxParams] = {
    "DPPC": LipidBoxParams(3.8, 64.0, 2.0, 2.5, 0.65),
    "POPC": LipidBoxParams(3.7, 65.0, 2.0, 2.5, 0.65),
    "POPE": LipidBoxParams(3.6, 60.0, 2.0, 2.5, 0.65),
    "DMPC": LipidBoxParams(3.5, 58.0, 2.0, 2.5, 0.65),
}

FALLBACK_PARAMS = LipidBoxParams(3.8, 64.0, 2.0, 2.5, 0.65)


def get_lipid_params(lipid: str) -> tuple[LipidBoxParams, bool]:
    """Return (LipidBoxParams, fallback_used).  fallback_used=True when lipid unknown."""
    key = lipid.upper()
    if key in LIPID_BOX_PARAMS:
        return LIPID_BOX_PARAMS[key], False
    return FALLBACK_PARAMS, True
