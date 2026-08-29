"""Residue-name classification for molecular-component inference.

SimForge has ~6 near-identical lipid/solvent resname sets scattered across
``core/`` and ``validators/``.  This module consolidates them for the campaign
layer and adds the amino-acid / nucleotide sets needed to tell polymers apart
from solvent.  It imports the canonical lipid/solvent sets from
``core.topology_models`` and only *extends* them, so the two stay compatible.
"""
from __future__ import annotations

from core.topology_models import _LIPID_RESNAMES, _SOLVENT_ION_RESNAMES

# Standard amino acids + common protonation/terminal variants seen in GRO/PDB.
AMINO_ACIDS: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation / tautomer / disulfide variants
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH",
    "LYN", "ARN", "TYM", "MSE", "SEC", "PYL",
    # terminal caps
    "ACE", "NME", "NMA", "NH2", "FOR",
})

NUCLEOTIDES: frozenset[str] = frozenset({
    "DA", "DT", "DG", "DC", "DU", "A", "U", "G", "C", "T", "I",
    "RA", "RU", "RG", "RC", "ADE", "THY", "GUA", "CYT", "URA",
    "DA5", "DA3", "DT5", "DT3", "DG5", "DG3", "DC5", "DC3",
})

WATER_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5", "TIP3P", "TIP4P", "TIP5P",
    "SPC", "SPCE", "T3P", "T4P", "OPC", "OPC3", "H2O",
})

ION_RESNAMES: frozenset[str] = frozenset({
    "NA", "CL", "K", "MG", "CA", "ZN", "SOD", "CLA", "POT", "MG2", "CAL",
    "LI", "RB", "CS", "BR", "IOD", "FE", "MN", "CU", "NI", "CO", "F",
    "NA+", "CL-", "K+", "IB",
})

# Lipids beyond the core set (Slipids / CHARMM-GUI names commonly encountered).
_EXTRA_LIPIDS: frozenset[str] = frozenset({
    "POP", "POE", "DPP", "DOPE", "DOPS", "DOPG", "DPPE", "DPPG", "DPPS",
    "PSM", "SSM", "DSPC", "DAPC", "DUPC", "SAPI", "SOPS", "PIP2", "PIP3",
    "CHL1", "CHOL", "ERG", "SITO", "STIG", "TRIO", "LPPC",
})
LIPID_RESNAMES: frozenset[str] = frozenset(_LIPID_RESNAMES) | _EXTRA_LIPIDS

# Everything the core module already treated as "not protein".
SOLVENT_LIKE: frozenset[str] = (
    frozenset(_SOLVENT_ION_RESNAMES) | WATER_RESNAMES | ION_RESNAMES
)

# Common small-molecule / cofactor codes that are neither polymer nor solvent.
_COMMON_COFACTORS: frozenset[str] = frozenset({
    "ATP", "ADP", "AMP", "GTP", "GDP", "GNP", "GSP", "NAD", "NAP", "NDP",
    "FAD", "FMN", "HEM", "HEC", "PLP", "SAM", "SAH", "COA", "ACO", "TPP",
    "RET", "BCL", "CLA",
})


def classify_resname(resname: str) -> str:
    """Return a coarse class for a residue name.

    One of: ``amino_acid``, ``nucleotide``, ``water``, ``ion``, ``lipid``,
    ``cofactor``, ``other``.
    """
    r = (resname or "").strip().upper()
    if r in AMINO_ACIDS:
        return "amino_acid"
    if r in NUCLEOTIDES:
        return "nucleotide"
    if r in WATER_RESNAMES:
        return "water"
    if r in ION_RESNAMES:
        return "ion"
    if r in LIPID_RESNAMES:
        return "lipid"
    if r in _COMMON_COFACTORS:
        return "cofactor"
    return "other"


def is_membrane_resname(resname: str) -> bool:
    return classify_resname(resname) == "lipid"


def is_solvent_resname(resname: str) -> bool:
    return classify_resname(resname) in ("water", "ion")


def is_polymer_resname(resname: str) -> bool:
    return classify_resname(resname) in ("amino_acid", "nucleotide")
