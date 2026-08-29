"""Residue-name classification."""
from __future__ import annotations

import pytest

from analysis.campaign.structure.residue_classes import (
    AMINO_ACIDS, classify_resname, is_membrane_resname, is_polymer_resname,
    is_solvent_resname,
)

_STANDARD_AA = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]


@pytest.mark.parametrize("aa", _STANDARD_AA)
def test_standard_amino_acids(aa):
    assert classify_resname(aa) == "amino_acid"


@pytest.mark.parametrize("aa", ["HID", "HIE", "HIP", "HSD", "HSE", "CYX", "CYM", "ACE", "NME"])
def test_amino_acid_variants(aa):
    assert classify_resname(aa) == "amino_acid"


@pytest.mark.parametrize("nt", ["DA", "DT", "DG", "DC", "RG", "RA", "A", "U", "G", "C"])
def test_nucleotides(nt):
    assert classify_resname(nt) == "nucleotide"


@pytest.mark.parametrize("w", ["SOL", "HOH", "TIP3", "TIP4P", "SPC", "SPCE", "WAT", "H2O"])
def test_water(w):
    assert classify_resname(w) == "water"


@pytest.mark.parametrize("ion", ["NA", "CL", "MG", "K", "CA", "ZN", "SOD", "CLA"])
def test_ions(ion):
    assert classify_resname(ion) == "ion"


@pytest.mark.parametrize("lip", ["DPPC", "POPC", "POP", "CHOL", "DOPE", "CHL1"])
def test_lipids(lip):
    assert classify_resname(lip) == "lipid"


@pytest.mark.parametrize("cof", ["ATP", "HEM", "NAD", "FAD", "GTP"])
def test_cofactors(cof):
    assert classify_resname(cof) == "cofactor"


@pytest.mark.parametrize("other", ["LIG", "XYZ", "UNK", "MOL", ""])
def test_unknown_is_other(other):
    assert classify_resname(other) == "other"


def test_case_and_whitespace_insensitive():
    assert classify_resname(" ala ") == "amino_acid"
    assert classify_resname("sol") == "water"


def test_helper_predicates_consistent():
    for name in list(AMINO_ACIDS)[:5]:
        assert is_polymer_resname(name)
        assert not is_solvent_resname(name)
        assert not is_membrane_resname(name)
    assert is_solvent_resname("SOL") and is_solvent_resname("NA")
    assert not is_solvent_resname("DPPC")
    assert is_membrane_resname("POPC")
    assert not is_membrane_resname("SOL")
    assert not is_polymer_resname("SOL")
    assert not is_polymer_resname("DPPC")
