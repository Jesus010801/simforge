"""Structure & topology parsers."""
from __future__ import annotations

from pathlib import Path

import pytest

from analysis.campaign.structure.parsers import (
    parse_gro_structure, parse_pdb_structure, parse_structure, parse_top_molecules,
)
from tests.analysis.campaign.conftest import pdb_atom, write_pdb, write_top


def test_parse_pdb_recovers_chains(minimal_pdb_receptor_peptide):
    model = parse_pdb_structure(minimal_pdb_receptor_peptide)
    assert model.source_format == "pdb"
    assert model.has_chain_ids is True
    by_id = {c.chain_id: c for c in model.chains}
    assert set(by_id) == {"A", "B"}
    assert by_id["A"].residue_count == 120
    assert by_id["B"].residue_count == 25
    assert by_id["A"].polymer_fraction() == pytest.approx(1.0)
    assert by_id["A"].dominant_class() == "amino_acid"
    hist = by_id["A"].class_histogram()
    assert hist["amino_acid"] == 120


def test_parse_pdb_multi_model_only_first(tmp_path):
    lines = ["MODEL        1"]
    lines += [pdb_atom(i, "CA", "ALA", "A", i, 1.0, 1.0, float(i)) for i in range(1, 11)]
    lines += ["ENDMDL", "MODEL        2"]
    lines += [pdb_atom(i, "CA", "GLY", "A", i, 1.0, 1.0, float(i)) for i in range(1, 11)]
    lines += ["ENDMDL", "END"]
    p = tmp_path / "multi.pdb"
    p.write_text("\n".join(lines) + "\n")

    model = parse_pdb_structure(p)
    assert model.chains[0].residue_count == 10
    # frame 2 (GLY) must not have been merged in
    assert model.chains[0].resname_set == {"ALA"}
    assert any("only MODEL 1" in w for w in model.warnings)


def test_parse_gro_reconstructs_single_polymer_segment(minimal_gro_solvated_protein):
    model = parse_gro_structure(minimal_gro_solvated_protein)
    assert model.source_format == "gro"
    assert model.has_chain_ids is False
    poly = [c for c in model.chains if c.chain_id != "_"]
    assert len(poly) == 1
    assert poly[0].residue_count == 50
    nonpoly = [c for c in model.chains if c.chain_id == "_"]
    assert len(nonpoly) == 1
    assert nonpoly[0].residue_count > 1000     # the SOL + ions bucket
    assert any("no chain IDs" in w for w in model.warnings)


def test_parse_structure_dispatch(minimal_pdb_single_chain, minimal_gro_solvated_protein):
    assert parse_structure(minimal_pdb_single_chain).source_format == "pdb"
    assert parse_structure(minimal_gro_solvated_protein).source_format == "gro"
    with pytest.raises(ValueError):
        parse_structure(Path("x.mol2"))


def test_parse_top_molecules_order_and_counts(minimal_top):
    tm = parse_top_molecules(minimal_top)
    assert tm.molecules == [
        ("Protein_chain_A", 1), ("Peptide", 1), ("SOL", 2000), ("NA", 10),
    ]
    assert tm.names() == ["Protein_chain_A", "Peptide", "SOL", "NA"]
    # the "[ moleculetype ]" nrexcl line ("Protein_chain_A   3") must be ignored
    assert ("Protein_chain_A", 3) not in tm.molecules


def test_parse_top_molecules_missing_section_warns(tmp_path):
    p = tmp_path / "nomol.top"
    p.write_text("[ moleculetype ]\nProtein 3\n\n[ atoms ]\n 1 N 1 ALA N 1\n")
    tm = parse_top_molecules(p)
    assert tm.molecules == []
    assert any("no [ molecules ]" in w for w in tm.warnings)


def test_parse_top_molecules_ignores_comments(tmp_path):
    p = write_top(tmp_path / "c.top", [("Protein", 2), ("SOL", 900)])
    tm = parse_top_molecules(p)
    assert tm.molecules == [("Protein", 2), ("SOL", 900)]
