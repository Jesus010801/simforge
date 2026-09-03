"""
Tests for ligand/chemical_perception.py.

Element normalization, the readiness gate, and residue-metadata rewriting
are pure Python (no RDKit) and always run. Geometry-based perception,
hydrogenation, and reference-guided mapping require RDKit >= 2022.09
(rdkit.Chem.rdDetermineBonds) and are skipped automatically when RDKit is
unavailable -- see ligand/test_pose_rewriter_heavy_atom_transfer.py and
ligand/test_preparation.py for the same convention elsewhere in this suite.

The A6 fixtures (tests/fixtures/ligands/a6/) are real data: a Gaussian-
derived, charge-annotated reference PDB (A6_con_cargas.pdb, ground truth:
C18H15ClO3, 22 heavy atoms, neutral) and two variants reproducing the exact
reported docking-tool defect (atom 16, named "CL", with its PDB element
column truncated to "C").
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ligand.chemical_perception import (
    ChemicalPerceptionResult,
    ElementCorrection,
    decide_ligpargen_readiness,
    normalize_ligand_elements,
    normalize_ligand_residue_metadata,
)

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligands" / "a6"
A6_REFERENCE = FIXTURE_DIR / "A6_reference_correct.pdb"
A6_DOCKED_BUGGY = FIXTURE_DIR / "A6_docked_buggy_element.pdb"
A6_HEAVY_ONLY_BUGGY = FIXTURE_DIR / "A6_docked_heavy_only_buggy_element.pdb"
A6_HEAVY_ONLY_CORRECT = FIXTURE_DIR / "A6_docked_heavy_only_correct_element.pdb"

_HAS_RDKIT_DETERMINE_BONDS = True
try:
    from rdkit.Chem import rdDetermineBonds  # noqa: F401
except ImportError:
    _HAS_RDKIT_DETERMINE_BONDS = False

requires_rdkit_determine_bonds = pytest.mark.skipif(
    not _HAS_RDKIT_DETERMINE_BONDS,
    reason="requires RDKit >= 2022.09 (rdkit.Chem.rdDetermineBonds)",
)


def _pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}"
    )


# ═══════════════════════════════════════════════════════════════════════════
# Element normalization (no RDKit)
# ═══════════════════════════════════════════════════════════════════════════

class TestNormalizeLigandElements:
    def test_cl_named_atom_with_wrong_c_element_is_corrected(self):
        line = _pdb_line("HETATM", 16, "CL", "LIG", "A", 1, 4.966, -0.995, 2.457, "C")
        result = normalize_ligand_elements(line + "\n")
        assert len(result.corrections) == 1
        c = result.corrections[0]
        assert c.atom_index == 16
        assert c.atom_name == "CL"
        assert c.original_element == "C"
        assert c.corrected_element == "Cl"
        assert result.pdb_text.splitlines()[0][76:78].strip().upper() == "CL"

    def test_br_named_atom_with_wrong_element_is_corrected(self):
        line = _pdb_line("HETATM", 5, "BR", "LIG", "A", 1, 1.0, 1.0, 1.0, "B")
        result = normalize_ligand_elements(line + "\n")
        assert len(result.corrections) == 1
        assert result.corrections[0].corrected_element == "Br"

    def test_f_named_atom_already_correct_no_correction(self):
        line = _pdb_line("HETATM", 3, "F", "LIG", "A", 1, 1.0, 1.0, 1.0, "F")
        result = normalize_ligand_elements(line + "\n")
        assert result.corrections == []

    def test_i_named_atom_already_correct_no_correction(self):
        line = _pdb_line("HETATM", 3, "I", "LIG", "A", 1, 1.0, 1.0, 1.0, "I")
        result = normalize_ligand_elements(line + "\n")
        assert result.corrections == []

    def test_protein_ca_atom_stays_carbon_not_calcium(self):
        """CA (alpha-carbon) is not in the halogen allowlist and must never
        be reinterpreted as calcium, whether in a protein ATOM record or a
        ligand HETATM record using generic carbon naming."""
        line = _pdb_line("ATOM", 2, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C")
        result = normalize_ligand_elements(line + "\n")
        assert result.corrections == []
        assert result.pdb_text.splitlines()[0][76:78].strip().upper() == "C"

    def test_ligand_ca_named_carbon_stays_carbon(self):
        line = _pdb_line("HETATM", 2, "CA", "LIG", "A", 1, 0.0, 0.0, 0.0, "C")
        result = normalize_ligand_elements(line + "\n")
        assert result.corrections == []

    def test_valid_explicit_elements_unchanged(self):
        lines = "\n".join([
            _pdb_line("HETATM", 1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, "C"),
            _pdb_line("HETATM", 2, "O1", "LIG", "A", 1, 1.0, 0.0, 0.0, "O"),
            _pdb_line("HETATM", 3, "N1", "LIG", "A", 1, 2.0, 0.0, 0.0, "N"),
        ])
        result = normalize_ligand_elements(lines + "\n")
        assert result.corrections == []
        assert result.pdb_text.rstrip("\n") == lines

    def test_confidence_is_high_when_corrections_applied(self):
        line = _pdb_line("HETATM", 16, "CL", "LIG", "A", 1, 0.0, 0.0, 0.0, "C")
        result = normalize_ligand_elements(line + "\n")
        assert result.confidence == "high"

    def test_multiple_corrections_all_recorded(self):
        lines = "\n".join([
            _pdb_line("HETATM", 1, "CL", "LIG", "A", 1, 0.0, 0.0, 0.0, "C"),
            _pdb_line("HETATM", 2, "BR", "LIG", "A", 1, 1.0, 0.0, 0.0, "N"),
        ])
        result = normalize_ligand_elements(lines + "\n")
        assert len(result.corrections) == 2
        assert {c.corrected_element for c in result.corrections} == {"Cl", "Br"}

    def test_real_a6_fixture_atom_16_corrected(self):
        result = normalize_ligand_elements(A6_DOCKED_BUGGY.read_text())
        assert len(result.corrections) == 1
        c = result.corrections[0]
        assert c.atom_index == 16
        assert c.atom_name == "CL"
        assert c.original_element == "C"
        assert c.corrected_element == "Cl"

    def test_real_a6_reference_fixture_has_no_corrections(self):
        """The ground-truth reference already has the correct element -- no
        correction should fire."""
        result = normalize_ligand_elements(A6_REFERENCE.read_text())
        assert result.corrections == []

    def test_coordinates_unchanged_by_normalization(self):
        original_coords = [
            (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            for l in A6_DOCKED_BUGGY.read_text().splitlines()
            if l.startswith("HETATM")
        ]
        result = normalize_ligand_elements(A6_DOCKED_BUGGY.read_text())
        corrected_coords = [
            (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            for l in result.pdb_text.splitlines()
            if l.startswith("HETATM")
        ]
        assert corrected_coords == original_coords

    # ── Short/loosely-formatted lines (real-world docking-tool output) ────────
    #
    # Many docking-tool "PDB" writers stop right after occupancy/temperature
    # factor and never write an element/segID column at all. The historical
    # bug: normalize_ligand_elements silently skipped any line shorter than
    # the full 78-column PDB width, so neither the CL/C conflict was detected
    # NOR was a valid element field ever produced downstream -- the eventual
    # ligand_for_ligpargen.pdb had a blank element column entirely.

    def test_short_line_with_halogen_conflict_is_corrected(self):
        short_line = "HETATM   16  CL  LIG A   1      30.302  28.869  42.389  0.00  0.00"
        assert len(short_line) < 78
        result = normalize_ligand_elements(short_line + "\n")
        assert len(result.corrections) == 1
        c = result.corrections[0]
        assert c.atom_name == "CL"
        assert c.corrected_element == "Cl"

    def test_short_line_output_is_properly_serialized_to_78_columns(self):
        short_line = "HETATM   16  CL  LIG A   1      30.302  28.869  42.389  0.00  0.00"
        result = normalize_ligand_elements(short_line + "\n")
        out = result.pdb_text.splitlines()[0]
        assert len(out) == 78
        assert out[76:78].strip() == "Cl"

    def test_short_line_coordinates_unchanged(self):
        short_line = "HETATM   16  CL  LIG A   1      30.302  28.869  42.389  0.00  0.00"
        result = normalize_ligand_elements(short_line + "\n")
        out = result.pdb_text.splitlines()[0]
        assert out[30:54] == short_line[30:54]

    def test_short_line_without_halogen_name_stays_uncorrected_but_gets_valid_element(self):
        """A plain carbon on a short line: no correction needed, but the
        output must still carry a valid, correctly positioned element field
        -- every atom in the file must be reliably parseable downstream,
        not only the ones that happened to need correction."""
        short_line = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  0.00  0.00"
        result = normalize_ligand_elements(short_line + "\n")
        assert result.corrections == []
        out = result.pdb_text.splitlines()[0]
        assert len(out) == 78
        assert out[76:78].strip() == "C"

    def test_entire_short_format_a6_fixture_normalizes_correctly(self):
        """The real regression: every line in the file uses the short
        (no element column) format, matching the exact reported production
        symptom."""
        short_lines = [
            l[:66] for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
        ]
        result = normalize_ligand_elements("\n".join(short_lines) + "\n")
        assert len(result.corrections) == 1
        assert result.corrections[0].atom_name == "CL"
        assert result.corrections[0].corrected_element == "Cl"

        formula = {}
        for l in result.pdb_text.splitlines():
            if not l.startswith("HETATM"):
                continue
            elem = l[76:78].strip().upper()
            assert elem, f"blank element field for line: {l!r}"
            formula[elem] = formula.get(elem, 0) + 1
        assert formula == {"C": 18, "O": 3, "CL": 1}

    def test_already_well_formed_line_byte_identical_after_normalization(self):
        """No regression for input that already follows strict PDB widths --
        rebuilding must reproduce it exactly when nothing needed correcting."""
        line = _pdb_line("HETATM", 1, "C1", "LIG", "A", 1, 1.234, 5.678, 9.012, "C")
        result = normalize_ligand_elements(line + "\n")
        assert result.corrections == []
        assert result.pdb_text.rstrip("\n") == line


# ═══════════════════════════════════════════════════════════════════════════
# Residue metadata normalization (no RDKit)
# ═══════════════════════════════════════════════════════════════════════════

class TestNormalizeLigandResidueMetadata:
    def test_unl_residue_rewritten_to_ligand_resname(self):
        pdb_text = (
            "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "HETATM    2  H1  UNL     1       0.500   0.500   0.500  1.00  0.00           H\n"
        )
        fixed = normalize_ligand_residue_metadata(pdb_text, resname="LIG", chain="A", resid=1)
        resnames = {l[17:20].strip() for l in fixed.splitlines() if l.startswith("HETATM")}
        assert resnames == {"LIG"}

    def test_blank_chain_rewritten(self):
        pdb_text = "HETATM    2  H1  UNL     1       0.500   0.500   0.500  1.00  0.00           H\n"
        fixed = normalize_ligand_residue_metadata(pdb_text, resname="LIG", chain="A", resid=1)
        chains = {l[21] for l in fixed.splitlines() if l.startswith("HETATM")}
        assert chains == {"A"}

    def test_resid_rewritten_consistently(self):
        pdb_text = (
            "HETATM    1  C1  LIG A   5       0.000   0.000   0.000  1.00  0.00           C\n"
            "HETATM    2  H1  UNL     1       0.500   0.500   0.500  1.00  0.00           H\n"
        )
        fixed = normalize_ligand_residue_metadata(pdb_text, resname="LIG", chain="A", resid=5)
        resids = {int(l[22:26]) for l in fixed.splitlines() if l.startswith("HETATM")}
        assert resids == {5}

    def test_coordinates_and_atom_names_untouched(self):
        pdb_text = "HETATM    2  H1  UNL     1       0.500   0.500   0.500  1.00  0.00           H\n"
        fixed = normalize_ligand_residue_metadata(pdb_text, resname="LIG")
        line = fixed.splitlines()[0]
        assert line[12:16] == " H1 "
        assert line[30:54] == "   0.500   0.500   0.500"

    def test_non_atom_lines_untouched(self):
        pdb_text = "REMARK some comment\nEND\n"
        fixed = normalize_ligand_residue_metadata(pdb_text, resname="LIG")
        assert fixed == pdb_text


# ═══════════════════════════════════════════════════════════════════════════
# Readiness gate (no RDKit)
# ═══════════════════════════════════════════════════════════════════════════

class TestDecideLigpargenReadiness:
    def test_all_high_and_complete_is_ready(self):
        status, reasons = decide_ligpargen_readiness(
            "high", "high", "high", "complete", "high", hydrogenation_performed=True,
        )
        assert status == "ready_for_ligpargen"
        assert reasons == []

    def test_addhs_success_alone_is_not_sufficient(self):
        """A successful hydrogenation call (hydrogenation_status='complete')
        with medium-confidence chemistry must NOT be ready when SimForge
        itself performed the hydrogenation -- this is the exact scenario
        that produced the historical mis-hydrogenation bug."""
        status, reasons = decide_ligpargen_readiness(
            element_assignment_confidence="high",
            connectivity_confidence="medium",
            bond_order_confidence="medium",
            hydrogenation_status="complete",
            hydrogenation_confidence="medium",
            hydrogenation_performed=True,
        )
        assert status == "chemical_validation_required"
        assert any("bond_order_confidence" in r for r in reasons)

    def test_incomplete_hydrogenation_status_blocks(self):
        status, reasons = decide_ligpargen_readiness(
            "high", "high", "high", "incomplete", "high", hydrogenation_performed=True,
        )
        assert status == "chemical_validation_required"
        assert any("hydrogenation_status" in r for r in reasons)

    def test_missing_hydrogenation_status_blocks(self):
        status, reasons = decide_ligpargen_readiness(
            "high", "low", "low", "missing", "unknown", hydrogenation_performed=False,
        )
        assert status == "chemical_validation_required"

    def test_element_confidence_always_requires_high(self):
        """Element normalization needs no RDKit and is unambiguous by
        construction -- the relaxed 'medium' bar for hydrogenation_performed
        =False must not apply to element assignment."""
        status, reasons = decide_ligpargen_readiness(
            element_assignment_confidence="medium",
            connectivity_confidence="medium",
            bond_order_confidence="medium",
            hydrogenation_status="complete",
            hydrogenation_confidence="medium",
            hydrogenation_performed=False,
        )
        assert status == "chemical_validation_required"
        assert any("element_assignment_confidence" in r for r in reasons)

    def test_medium_sufficient_when_hydrogenation_not_performed(self):
        """SimForge did not compute or add anything -- a heuristic-level
        assessment of an already-complete input is an acceptable bar."""
        status, reasons = decide_ligpargen_readiness(
            element_assignment_confidence="high",
            connectivity_confidence="medium",
            bond_order_confidence="medium",
            hydrogenation_status="complete",
            hydrogenation_confidence="medium",
            hydrogenation_performed=False,
        )
        assert status == "ready_for_ligpargen"

    def test_low_confidence_blocks_even_when_not_performed(self):
        status, reasons = decide_ligpargen_readiness(
            "high", "low", "low", "complete", "low", hydrogenation_performed=False,
        )
        assert status == "chemical_validation_required"

    def test_reference_guided_high_confidence_is_ready(self):
        status, reasons = decide_ligpargen_readiness(
            "high", "high", "high", "complete", "high", hydrogenation_performed=True,
        )
        assert status == "ready_for_ligpargen"


# ═══════════════════════════════════════════════════════════════════════════
# Geometry-based perception / hydrogenation / reference mapping (RDKit)
# ═══════════════════════════════════════════════════════════════════════════

@requires_rdkit_determine_bonds
class TestPerceiveChemistryA6Regression:
    """The real A6 regression: an 18-carbon, 3-oxygen, 1-chlorine, neutral
    ligand (C18H15ClO3, 22 heavy atoms) whose docking-tool output mistypes
    atom 16 ("CL") with PDB element column "C". The historical bug produced
    an approximately C19H34O3 molecule and reported it as
    hydrogenation-complete and ready_for_ligpargen."""

    def test_buggy_element_uncorrected_does_not_silently_converge_wrong(self):
        """Even feeding perceive_chemistry the UNCORRECTED buggy file (i.e.
        without running normalize_ligand_elements first) must never produce
        the wrong ~C19H34O3 structure -- either it fails to converge, or it
        converges to something that is not silently treated as correct.
        This guards the underlying mechanism, independent of the element
        normalization safety net."""
        from ligand.chemical_perception import perceive_chemistry
        from rdkit import Chem

        result = perceive_chemistry(A6_DOCKED_BUGGY)
        if result.success:
            n_c = sum(1 for a in result.mol.GetAtoms() if a.GetSymbol() == "C")
            assert n_c != 19  # never the historically-wrong carbon count

    def test_normalized_buggy_fixture_converges_to_correct_formula(self):
        from collections import Counter

        from ligand.chemical_perception import normalize_ligand_elements, perceive_chemistry
        from rdkit import Chem

        norm = normalize_ligand_elements(A6_DOCKED_BUGGY.read_text())
        tmp = FIXTURE_DIR  # read-only fixtures dir; write to a temp copy instead
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "normalized.pdb"
            p.write_text(norm.pdb_text)
            result = perceive_chemistry(p)

        assert result.success, result.warning
        assert result.connectivity_confidence == "high"
        assert result.bond_order_confidence == "high"
        assert result.used_charge == 0

        mol_h = Chem.AddHs(result.mol, addCoords=True)
        formula = Counter(a.GetSymbol() for a in mol_h.GetAtoms())
        assert dict(formula) == {"C": 18, "O": 3, "Cl": 1, "H": 15}

    def test_heavy_atom_only_no_reference_fails_to_converge(self):
        """A heavy-atom-only pose (no explicit hydrogens, the realistic
        campaign-discovery scenario) with no chemical reference: geometry
        alone cannot resolve valence completion. Must fail explicitly, not
        guess."""
        from ligand.chemical_perception import perceive_chemistry

        result = perceive_chemistry(A6_HEAVY_ONLY_CORRECT)
        assert result.success is False
        assert result.connectivity_confidence == "low"
        assert result.bond_order_confidence == "low"


@requires_rdkit_determine_bonds
class TestAddHydrogensWithConfidenceA6:
    def test_buggy_fixture_normalized_hydrogenates_correctly(self, tmp_path):
        from collections import Counter

        from ligand.chemical_perception import add_hydrogens_with_confidence, normalize_ligand_elements

        norm = normalize_ligand_elements(A6_DOCKED_BUGGY.read_text())
        input_pdb = tmp_path / "normalized.pdb"
        input_pdb.write_text(norm.pdb_text)

        result = add_hydrogens_with_confidence(input_pdb, tmp_path / "out.pdb")
        assert result.success, result.warning
        assert result.n_hydrogen_atoms == 15
        assert result.formal_charge == 0
        assert result.connectivity_confidence == "high"
        assert result.bond_order_confidence == "high"
        assert result.hydrogenation_confidence == "high"

        formula = Counter(
            l[76:78].strip().upper() for l in result.output_path.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))
        )
        assert dict(formula) == {"C": 18, "O": 3, "CL": 1, "H": 15}

    def test_heavy_atom_only_input_fails_explicitly(self, tmp_path):
        """No 34-hydrogen (or any other guessed) output is accepted for a
        heavy-atom-only pose -- hydrogenation must fail outright."""
        from ligand.chemical_perception import add_hydrogens_with_confidence

        result = add_hydrogens_with_confidence(A6_HEAVY_ONLY_CORRECT, tmp_path / "out.pdb")
        assert result.success is False
        assert not (tmp_path / "out.pdb").exists()

    def test_heavy_atom_coordinates_preserved_exactly(self, tmp_path):
        """Chemical normalization must never modify receptor-specific
        heavy-atom coordinates -- only metadata/connectivity/H positions
        may change."""
        from ligand.chemical_perception import add_hydrogens_with_confidence, normalize_ligand_elements

        original_heavy_coords = {}
        for l in A6_DOCKED_BUGGY.read_text().splitlines():
            if l.startswith("HETATM") and l[76:78].strip().upper() != "H":
                serial = int(l[6:11])
                original_heavy_coords[serial] = (
                    round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3)
                )

        norm = normalize_ligand_elements(A6_DOCKED_BUGGY.read_text())
        input_pdb = tmp_path / "normalized.pdb"
        input_pdb.write_text(norm.pdb_text)
        result = add_hydrogens_with_confidence(input_pdb, tmp_path / "out.pdb")
        assert result.success, result.warning

        out_heavy_coords = []
        for l in result.output_path.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")) and l[76:78].strip().upper() != "H":
                out_heavy_coords.append(
                    (round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
                )

        assert sorted(out_heavy_coords) == sorted(original_heavy_coords.values())


@requires_rdkit_determine_bonds
class TestMapReferenceOntoPoseA6:
    def test_chlorine_maps_by_name_despite_docked_element_bug(self):
        """For the A6 reference, the chlorine must map to the docked atom
        named 'CL' even when the docking PDB element field incorrectly
        says carbon (fed here WITHOUT prior element normalization, to
        isolate the mapping algorithm's own robustness)."""
        from ligand.chemical_perception import map_reference_onto_pose

        result = map_reference_onto_pose(A6_REFERENCE, A6_HEAVY_ONLY_BUGGY)
        assert result.success, result.error
        assert any("CL" in w for w in result.warnings)

        cl_atoms = [a for a in result.mol.GetAtoms() if a.GetSymbol() == "Cl"]
        assert len(cl_atoms) == 1

        conf = result.mol.GetConformer()
        pos = conf.GetAtomPosition(cl_atoms[0].GetIdx())

        docked_cl_line = next(
            l for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l[12:16].strip() == "CL"
        )
        expected = (
            float(docked_cl_line[30:38]), float(docked_cl_line[38:46]), float(docked_cl_line[46:54])
        )
        assert (round(pos.x, 3), round(pos.y, 3), round(pos.z, 3)) == pytest.approx(expected, abs=1e-3)

    def test_reference_guided_produces_correct_formula(self):
        from collections import Counter

        from ligand.chemical_perception import map_reference_onto_pose

        result = map_reference_onto_pose(A6_REFERENCE, A6_HEAVY_ONLY_BUGGY)
        assert result.success, result.error
        formula = Counter(a.GetSymbol() for a in result.mol.GetAtoms())
        assert dict(formula) == {"C": 18, "O": 3, "Cl": 1, "H": 15}
        assert result.formal_charge == 0

    def test_heavy_atom_count_mismatch_fails_explicitly(self, tmp_path):
        from ligand.chemical_perception import map_reference_onto_pose

        lines = [
            l for l in A6_HEAVY_ONLY_CORRECT.read_text().splitlines()
            if l.startswith("HETATM")
        ][:-1] + ["END"]
        truncated = tmp_path / "truncated.pdb"
        truncated.write_text("\n".join(lines) + "\n")

        result = map_reference_onto_pose(A6_REFERENCE, truncated)
        assert result.success is False
        assert "mismatch" in result.error.lower()

    def test_missing_reference_fails_explicitly(self, tmp_path):
        from ligand.chemical_perception import map_reference_onto_pose

        result = map_reference_onto_pose(tmp_path / "nonexistent.pdb", A6_HEAVY_ONLY_CORRECT)
        assert result.success is False
        assert result.error is not None
