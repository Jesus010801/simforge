"""
Regression tests for the element-normalization integration bug: a real
campaign run against the A6 fixture (atom 16 named "CL", PDB element column
"C") produced a final ligand_for_ligpargen.pdb that STILL showed 19
carbons/0 chlorine and a manifest claiming
``element_assignment_confidence: high`` / ``element_corrections: []`` --
i.e. the campaign path silently skipped normalization even though the
underlying ``ligand.chemical_perception.normalize_ligand_elements`` helper
worked correctly in isolation.

Root cause: ``CampaignManifest.validate_paths()`` (the check that decides
whether ``prepare_campaign()`` may resume from an existing manifest rather
than regenerate) only checked that referenced files existed, never that
their content still reflected current element-normalization logic. Once a
campaign directory had been prepared once -- including, critically, by an
earlier code revision, or via any process that reintroduced an
un-normalized element into ``ligand_for_ligpargen.pdb`` -- every subsequent
``batch-prepare`` silently trusted the stale artifacts and manifest
forever.

These tests exercise the REAL high-level entrypoints
(``ligand.campaign.prepare_campaign`` and
``ligand.prepare.extract_ligand_from_complex``), not just the
``normalize_ligand_elements`` helper in isolation, and cover both the
Open Babel and RDKit hydrogenation backends (mocked, per this repo's
existing convention -- see ligand/test_integrate.py) so the fix cannot
depend on which backend a given environment happens to select.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from unittest.mock import patch

import pytest

from ligand.campaign import CampaignManifest, prepare_campaign
from ligand.chemical_perception import ChemicalPreparationResult
from ligand.prepare import extract_ligand_from_complex

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligands" / "a6"
A6_DOCKED_BUGGY = FIXTURE_DIR / "A6_docked_buggy_element.pdb"
A6_HEAVY_ONLY_BUGGY = FIXTURE_DIR / "A6_docked_heavy_only_buggy_element.pdb"

_RDKIT_AVAIL = "ligand.hydrogenation._rdkit_available"
_OBABEL_AVAIL = "ligand.hydrogenation._obabel_available"
_SUBPROC_RUN = "ligand.hydrogenation.subprocess.run"

_EXPECTED_FORMULA_NO_H = {"C": 18, "O": 3, "CL": 1}
_EXPECTED_FORMULA_WITH_H = {"C": 18, "O": 3, "CL": 1, "H": 15}


def _pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}"
    )


def _protein_lines() -> list[str]:
    return [_pdb_line("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")]


def _formula(pdb_path: Path) -> dict:
    lines = [l for l in pdb_path.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    return dict(Counter(l[76:78].strip().upper() for l in lines))


def _mock_obabel_success(cmd, *args, **kwargs):
    """side_effect for subprocess.run: copies input PDB to output, adding
    nothing -- exercises whatever the RDKit-side normalization already
    wrote, so this test genuinely observes what obabel *received*."""
    from unittest.mock import MagicMock

    if "-O" in cmd:
        out_path = Path(cmd[cmd.index("-O") + 1])
        in_path = Path(cmd[cmd.index("-ipdb") + 1]) if "-ipdb" in cmd else None
        # obabel adds explicit hydrogens; approximate by copying input as-is
        # plus 15 synthetic H atoms so downstream H-count logic is exercised
        # without needing a real obabel binary.
        src_lines = [l for l in in_path.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
        h_lines = [
            _pdb_line("HETATM", 100 + i, f"H{i}", "LIG", "A", 1, float(i), 0.0, 0.0, "H")
            for i in range(15)
        ]
        out_path.write_text("\n".join(src_lines + h_lines) + "\nEND\n")
        return MagicMock(returncode=0, stderr="")
    return MagicMock(returncode=0, stderr="")


# ═══════════════════════════════════════════════════════════════════════════
# Single-complex ligand prepare -- real execution path
# ═══════════════════════════════════════════════════════════════════════════

class TestSingleComplexPrepareElementOrdering:
    def test_openbabel_backend_receives_corrected_element(self, tmp_path):
        """Requirement: hydrogenation must never receive the
        pre-normalization PDB when corrections were identified. Verified by
        capturing exactly what obabel's subprocess call was invoked with.

        Uses the heavy-atom-only (no explicit H) fixture so hydrogenation is
        actually triggered (a ligand that already has a plausible H count
        classifies as "unknown"/no backend call needed -- see
        ligand.prepare._determine_hydrogenation_status)."""
        complex_pdb = tmp_path / "complex.pdb"
        ligand_lines = [l for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")]
        complex_pdb.write_text("\n".join(_protein_lines() + ["TER"] + ligand_lines + ["END"]) + "\n")

        captured_input_text = {}

        def _capture_obabel(cmd, *args, **kwargs):
            in_path = Path(cmd[cmd.index("-ipdb") + 1])
            captured_input_text["text"] = in_path.read_text()
            return _mock_obabel_success(cmd, *args, **kwargs)

        with (
            patch(_RDKIT_AVAIL, return_value=False),
            patch(_OBABEL_AVAIL, return_value=True),
            patch(_SUBPROC_RUN, side_effect=_capture_obabel),
        ):
            result = extract_ligand_from_complex(complex_pdb, "LIG", tmp_path / "out")

        assert result.success, result.error
        assert result.element_corrections, "normalization must have run and found the CL/C conflict"

        # The exact bytes obabel's subprocess received must already show Cl.
        obabel_input_lines = [
            l for l in captured_input_text["text"].splitlines() if l.startswith("HETATM")
        ]
        cl_line = next(l for l in obabel_input_lines if l[12:16].strip() == "CL")
        assert cl_line[76:78].strip().upper() == "CL", (
            "Open Babel received the UN-corrected element -- hydrogenation ran "
            "before (or without) normalization"
        )

    def test_rdkit_backend_receives_corrected_element(self, tmp_path):
        """Same guarantee for the RDKit backend, via a mocked
        add_hydrogens_with_confidence so this runs without RDKit installed --
        captures the exact file content the backend was handed."""
        complex_pdb = tmp_path / "complex.pdb"
        ligand_lines = [l for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")]
        complex_pdb.write_text("\n".join(_protein_lines() + ["TER"] + ligand_lines + ["END"]) + "\n")

        captured_input_text = {}

        def _fake_add_hydrogens(input_pdb, output_pdb):
            captured_input_text["text"] = Path(input_pdb).read_text()
            heavy = [l for l in captured_input_text["text"].splitlines() if l.startswith("HETATM")]
            h = [
                _pdb_line("HETATM", 100 + i, f"H{i}", "LIG", "A", 1, float(i), 0.0, 0.0, "H")
                for i in range(15)
            ]
            Path(output_pdb).write_text("\n".join(heavy + h) + "\nEND\n")
            return ChemicalPreparationResult(
                success=True, output_path=Path(output_pdb), n_hydrogen_atoms=15,
                formal_charge=0, method="geometry_valence_optimization",
                connectivity_confidence="high", bond_order_confidence="high",
                hydrogenation_confidence="high",
            )

        with (
            patch(_RDKIT_AVAIL, return_value=True),
            patch("ligand.chemical_perception.add_hydrogens_with_confidence", side_effect=_fake_add_hydrogens),
        ):
            result = extract_ligand_from_complex(complex_pdb, "LIG", tmp_path / "out")

        assert result.success, result.error
        assert result.element_corrections
        assert result.hydrogenation_backend == "rdkit"

        rdkit_input_lines = [l for l in captured_input_text["text"].splitlines() if l.startswith("HETATM")]
        cl_line = next(l for l in rdkit_input_lines if l[12:16].strip() == "CL")
        assert cl_line[76:78].strip().upper() == "CL", (
            "RDKit backend received the UN-corrected element -- hydrogenation ran "
            "before (or without) normalization"
        )

    def test_final_ligpargen_file_and_report_agree_with_corrected_element(self, tmp_path):
        complex_pdb = tmp_path / "complex.pdb"
        ligand_lines = [l for l in A6_DOCKED_BUGGY.read_text().splitlines() if l.startswith("HETATM")]
        complex_pdb.write_text("\n".join(_protein_lines() + ["TER"] + ligand_lines + ["END"]) + "\n")

        with (
            patch(_RDKIT_AVAIL, return_value=False),
            patch(_OBABEL_AVAIL, return_value=True),
            patch(_SUBPROC_RUN, side_effect=_mock_obabel_success),
        ):
            result = extract_ligand_from_complex(complex_pdb, "LIG", tmp_path / "out")

        assert result.success, result.error
        submitted = result.recommended_ligpargen_input or result.ligand_pdb
        assert _formula(submitted) == _EXPECTED_FORMULA_WITH_H

        import yaml
        report = yaml.safe_load(result.report_path.read_text())
        assert report["element_corrections"]
        assert report["element_corrections"][0]["atom_name"] == "CL"
        assert report["element_corrections"][0]["corrected_element"] == "Cl"
        assert report["element_assignment_confidence"] == "high"


# ═══════════════════════════════════════════════════════════════════════════
# Campaign prepare -- real execution path (the exact command from the bug)
# ═══════════════════════════════════════════════════════════════════════════

def _campaign_root(tmp_path) -> Path:
    lines = _protein_lines() + ["TER"] + [
        l for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
    ] + ["END"]
    (tmp_path / "AA").mkdir(parents=True, exist_ok=True)
    (tmp_path / "AA" / "AA.pdb").write_text("\n".join(lines) + "\n")
    return tmp_path


class TestCampaignPrepareElementOrdering:
    """Exercises the actual command from the bug report:
    `simforge ligand batch-prepare . --ligand-resname LIG`."""

    def test_campaign_invokes_element_normalization(self, tmp_path):
        root = _campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")

        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections, "campaign path must invoke normalize_ligand_elements"
        assert lig.element_corrections[0]["atom_name"] == "CL"
        assert lig.element_corrections[0]["corrected_element"] == "Cl"
        assert lig.element_assignment_confidence == "high"

    def test_generated_ligpargen_file_contains_chlorine(self, tmp_path):
        root = _campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")

        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        param_pdb = root / lig.ligand_for_ligpargen
        formula = _formula(param_pdb)
        assert formula.get("CL") == 1, f"ligand_for_ligpargen.pdb still missing Cl: {formula}"
        assert formula == _EXPECTED_FORMULA_NO_H

    def test_correction_propagates_into_campaign_manifest_provenance(self, tmp_path):
        root = _campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")

        assert result.success, result.error
        loaded = CampaignManifest.load(result.manifest_path)
        lig = loaded.ligands["LIG"]
        assert lig.element_corrections
        assert lig.element_corrections[0] == {
            "atom_index": 16,
            "atom_name": "CL",
            "original_element": "C",
            "corrected_element": "Cl",
            "reason": "ligand_atom_name_element_conflict",
        }

    def test_heavy_atom_coordinates_unchanged_by_normalization(self, tmp_path):
        root = _campaign_root(tmp_path)
        original_coords = {
            int(l[6:11]): (round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
            for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
        }

        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error

        pose_path = root / "AA" / "simforge" / "ligand_pose.pdb"
        pose_coords = {
            int(l[6:11]): (round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
            for l in pose_path.read_text().splitlines() if l.startswith("HETATM")
        }
        assert pose_coords == original_coords

    def test_openbabel_backend_campaign_path(self, tmp_path):
        root = _campaign_root(tmp_path)
        with (
            patch(_RDKIT_AVAIL, return_value=False),
            patch(_OBABEL_AVAIL, return_value=True),
            patch(_SUBPROC_RUN, side_effect=_mock_obabel_success),
        ):
            result = prepare_campaign(root, "LIG")

        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections
        param_pdb = root / lig.ligand_for_ligpargen
        assert _formula(param_pdb) == _EXPECTED_FORMULA_WITH_H
        assert lig.hydrogenation_backend == "openbabel"

    def test_rdkit_backend_campaign_path(self, tmp_path):
        root = _campaign_root(tmp_path)

        def _fake_add_hydrogens(input_pdb, output_pdb):
            heavy = [l for l in Path(input_pdb).read_text().splitlines() if l.startswith("HETATM")]
            h = [
                _pdb_line("HETATM", 100 + i, f"H{i}", "LIG", "A", 1, float(i), 0.0, 0.0, "H")
                for i in range(15)
            ]
            Path(output_pdb).write_text("\n".join(heavy + h) + "\nEND\n")
            return ChemicalPreparationResult(
                success=True, output_path=Path(output_pdb), n_hydrogen_atoms=15,
                formal_charge=0, method="geometry_valence_optimization",
                connectivity_confidence="high", bond_order_confidence="high",
                hydrogenation_confidence="high",
            )

        with (
            patch(_RDKIT_AVAIL, return_value=True),
            patch("ligand.chemical_perception.add_hydrogens_with_confidence", side_effect=_fake_add_hydrogens),
        ):
            result = prepare_campaign(root, "LIG")

        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections
        param_pdb = root / lig.ligand_for_ligpargen
        assert _formula(param_pdb) == _EXPECTED_FORMULA_WITH_H
        assert lig.hydrogenation_backend == "rdkit"


# ═══════════════════════════════════════════════════════════════════════════
# Resume/staleness detection -- the actual root cause
# ═══════════════════════════════════════════════════════════════════════════

class TestCampaignResumeStalenessDetection:
    def test_resume_rejects_manifest_with_stale_unnormalized_artifact(self, tmp_path):
        """Reproduces the exact reported bug: a previously-prepared campaign
        whose on-disk ligand_for_ligpargen.pdb still contains an
        un-normalized element must NOT be silently resumed -- the whole
        point of resume is to avoid reprocessing valid state, not to trust
        stale/pre-fix state forever."""
        root = _campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            first = prepare_campaign(root, "LIG")
        assert first.success, first.error

        # Simulate a stale artifact -- e.g. written by a prior code revision,
        # or externally modified -- by corrupting the element back.
        lig = first.manifest.ligands["LIG"]
        param_pdb = root / lig.ligand_for_ligpargen
        lines = param_pdb.read_text().splitlines()
        for i, l in enumerate(lines):
            if l.startswith("HETATM") and l[12:16].strip() == "CL":
                lines[i] = l[:76] + " C"
        param_pdb.write_text("\n".join(lines) + "\n")

        second = prepare_campaign(root, "LIG")
        assert second.success is False
        assert "un-normalized" in second.error.lower()
        assert "CL" in second.error

    def test_resume_succeeds_when_artifact_is_already_normalized(self, tmp_path):
        """Sanity check: the new validation does not false-positive on a
        genuinely already-normalized campaign."""
        root = _campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            first = prepare_campaign(root, "LIG")
        assert first.success, first.error

        second = prepare_campaign(root, "LIG")
        assert second.success, second.error
        assert second.resumed is True


# ═══════════════════════════════════════════════════════════════════════════
# Short/loosely-formatted line serialization -- the actual reported bug
# ═══════════════════════════════════════════════════════════════════════════
#
# The real production docking PDB stops right after occupancy/temperature
# factor (no element/segID columns at all -- a very common real-world
# "PDB-like" output convention). normalize_ligand_elements previously
# skipped any such line outright (len(line) < 78), so the CL/C conflict was
# never detected NOR was a valid element field ever written downstream: the
# final ligand_for_ligpargen.pdb had a blank element column for every atom
# on such a line, and the report showed element_corrections: [] despite the
# contradiction being present in the input.

def _short_format_campaign_root(tmp_path) -> Path:
    """Same A6 heavy-only, element-buggy pose as _campaign_root, but every
    HETATM line truncated to end right after occupancy/temperature-factor
    (66 columns) -- matching the exact format reported in production, with
    NO element/segID columns present at all."""
    short_lines = [
        l[:66] for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
    ]
    lines = _protein_lines() + ["TER"] + short_lines + ["END"]
    (tmp_path / "AA").mkdir(parents=True, exist_ok=True)
    (tmp_path / "AA" / "AA.pdb").write_text("\n".join(lines) + "\n")
    return tmp_path


class TestCampaignShortLinePdbElementSerialization:
    """End-to-end regression through the actual `prepare_campaign()` path
    (the same one `simforge ligand batch-prepare` calls), starting from a
    docking ligand whose atom 16 has name=CL with NO element field at all
    (not even the wrong "C") -- the exact production symptom reported."""

    def test_internal_element_normalized_to_cl(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections
        assert lig.element_corrections[0]["atom_name"] == "CL"
        assert lig.element_corrections[0]["corrected_element"] == "Cl"

    def test_output_pdb_element_column_is_populated(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        param_pdb = root / lig.ligand_for_ligpargen
        for line in param_pdb.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                assert len(line) >= 78 and line[76:78].strip(), f"blank element column: {line!r}"

    def test_heavy_atom_composition_matches_ground_truth(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        param_pdb = root / lig.ligand_for_ligpargen
        assert _formula(param_pdb) == _EXPECTED_FORMULA_NO_H  # {"C": 18, "O": 3, "CL": 1}

    def test_element_corrections_records_c_to_cl(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections[0]["reason"] == "ligand_atom_name_element_conflict"
        assert lig.element_corrections[0]["corrected_element"] == "Cl"

    def test_element_assignment_confidence_is_high(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        assert result.manifest.ligands["LIG"].element_assignment_confidence == "high"

    def test_heavy_atom_coordinates_unchanged(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        original_coords = {
            int(l[6:11]): (round(float(l[30:38]), 3), round(float(l[38:46]), 3))
            for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
        }
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        pose_path = root / "AA" / "simforge" / "ligand_pose.pdb"
        pose_coords = {
            int(l[6:11]): (round(float(l[30:38]), 3), round(float(l[38:46]), 3))
            for l in pose_path.read_text().splitlines() if l.startswith("HETATM")
        }
        assert pose_coords == original_coords

    def test_same_behavior_with_openbabel_backend(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)
        with (
            patch(_RDKIT_AVAIL, return_value=False),
            patch(_OBABEL_AVAIL, return_value=True),
            patch(_SUBPROC_RUN, side_effect=_mock_obabel_success),
        ):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections
        param_pdb = root / lig.ligand_for_ligpargen
        formula = _formula(param_pdb)
        assert formula["CL"] == 1 and formula["C"] == 18 and formula["O"] == 3

    def test_same_behavior_with_rdkit_backend(self, tmp_path):
        root = _short_format_campaign_root(tmp_path)

        def _fake_add_hydrogens(input_pdb, output_pdb):
            heavy = [l for l in Path(input_pdb).read_text().splitlines() if l.startswith("HETATM")]
            h = [
                _pdb_line("HETATM", 100 + i, f"H{i}", "LIG", "A", 1, float(i), 0.0, 0.0, "H")
                for i in range(15)
            ]
            Path(output_pdb).write_text("\n".join(heavy + h) + "\nEND\n")
            return ChemicalPreparationResult(
                success=True, output_path=Path(output_pdb), n_hydrogen_atoms=15,
                formal_charge=0, method="geometry_valence_optimization",
                connectivity_confidence="high", bond_order_confidence="high",
                hydrogenation_confidence="high",
            )

        with (
            patch(_RDKIT_AVAIL, return_value=True),
            patch("ligand.chemical_perception.add_hydrogens_with_confidence", side_effect=_fake_add_hydrogens),
        ):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_corrections
        param_pdb = root / lig.ligand_for_ligpargen
        formula = _formula(param_pdb)
        assert formula["CL"] == 1 and formula["C"] == 18 and formula["O"] == 3

    def test_single_complex_prepare_path_matches_campaign_path(self, tmp_path):
        """Campaign and non-campaign preparation must not diverge: same
        short-line-format input through extract_ligand_from_complex()
        directly produces the same corrected composition."""
        complex_pdb = tmp_path / "complex.pdb"
        short_lines = [
            l[:66] for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")
        ]
        complex_pdb.write_text("\n".join(_protein_lines() + ["TER"] + short_lines + ["END"]) + "\n")

        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = extract_ligand_from_complex(complex_pdb, "LIG", tmp_path / "out")

        assert result.success, result.error
        assert result.element_corrections
        assert result.element_corrections[0]["corrected_element"] == "Cl"
        assert _formula(result.ligand_pdb) == _EXPECTED_FORMULA_NO_H

    def test_remains_blocked_without_reference_because_bond_order_confidence_insufficient(self, tmp_path):
        """The scientific safety behavior must be untouched by this fix:
        heavy-atom-only input with no chemical reference still cannot
        establish bond orders and must remain chemical_validation_required,
        even though elements are now correctly serialized."""
        root = _short_format_campaign_root(tmp_path)
        with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
            result = prepare_campaign(root, "LIG")
        assert result.success, result.error
        lig = result.manifest.ligands["LIG"]
        assert lig.element_assignment_confidence == "high"
        assert lig.parameterization_status == "chemical_validation_required"
        assert any("bond_order_confidence" in r for r in lig.readiness_block_reasons)
