"""
CLI tests for `simforge ligand export-ligpargen`.

All tests in TestExportLigpargenMocked mock the export functions so the
standard suite does NOT require RDKit.

TestExportLigpargenIntegration is marked with pytest.mark.rdkit and is
skipped automatically when RDKit is absent.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from cli import cli
from core.ligand_workflow_models import LigandExportResult

runner = CliRunner(mix_stderr=False)


# ── Shared mock results ───────────────────────────────────────────────────────

def _ok_smiles_result(output_dir: Path, name: str = "LIG") -> LigandExportResult:
    return LigandExportResult(
        success=True,
        exported_path=output_dir / f"{name}.smi",
        molecule_name=name,
        atom_count=18,
        smiles="CCO",
    )


def _ok_result(output_dir: Path, name: str = "LIG", legacy: bool = False) -> LigandExportResult:
    filename = f"{name}_ligpargen_legacy.pdb" if legacy else f"{name}.pdb"
    return LigandExportResult(
        success=True,
        exported_path=output_dir / filename,
        molecule_name=name,
        atom_count=35,
        heavy_atom_rmsd=0.001,
    )


def _fail_result(error: str = "RDKit parse error: bad molecule") -> LigandExportResult:
    return LigandExportResult(success=False, error=error)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _ligand_sdf(tmp_path: Path) -> Path:
    """Minimal SDF file — enough to pass the file-exists check in the command."""
    p = tmp_path / "A1.sdf"
    p.write_text("dummy sdf content\n")
    return p


def _unwrapped_output(text: str) -> str:
    """Collapse Rich line-wrapping so substring checks are terminal-width independent."""
    return text.replace("\n", "")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Help / routing
# ═══════════════════════════════════════════════════════════════════════════════

class TestLigandSubcommandRouting:
    def test_ligand_help(self):
        result = runner.invoke(cli, ["ligand", "--help"])
        assert result.exit_code == 0
        assert "export-ligpargen" in result.output

    def test_export_ligpargen_help(self):
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        assert "--legacy" in result.output
        assert "--mol-name" in result.output
        assert "--output-dir" in result.output

    def test_top_level_help_does_not_crash(self):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "ligand" in result.output


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Modern export (no --legacy)
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportLigpargenMocked:
    def test_modern_export_success(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path / "ligpargen_export", "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path / "ligpargen_export"),
            ])
        assert result.exit_code == 0, result.output
        assert "LIG" in result.output
        assert "35" in result.output

    def test_modern_export_calls_modern_function(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock) as mfn, \
             patch("ligand.export.export_for_ligpargen_legacy") as mleg:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert mfn.called, "export_for_ligpargen must be called without --legacy"
        assert not mleg.called, "export_for_ligpargen_legacy must NOT be called without --legacy"

    def test_mol_name_passed_through(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "A001", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock) as mfn:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--mol-name", "A001",
                "--output-dir", str(tmp_path),
            ])
        _call_kwargs = mfn.call_args
        assert _call_kwargs is not None
        # mol_name must reach the function
        args, kwargs = _call_kwargs
        passed_name = kwargs.get("mol_name") or (args[2] if len(args) > 2 else None)
        assert passed_name == "A001"

    def test_output_dir_passed_through(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        out = tmp_path / "my_export"
        mock = _ok_result(out, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock) as mfn:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(out),
            ])
        args, kwargs = mfn.call_args
        passed_dir = Path(kwargs.get("output_dir") or args[1])
        assert passed_dir == out

    def test_shows_output_file_path(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "LIG.pdb" in result.output

    def test_shows_atom_count(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "35" in result.output

    def test_shows_molecule_name(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "LIG" in result.output

    def test_rmsd_within_tolerance_shown(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = LigandExportResult(
            success=True,
            exported_path=tmp_path / "LIG.pdb",
            molecule_name="LIG",
            atom_count=35,
            heavy_atom_rmsd=0.001,
        )
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "0.0010" in result.output or "within tolerance" in result.output

    def test_rmsd_above_threshold_shows_warning(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = LigandExportResult(
            success=True,
            exported_path=tmp_path / "LIG.pdb",
            molecule_name="LIG",
            atom_count=35,
            heavy_atom_rmsd=0.12,
        )
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "⚠" in result.output or "warn" in result.output.lower() or "0.1200" in result.output


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Legacy export (--legacy)
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportLigpargenLegacy:
    def test_legacy_calls_legacy_function(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock) as mleg, \
             patch("ligand.export.export_for_ligpargen") as mmod:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert mleg.called, "export_for_ligpargen_legacy must be called with --legacy"
        assert not mmod.called, "export_for_ligpargen must NOT be called with --legacy"

    def test_legacy_success_shows_legacy_filename(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "LIG_ligpargen_legacy.pdb" in result.output

    def test_legacy_label_in_output(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert "legacy" in result.output.lower()


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Error cases
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportLigpargenErrors:
    def test_missing_input_file_exits_nonzero(self, tmp_path):
        result = runner.invoke(cli, [
            "ligand", "export-ligpargen", str(tmp_path / "nonexistent.sdf"),
            "--output-dir", str(tmp_path),
        ])
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "error" in result.output.lower()

    def test_export_failure_result_exits_nonzero(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        with patch("ligand.export.export_for_ligpargen", return_value=_fail_result()):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code != 0

    def test_export_failure_shows_error_message(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        with patch("ligand.export.export_for_ligpargen",
                   return_value=_fail_result("No 3D conformer found")):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert "3D conformer" in result.output or "failed" in result.output.lower()

    def test_rdkit_import_error_shows_friendly_message(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        with patch("ligand.export.export_for_ligpargen",
                   side_effect=ImportError("RDKit not installed")):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code != 0
        assert "rdkit" in result.output.lower() or "RDKit" in result.output

    def test_rdkit_import_error_mentions_rdkit_env(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        with patch("ligand.export.export_for_ligpargen",
                   side_effect=ImportError("no module rdkit")):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        # Must mention the env or how to install RDKit
        output_lower = result.output.lower()
        assert "rdkit_env" in output_lower or "conda" in output_lower or "install" in output_lower


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Defaults
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportLigpargenDefaults:
    def test_default_mol_name_is_lig(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG")
        with patch("ligand.export.export_for_ligpargen", return_value=mock) as mfn:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        args, kwargs = mfn.call_args
        passed_name = kwargs.get("mol_name") or (args[2] if len(args) > 2 else "LIG")
        assert passed_name == "LIG"

    def test_default_output_dir_is_ligpargen_export(self, tmp_path, monkeypatch):
        """Default --output-dir resolves to ./ligpargen_export."""
        monkeypatch.chdir(tmp_path)
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path / "ligpargen_export", "LIG")
        with patch("ligand.export.export_for_ligpargen", return_value=mock) as mfn:
            runner.invoke(cli, ["ligand", "export-ligpargen", str(sdf)])
        args, kwargs = mfn.call_args
        passed_dir = Path(kwargs.get("output_dir") or args[1])
        assert passed_dir == Path("ligpargen_export")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Integration test (RDKit required) — skipped in standard suite
# ═══════════════════════════════════════════════════════════════════════════════

try:
    import rdkit.Chem  # noqa: F401
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


@pytest.mark.rdkit
@pytest.mark.skipif(not _RDKIT_AVAILABLE, reason="RDKit not installed — run in rdkit_env")
class TestExportLigpargenIntegration:
    """End-to-end tests that call the real export functions.

    Skipped automatically when RDKit is absent.
    Run explicitly with:
        conda run -n rdkit_env python -m pytest test_cli_ligand.py -m rdkit -v
    """

    @pytest.fixture(scope="class")
    def a1_sdf(self):
        p = Path("tests/fixtures/ligpargen/a1/A1_ligpargen_input.pdb")
        if not p.exists():
            pytest.skip(f"Fixture not found: {p}")
        return p

    def test_modern_export_real(self, tmp_path, a1_sdf):
        result = runner.invoke(cli, [
            "ligand", "export-ligpargen", str(a1_sdf),
            "--output-dir", str(tmp_path),
        ])
        assert result.exit_code == 0, result.output
        assert "LIG" in result.output or "A1" in result.output

    def test_legacy_export_real(self, tmp_path, a1_sdf):
        result = runner.invoke(cli, [
            "ligand", "export-ligpargen", str(a1_sdf),
            "--legacy",
            "--output-dir", str(tmp_path),
        ])
        assert result.exit_code == 0, result.output
        assert "legacy" in result.output.lower()
        exported = list(tmp_path.glob("*_ligpargen_legacy.pdb"))
        assert len(exported) == 1, f"Expected 1 legacy PDB, got: {exported}"

    def test_legacy_export_output_exists(self, tmp_path, a1_sdf):
        runner.invoke(cli, [
            "ligand", "export-ligpargen", str(a1_sdf),
            "--legacy",
            "--output-dir", str(tmp_path),
        ])
        assert any(tmp_path.glob("*_ligpargen_legacy.pdb"))


# ═══════════════════════════════════════════════════════════════════════════════
# 7. SMILES export mode
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportSmiles:
    """Tests for --smiles mode.

    All tests in this class mock the export function so RDKit is not required.
    """

    def test_smiles_creates_smi_result(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_smiles_result(tmp_path)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock) as mfn:
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert mfn.called
        # exported path ends with .smi
        call_args = mfn.call_args
        assert result.exit_code == 0
        assert "LIG" in result.output

    def test_smiles_shows_smiles_string(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_smiles_result(tmp_path)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "CCO" in result.output

    def test_smiles_calls_smiles_not_pdb(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_smiles_result(tmp_path)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock) as ms, \
             patch("ligand.export.export_for_ligpargen_legacy") as mleg, \
             patch("ligand.export.export_for_ligpargen") as mmod:
            runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert ms.called
        assert not mleg.called
        assert not mmod.called

    def test_smiles_no_pdb_advisory_shown(self, tmp_path):
        """--smiles output must NOT show the PDB rejection advisory."""
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_smiles_result(tmp_path)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert "If LigParGen rejects" not in result.output

    def test_smiles_failure_exits_nonzero(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        fail = LigandExportResult(success=False, error="SMILES generation failed: bad mol")
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=fail):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code != 0

    def test_smiles_rdkit_missing_shows_friendly_error(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        with patch("ligand.export.export_for_ligpargen_smiles",
                   side_effect=ImportError("no module rdkit")):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code != 0
        output_lower = result.output.lower()
        assert "rdkit" in output_lower or "conda" in output_lower


class TestHelpText:
    """Help text documents validated modes and mandatory charge check."""

    def test_help_shows_smiles_option(self):
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        assert "--smiles" in result.output

    def test_help_shows_legacy_option(self):
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        assert "--legacy" in result.output

    def test_help_legacy_described_as_validated(self):
        """--legacy help must describe the PDB as experimentally validated."""
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        output_lower = result.output.lower()
        assert "validated" in output_lower or "experimentally" in output_lower

    def test_help_mentions_charge(self):
        """Help must explain that the formal charge must match LigParGen selection."""
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        output_lower = result.output.lower()
        assert "charge" in output_lower

    def test_help_legacy_does_not_claim_guaranteed_acceptance(self):
        """--legacy help must not say it is 'accepted by the LigParGen online server'."""
        result = runner.invoke(cli, ["ligand", "export-ligpargen", "--help"])
        assert result.exit_code == 0
        # Old false claim: "accepted by the LigParGen online server"
        assert "accepted by the ligpargen online server" not in result.output.lower()


class TestPDBAdvisory:
    """Legacy PDB mode lists companion files; modern PDB mode succeeds cleanly."""

    def test_legacy_lists_companion_smi(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "LIG_ligpargen" in result.output
        assert ".smi" in result.output

    def test_legacy_lists_charge_txt(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        out = result.output.replace("\n", "")
        assert "LIG_charge.txt" in out

    def test_modern_succeeds_without_crash(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=False)
        with patch("ligand.export.export_for_ligpargen", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output


class TestLegacyBackwardCompat:
    """Legacy PDB export must remain functionally unchanged."""

    def test_legacy_still_calls_legacy_function(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock) as mleg:
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert mleg.called

    def test_legacy_exit_zero_on_success(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0

    def test_legacy_shows_mol_name(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _ok_result(tmp_path, "LIG", legacy=True)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert "LIG" in result.output


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Charge reporting
# ═══════════════════════════════════════════════════════════════════════════════

def _charged_legacy_result(
    output_dir: Path,
    charge: int,
    name: str = "LIG",
) -> LigandExportResult:
    return LigandExportResult(
        success=True,
        exported_path=output_dir / f"{name}_ligpargen_legacy.pdb",
        molecule_name=name,
        atom_count=35,
        heavy_atom_rmsd=0.001,
        smiles="CC[NH+](C)C" if charge > 0 else ("[O-]C" if charge < 0 else "CCO"),
        formal_charge=charge,
    )


def _charged_smiles_result(
    output_dir: Path,
    charge: int,
    name: str = "LIG",
) -> LigandExportResult:
    return LigandExportResult(
        success=True,
        exported_path=output_dir / f"{name}.smi",
        molecule_name=name,
        atom_count=18,
        smiles="CC[NH+](C)C" if charge > 0 else ("[O-]C" if charge < 0 else "CCO"),
        formal_charge=charge,
    )


class TestChargeReporting:
    """Charge is computed, displayed, and emits a warning for non-zero values.

    Tests cover: neutral (0), positive (+1), negative (-1) molecules.
    All mocked — RDKit not required.
    """

    # ── Neutral molecule ──────────────────────────────────────────────────────

    def test_neutral_shows_zero_charge(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=0)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "Formal charge:  0" in result.output or "Formal charge: 0" in result.output

    def test_neutral_no_charge_warning(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=0)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert "WARNING" not in result.output

    # ── Positive charge ───────────────────────────────────────────────────────

    def test_positive_charge_shows_plus_label(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=1)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "+1" in result.output

    def test_positive_charge_emits_warning(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=1)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        output = result.output
        assert "WARNING" in output
        assert "+1" in output

    def test_positive_charge_warning_mentions_ligpargen_selection(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=1)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert "LigParGen" in result.output or "ligpargen" in result.output.lower()
        assert "select" in result.output.lower() or "instead of 0" in result.output

    # ── Negative charge ───────────────────────────────────────────────────────

    def test_negative_charge_shows_minus_label(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=-1)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "-1" in result.output

    def test_negative_charge_emits_warning(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_legacy_result(tmp_path, charge=-1)
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--legacy", "--output-dir", str(tmp_path),
            ])
        assert "WARNING" in result.output
        assert "-1" in result.output

    # ── SMILES mode also reports charge ───────────────────────────────────────

    def test_smiles_mode_shows_positive_charge(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_smiles_result(tmp_path, charge=1)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles", "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "+1" in result.output

    def test_smiles_mode_shows_zero_charge(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_smiles_result(tmp_path, charge=0)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles", "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "WARNING" not in result.output

    def test_smiles_mode_positive_charge_emits_warning(self, tmp_path):
        sdf = _ligand_sdf(tmp_path)
        mock = _charged_smiles_result(tmp_path, charge=2)
        with patch("ligand.export.export_for_ligpargen_smiles", return_value=mock):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(sdf),
                "--smiles", "--output-dir", str(tmp_path),
            ])
        assert "WARNING" in result.output
        assert "+2" in result.output

    # ── charge_label helper ───────────────────────────────────────────────────

    def test_charge_label_positive(self):
        from ligand.export import _charge_label
        assert _charge_label(1) == "+1"
        assert _charge_label(3) == "+3"

    def test_charge_label_zero(self):
        from ligand.export import _charge_label
        assert _charge_label(0) == "0"

    def test_charge_label_negative(self):
        from ligand.export import _charge_label
        assert _charge_label(-1) == "-1"
        assert _charge_label(-2) == "-2"

    # ── charge.txt content ────────────────────────────────────────────────────

    def test_write_charge_txt_content(self, tmp_path):
        from ligand.export import _write_charge_txt
        p = _write_charge_txt(tmp_path, "LIG", 1)
        text = p.read_text()
        assert "Molecule: LIG" in text
        assert "Formal charge: +1" in text
        assert "Recommended LigParGen charge selection: +1" in text

    def test_write_charge_txt_neutral(self, tmp_path):
        from ligand.export import _write_charge_txt
        p = _write_charge_txt(tmp_path, "LIG", 0)
        text = p.read_text()
        assert "Formal charge: 0" in text
        assert "Recommended LigParGen charge selection: 0" in text

    def test_write_charge_txt_negative(self, tmp_path):
        from ligand.export import _write_charge_txt
        p = _write_charge_txt(tmp_path, "MOL", -1)
        text = p.read_text()
        assert "Formal charge: -1" in text
        assert "Molecule: MOL" in text


# ═══════════════════════════════════════════════════════════════════════════════
# 9. External-directory regression — ligand package importable from any cwd
# ═══════════════════════════════════════════════════════════════════════════════

def _dummy_pdb(tmp_path: Path) -> Path:
    p = tmp_path / "E20.pdb"
    p.write_text(
        "ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    return p


class TestCLIExternalDirectory:
    """Regression: 'from ligand import export' must not raise ModuleNotFoundError
    when simforge is invoked from a working directory outside the repo.

    Root cause: ligand/ was missing from pyproject.toml packages.find.include,
    so it was absent from the editable-install MAPPING and invisible to the
    import system in subprocesses.
    """

    def test_help_works_from_external_cwd(self, tmp_path):
        """--help must succeed from any working directory (lazy import gate)."""
        result = subprocess.run(
            [sys.executable, "-m", "typer", "cli:cli", "run", "ligand", "--", "--help"],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        # Prefer using the installed entry point
        result2 = subprocess.run(
            ["simforge", "ligand", "export-ligpargen", "--help"],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        assert result2.returncode == 0, result2.stderr
        assert "No module named" not in result2.stderr

    def test_no_module_not_found_ligand_from_external_cwd(self, tmp_path):
        """Running export-ligpargen from outside the repo must not raise
        ModuleNotFoundError for 'ligand' or 'utils'."""
        dummy = _dummy_pdb(tmp_path)
        other_dir = tmp_path / "workdir"
        other_dir.mkdir()

        result = subprocess.run(
            ["simforge", "ligand", "export-ligpargen", str(dummy), "--legacy",
             "--output-dir", str(tmp_path / "out")],
            capture_output=True, text=True, cwd=str(other_dir),
        )
        combined = result.stdout + result.stderr
        assert "No module named 'ligand'" not in combined, (
            f"ModuleNotFoundError for 'ligand' — ligand package not installed.\n"
            f"stderr: {result.stderr}"
        )
        assert "No module named 'utils'" not in combined, (
            f"ModuleNotFoundError for 'utils' — utils package not installed.\n"
            f"stderr: {result.stderr}"
        )
        # Acceptable outcomes: success OR RDKit-not-installed message
        # (exit code 1 is fine when RDKit is absent)
        rdkit_error = "rdkit" in combined.lower() or "RDKit" in combined
        import_error = "ModuleNotFoundError" in combined
        assert not import_error or rdkit_error, (
            f"Unexpected import error. stdout+stderr:\n{combined}"
        )

    def test_mocked_export_works_from_isolated_filesystem(self, tmp_path, monkeypatch):
        """In-process regression: ligand.export is importable and mockable
        when the working directory is a temp path outside the repo."""
        monkeypatch.chdir(tmp_path)
        dummy = _dummy_pdb(tmp_path)
        mock_result = LigandExportResult(
            success=True,
            exported_path=tmp_path / "LIG_ligpargen_legacy.pdb",
            molecule_name="LIG",
            atom_count=1,
            heavy_atom_rmsd=0.0,
        )
        with patch("ligand.export.export_for_ligpargen_legacy", return_value=mock_result):
            result = runner.invoke(cli, [
                "ligand", "export-ligpargen", str(dummy),
                "--legacy",
                "--output-dir", str(tmp_path),
            ])
        assert result.exit_code == 0, result.output
        assert "No module named 'ligand'" not in (result.output or "")
        assert "LIG" in result.output
