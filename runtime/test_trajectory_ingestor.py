"""Phase A — unit tests for runtime/trajectory_ingestor.py."""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from runtime.trajectory_ingestor import (
    TrajectoryManifest,
    discover_trajectory,
    load_xvg_files,
)


# ─────────────────────────────────────────────────────────────────────────────
# Minimal XVG content factories
# ─────────────────────────────────────────────────────────────────────────────

def _write_xvg(path: Path, title: str = "", n_points: int = 20, value: float = 0.2) -> None:
    """Write a minimal valid XVG file."""
    lines = []
    if title:
        lines.append(f'@ title "{title}"')
    lines.append('@ xaxis label "Time (ps)"')
    lines.append('@ yaxis label "Value"')
    lines.append('@ s0 legend "Series"')
    for i in range(n_points):
        lines.append(f"{float(i * 100):.1f}  {value:.4f}")
    path.write_text("\n".join(lines) + "\n")


def _write_xvg_rmsd(path: Path, n_points: int = 100) -> None:
    _write_xvg(path, title="RMSD", n_points=n_points, value=0.25)


def _write_xvg_energy(path: Path, n_points: int = 100) -> None:
    _write_xvg(path, title="Potential Energy", n_points=n_points, value=-50000.0)


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestDiscoverFlat:
    """Flat directory with named XVG files."""

    def test_rmsd_and_energy_labeled(self, tmp_path):
        _write_xvg_rmsd(tmp_path / "md_rmsd.xvg")
        _write_xvg_energy(tmp_path / "energy.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files
        assert "potential_energy" in manifest.xvg_files

    def test_rmsd_named_file(self, tmp_path):
        _write_xvg_rmsd(tmp_path / "md_rmsd.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files
        assert manifest.xvg_files["rmsd"].name == "md_rmsd.xvg"

    def test_energy_named_file(self, tmp_path):
        _write_xvg_energy(tmp_path / "energy.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "potential_energy" in manifest.xvg_files

    def test_non_xvg_files_ignored(self, tmp_path):
        (tmp_path / "some.log").write_text("log content")
        (tmp_path / "traj.xtc").write_bytes(b"\x00" * 16)
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files == {}
        assert len(manifest.log_files) == 1
        assert len(manifest.xtc_files) == 1

    def test_xtc_edr_tpr_gro_collected(self, tmp_path):
        (tmp_path / "md.xtc").write_bytes(b"\x00")
        (tmp_path / "md.edr").write_bytes(b"\x00")
        (tmp_path / "md.tpr").write_bytes(b"\x00")
        (tmp_path / "md.gro").write_text("GRO\n0\n0 0 0\n")
        manifest = discover_trajectory(tmp_path)
        assert len(manifest.xtc_files) == 1
        assert len(manifest.edr_files) == 1
        assert len(manifest.tpr_files) == 1
        assert len(manifest.gro_files) == 1


class TestDiscoverSimForgeWorkspace:
    """SimForge workspace detection (steps/ + metadata/)."""

    def test_simforge_workspace_detected(self, tmp_path):
        (tmp_path / "steps").mkdir()
        (tmp_path / "metadata").mkdir()
        # Put an XVG in a step subdir
        step_dir = tmp_path / "steps" / "10_production_1"
        step_dir.mkdir(parents=True)
        _write_xvg_rmsd(step_dir / "rmsd.xvg")
        manifest = discover_trajectory(tmp_path)
        assert manifest.simforge_workspace is True
        assert "rmsd" in manifest.xvg_files

    def test_non_simforge_not_flagged(self, tmp_path):
        _write_xvg_rmsd(tmp_path / "rmsd.xvg")
        manifest = discover_trajectory(tmp_path)
        assert manifest.simforge_workspace is False

    def test_simforge_root_path_stored(self, tmp_path):
        (tmp_path / "steps").mkdir()
        (tmp_path / "metadata").mkdir()
        manifest = discover_trajectory(tmp_path)
        assert manifest.root == tmp_path


class TestEmptyDirectory:
    def test_empty_dir_returns_empty_manifest(self, tmp_path):
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files == {}
        assert manifest.xtc_files == []
        assert manifest.edr_files == []
        assert manifest.simforge_workspace is False

    def test_empty_dir_no_error(self, tmp_path):
        # Must not raise
        manifest = discover_trajectory(tmp_path)
        assert isinstance(manifest, TrajectoryManifest)

    def test_nonexistent_dir_returns_empty(self, tmp_path):
        fake = tmp_path / "does_not_exist"
        manifest = discover_trajectory(fake)
        assert manifest.xvg_files == {}
        assert manifest.simforge_workspace is False


class TestMixedLayout:
    """Mixed flat and nested directories."""

    def test_flat_and_nested_xvg_discovered(self, tmp_path):
        # Flat RMSD
        _write_xvg_rmsd(tmp_path / "rmsd.xvg")
        # Nested energy
        nested = tmp_path / "analysis" / "energy"
        nested.mkdir(parents=True)
        _write_xvg_energy(nested / "energy.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files
        assert "potential_energy" in manifest.xvg_files

    def test_multiple_xtc_collected(self, tmp_path):
        for i in range(3):
            (tmp_path / f"part{i}.xtc").write_bytes(b"\x00")
        manifest = discover_trajectory(tmp_path)
        assert len(manifest.xtc_files) == 3

    def test_rmsf_label(self, tmp_path):
        _write_xvg(tmp_path / "rmsf.xvg", title="RMSF")
        manifest = discover_trajectory(tmp_path)
        assert "rmsf" in manifest.xvg_files

    def test_temperature_label(self, tmp_path):
        _write_xvg(tmp_path / "temperature.xvg", title="Temperature")
        manifest = discover_trajectory(tmp_path)
        assert "temperature" in manifest.xvg_files

    def test_pressure_label(self, tmp_path):
        _write_xvg(tmp_path / "pressure.xvg", title="Pressure")
        manifest = discover_trajectory(tmp_path)
        assert "pressure" in manifest.xvg_files

    def test_gyration_label(self, tmp_path):
        _write_xvg(tmp_path / "gyrate.xvg", title="Radius of gyration")
        manifest = discover_trajectory(tmp_path)
        assert "gyration" in manifest.xvg_files

    def test_unknown_xvg_uses_stem(self, tmp_path):
        _write_xvg(tmp_path / "my_weird_data.xvg", title="Some custom data")
        manifest = discover_trajectory(tmp_path)
        assert "my_weird_data" in manifest.xvg_files


class TestLoadXvgFiles:
    def test_loads_parseable_xvg(self, tmp_path):
        _write_xvg_rmsd(tmp_path / "rmsd.xvg")
        manifest = discover_trajectory(tmp_path)
        data = load_xvg_files(manifest)
        assert "rmsd" in data
        assert len(data["rmsd"].time_ps) > 0

    def test_skips_corrupted_xvg(self, tmp_path):
        # Valid file
        _write_xvg_rmsd(tmp_path / "rmsd.xvg")
        # Corrupted file (but parse_xvg is very lenient — it returns empty, not error)
        bad = tmp_path / "bad_file.xvg"
        bad.write_bytes(b"\xff\xfe invalid binary")
        # discover will label bad_file.xvg → stem "bad_file"
        manifest = discover_trajectory(tmp_path)
        data = load_xvg_files(manifest)
        # Must not raise regardless of what bad_file.xvg contains
        assert isinstance(data, dict)

    def test_empty_manifest_returns_empty_dict(self, tmp_path):
        manifest = discover_trajectory(tmp_path)
        data = load_xvg_files(manifest)
        assert data == {}

    def test_returns_xvgdata_objects(self, tmp_path):
        from runtime.xvg_parser import XVGData
        _write_xvg_rmsd(tmp_path / "rmsd.xvg")
        manifest = discover_trajectory(tmp_path)
        data = load_xvg_files(manifest)
        for v in data.values():
            assert isinstance(v, XVGData)


# ─────────────────────────────────────────────────────────────────────────────
# TestNameTolerantDiscovery — real-world GROMACS naming patterns
# ─────────────────────────────────────────────────────────────────────────────

class TestNameTolerantDiscovery:
    """Filenames that are common in practice but were not matched by the
    original _LABEL_RULES (hb_, rg, sasa) or need regression coverage
    (rmsd-aldosa, rmsd-ligand, rmsf)."""

    # ── RMSD variants ─────────────────────────────────────────────────────────

    def test_rmsd_aldosa_labeled_rmsd(self, tmp_path):
        _write_xvg(tmp_path / "rmsd-aldosa.xvg", title="RMSD")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files
        assert manifest.xvg_files["rmsd"].name == "rmsd-aldosa.xvg"

    def test_rmsd_ligand_labeled_rmsd(self, tmp_path):
        _write_xvg(tmp_path / "rmsd-ligand.xvg", title="RMSD")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files

    def test_two_rmsd_files_both_handled(self, tmp_path):
        """Both rmsd-*.xvg files map to 'rmsd'; dedup keeps the larger one."""
        _write_xvg(tmp_path / "rmsd-aldosa.xvg", title="RMSD", n_points=50)
        _write_xvg(tmp_path / "rmsd-ligand.xvg", title="RMSD", n_points=100)
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files
        # Larger file (more points) wins
        assert manifest.xvg_files["rmsd"].name == "rmsd-ligand.xvg"

    def test_rmsd_protein_labeled_rmsd(self, tmp_path):
        _write_xvg(tmp_path / "rmsd-protein.xvg", title="RMSD")
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files.get("rmsd") is not None

    def test_rmsd_underscore_labeled_rmsd(self, tmp_path):
        _write_xvg(tmp_path / "rmsd_backbone.xvg", title="RMSD")
        manifest = discover_trajectory(tmp_path)
        assert "rmsd" in manifest.xvg_files

    # ── Hydrogen bonds — hb prefix ────────────────────────────────────────────

    def test_hb_underscore_labeled_hydrogen_bonds(self, tmp_path):
        """hb_.xvg is the real-world filename that triggered this fix."""
        _write_xvg(tmp_path / "hb_.xvg", title="Number of Hydrogen Bonds")
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files.get("hydrogen_bonds") is not None
        assert manifest.xvg_files["hydrogen_bonds"].name == "hb_.xvg"

    def test_hbond_labeled_hydrogen_bonds(self, tmp_path):
        _write_xvg(tmp_path / "hbond.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "hydrogen_bonds" in manifest.xvg_files

    def test_hbonds_labeled_hydrogen_bonds(self, tmp_path):
        _write_xvg(tmp_path / "hbonds.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "hydrogen_bonds" in manifest.xvg_files

    def test_hb_dash_labeled_hydrogen_bonds(self, tmp_path):
        _write_xvg(tmp_path / "hb-prot-lig.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "hydrogen_bonds" in manifest.xvg_files

    # ── Radius of gyration — rg prefix ───────────────────────────────────────

    def test_rg_labeled_gyration(self, tmp_path):
        """rg.xvg is the bare GROMACS default output name."""
        _write_xvg(tmp_path / "rg.xvg", title="Radius of gyration")
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files.get("gyration") is not None
        assert manifest.xvg_files["gyration"].name == "rg.xvg"

    def test_rg_protein_labeled_gyration(self, tmp_path):
        _write_xvg(tmp_path / "rg_protein.xvg")
        manifest = discover_trajectory(tmp_path)
        assert "gyration" in manifest.xvg_files

    def test_rg_does_not_match_energy_file(self, tmp_path):
        """'rg' is a substring of 'energy' — must NOT mislabel energy files."""
        _write_xvg(tmp_path / "energy.xvg", title="Potential Energy")
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files.get("potential_energy") is not None
        assert "gyration" not in manifest.xvg_files

    # ── SASA ─────────────────────────────────────────────────────────────────

    def test_sasa_labeled_sasa(self, tmp_path):
        _write_xvg(tmp_path / "sasa.xvg", title="Solvent Accessible Surface")
        manifest = discover_trajectory(tmp_path)
        assert manifest.xvg_files.get("sasa") is not None
        assert manifest.xvg_files["sasa"].name == "sasa.xvg"

    def test_sasa_no_title_labeled_sasa(self, tmp_path):
        _write_xvg(tmp_path / "sasa.xvg", title="")
        manifest = discover_trajectory(tmp_path)
        assert "sasa" in manifest.xvg_files

    # ── RMSF (regression) ────────────────────────────────────────────────────

    def test_rmsf_labeled_rmsf(self, tmp_path):
        _write_xvg(tmp_path / "rmsf.xvg", title="RMSF")
        manifest = discover_trajectory(tmp_path)
        assert "rmsf" in manifest.xvg_files

    # ── Full set from the real-world report ──────────────────────────────────

    def test_full_set_all_mapped(self, tmp_path):
        """Exact filenames from the user-reported failure case."""
        files = {
            "rmsd-aldosa.xvg":  "RMSD",
            "rmsd-ligand.xvg":  "RMSD",
            "rmsf.xvg":         "RMSF",
            "rg.xvg":           "Radius of gyration",
            "sasa.xvg":         "Solvent Accessible Surface",
            "hb_.xvg":          "Number of Hydrogen Bonds",
        }
        for fname, title in files.items():
            _write_xvg(tmp_path / fname, title=title)
        manifest = discover_trajectory(tmp_path)
        assert "rmsd"            in manifest.xvg_files
        assert "rmsf"            in manifest.xvg_files
        assert "gyration"        in manifest.xvg_files
        assert "sasa"            in manifest.xvg_files
        assert "hydrogen_bonds"  in manifest.xvg_files
        # Both rmsd files present — one kept, one deduped
        assert len(manifest.xvg_files) == 5  # 6 files → 5 labels (rmsd deduped)

    # ── _label_xvg unit tests ────────────────────────────────────────────────

    def test_label_xvg_hb_underscore(self):
        from runtime.trajectory_ingestor import _label_xvg
        assert _label_xvg(Path("hb_.xvg"), "") == "hydrogen_bonds"

    def test_label_xvg_rg(self):
        from runtime.trajectory_ingestor import _label_xvg
        assert _label_xvg(Path("rg.xvg"), "") == "gyration"

    def test_label_xvg_sasa(self):
        from runtime.trajectory_ingestor import _label_xvg
        assert _label_xvg(Path("sasa.xvg"), "") == "sasa"

    def test_label_xvg_rmsd_aldosa(self):
        from runtime.trajectory_ingestor import _label_xvg
        assert _label_xvg(Path("rmsd-aldosa.xvg"), "") == "rmsd"

    def test_label_xvg_energy_not_gyration(self):
        from runtime.trajectory_ingestor import _label_xvg
        assert _label_xvg(Path("energy.xvg"), "Potential Energy") == "potential_energy"
