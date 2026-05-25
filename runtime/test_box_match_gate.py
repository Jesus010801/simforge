"""
Tests for runtime/box_match_gate.py.

Gate rules:
  - errors > 0  → blocked
  - warnings > 0 → advisory (not blocked)
  - clean       → pass
"""
import json
import pytest
from pathlib import Path

from runtime.box_match_gate import evaluate_box_match_gate, read_box_match_report


# ─── helpers ──────────────────────────────────────────────────────────────────

def _write_report(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data))


def _base_report(**kwargs) -> dict:
    base = {
        "lipid_type": "DPPC",
        "fallback_used": False,
        "confidence": 0.85,
        "passed": True,
        "protein_geometry": {
            "x_extent_nm": 3.0, "y_extent_nm": 3.5, "z_extent_nm": 5.0,
        },
        "recommended_box": {
            "box_x_nm": 7.0, "box_y_nm": 7.5, "box_z_nm": 14.3,
            "bilayer_thickness_nm": 3.8,
            "lateral_padding_nm": 2.0, "solvent_z_padding_nm": 2.5,
        },
        "estimates": {
            "n_lipids_estimate": 252, "box_area_nm2": 52.5,
            "box_volume_nm3": 750.3, "solvent_volume_nm3": 220.0,
            "protein_xy_coverage": 0.27,
        },
        "checks": [],
        "warnings": [],
        "errors": [],
    }
    base.update(kwargs)
    return base


# ─── evaluate_box_match_gate ───────────────────────────────────────────────────

class TestEvaluateBoxMatchGate:

    def test_missing_report_returns_none(self, tmp_path):
        assert evaluate_box_match_gate(tmp_path) is None

    def test_corrupted_json_returns_none(self, tmp_path):
        (tmp_path / "box_match_report.json").write_text("not json")
        assert evaluate_box_match_gate(tmp_path) is None

    def test_clean_pass(self, tmp_path):
        _write_report(tmp_path / "box_match_report.json", _base_report())
        r = evaluate_box_match_gate(tmp_path)
        assert r is not None
        assert r.passed is True
        assert r.blocked is False
        assert r.errors == []
        assert r.warnings == []
        assert r.confidence == pytest.approx(0.85)

    def test_box_dimensions_exposed(self, tmp_path):
        _write_report(tmp_path / "box_match_report.json", _base_report())
        r = evaluate_box_match_gate(tmp_path)
        assert r.box_x_nm == pytest.approx(7.0)
        assert r.box_y_nm == pytest.approx(7.5)
        assert r.box_z_nm == pytest.approx(14.3)
        assert r.n_lipids_estimate == 252
        assert r.protein_xy_coverage == pytest.approx(0.27)

    def test_warnings_only_not_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.70, passed=True,
            warnings=["Protein XY footprint covers 68% of box XY."],
        )
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.blocked is False
        assert r.passed is True
        assert len(r.warnings) == 1

    def test_errors_trigger_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.5, passed=False,
            errors=["Protein Z extent (0.3 nm < 1.0 nm). Structure may not be membrane-oriented."],
        )
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.blocked is True
        assert r.passed is False

    def test_multiple_errors_all_captured(self, tmp_path):
        data = _base_report(
            confidence=0.3, passed=False,
            errors=["Protein Z too small.", "Protein XY too large."],
        )
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.blocked is True
        assert len(r.errors) == 2

    def test_errors_and_warnings_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.35, passed=False,
            errors=["Protein Z too small."],
            warnings=["Protein very tall."],
        )
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.blocked is True
        assert len(r.warnings) == 1

    def test_status_str_blocked(self, tmp_path):
        data = _base_report(passed=False, errors=["Z too small."])
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.status_str == "BLOCKED"

    def test_status_str_advisory(self, tmp_path):
        data = _base_report(warnings=["Large protein."])
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.status_str == "PASS (advisory)"

    def test_status_str_pass(self, tmp_path):
        _write_report(tmp_path / "box_match_report.json", _base_report())
        r = evaluate_box_match_gate(tmp_path)
        assert r.status_str == "PASS"

    def test_fallback_lipid_exposed(self, tmp_path):
        data = _base_report(lipid_type="UNKNOWN_LIPID", fallback_used=True,
                            warnings=["Lipid 'UNKNOWN_LIPID' not in parameter table."])
        _write_report(tmp_path / "box_match_report.json", data)
        r = evaluate_box_match_gate(tmp_path)
        assert r.blocked is False
        assert len(r.warnings) == 1


# ─── read_box_match_report ────────────────────────────────────────────────────

class TestReadBoxMatchReport:

    def test_no_steps_dir(self, tmp_path):
        assert read_box_match_report(tmp_path) is None

    def test_steps_dir_empty(self, tmp_path):
        (tmp_path / "steps").mkdir()
        assert read_box_match_report(tmp_path) is None

    def test_no_box_match_in_steps(self, tmp_path):
        step = tmp_path / "steps" / "01_orient_protein_1"
        step.mkdir(parents=True)
        (step / "orientation_report.json").write_text("{}")
        assert read_box_match_report(tmp_path) is None

    def test_finds_report_in_step(self, tmp_path):
        step = tmp_path / "steps" / "02_match_box_to_bilayer_1"
        step.mkdir(parents=True)
        data = _base_report(confidence=0.92)
        _write_report(step / "box_match_report.json", data)
        r = read_box_match_report(tmp_path)
        assert r is not None
        assert r["confidence"] == pytest.approx(0.92)

    def test_corrupted_report_returns_none(self, tmp_path):
        step = tmp_path / "steps" / "02_match_box_to_bilayer_1"
        step.mkdir(parents=True)
        (step / "box_match_report.json").write_text("{bad json")
        assert read_box_match_report(tmp_path) is None
