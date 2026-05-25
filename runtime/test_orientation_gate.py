"""
Tests for runtime/orientation_gate.py.

Covers all three gate outcomes (blocked / advisory / pass) plus
the read_orientation_report workspace scanner.
"""
import json
import pytest
from pathlib import Path

from runtime.orientation_gate import evaluate_orientation_gate, read_orientation_report


# ─── helpers ──────────────────────────────────────────────────────────────────

def _write_report(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data))


def _base_report(**kwargs) -> dict:
    base = {
        "confidence": 0.9,
        "passed": True,
        "geometry": {"ec_com_z_nm": 1.5, "ic_com_z_nm": -1.5, "tilt_angle_deg": 5.0},
        "checks": [],
        "warnings": [],
        "errors": [],
    }
    base.update(kwargs)
    return base


# ─── evaluate_orientation_gate ─────────────────────────────────────────────────

class TestEvaluateOrientationGate:

    def test_missing_report_returns_none(self, tmp_path):
        assert evaluate_orientation_gate(tmp_path) is None

    def test_corrupted_json_returns_none(self, tmp_path):
        (tmp_path / "orientation_report.json").write_text("not json {{")
        assert evaluate_orientation_gate(tmp_path) is None

    def test_clean_pass(self, tmp_path):
        _write_report(tmp_path / "orientation_report.json", _base_report())
        result = evaluate_orientation_gate(tmp_path)
        assert result is not None
        assert result.passed is True
        assert result.blocked is False
        assert result.errors == []
        assert result.warnings == []
        assert result.confidence == pytest.approx(0.9)

    def test_warnings_only_not_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.75,
            passed=True,
            warnings=["TM COM slightly off-center (z=+2.1 nm)"],
        )
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.blocked is False
        assert result.passed is True
        assert len(result.warnings) == 1
        assert "TM COM" in result.warnings[0]

    def test_errors_trigger_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.5,
            passed=False,
            errors=["Orientation inverted: EC COM_z=-1.2 nm, IC COM_z=+1.4 nm."],
        )
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.blocked is True
        assert result.passed is False
        assert len(result.errors) == 1

    def test_multiple_errors_all_captured(self, tmp_path):
        data = _base_report(
            confidence=0.3,
            passed=False,
            errors=[
                "No CA atoms found for EC residues.",
                "No CA atoms found for IC residues.",
            ],
        )
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.blocked is True
        assert len(result.errors) == 2

    def test_errors_and_warnings_still_blocked(self, tmp_path):
        data = _base_report(
            confidence=0.4,
            passed=False,
            errors=["Orientation inverted."],
            warnings=["TM COM off-center."],
        )
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.blocked is True
        assert len(result.warnings) == 1

    def test_status_str_blocked(self, tmp_path):
        data = _base_report(passed=False, errors=["inverted"])
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.status_str == "BLOCKED"

    def test_status_str_advisory(self, tmp_path):
        data = _base_report(warnings=["small separation"])
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.status_str == "PASS (advisory)"

    def test_status_str_pass(self, tmp_path):
        _write_report(tmp_path / "orientation_report.json", _base_report())
        result = evaluate_orientation_gate(tmp_path)
        assert result.status_str == "PASS"

    def test_confidence_zero_with_errors(self, tmp_path):
        data = _base_report(confidence=0.0, passed=False,
                            errors=["No CA atoms for EC."])
        _write_report(tmp_path / "orientation_report.json", data)
        result = evaluate_orientation_gate(tmp_path)
        assert result.confidence == pytest.approx(0.0)
        assert result.blocked is True


# ─── read_orientation_report ──────────────────────────────────────────────────

class TestReadOrientationReport:

    def test_no_steps_dir(self, tmp_path):
        assert read_orientation_report(tmp_path) is None

    def test_steps_dir_empty(self, tmp_path):
        (tmp_path / "steps").mkdir()
        assert read_orientation_report(tmp_path) is None

    def test_no_orientation_report_in_steps(self, tmp_path):
        step = tmp_path / "steps" / "01_prepare_protein_1"
        step.mkdir(parents=True)
        (step / "metadata.json").write_text("{}")
        assert read_orientation_report(tmp_path) is None

    def test_finds_report_in_step(self, tmp_path):
        step = tmp_path / "steps" / "01_orient_protein_1"
        step.mkdir(parents=True)
        data = _base_report(confidence=0.88)
        _write_report(step / "orientation_report.json", data)
        result = read_orientation_report(tmp_path)
        assert result is not None
        assert result["confidence"] == pytest.approx(0.88)

    def test_returns_first_alphabetically(self, tmp_path):
        # Creates two steps with orientation reports — should pick the first sorted one.
        for name, conf in [("01_orient_protein_1", 0.91), ("02_some_other_1", 0.50)]:
            step = tmp_path / "steps" / name
            step.mkdir(parents=True)
            _write_report(step / "orientation_report.json", _base_report(confidence=conf))
        result = read_orientation_report(tmp_path)
        assert result["confidence"] == pytest.approx(0.91)

    def test_corrupted_report_returns_none(self, tmp_path):
        step = tmp_path / "steps" / "01_orient_protein_1"
        step.mkdir(parents=True)
        (step / "orientation_report.json").write_text("bad json {{{")
        assert read_orientation_report(tmp_path) is None
