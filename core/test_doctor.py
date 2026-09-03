"""Tests for core.doctor — environment inspection."""

from __future__ import annotations

import core.doctor as doctor_mod
from core.doctor import (
    AVAILABLE, MISSING, PASS, REQUIRED, WRONG_ENVIRONMENT, DoctorCheck,
    DoctorReport, run_doctor,
)


def test_report_ready_only_depends_on_required_checks():
    checks = [
        DoctorCheck("SimForge", REQUIRED, PASS),
        DoctorCheck("RDKit", "OPTIONAL", MISSING),
    ]
    assert DoctorReport(checks).ready is True
    assert DoctorReport(checks).overall == "READY"


def test_report_not_ready_when_required_missing():
    checks = [DoctorCheck("GROMACS", REQUIRED, MISSING)]
    assert DoctorReport(checks).ready is False
    assert DoctorReport(checks).overall == "NOT READY"


def test_run_doctor_reports_gromacs(monkeypatch):
    monkeypatch.setattr(doctor_mod.shutil, "which", lambda name: "/usr/bin/gmx" if name == "gmx" else None)
    monkeypatch.setattr(doctor_mod, "gromacs_version", lambda: "2025.2")
    report = run_doctor()
    gmx = next(c for c in report.checks if c.name == "GROMACS")
    assert gmx.status == PASS
    assert gmx.version == "2025.2"


def test_run_doctor_flags_missing_gromacs(monkeypatch):
    monkeypatch.setattr(doctor_mod.shutil, "which", lambda name: None)
    monkeypatch.setattr(doctor_mod, "gromacs_version", lambda: None)
    report = run_doctor()
    gmx = next(c for c in report.checks if c.name == "GROMACS")
    assert gmx.status == MISSING
    assert report.ready is False


def test_run_doctor_detects_wrong_environment(monkeypatch):
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "base")
    monkeypatch.setattr(doctor_mod, "conda_env", lambda: "/opt/conda")
    monkeypatch.setattr(doctor_mod, "_module_version", lambda mod: None)
    report = run_doctor()
    env = next(c for c in report.checks if c.name == "Environment")
    assert env.status == WRONG_ENVIRONMENT
    # optional check -> does not by itself make the report NOT READY
    assert env.level == "OPTIONAL"


def test_run_doctor_rdkit_available(monkeypatch):
    monkeypatch.setattr(doctor_mod, "_module_version", lambda mod: "2026.03.3" if mod == "rdkit" else None)
    report = run_doctor()
    rd = next(c for c in report.checks if c.name == "RDKit")
    assert rd.status == AVAILABLE
    assert rd.version == "2026.03.3"


def test_doctor_report_json_roundtrip():
    report = run_doctor()
    d = report.to_dict()
    assert set(d) == {"overall", "ready", "checks"}
    assert all({"name", "level", "status"} <= set(c) for c in d["checks"])
