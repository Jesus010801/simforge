"""
``simforge doctor`` — inspect the runtime environment and report the state of
the scientific software SimForge relies on.

Each check yields a ``DoctorCheck`` with a requirement level and a status:

    level  : REQUIRED | OPTIONAL
    status : PASS | AVAILABLE | MISSING | WARN | WRONG_ENVIRONMENT

Overall status is READY unless a REQUIRED check is MISSING/WARN, in which case
it is NOT READY. Missing OPTIONAL software never fails the report.
"""

from __future__ import annotations

import os
import platform
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from core.provenance import (
    conda_env,
    gromacs_data_prefix,
    gromacs_version,
    openbabel_version,
    simforge_version,
    _module_version,
)

REQUIRED = "REQUIRED"
OPTIONAL = "OPTIONAL"

PASS = "PASS"
AVAILABLE = "AVAILABLE"
MISSING = "MISSING"
WARN = "WARN"
WRONG_ENVIRONMENT = "WRONG_ENVIRONMENT"

_OK_STATUSES = {PASS, AVAILABLE}


@dataclass
class DoctorCheck:
    name: str
    level: str
    status: str
    detail: str = ""
    version: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "level": self.level,
            "status": self.status,
            "version": self.version,
            "detail": self.detail,
        }


@dataclass
class DoctorReport:
    checks: list[DoctorCheck] = field(default_factory=list)

    @property
    def ready(self) -> bool:
        return all(
            c.status in _OK_STATUSES
            for c in self.checks
            if c.level == REQUIRED
        )

    @property
    def overall(self) -> str:
        return "READY" if self.ready else "NOT READY"

    def to_dict(self) -> dict[str, Any]:
        return {
            "overall": self.overall,
            "ready": self.ready,
            "checks": [c.to_dict() for c in self.checks],
        }


def _forcefield_available(ff: str = "oplsaa") -> tuple[str, str]:
    """Return (status, detail) for a base force field bundled with GROMACS."""
    top = gromacs_data_prefix()
    if top is None:
        return (WARN, "GROMACS data directory not found — cannot verify force fields")
    ff_dir = Path(top) / f"{ff}.ff"
    if ff_dir.is_dir() and (ff_dir / "forcefield.itp").exists():
        return (PASS, str(ff_dir))
    return (MISSING, f"{ff}.ff not found under {top}")


def run_doctor(forcefield: str = "oplsaa") -> DoctorReport:
    checks: list[DoctorCheck] = []

    # ── SimForge ────────────────────────────────────────────────────────────
    checks.append(DoctorCheck(
        "SimForge", REQUIRED, PASS, version=simforge_version(),
        detail=str(Path(__file__).resolve().parent.parent),
    ))

    # ── Python ──────────────────────────────────────────────────────────────
    py_ok = sys.version_info >= (3, 11)
    checks.append(DoctorCheck(
        "Python", REQUIRED, PASS if py_ok else WARN,
        version=platform.python_version(),
        detail=sys.executable + ("" if py_ok else "  (SimForge requires >= 3.11)"),
    ))

    # ── Environment (conda) ─────────────────────────────────────────────────
    env = conda_env()
    env_name = os.environ.get("CONDA_DEFAULT_ENV")
    rdkit_v = _module_version("rdkit")
    if env is None:
        checks.append(DoctorCheck(
            "Environment", OPTIONAL, AVAILABLE,
            detail="no conda environment detected (system Python)",
        ))
    else:
        wrong = env_name in (None, "base") and rdkit_v is None
        checks.append(DoctorCheck(
            "Environment", OPTIONAL,
            WRONG_ENVIRONMENT if wrong else AVAILABLE,
            detail=(env_name or env) + (
                "  — RDKit-dependent steps need a dedicated env (e.g. rdkit_env)"
                if wrong else ""
            ),
        ))

    # ── GROMACS ─────────────────────────────────────────────────────────────
    gmx = shutil.which("gmx")
    gv = gromacs_version()
    checks.append(DoctorCheck(
        "GROMACS", REQUIRED,
        PASS if gmx else MISSING,
        version=gv,
        detail=gmx or "gmx not found in PATH — required for pdb2gmx and grompp",
    ))

    # ── Force field ─────────────────────────────────────────────────────────
    ff_status, ff_detail = _forcefield_available(forcefield)
    # Only block when we could actually look and the FF is absent; if the
    # GROMACS data directory itself is not locatable, this is informational.
    ff_level = OPTIONAL if ff_status == WARN else REQUIRED
    checks.append(DoctorCheck(
        forcefield.upper() if forcefield.islower() else forcefield,
        ff_level, ff_status, detail=ff_detail,
    ))

    # ── RDKit ───────────────────────────────────────────────────────────────
    checks.append(DoctorCheck(
        "RDKit", OPTIONAL,
        AVAILABLE if rdkit_v else MISSING,
        version=rdkit_v,
        detail="chemical perception / hydrogenation / pose reconstruction"
        if rdkit_v else "not importable — ligand chemistry features degrade to heuristics",
    ))

    # ── Open Babel ──────────────────────────────────────────────────────────
    ob_v = openbabel_version()
    checks.append(DoctorCheck(
        "Open Babel", OPTIONAL,
        AVAILABLE if ob_v else MISSING,
        version=ob_v,
        detail="fallback hydrogenation backend" if ob_v else "not found (optional)",
    ))

    return DoctorReport(checks=checks)
