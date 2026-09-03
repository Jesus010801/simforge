"""
Formal provenance / reproducibility layer for SimForge system builds.

``collect_environment()`` snapshots the toolchain (SimForge, Python, git,
GROMACS, RDKit, Open Babel). ``build_provenance()`` combines that with the
normalized system spec, input/output file hashes, and the assembly decisions
(protein mode, pdb2gmx used, discovered topology files, atomtype rename maps,
warnings, validation status) into a single ``provenance.json``-shaped dict.

Everything here is best-effort and side-effect free: a missing tool degrades
to ``None``/``"unknown"`` rather than raising.
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


def _run(cmd: list[str], timeout: int = 10) -> Optional[str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if out.returncode != 0:
            return None
        return (out.stdout or out.stderr).strip()
    except Exception:
        return None


def simforge_version() -> str:
    try:
        from importlib.metadata import version

        return version("simforge")
    except Exception:
        pass
    try:
        pp = Path(__file__).resolve().parent.parent / "pyproject.toml"
        m = re.search(r'^version\s*=\s*"([^"]+)"', pp.read_text(), re.M)
        if m:
            return m.group(1)
    except Exception:
        pass
    return "0.1.0"


def git_info(repo_dir: Optional[Path] = None) -> dict[str, Any]:
    repo_dir = Path(repo_dir or Path(__file__).resolve().parent.parent)
    info: dict[str, Any] = {"commit": None, "branch": None, "dirty": None}
    if not shutil.which("git"):
        return info

    def g(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(
                ["git", "-C", str(repo_dir), *args],
                capture_output=True, text=True, timeout=10,
            )
            return out.stdout.strip() if out.returncode == 0 else None
        except Exception:
            return None

    info["commit"] = g("rev-parse", "HEAD")
    info["branch"] = g("rev-parse", "--abbrev-ref", "HEAD")
    status = g("status", "--porcelain")
    if status is not None:
        info["dirty"] = bool(status.strip())
    return info


def gromacs_version() -> Optional[str]:
    if not shutil.which("gmx"):
        return None
    out = _run(["gmx", "--version"])
    if not out:
        return None
    m = re.search(r"GROMACS version:?\s*(\S+)", out) or re.search(r"gmx,\s*([\d.]+)", out)
    return m.group(1) if m else None


def _normalize_top_dir(base: Path) -> Optional[Path]:
    """Given a candidate dir, return the one that actually holds ``*.ff/``."""
    for cand in (base, base / "top", base / "share" / "gromacs" / "top"):
        if cand.is_dir() and any(cand.glob("*.ff")):
            return cand
    return None


def gromacs_data_prefix() -> Optional[Path]:
    for env in ("GMXDATA", "GMXLIB"):
        v = os.environ.get(env)
        if v and Path(v).exists():
            found = _normalize_top_dir(Path(v))
            if found:
                return found
    out = _run(["gmx", "--version"]) if shutil.which("gmx") else None
    if out:
        m = re.search(r"Data prefix:\s*(\S+)", out)
        if m:
            found = _normalize_top_dir(Path(m.group(1)))
            if found:
                return found
    return None


def _module_version(mod: str) -> Optional[str]:
    try:
        m = __import__(mod)
        return getattr(m, "__version__", "unknown")
    except Exception:
        return None


def openbabel_version() -> Optional[str]:
    v = _module_version("openbabel")
    if v:
        return v
    if shutil.which("obabel"):
        out = _run(["obabel", "-V"])
        if out:
            m = re.search(r"([\d.]+)", out)
            return m.group(1) if m else out
    return None


def conda_env() -> Optional[str]:
    return os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("CONDA_PREFIX")


def collect_environment() -> dict[str, Any]:
    """Snapshot the tool + interpreter environment."""
    return {
        "simforge_version": simforge_version(),
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "conda_env": conda_env(),
        "git": git_info(),
        "gromacs_version": gromacs_version(),
        "rdkit_version": _module_version("rdkit"),
        "openbabel_version": openbabel_version(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def sha256_file(path: str | Path) -> Optional[str]:
    p = Path(path)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_files(paths: dict[str, Any]) -> dict[str, Any]:
    """``{label: path}`` -> ``{label: {path, sha256, bytes}}`` (skips ``None``)."""
    out: dict[str, Any] = {}
    for label, path in paths.items():
        if path is None:
            continue
        p = Path(path)
        out[label] = {
            "path": str(p),
            "exists": p.exists(),
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size if p.is_file() else None,
        }
    return out


def build_provenance(
    *,
    spec: Any,
    spec_path: Optional[Path],
    command: Optional[list[str]],
    inputs: dict[str, Any],
    outputs: dict[str, Any],
    assembly: dict[str, Any],
    validation: Optional[dict[str, Any]] = None,
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Assemble the full provenance record for one system build."""
    env = collect_environment()
    record: dict[str, Any] = {
        "schema": "simforge/provenance/v1",
        "generated_utc": env["timestamp_utc"],
        "environment": env,
        "command": command or [],
        "spec_path": str(spec_path) if spec_path else None,
        "system_spec": spec.normalized_dict() if hasattr(spec, "normalized_dict") else spec,
        "forcefield": getattr(spec, "forcefield", assembly.get("forcefield")),
        "water_model": getattr(spec, "water_model", assembly.get("water_model")),
        "protein_mode": assembly.get("protein_mode"),
        "pdb2gmx_used": assembly.get("pdb2gmx_used"),
        "pdb2gmx_ignh_fallback": assembly.get("pdb2gmx_ignh_fallback", False),
        "discovered_protein_topology": assembly.get("discovered_protein_topology", []),
        "discovery_notes": assembly.get("discovery_notes", []),
        "protein_molecules": assembly.get("protein_molecules", []),
        "component_molecule_names": assembly.get("component_molecule_names", {}),
        "atomtype_renames": assembly.get("atomtype_renames", {}),
        "total_atom_count": assembly.get("total_atom_count"),
        "warnings": assembly.get("warnings", []),
        "inputs": hash_files(inputs),
        "outputs": hash_files(outputs),
        "validation": validation or {},
    }
    if extra:
        record.update(extra)
    return record
