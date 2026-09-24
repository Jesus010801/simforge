"""Execution context and read-only probes for trajectory diagnostics.

The context lazily provides — once, shared by all detectors:

* the Phase 1 :class:`FrameTimeIndex` (timestamps, MD steps, box per frame);
* per-component centre series (:class:`ComponentMotionSeries`) from one
  read-only ``gmx trajectory -nopbc`` pass per set of groups, cached on disk
  by trajectory content + group atom-set evidence + convention.

Nothing here writes next to, or transforms, the trajectory.  Probe outputs are
diagnostic intermediates in a temporary directory / the diagnostics cache,
never :class:`TrajectoryView` entries.
"""
from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from analysis.campaign.diagnostics.params import DiagnosticParams
from analysis.campaign.fingerprint import fingerprint_file
from analysis.campaign.gmx import gmx_available, gmx_version, run_gmx
from analysis.campaign.models import FrameTimeIndex, SourceFingerprint, ValidationState

MOTION_PROBE_SCHEMA = "simforge/diagnostic-probe/component-motion/v1"

#: roles whose centres the generic motion detectors consider (in this order)
MOTION_ROLES = ("system", "receptor", "peptide", "ligand", "protein_partner",
                "cofactor", "nucleic_acid", "complex", "membrane")


class ProbeError(RuntimeError):
    """A read-only probe could not produce trustworthy evidence."""


class CenterConvention:
    COM = "center_of_mass"               # masses from a .tpr
    COG = "center_of_geometry"           # no reliable masses (structure file only)


@dataclass
class ComponentMotionSeries:
    role: str
    group: str
    convention: str                      # CenterConvention.*
    positions: np.ndarray                # (n_frames, 3) nm, raw stored coordinates
    coordinates: str = "raw_stored_no_pbc_unwrapping"
    mass_source: Optional[str] = None
    group_evidence: dict = field(default_factory=dict)
    cache_key: str = ""
    command: list[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════════════
# Box geometry
# ═══════════════════════════════════════════════════════════════════════════════

def box_matrix(box9) -> np.ndarray:
    """GRO box order (v1x v2y v3z v1y v1z v2x v2z v3x v3y) → rows a, b, c (nm)."""
    b = [float(v) for v in box9]
    return np.array([[b[0], b[3], b[4]],
                     [b[5], b[1], b[6]],
                     [b[7], b[8], b[2]]])


def box_valid(B: np.ndarray) -> bool:
    return bool(np.isfinite(B).all()) and abs(np.linalg.det(B)) > 1e-9


def minimum_image(d: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(d_min_image, fractional_d)`` for displacement ``d`` in box ``B``.

    Rounds fractional coordinates (exact for rectangular boxes; the standard
    approximation for moderately skewed triclinic boxes).
    """
    f = d @ np.linalg.inv(B)
    return d - np.round(f) @ B, f


def min_edge(B: np.ndarray) -> float:
    return float(np.linalg.norm(B, axis=1).min())


# ═══════════════════════════════════════════════════════════════════════════════
# Context
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class DiagnosticContext:
    trajectory_path: str
    structure_path: Optional[str] = None          # -s for probes (.tpr gives real masses)
    index_path: Optional[str] = None              # semantic index
    groups: dict[str, str] = field(default_factory=dict)    # role -> index group name
    view_ref: Optional[str] = None
    membrane_present: Optional[bool] = None
    params: DiagnosticParams = field(default_factory=DiagnosticParams)
    gmx: str = "gmx"
    cache_dir: Optional[str] = None               # time index + probe cache (None = no disk cache)
    # injectable evidence (tests / callers that already hold it)
    time_index: Optional[FrameTimeIndex] = None
    motion_series: dict[str, ComponentMotionSeries] = field(default_factory=dict)
    probe_records: list[dict] = field(default_factory=list)
    _fingerprint: Optional[SourceFingerprint] = None
    _group_ev: dict = field(default_factory=dict)

    # ── identity ───────────────────────────────────────────────────────────
    def trajectory_fingerprint(self) -> Optional[SourceFingerprint]:
        if self._fingerprint is None and Path(self.trajectory_path).is_file():
            self._fingerprint = fingerprint_file(Path(self.trajectory_path))
        return self._fingerprint

    def group_evidence(self, role: str) -> Optional[dict]:
        if role in self._group_ev:
            return self._group_ev[role]
        ev = None
        name = self.groups.get(role)
        if name and self.index_path:
            from analysis.campaign.results import index_group_evidence
            ev = index_group_evidence(self.index_path, name)
        if ev is None and role in self.motion_series:
            ev = self.motion_series[role].group_evidence or None
        self._group_ev[role] = ev
        return ev

    def has_role(self, role: str) -> bool:
        return role in self.motion_series or (
            role in self.groups and self.group_evidence(role) is not None)

    # ── Phase 1 timeline ───────────────────────────────────────────────────
    def get_time_index(self) -> FrameTimeIndex:
        if self.time_index is None:
            from analysis.campaign.trajectory.time_index import (
                GromacsTimeIndexBackend, index_trajectory,
            )
            ti_cache = Path(self.cache_dir) / "time_index" if self.cache_dir else None
            self.time_index, cached = index_trajectory(
                self.trajectory_path, backend=GromacsTimeIndexBackend(self.gmx),
                cache_dir=ti_cache)
            self.probe_records.append({
                "probe": "frame_time_index", "backend": self.time_index.backend,
                "backend_version": self.time_index.backend_version,
                "from_cache": cached, "command": self.time_index.command})
        return self.time_index

    def boxes(self) -> Optional[np.ndarray]:
        ti = self.get_time_index()
        if not ti.boxes:
            return None
        return np.asarray(ti.boxes, dtype=float).reshape(-1, 9)

    # ── component motion probe ─────────────────────────────────────────────
    def convention(self) -> tuple[str, Optional[str]]:
        if self.structure_path and self.structure_path.lower().endswith(".tpr"):
            return CenterConvention.COM, "tpr"
        return CenterConvention.COG, None

    def motion(self, roles: list[str]) -> dict[str, ComponentMotionSeries]:
        """Centre series for ``roles`` (one gmx pass for all uncached roles)."""
        missing = [r for r in roles if r not in self.motion_series]
        if missing:
            self._probe_motion(missing)
        return {r: self.motion_series[r] for r in roles}

    def _motion_key(self, role: str, convention: str) -> str:
        from analysis.campaign.results import definition_token
        fp = self.trajectory_fingerprint()
        sfp = fingerprint_file(Path(self.structure_path)) if self.structure_path else None
        return definition_token({
            "probe": MOTION_PROBE_SCHEMA,
            "trajectory": fp.digest if fp else None,
            "group": self.group_evidence(role),
            "convention": convention,
            "structure": sfp.digest if sfp else None,
            "gmx": gmx_version(self.gmx),
        })

    def _cache_file(self, key: str) -> Optional[Path]:
        return Path(self.cache_dir) / "motion" / f"{key}.npz" if self.cache_dir else None

    def _probe_motion(self, roles: list[str]) -> None:
        n_frames = self.get_time_index().n_frames
        convention, mass_source = self.convention()
        todo: list[tuple[str, str, str]] = []
        for role in roles:
            group = self.groups.get(role)
            if not group or self.group_evidence(role) is None:
                raise ProbeError(f"no semantic index group for role '{role}'")
            key = self._motion_key(role, convention)
            cf = self._cache_file(key)
            if cf is not None and cf.is_file():
                with np.load(cf) as z:
                    pos = z["positions"]
                if pos.shape == (n_frames, 3):
                    self.motion_series[role] = ComponentMotionSeries(
                        role=role, group=group, convention=convention, positions=pos,
                        mass_source=mass_source, group_evidence=self.group_evidence(role),
                        cache_key=key)
                    self.probe_records.append({"probe": "component_motion", "role": role,
                                               "cache_key": key, "from_cache": True})
                    continue
            todo.append((role, group, key))
        if not todo:
            return
        if not self.structure_path or not self.index_path:
            raise ProbeError("component motion probe needs a structure (-s) and a semantic index")
        if not gmx_available(self.gmx):
            raise ProbeError(f"'{self.gmx}' not available")
        kw = "com" if convention == CenterConvention.COM else "cog"
        selections = [f'{kw} of group "{g}"' for _r, g, _k in todo]
        with tempfile.TemporaryDirectory(prefix="simforge_diag_") as tmp:
            out = Path(tmp) / "centres.xvg"
            args = ["trajectory", "-f", str(Path(self.trajectory_path).resolve()),
                    "-s", str(self.structure_path), "-n", str(self.index_path),
                    "-nopbc", "-ox", str(out), "-select", *selections]
            res = run_gmx(args, gmx=self.gmx, timeout=4 * 3600)
            if not res.ok or not out.is_file():
                raise ProbeError(f"gmx trajectory failed (rc={res.returncode}): "
                                 f"{res.stderr.strip()[-300:]}")
            rows = _read_xvg_rows(out)
        expected_cols = 1 + 3 * len(todo)
        if len(rows) != n_frames:
            raise ProbeError(f"motion probe returned {len(rows)} rows for {n_frames} indexed frames")
        if any(len(r) != expected_cols for r in rows):
            raise ProbeError("motion probe returned an unexpected number of columns")
        data = np.asarray(rows, dtype=float)
        for j, (role, group, key) in enumerate(todo):
            pos = data[:, 1 + 3 * j: 4 + 3 * j].copy()
            series = ComponentMotionSeries(
                role=role, group=group, convention=convention, positions=pos,
                mass_source=mass_source, group_evidence=self.group_evidence(role),
                cache_key=key, command=list(res.argv))
            self.motion_series[role] = series
            cf = self._cache_file(key)
            if cf is not None:
                cf.parent.mkdir(parents=True, exist_ok=True)
                tmpf = cf.with_suffix(".tmp.npz")
                np.savez(tmpf, positions=pos)
                os.replace(tmpf, cf)
                meta = {"schema": MOTION_PROBE_SCHEMA, "role": role, "group": group,
                        "convention": convention, "group_evidence": series.group_evidence,
                        "command": list(res.argv)}
                cf.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
            self.probe_records.append({"probe": "component_motion", "role": role,
                                       "group": group, "convention": convention,
                                       "cache_key": key, "from_cache": False,
                                       "command": list(res.argv)})


def _read_xvg_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        rows.append([float(v) for v in s.split()])
    return rows


def groups_from_semantic_index(index_path: str, components=None) -> dict[str, str]:
    """Role → group mapping for the standard semantic groups present in the index.

    Only groups that exist are mapped; a role whose component is not RESOLVED
    (when ``components`` is given) is left out.
    """
    from analysis.campaign.models import ClassificationState
    from analysis.campaign.structure.index_groups import parse_index_groups
    names = {g.name for g in parse_index_groups(index_path)}
    candidates = {"system": "System", "receptor": "Receptor", "peptide": "Peptide",
                  "protein_partner": "ProteinPartner", "complex": "Complex",
                  "membrane": "Membrane", "ligand": "Ligand", "cofactor": "Cofactor",
                  "nucleic_acid": "NucleicAcid"}
    resolved = None
    if components is not None:
        resolved = {c.component_type for c in components
                    if c.classification_state == ClassificationState.RESOLVED}
    out: dict[str, str] = {}
    for role, name in candidates.items():
        if name not in names:
            continue
        if resolved is not None and role not in ("system",) and role not in resolved:
            continue
        out[role] = name
    return out
