"""Controlled execution of the standardized mechanistic descriptor suite.

Companion to :mod:`analysis.campaign.mechanistic` (the dry-run planner).  This
module *executes* GROMACS against real mechanistic trajectories, with the same
guarantees as :mod:`analysis.campaign.xanthone_short`:

* the analysis trajectory is the PBC-corrected rot+trans fit **``mdfit.xtc``**
  only — raw ``md.xtc`` is never an analysis input.  References that lack
  ``mdfit.xtc`` are reprocessed **into the managed output tree**
  (``md.xtc`` -> ``trjconv -pbc res -ur compact -center`` -> ``mdcenter.xtc``
  -> ``trjconv -fit rot+trans`` -> ``mdfit.xtc``); the source directory is only
  ever read.
* every analysis writes a provenance JSON (source trajectory + sha, tpr, index,
  semantic groups by name, argv, gmx return code, output sha).
* the scientific source trees are snapshotted before/after and verified
  byte-unchanged.
* semantic groups are resolved by *name* from each system's ``index.ndx``.

The protocol is `reports/trajectory_protocol_correction_2026-09-09/
mechanistic_standardized_protocol.yaml`.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.gmx import run_gmx, gmx_available, gmx_version
from analysis.campaign.xanthone_short import snapshot_tree, verify_integrity, _sha256
from analysis.campaign.trajectory.inspector import _parse_gmx_check

PROFILE = "xanthone_mechanistic"
PROTOCOL_VERSION = "mechanistic/1.0"
PROVENANCE_SCHEMA = "simforge/study-campaign/mechanistic/provenance/v1"
FIT_TRAJECTORY_NAME = "mdfit.xtc"

#: common analysis window and resampling stride for cross-system comparability
WINDOW_PS = (0.0, 200000.0)
COMMON_STRIDE_PS = 200

MECHANISTIC_SOURCES = {
    "Mecanismo_inhibitorio": "/home/jesusxd/Escritorio/Mecanismo_inhibitorio",
    "Nuevos_sistemas": "/home/jesusxd/Escritorio/Nuevos_sistemas",
}


@dataclass(frozen=True)
class SystemDef:
    name: str
    directory: str
    role: str                       # "reference" | "new"
    receptor: str
    ligand_group: Optional[str]     # xanthone / substrate group name in index.ndx, or None
    ligand_kind: Optional[str]      # "xanthone" | "substrate" | None
    coa_group: Optional[str]        # "COA" when a CoA molecule is present
    active_site_group: Optional[str]
    catalytic_group: Optional[str]
    comparison_label: Optional[str]  # "A1" / "A3" / "A6" for the shared-basis trio
    needs_reprocess: bool            # True when the source has no mdfit.xtc


SYSTEMS: tuple[SystemDef, ...] = (
    SystemDef("A1", f"{MECHANISTIC_SOURCES['Mecanismo_inhibitorio']}/A1", "reference",
              "HMG-CoA reductase", "LIG", "xanthone", None,
              "ActiveSite_HMG", "Catalytic_HMG", "A1", True),
    SystemDef("APO", f"{MECHANISTIC_SOURCES['Mecanismo_inhibitorio']}/APO", "reference",
              "HMG-CoA reductase", None, None, None, None, None, None, True),
    SystemDef("COA", f"{MECHANISTIC_SOURCES['Mecanismo_inhibitorio']}/COA", "reference",
              "HMG-CoA reductase", "COA", "substrate", "COA",
              "ActiveSite_HMG", "Catalytic_HMG", None, True),
    SystemDef("COMP", f"{MECHANISTIC_SOURCES['Mecanismo_inhibitorio']}/COMP", "reference",
              "HMG-CoA reductase", "LIG", "xanthone", "COA",
              "ActiveSite_HMG", "Catalytic_HMG", None, True),
    SystemDef("A3-HMG-R", f"{MECHANISTIC_SOURCES['Nuevos_sistemas']}/A3-HMG-R", "new",
              "HMG-CoA reductase", "LIG", "xanthone", None,
              "ActiveSite_HMG", "Catalytic_HMG", "A3", False),
    SystemDef("HMG-R-200ns-A6", f"{MECHANISTIC_SOURCES['Nuevos_sistemas']}/HMG-R-200ns-A6", "new",
              "HMG-CoA reductase", "LIG", "xanthone", None,
              "ActiveSite_HMG", "Catalytic_HMG", "A6", False),
    SystemDef("system_A3_COA", f"{MECHANISTIC_SOURCES['Nuevos_sistemas']}/system_A3_COA", "new",
              "HMG-CoA reductase", "LIG", "xanthone", "COA",
              None, None, None, False),
    SystemDef("system_A6_COA", f"{MECHANISTIC_SOURCES['Nuevos_sistemas']}/system_A6_COA", "new",
              "HMG-CoA reductase", "LIG", "xanthone", "COA",
              None, None, None, False),
)

SYSTEM_BY_NAME = {s.name: s for s in SYSTEMS}
SHARED_BASIS_TRIO = ("A1", "A3-HMG-R", "HMG-R-200ns-A6")   # labels A1 / A3 / A6


# ═══════════════════════════════════════════════════════════════════════════════
# Compatibility pre-flight
# ═══════════════════════════════════════════════════════════════════════════════

def _index_group_atoms(index_path: Path) -> dict[str, int]:
    """{group name: atom count} for every group in an index.ndx."""
    out: dict[str, int] = {}
    if not index_path.is_file():
        return out
    name = None
    n = 0
    for line in index_path.read_text(errors="replace").splitlines():
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            if name is not None:
                out[name] = n
            name = s[1:-1].strip()
            n = 0
        elif s:
            n += len(s.split())
    if name is not None:
        out[name] = n
    return out


def _check_span(traj: Path, gmx: str) -> dict:
    res = run_gmx(["check", "-f", str(traj)], gmx=gmx, timeout=7200)
    if not res.ok:
        return {"ok": False, "error": res.stderr.strip()[-300:]}
    p = _parse_gmx_check(res.stderr + "\n" + res.stdout)
    end = p.get("end_time_ps")
    dt = p.get("dt_ps")
    n = p.get("n_frames")
    if end is None and n and dt is not None:
        end = p.get("start_time_ps", 0.0) + (n - 1) * dt
    return {"ok": True, "n_frames": n, "dt_ps": dt,
            "start_ps": p.get("start_time_ps"), "end_ps": end,
            "atoms": p.get("atom_count")}


@dataclass
class Compat:
    system: str
    role: str
    status: str = "OK"                       # OK | REVIEW_REQUIRED
    reasons: list[str] = field(default_factory=list)
    fit_source: str = ""                     # source mdfit.xtc, or "reprocess"
    raw_md_xtc: Optional[str] = None
    md_tpr: Optional[str] = None
    index_ndx: Optional[str] = None
    n_frames: Optional[int] = None
    dt_ps: Optional[float] = None
    span_ps: Optional[float] = None
    traj_atoms: Optional[int] = None
    protein_atoms: Optional[int] = None
    c_alpha: Optional[int] = None
    ligand_group: Optional[str] = None
    ligand_atoms: Optional[int] = None
    coa_group: Optional[str] = None
    coa_atoms: Optional[int] = None
    active_site_group: Optional[str] = None
    catalytic_group: Optional[str] = None

    def to_row(self) -> dict:
        return {k: ("" if v is None else v) for k, v in self.__dict__.items()
                if k != "reasons"} | {"reasons": "; ".join(self.reasons)}


def preflight(system: SystemDef, *, gmx: str = "gmx") -> Compat:
    d = Path(system.directory)
    c = Compat(system=system.name, role=system.role)
    c.active_site_group = system.active_site_group
    c.catalytic_group = system.catalytic_group
    c.ligand_group = system.ligand_group
    c.coa_group = system.coa_group

    if not d.is_dir():
        c.status, c.reasons = "REVIEW_REQUIRED", [f"source dir missing: {d}"]
        return c

    raw = d / "md.xtc"
    mdc = d / "mdcenter.xtc"
    fit = d / FIT_TRAJECTORY_NAME
    tpr = d / "md.tpr"
    ndx = d / "index.ndx"
    c.raw_md_xtc = str(raw) if raw.is_file() else None
    c.md_tpr = str(tpr) if tpr.is_file() else None
    c.index_ndx = str(ndx) if ndx.is_file() else None

    if not tpr.is_file():
        c.reasons.append("md.tpr missing")
    if not ndx.is_file():
        c.reasons.append("index.ndx missing")

    # which trajectory will be the analysis input?
    if fit.is_file():
        c.fit_source = str(fit)
        probe = fit
    elif raw.is_file() or mdc.is_file():
        c.fit_source = "reprocess"
        probe = mdc if mdc.is_file() else raw
    else:
        c.reasons.append("no mdfit.xtc / mdcenter.xtc / md.xtc — cannot analyse")
        c.status = "REVIEW_REQUIRED"
        return c

    span = _check_span(probe, gmx) if gmx_available(gmx) else {"ok": False}
    if span.get("ok"):
        c.n_frames = span["n_frames"]
        c.dt_ps = span["dt_ps"]
        c.span_ps = span["end_ps"]
        c.traj_atoms = span["atoms"]
        if span["end_ps"] is not None and span["end_ps"] + (span["dt_ps"] or 0) < 190000:
            c.reasons.append(
                f"trajectory span {span['end_ps']:g} ps < ~200 ns — analysed over "
                f"its full length, flagged")

    groups = _index_group_atoms(ndx)
    c.protein_atoms = groups.get("Protein")
    c.c_alpha = groups.get("C-alpha")
    if system.ligand_group:
        c.ligand_atoms = groups.get(system.ligand_group)
        if c.ligand_atoms is None:
            c.reasons.append(f"ligand group '{system.ligand_group}' not in index.ndx")
    if system.coa_group:
        c.coa_atoms = groups.get(system.coa_group)
        if c.coa_atoms is None:
            c.reasons.append(f"CoA group '{system.coa_group}' not in index.ndx")
    for g in (system.active_site_group, system.catalytic_group):
        if g and g not in groups:
            c.reasons.append(f"site group '{g}' not in index.ndx "
                             f"(active-site/catalytic observables -> REVIEW_REQUIRED)")

    if not groups.get("C-alpha"):
        c.reasons.append("no C-alpha index group — PCA/DCCM blocked")
        c.status = "REVIEW_REQUIRED"

    # hard blockers
    if not tpr.is_file() or (not ndx.is_file() and system.name != "APO"):
        c.status = "REVIEW_REQUIRED"
    return c


# ═══════════════════════════════════════════════════════════════════════════════
# Provenance
# ═══════════════════════════════════════════════════════════════════════════════

def _fp(path: Optional[str], limit: int = 2 * 1024**3) -> dict:
    if not path:
        return {"present": False}
    p = Path(path)
    if not p.is_file():
        return {"path": path, "present": False}
    st = p.stat()
    return {"path": str(p.resolve()), "bytes": st.st_size, "mtime_ns": st.st_mtime_ns,
            "sha256": _sha256(p, limit=limit)}


def write_provenance(out: Path, *, analysis: str, system: str, argv: list[str],
                     stdin: str, source_trajectory: str, tpr: str, index: str,
                     semantic_groups: dict, gmx: str, returncode: int,
                     outputs: list[Path], extra: Optional[dict] = None,
                     raw_md_xtc: Optional[str] = None) -> dict:
    prov = {
        "schema": PROVENANCE_SCHEMA,
        "profile": PROFILE,
        "protocol_version": PROTOCOL_VERSION,
        "analysis": analysis,
        "system": system,
        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "gromacs_version": gmx_version(gmx),
        "source_trajectory": source_trajectory,
        "analysis_trajectory_role": "pbc-corrected rot+trans fit (mdfit.xtc), "
                                    "resampled to a common 200 ps stride for cross-system comparison",
        "raw_production_trajectory": raw_md_xtc,
        "preprocessing_chain": ["md.xtc",
                                "trjconv -pbc res -ur compact -center -> mdcenter.xtc",
                                "trjconv -fit rot+trans -> mdfit.xtc",
                                f"trjconv -dt {COMMON_STRIDE_PS} -o <analysis trajectory>"],
        "source_tpr": tpr,
        "index_file": index,
        "semantic_groups": semantic_groups,
        "analysis_window_ps": list(WINDOW_PS),
        "command_argv": argv,
        "command_stdin": stdin,
        "gmx_returncode": returncode,
        "source_fingerprints": {
            "analysis_trajectory": _fp(source_trajectory),
            "tpr": _fp(tpr),
            "index": _fp(index),
            "raw_production_trajectory": _fp(raw_md_xtc),
        },
        "outputs": [{"path": str(o.resolve()), "bytes": o.stat().st_size,
                     "sha256": _sha256(o)} for o in outputs if o.is_file()],
    }
    if extra:
        prov.update(extra)
    out.write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return prov


def shell_preview(argv: list[str], stdin: str = "") -> str:
    pre = f"printf %s {shlex.quote(stdin)} | " if stdin else ""
    return pre + " ".join(shlex.quote(a) for a in argv)
