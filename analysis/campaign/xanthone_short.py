"""Controlled, opt-in execution of the validated ``xanthone_short`` thesis core.

This module is the *only* place in the study-campaign layer that invokes
GROMACS analysis tools against real scientific data.  It is deliberately narrow:

* Exactly five analyses are executed, and no others:
    - protein RMSD                    (``gmx rms``)
    - ligand RMSD                     (``gmx rms``)
    - ligand-active-site min distance (``gmx mindist``)
    - ligand-active-site contacts     (``gmx mindist``, same invocation)
    - ligand-catalytic-site COM dist  (``gmx distance``)
* The analysis window is always ``-b 0 -e 25000`` (picoseconds); a production
  trajectory longer than 25 ns is *analysed* over 0-25 ns, never truncated.
* The historical contact cutoff is ``-d 0.6`` (nm).
* Catalytic distance is ``COM(LIG)`` vs ``COM(Catalytic_*)`` via ``gmx distance``.
* Semantic groups are resolved by *name* from each system's ``index.ndx``
  (``Protein`` / ``LIG`` / ``ActiveSite_*`` / ``Catalytic_*``).  Numeric group
  ids are recorded for provenance but never fed to GROMACS.

Nothing here transforms coordinates and nothing is ever written inside a source
directory.  Original trajectories are only read.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.gmx import gmx_available, gmx_version, run_gmx
from analysis.campaign.legacy import PATTERNS, _valid_xvg, existing_results
from analysis.campaign.manifest import build_manifest
from analysis.campaign.models import SystemRecord
from analysis.campaign.structure.index_groups import parse_index_groups

PROFILE = "xanthone_short"
PROTOCOL_VERSION = "xanthone_short/1.0"
WINDOW_PS: tuple[float, float] = (0.0, 25000.0)
CONTACT_CUTOFF_NM = 0.6
PROVENANCE_SCHEMA = "simforge/study-campaign/xanthone-short/provenance/v1"

#: Trajectories larger than this are fingerprinted by (size, mtime) only — a
#: full SHA-256 of a multi-GB trajectory on every run is not worth the I/O.
_HASH_SIZE_LIMIT = 2 * 1024**3

# Status vocabulary (spec §5)
SKIP_EXISTING = "SKIP_EXISTING"
EXECUTED = "EXECUTED"
REVIEW_REQUIRED = "REVIEW_REQUIRED"
FAILED = "FAILED"

#: Every semantic role the profile can use.
REQUIRED_ROLES = ("protein", "ligand", "active_site", "catalytic_site")
#: Roles every task needs — gated at the system level.  ``active_site`` /
#: ``catalytic_site`` are gated per-task instead, so a system missing only a
#: site group can still run protein/ligand RMSD (spec §7: "do not run *that*
#: analysis").
IDENTITY_ROLES = ("protein", "ligand")

# Issues the legacy annotator may raise that an *explicit* ``--profile
# xanthone_short`` request is allowed to waive (spec §2/§3): they are about the
# folder name / intended duration, not identity or data integrity.
_WAIVED_ISSUES = (
    "folder duration conflicts with production metadata",
    "production duration requires verification",
    "insufficient evidence for study profile",
)
# Substrings that always block execution regardless of the requested profile.
_BLOCKING_ISSUES = (
    "conflicting compound identity",
    "conflicting target identity",
    "unresolved compound or target",
    "relies on folder hints",
    "missing production trajectory",
    "production trajectory path does not exist",
    "multiple index candidates",
    "missing index",
    "invalid index",
    "empty or out-of-range atoms",
)


# ═══════════════════════════════════════════════════════════════════════════════
# Task model
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Task:
    """One GROMACS invocation and the logical analyses it satisfies."""
    key: str
    tool: str                       # rms | mindist | distance
    analyses: tuple[str, ...]
    roles: tuple[str, ...]          # semantic roles, in GROMACS prompt order
    outputs: dict[str, str]         # logical name -> output filename


TASKS: tuple[Task, ...] = (
    Task("protein_rmsd", "rms", ("protein_rmsd",), ("protein", "protein"),
         {"rmsd": "rmsd_protein.xvg"}),
    Task("ligand_rmsd", "rms", ("ligand_rmsd",), ("protein", "ligand"),
         {"rmsd": "rmsd_ligand.xvg"}),
    Task("active_site", "mindist", ("active_site_mindist", "active_site_contacts"),
         ("ligand", "active_site"),
         {"mindist": "mindist_lig_active.xvg", "contacts": "contacts_lig_active.xvg"}),
    Task("catalytic_com_distance", "distance", ("catalytic_com_distance",),
         ("ligand", "catalytic_site"),
         {"dist": "dist_lig_catalytic.xvg"}),
)

ALL_ANALYSES: tuple[str, ...] = tuple(a for t in TASKS for a in t.analyses)


# ═══════════════════════════════════════════════════════════════════════════════
# Eligibility (identity / selection gating — spec §3, §7)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Eligibility:
    system: str
    eligible: bool
    reasons: list[str] = field(default_factory=list)          # why blocked
    missing_evidence: list[str] = field(default_factory=list)
    waived: list[str] = field(default_factory=list)
    override_applied: bool = False

    def to_dict(self) -> dict:
        return {
            "system": self.system, "eligible": self.eligible,
            "reasons": self.reasons, "missing_evidence": self.missing_evidence,
            "waived": self.waived, "override_applied": self.override_applied,
        }


def _folder_short_scope(name: str) -> bool:
    """The folder name itself claims a <=50 ns scope (e.g. ``HMG-R-25ns-A6``)."""
    m = re.search(r"(\d+(?:\.\d+)?)\s*ns", name, re.I)
    return bool(m) and float(m.group(1)) <= 50.0


def assess_eligibility(legacy: dict, rec: SystemRecord) -> Eligibility:
    """Decide whether a system may run the explicit ``xanthone_short`` core."""
    name = legacy.get("system", "?")
    el = Eligibility(system=name, eligible=True)

    compound = legacy.get("compound", {})
    target = legacy.get("target", {})
    selections = legacy.get("selections", {})
    profile = (legacy.get("profile") or {}).get("value")
    issues = list(legacy.get("issues", []))

    if legacy.get("cofactor", {}).get("value"):
        el.eligible = False
        el.reasons.append("cofactor system requires a dedicated mechanistic mapping")

    if not compound.get("value") or compound.get("confidence", 0) < 0.9:
        el.eligible = False
        el.missing_evidence.append(
            f"compound identity not structurally established "
            f"(value={compound.get('value')!r}, confidence={compound.get('confidence', 0)})")
    if not target.get("value") or target.get("confidence", 0) < 0.85:
        el.eligible = False
        el.missing_evidence.append(
            f"biological target not established from index headers "
            f"(value={target.get('value')!r}, confidence={target.get('confidence', 0)})")

    if not legacy.get("index"):
        el.eligible = False
        el.missing_evidence.append("no single canonical index.ndx resolved")

    for role in IDENTITY_ROLES:
        sel = selections.get(role, {})
        if sel.get("group_id") is None or not sel.get("confidence"):
            el.eligible = False
            el.missing_evidence.append(f"semantic group for role '{role}' not resolved")

    prod = [p for p in rec.production_trajectory_paths if Path(p).is_file()]
    if not prod:
        el.eligible = False
        el.missing_evidence.append("no readable production trajectory")

    for w in rec.warnings:
        if getattr(w, "severity", "") in ("review_required", "error"):
            el.eligible = False
            el.reasons.append(f"discovery warning: {w.message}")

    # Issue triage: waivable vs blocking.
    unresolved = []
    for issue in issues:
        low = issue.lower()
        if any(b in low for b in _BLOCKING_ISSUES):
            el.eligible = False
            el.reasons.append(issue)
        elif issue in _WAIVED_ISSUES or low.startswith("missing required groups"):
            # site-group gaps are handled per-task below, not at system level.
            el.waived.append(issue)
        else:
            unresolved.append(issue)

    # Profile gate: accept a native xanthone_short classification, or an explicit
    # override for a folder that itself claims short scope and whose only
    # discrepancies are the waived (naming / intended-duration) ones.
    if profile != PROFILE:
        if (el.eligible and _folder_short_scope(name)
                and not unresolved and el.waived):
            el.override_applied = True
        else:
            el.eligible = False
            el.reasons.append(
                f"classified profile is {profile!r}, not {PROFILE!r}, and the "
                f"folder does not unambiguously claim a <=50 ns scope"
                if not _folder_short_scope(name) else
                f"classified profile is {profile!r}; unresolved issues block the "
                f"explicit-{PROFILE} override: {unresolved}")

    for issue in unresolved:
        el.reasons.append(f"unresolved issue: {issue}")
        el.eligible = False

    return el


# ═══════════════════════════════════════════════════════════════════════════════
# Command construction
# ═══════════════════════════════════════════════════════════════════════════════

def _group_name(legacy: dict, role: str) -> Optional[str]:
    return legacy.get("selections", {}).get(role, {}).get("value")


def _group_id(legacy: dict, role: str):
    return legacy.get("selections", {}).get(role, {}).get("group_id")


def build_command(task: Task, *, tpr: str, traj: str, ndx: str,
                  legacy: dict, out_dir: Path, gmx: str = "gmx"
                  ) -> tuple[list[str], str]:
    """Return ``(argv, stdin)`` for *task*.  GROMACS is not invoked here."""
    b, e = WINDOW_PS
    base = [gmx, task.tool, "-s", tpr, "-f", traj, "-n", ndx,
            "-b", _fmt(b), "-e", _fmt(e)]
    names = [_group_name(legacy, r) for r in task.roles]
    if any(n is None for n in names):
        missing = [r for r, n in zip(task.roles, names) if n is None]
        raise ValueError(
            f"cannot build '{task.key}' command: unresolved semantic group(s) "
            f"for role(s) {missing} — a numeric/None id must never reach GROMACS")

    if task.tool == "rms":
        argv = base + ["-o", str(out_dir / task.outputs["rmsd"])]
        stdin = "".join(f"{n}\n" for n in names)
        return argv, stdin

    if task.tool == "mindist":
        argv = base + [
            "-d", _fmt(CONTACT_CUTOFF_NM),
            "-od", str(out_dir / task.outputs["mindist"]),
            "-on", str(out_dir / task.outputs["contacts"]),
        ]
        stdin = "".join(f"{n}\n" for n in names)
        return argv, stdin

    if task.tool == "distance":
        lig, cat = names
        select = f'com of group "{lig}" plus com of group "{cat}"'
        argv = base + ["-select", select, "-oall", str(out_dir / task.outputs["dist"])]
        return argv, ""

    raise ValueError(f"unknown tool {task.tool!r}")


def _fmt(x: float) -> str:
    return str(int(x)) if float(x).is_integer() else repr(x)


def shell_preview(argv: list[str], stdin: str) -> str:
    prefix = f"printf %s {shlex.quote(stdin)} | " if stdin else ""
    return prefix + " ".join(shlex.quote(a) for a in argv)


# ═══════════════════════════════════════════════════════════════════════════════
# Existing-result detection (missing-only semantics — spec §5)
# ═══════════════════════════════════════════════════════════════════════════════

def _existing_analysis(analysis: str, source_dir: Path, obs_dir: Path,
                       filename: str, source_families: dict) -> Optional[str]:
    """Return the path of a pre-existing valid result for *analysis*, or None."""
    if analysis in source_families and source_families[analysis]:
        return source_families[analysis][0]
    produced = obs_dir / filename
    header = PATTERNS.get(analysis, [(None, r".")])[0][1]
    if produced.is_file() and _valid_xvg(produced, header):
        return str(produced)
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Source-integrity verification (spec §11)
# ═══════════════════════════════════════════════════════════════════════════════

def _sha256(path: Path, *, limit: Optional[int] = None) -> Optional[str]:
    size = path.stat().st_size
    if limit is not None and size > limit:
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def snapshot_tree(root: Path) -> dict:
    """File inventory of *root*: relative path -> (size, mtime_ns).

    Plus a SHA-256 for every canonical ``index.ndx`` (GROMACS ``#...#`` backups
    are excluded — they are not canonical inputs)."""
    root = Path(root)
    files: dict[str, list] = {}
    index_hashes: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = str(p.relative_to(root))
        st = p.stat()
        files[rel] = [st.st_size, st.st_mtime_ns]
        if p.name == "index.ndx":
            index_hashes[rel] = _sha256(p)
    return {"root": str(root), "n_files": len(files),
            "files": files, "index_sha256": index_hashes}


@dataclass
class IntegrityReport:
    ok: bool
    added: list[str] = field(default_factory=list)
    removed: list[str] = field(default_factory=list)
    modified: list[str] = field(default_factory=list)
    index_changed: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "ok": self.ok, "added": self.added, "removed": self.removed,
            "modified": self.modified, "index_changed": self.index_changed,
        }


def verify_integrity(before: dict, after: dict) -> IntegrityReport:
    b, a = before["files"], after["files"]
    added = sorted(set(a) - set(b))
    removed = sorted(set(b) - set(a))
    modified = sorted(k for k in set(a) & set(b) if a[k] != b[k])
    index_changed = sorted(
        k for k in before["index_sha256"]
        if before["index_sha256"].get(k) != after["index_sha256"].get(k))
    ok = not (added or removed or modified or index_changed)
    return IntegrityReport(ok, added, removed, modified, index_changed)


# ═══════════════════════════════════════════════════════════════════════════════
# Output-path safety (spec §7)
# ═══════════════════════════════════════════════════════════════════════════════

def assert_safe_output(out_root: Path, study_root: Path,
                       system_dirs: list[Path]) -> None:
    out_root = Path(out_root).resolve()
    guarded = [Path(study_root).resolve(), *(Path(d).resolve() for d in system_dirs)]
    for g in guarded:
        if out_root == g or g in out_root.parents or out_root in g.parents:
            raise ValueError(
                f"unsafe output path: {out_root} overlaps the scientific source "
                f"tree at {g}")


# ═══════════════════════════════════════════════════════════════════════════════
# Execution
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AnalysisRun:
    analysis: str
    status: str
    task_key: str
    message: str = ""
    argv: list[str] = field(default_factory=list)
    stdin: str = ""
    output_files: list[str] = field(default_factory=list)
    checksums: dict[str, str] = field(default_factory=dict)
    provenance_path: Optional[str] = None
    xvg_check: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "analysis": self.analysis, "status": self.status,
            "task": self.task_key, "message": self.message,
            "command": shell_preview(self.argv, self.stdin) if self.argv else "",
            "argv": self.argv, "stdin": self.stdin,
            "output_files": self.output_files, "checksums": self.checksums,
            "provenance": self.provenance_path, "xvg_check": self.xvg_check,
        }


@dataclass
class SystemRun:
    system: str
    directory: str
    output_dir: str
    eligibility: dict
    window_ps: list[float]
    trajectory: Optional[str] = None
    tpr: Optional[str] = None
    index: Optional[str] = None
    measured_end_ps: Optional[float] = None
    analyses: list[AnalysisRun] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "system": self.system, "directory": self.directory,
            "output_dir": self.output_dir, "eligibility": self.eligibility,
            "analysis_window_ps": self.window_ps, "trajectory": self.trajectory,
            "tpr": self.tpr, "index": self.index,
            "measured_end_ps": self.measured_end_ps,
            "analyses": [a.to_dict() for a in self.analyses],
        }


def _coverage_ok(rec: SystemRecord) -> tuple[bool, Optional[float], str]:
    """Confirm the production trajectory spans at least 0-25 000 ps."""
    insp = rec.trajectory_inspection
    need = WINDOW_PS[1]
    if insp and insp.inspected and insp.end_time_ps is not None:
        dt = insp.dt_ps or 0.0
        if insp.end_time_ps + max(dt, 1.0) >= need:
            return True, insp.end_time_ps, "measured trajectory span"
        return False, insp.end_time_ps, (
            f"measured trajectory ends at {insp.end_time_ps:g} ps (< {need:g} ps)")
    if insp and insp.tpr_duration_ps and insp.tpr_duration_ps >= need:
        return True, None, "intended TPR duration (trajectory span not measured)"
    return False, None, "could not confirm 0-25 000 ps coverage (gmx check unavailable)"


def _resolve_index_names(index_path: str) -> set[str]:
    try:
        return {g.name for g in parse_index_groups(index_path)}
    except (OSError, ValueError):
        return set()


def _write_provenance(path: Path, *, task: Task, analysis: str, system: str,
                      argv: list[str], stdin: str, legacy: dict, rec: SystemRecord,
                      tpr: str, traj: str, ndx: str, gmx: str, returncode: int,
                      outputs: list[Path], measured_end_ps, source_fp: dict) -> dict:
    roles = list(dict.fromkeys(task.roles))
    prov = {
        "schema": PROVENANCE_SCHEMA,
        "profile": PROFILE,
        "protocol_version": PROTOCOL_VERSION,
        "analysis": analysis,
        "task": task.key,
        "system": system,
        "timestamp": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "gromacs_version": gmx_version(gmx),
        "source_trajectory": traj,
        "source_tpr": tpr,
        "index_file": ndx,
        "semantic_groups": {r: _group_name(legacy, r) for r in roles},
        "numeric_group_ids": {
            _group_name(legacy, r): _group_id(legacy, r) for r in roles},
        "analysis_window_ps": list(WINDOW_PS),
        "contact_cutoff_nm": CONTACT_CUTOFF_NM if task.tool == "mindist" else None,
        "catalytic_semantics": (
            'COM(LIG) vs COM(Catalytic_*) via gmx distance -select'
            if task.tool == "distance" else None),
        "command_argv": argv,
        "command_stdin": stdin,
        "gmx_returncode": returncode,
        "source_fingerprints": source_fp,
        "trajectory_measured_end_ps": measured_end_ps,
        "production_intended_duration_ns": (legacy.get("production_duration") or {}),
        "compound": legacy.get("compound"),
        "target": legacy.get("target"),
        "outputs": [
            {"path": str(o.resolve()), "bytes": o.stat().st_size,
             "sha256": _sha256(o)} for o in outputs if o.is_file()],
    }
    path.write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return prov


def _source_fingerprints(traj: str, tpr: str, ndx: str) -> dict:
    out = {}
    for label, p in (("trajectory", traj), ("tpr", tpr), ("index", ndx)):
        fp = Path(p)
        if not fp.is_file():
            out[label] = {"path": p, "present": False}
            continue
        st = fp.stat()
        out[label] = {
            "path": str(fp.resolve()), "bytes": st.st_size,
            "mtime_ns": st.st_mtime_ns,
            "sha256": _sha256(fp, limit=_HASH_SIZE_LIMIT),
        }
    return out


def _post_validate_xvg(path: Path, dt_hint: Optional[float]) -> dict:
    """Structural / range / finiteness checks on a produced XVG (spec §10)."""
    info: dict = {"path": str(path), "readable": False}
    try:
        from runtime.xvg_parser import parse_xvg
        data = parse_xvg(path)
    except Exception as exc:  # noqa: BLE001
        info["error"] = f"parse failed: {exc}"
        return info
    ts = [float(t) for t in getattr(data, "time_ps", [])]
    series = data.series[0].values if data.series else []
    ys = [float(v) for v in series]
    info["readable"] = True
    info["n_points"] = len(ts)
    info["time_start_ps"] = ts[0] if ts else None
    info["time_end_ps"] = ts[-1] if ts else None
    info["all_finite"] = bool(ys) and all(v == v and abs(v) != float("inf") for v in ys)
    info["xlabel"] = data.xlabel
    info["ylabel"] = data.ylabel
    info["time_range_ok"] = bool(ts) and ts[0] <= 1.0 and abs(ts[-1] - WINDOW_PS[1]) <= max(dt_hint or 0.0, 10.0)
    if ts and len(ts) > 1:
        step = round((ts[-1] - ts[0]) / (len(ts) - 1), 6)
        info["dt_ps"] = step
        info["expected_points_if_10ps"] = int(round(WINDOW_PS[1] / 10.0)) + 1
        info["points_match_10ps_sampling"] = len(ts) == info["expected_points_if_10ps"]
    return info


def run_system(rec: SystemRecord, *, out_root: Path, gmx: str,
               force: bool, dry_run: bool) -> SystemRun:
    legacy = rec.legacy_study or {}
    name = legacy.get("system") or rec.system_id
    directory = legacy.get("directory")
    source_dir = Path(directory).resolve() if directory else None
    obs_root = Path(out_root) / name / "observables"

    el = assess_eligibility(legacy, rec)
    if not source_dir or not source_dir.is_dir():
        el.eligible = False
        el.reasons.append("legacy annotation has no resolved source directory")
    run = SystemRun(system=name, directory=str(source_dir or ""),
                    output_dir=str(Path(out_root) / name),
                    eligibility=el.to_dict(), window_ps=list(WINDOW_PS))

    traj = next((p for p in rec.production_trajectory_paths if Path(p).is_file()), None)
    # Reference is the *production* topology only — the historical protocol used
    # ``md.tpr``.  Never silently fall back to an equilibration/EM ``.tpr``.
    tpr = (str(source_dir / "md.tpr")
           if source_dir and (source_dir / "md.tpr").is_file() else None)
    ndx = legacy.get("index")
    run.trajectory, run.tpr, run.index = traj, tpr, ndx

    source_families, _ = (existing_results(source_dir)
                          if source_dir and source_dir.is_dir() else ({}, []))

    # ---- system-wide safety gates (spec §7) --------------------------------
    blockers: list[str] = []
    if not el.eligible:
        blockers.append("system not eligible: " + "; ".join(
            el.reasons + el.missing_evidence))
    if not traj:
        blockers.append("no production trajectory file")
    if not tpr or not Path(tpr).is_file():
        blockers.append("md.tpr not found")
    if not ndx or not Path(ndx).is_file():
        blockers.append("index.ndx not found")

    cov_ok, measured_end, cov_msg = _coverage_ok(rec)
    run.measured_end_ps = measured_end
    if not cov_ok:
        blockers.append(cov_msg)

    if tpr and Path(tpr).is_file() and gmx_available(gmx) and not dry_run:
        insp = rec.trajectory_inspection
        if not (insp and insp.tpr_duration_ps):
            # inspector did not recover the intended duration — probe directly.
            probe = run_gmx(["dump", "-s", tpr], gmx=gmx, timeout=300)
            if not probe.ok or "inputrec" not in probe.stdout:
                blockers.append(
                    f"gmx dump could not read the TPR (rc={probe.returncode})")

    index_names = _resolve_index_names(ndx) if ndx else set()

    dt_hint = rec.trajectory_inspection.dt_ps if rec.trajectory_inspection else None

    for task in TASKS:
        obs_dir = obs_root / task.key
        needed_groups = {_group_name(legacy, r) for r in task.roles}
        task_blockers = list(blockers)
        missing_groups = sorted(g for g in needed_groups if g and g not in index_names)
        if missing_groups:
            task_blockers.append(f"index.ndx missing group(s): {missing_groups}")
        if any(_group_name(legacy, r) is None for r in task.roles):
            task_blockers.append("unresolved semantic group for this task")

        # existing-result check per logical analysis
        existing = {
            a: _existing_analysis(a, source_dir, obs_dir,
                                  _analysis_filename(task, a), source_families)
            for a in task.analyses
        }
        all_existing = all(existing.values())

        if all_existing and not force:
            for a in task.analyses:
                run.analyses.append(AnalysisRun(
                    a, SKIP_EXISTING, task.key,
                    message=f"validated result already present: {existing[a]}"))
            continue

        if task_blockers:
            for a in task.analyses:
                run.analyses.append(AnalysisRun(
                    a, REVIEW_REQUIRED, task.key,
                    message="; ".join(dict.fromkeys(task_blockers))))
            continue

        argv, stdin = build_command(task, tpr=tpr, traj=traj, ndx=ndx,
                                    legacy=legacy, out_dir=obs_dir, gmx=gmx)

        if dry_run:
            for a in task.analyses:
                run.analyses.append(AnalysisRun(
                    a, "PLANNED", task.key, argv=argv, stdin=stdin,
                    message="dry-run: command constructed, not executed"))
            continue

        obs_dir.mkdir(parents=True, exist_ok=True)
        res = run_gmx(argv[1:], stdin=stdin or None, gmx=gmx, timeout=14400)
        (obs_dir / "gmx.log").write_text(
            f"$ {shell_preview(argv, stdin)}\n\n"
            f"return code: {res.returncode}\n\nSTDERR:\n{res.stderr}\n\nSTDOUT:\n{res.stdout}\n")

        produced = [obs_dir / fn for fn in task.outputs.values()]
        present = [p for p in produced if p.is_file()]
        source_fp = _source_fingerprints(traj, tpr, ndx)

        for a in task.analyses:
            fn = _analysis_filename(task, a)
            out_path = obs_dir / fn
            ar = AnalysisRun(a, EXECUTED, task.key, argv=argv, stdin=stdin)
            if not res.ok or not out_path.is_file():
                ar.status = FAILED
                ar.message = (f"gmx {task.tool} failed (rc={res.returncode}): "
                              f"{res.stderr.strip()[-400:]}")
            else:
                ar.output_files = [str(out_path.resolve())]
                ar.checksums = {out_path.name: _sha256(out_path)}
                ar.xvg_check = _post_validate_xvg(out_path, dt_hint)
                ar.message = f"executed: {shell_preview(argv, stdin)}"
            prov_path = obs_dir / f"provenance.{a}.json"
            _write_provenance(
                prov_path, task=task, analysis=a, system=name, argv=argv,
                stdin=stdin, legacy=legacy, rec=rec, tpr=tpr, traj=traj, ndx=ndx,
                gmx=gmx, returncode=res.returncode, outputs=present,
                measured_end_ps=measured_end, source_fp=source_fp)
            ar.provenance_path = str(prov_path.resolve())
            run.analyses.append(ar)

    # per-system run record
    sys_dir = Path(out_root) / name
    sys_dir.mkdir(parents=True, exist_ok=True)
    (sys_dir / "system_run.json").write_text(
        json.dumps(run.to_dict(), indent=2, sort_keys=True) + "\n")
    return run


def _analysis_filename(task: Task, analysis: str) -> str:
    if task.key == "active_site":
        return task.outputs["contacts" if analysis.endswith("contacts") else "mindist"]
    return next(iter(task.outputs.values()))


# ═══════════════════════════════════════════════════════════════════════════════
# Study-level driver
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class StudyRunReport:
    study_root: str
    output_dir: str
    profile: str
    gromacs_version: Optional[str]
    dry_run: bool
    systems: list[SystemRun] = field(default_factory=list)
    integrity: Optional[dict] = None
    skipped_systems: list[dict] = field(default_factory=list)

    def counts(self) -> dict:
        c = {EXECUTED: 0, SKIP_EXISTING: 0, REVIEW_REQUIRED: 0, FAILED: 0, "PLANNED": 0}
        for s in self.systems:
            for a in s.analyses:
                c[a.status] = c.get(a.status, 0) + 1
        return c

    def to_dict(self) -> dict:
        return {
            "study_root": self.study_root, "output_dir": self.output_dir,
            "profile": self.profile, "protocol_version": PROTOCOL_VERSION,
            "gromacs_version": self.gromacs_version, "dry_run": self.dry_run,
            "analysis_window_ps": list(WINDOW_PS),
            "contact_cutoff_nm": CONTACT_CUTOFF_NM,
            "counts": self.counts(),
            "systems": [s.to_dict() for s in self.systems],
            "skipped_systems": self.skipped_systems,
            "source_integrity": self.integrity,
        }


def run_study(study_root: str | Path, *, output_dir: Optional[str | Path] = None,
              gmx: str = "gmx", only: Optional[list[str]] = None,
              force: bool = False, dry_run: bool = False,
              profile: str = PROFILE) -> StudyRunReport:
    if profile != PROFILE:
        raise ValueError(f"only the {PROFILE!r} profile is implemented")
    root = Path(study_root).resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)

    out_root = (Path(output_dir).resolve() if output_dir
                else (Path.cwd() / "analysis_outputs" / PROFILE).resolve())

    before = snapshot_tree(root)

    manifest = build_manifest(root, inspect_trajectories=True, gmx=gmx)
    system_dirs = [Path((s.legacy_study or {}).get("directory") or s.system_id)
                   for s in manifest.systems]
    assert_safe_output(out_root, root, system_dirs)
    out_root.mkdir(parents=True, exist_ok=True)

    report = StudyRunReport(
        study_root=str(root), output_dir=str(out_root), profile=profile,
        gromacs_version=gmx_version(gmx), dry_run=dry_run)

    for rec in manifest.systems:
        legacy = rec.legacy_study or {}
        name = legacy.get("system") or rec.system_id
        if only and name not in only:
            report.skipped_systems.append(
                {"system": name, "reason": "not in --system filter"})
            continue
        run = run_system(rec, out_root=out_root, gmx=gmx, force=force, dry_run=dry_run)
        report.systems.append(run)

    after = snapshot_tree(root)
    integ = verify_integrity(before, after)
    report.integrity = {
        "before": before, "after": after, "verdict": integ.to_dict()}
    (out_root / "source_integrity.json").write_text(
        json.dumps(report.integrity, indent=2, sort_keys=True) + "\n")
    (out_root / "run_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n")
    (out_root / "run_report.md").write_text(render_markdown(report))
    return report


def render_markdown(report: StudyRunReport) -> str:
    lines = [
        f"# xanthone_short controlled execution", "",
        f"- study root: `{report.study_root}`",
        f"- output dir: `{report.output_dir}`",
        f"- profile / protocol: `{report.profile}` / `{PROTOCOL_VERSION}`",
        f"- GROMACS: `{report.gromacs_version}`",
        f"- analysis window: `-b 0 -e 25000` ps  |  contact cutoff: `{CONTACT_CUTOFF_NM}` nm",
        f"- dry-run: {report.dry_run}", "",
        f"## Totals", "",
        "| status | count |", "| --- | --- |",
    ]
    for k, v in report.counts().items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    for s in report.systems:
        lines.append(f"## {s.system}")
        lines.append("")
        el = s.eligibility
        lines.append(f"- eligible: **{el['eligible']}**"
                     + (f" (override applied)" if el.get("override_applied") else ""))
        if el.get("waived"):
            lines.append(f"- waived issues: {el['waived']}")
        if el.get("reasons"):
            lines.append(f"- blocking: {el['reasons']}")
        if el.get("missing_evidence"):
            lines.append(f"- missing evidence: {el['missing_evidence']}")
        lines.append(f"- trajectory: `{s.trajectory}`")
        lines.append(f"- tpr: `{s.tpr}`")
        lines.append(f"- index: `{s.index}`")
        lines.append(f"- measured trajectory end: {s.measured_end_ps} ps")
        lines.append("")
        lines.append("| analysis | status | detail |")
        lines.append("| --- | --- | --- |")
        for a in s.analyses:
            detail = a.message.replace("\n", " ")[:160]
            lines.append(f"| {a.analysis} | {a.status} | {detail} |")
        lines.append("")
    integ = (report.integrity or {}).get("verdict", {})
    lines.append("## Source integrity")
    lines.append("")
    lines.append(f"- verdict: **{'UNCHANGED' if integ.get('ok') else 'CHANGED'}**")
    for key in ("added", "removed", "modified", "index_changed"):
        if integ.get(key):
            lines.append(f"- {key}: {integ[key]}")
    lines.append("")
    return "\n".join(lines)
