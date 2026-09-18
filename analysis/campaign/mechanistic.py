"""Dry-run planning for the standardized mechanistic descriptor profile.

This module is **planning only** — it never invokes GROMACS and never writes
inside a scientific source tree.  It answers one question for every
(system, observable) pair:

* what is the *source trajectory* (always the PBC-corrected rot+trans fit
  ``mdfit.xtc`` — raw ``md.xtc`` is never an analysis input), and
* what is its *status*:

  ``AVAILABLE``           a valid output already exists **and** was produced on
                          ``mdfit.xtc`` under this protocol — reuse it.
  ``SKIP_EXISTING``       synonym of AVAILABLE for callers that prefer it.
  ``PLANNED``             ``mdfit.xtc`` is present; the command can be built and
                          run when execution is authorised.
  ``REVIEW_REQUIRED``     ``mdfit.xtc`` (or the TPR / index) is missing, or the
                          system identity is unresolved — blocked, no fallback.
  ``SUPERSEDED_PENDING``  a historical output exists but was computed from raw
                          ``md.xtc`` — it must be recomputed on ``mdfit.xtc``.
  ``NEW_EXTENSION``       an analysis with no historical precedent for A1
                          (clustering / conformational-state populations).

The standardized protocol itself is documented in
``reports/trajectory_protocol_correction_2026-09-09/mechanistic_standardized_protocol.yaml``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.manifest import build_manifest
from analysis.campaign.xanthone_short import _fit_trajectory, FIT_TRAJECTORY_NAME

PROFILE = "mechanistic"
PROTOCOL_VERSION = "mechanistic/1.0-draft"

AVAILABLE = "AVAILABLE"
SKIP_EXISTING = "SKIP_EXISTING"
PLANNED = "PLANNED"
REVIEW_REQUIRED = "REVIEW_REQUIRED"
SUPERSEDED_PENDING = "SUPERSEDED_PENDING"
NEW_EXTENSION = "NEW_EXTENSION"


@dataclass(frozen=True)
class Observable:
    key: str
    tool: str
    #: glob patterns (relative to the system dir) that a *historical* output for
    #: this observable would match — used to detect raw-``md.xtc`` results that
    #: must be superseded.
    historical_globs: tuple[str, ...]
    needs_ligand: bool = False
    new_extension: bool = False
    note: str = ""


OBSERVABLES: tuple[Observable, ...] = (
    Observable("protein_rmsd", "gmx rms", ("rmsd*.xvg", "*rmsd*.xvg")),
    Observable("ligand_rmsd", "gmx rms", ("rmsd_*LIG*.xvg", "*lig*rmsd*.xvg",
                                          "rmsd_*-HMG.xvg"), needs_ligand=True),
    Observable("rmsf", "gmx rmsf", ("rmsf*.xvg",)),
    Observable("radius_of_gyration", "gmx gyrate", ("rg*.xvg", "gyrate*.xvg")),
    Observable("sasa", "gmx sasa", ("sasa*.xvg", "area*.xvg")),
    Observable("hbonds", "gmx hbond", ("hb*.xvg", "hbond*.xvg", "hbnum*.xvg")),
    Observable("active_site_mindist", "gmx mindist",
               ("mindist*active*.xvg", "mindist_lig_active.xvg"), needs_ligand=True),
    Observable("active_site_contacts", "gmx mindist",
               ("contacts*active*.xvg", "contacts_lig_active.xvg"), needs_ligand=True),
    Observable("catalytic_com_distance", "gmx distance",
               ("dist*catalytic*.xvg", "dist_lig_catalytic.xvg"), needs_ligand=True),
    Observable("pca", "gmx covar + gmx anaeig",
               ("eigenval*.xvg", "pc1*.xvg", "pc2*.xvg", "proj*.xvg", "PCA*.pdf")),
    Observable("fel", "gmx sham / free-energy landscape",
               ("FEL*.pdf", "fel*.xvg", "sham*.xpm")),
    Observable("dccm", "dynamic cross-correlation (CA)",
               ("DCCM*.pdf", "dccm*.dat", "covar*.dat")),
    Observable("mmpbsa", "g_mmpbsa run -pbsa",
               ("binding_energy.xvg", "energy_MM.xvg", "energy_summary.csv")),
    Observable("per_residue_decomposition", "g_mmpbsa run -decomp",
               ("residues_energy_summary.csv", "contrib_MM.dat")),
    # ── NEW extensions — no historical A1 precedent ────────────────────────────
    Observable("clustering", "gmx cluster (gromos)", (), new_extension=True,
               note="NEW: not historically present for A1"),
    Observable("conformational_state_populations",
               "state assignment from clustering + shared PCA basin occupancy",
               (), new_extension=True,
               note="NEW: not historically present for A1"),
)


@dataclass
class ObservablePlan:
    observable: str
    tool: str
    status: str
    source_trajectory: Optional[str]
    reason: str
    existing_outputs: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "observable": self.observable, "tool": self.tool,
            "status": self.status, "source_trajectory": self.source_trajectory,
            "reason": self.reason, "existing_outputs": self.existing_outputs,
        }


@dataclass
class SystemPlan:
    system: str
    directory: str
    role: str                       # "reference" | "new"
    fit_trajectory: Optional[str]
    raw_trajectory: Optional[str]
    tpr: Optional[str]
    index: Optional[str]
    observables: list[ObservablePlan] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "system": self.system, "directory": self.directory, "role": self.role,
            "analysis_trajectory": self.fit_trajectory,
            "analysis_trajectory_name": FIT_TRAJECTORY_NAME,
            "raw_production_trajectory": self.raw_trajectory,
            "tpr": self.tpr, "index": self.index,
            "observables": [o.to_dict() for o in self.observables],
        }


def _cmdline_trajectory(xvg: Path) -> Optional[str]:
    """The ``-f`` argument recorded in a GROMACS/g_mmpbsa XVG header, if any."""
    if not xvg.is_file():
        return None
    try:
        head = xvg.read_text(errors="replace")[:4000]
    except OSError:
        return None
    for line in head.splitlines():
        if "-f " in line and ("Command line" in line or line.strip().startswith("#")):
            toks = line.replace("'", " ").split()
            if "-f" in toks:
                return toks[toks.index("-f") + 1]
    return None


def _has_fit_provenance(system_dir: Path, obs_key: str) -> bool:
    """Evidence that *obs_key* was computed on mdfit.xtc.

    A SimForge/manual sidecar, or — for the g_mmpbsa observables — the ``-f``
    argument captured in the ``energy_MM.xvg`` / ``binding_energy.xvg`` header.
    """
    for name in (f"provenance.{obs_key}.json", f"{obs_key}.mdfit.json"):
        if (system_dir / name).is_file():
            return True
    if obs_key in ("mmpbsa", "per_residue_decomposition"):
        for xvg in ("energy_MM.xvg", "binding_energy.xvg"):
            f = _cmdline_trajectory(system_dir / xvg)
            if f and Path(f).name == FIT_TRAJECTORY_NAME:
                return True
    return False


def _matches(system_dir: Path, globs: tuple[str, ...]) -> list[str]:
    hits: list[str] = []
    for g in globs:
        hits += [p.name for p in sorted(system_dir.glob(g)) if p.is_file()]
    return sorted(set(hits))


def plan_system(system_dir: Path, *, role: str, gmx: str = "gmx",
                system_name: Optional[str] = None,
                fit_trajectory: Optional[str] = None,
                raw_trajectory: Optional[str] = None,
                has_ligand: bool = True) -> SystemPlan:
    system_dir = Path(system_dir)
    name = system_name or system_dir.name
    fit = fit_trajectory or (
        str((system_dir / FIT_TRAJECTORY_NAME).resolve())
        if (system_dir / FIT_TRAJECTORY_NAME).is_file() else None)
    raw = raw_trajectory or (
        str((system_dir / "md.xtc").resolve())
        if (system_dir / "md.xtc").is_file() else None)
    tpr = str((system_dir / "md.tpr").resolve()) if (system_dir / "md.tpr").is_file() else None
    ndx = str((system_dir / "index.ndx").resolve()) if (system_dir / "index.ndx").is_file() else None

    plan = SystemPlan(system=name, directory=str(system_dir.resolve()), role=role,
                      fit_trajectory=fit, raw_trajectory=raw, tpr=tpr, index=ndx)

    base_blocked = None
    if fit is None:
        base_blocked = (f"{FIT_TRAJECTORY_NAME} not present — canonical PBC-corrected "
                        f"rot+trans fit trajectory required; raw md.xtc is not accepted")
    elif tpr is None:
        base_blocked = "md.tpr not present"
    elif ndx is None:
        base_blocked = "index.ndx not present"

    for obs in OBSERVABLES:
        existing = _matches(system_dir, obs.historical_globs)
        if obs.new_extension:
            plan.observables.append(ObservablePlan(
                obs.key, obs.tool, NEW_EXTENSION, fit,
                obs.note + "; source trajectory = mdfit.xtc when executed",
                existing))
            continue
        if obs.needs_ligand and not has_ligand:
            plan.observables.append(ObservablePlan(
                obs.key, obs.tool, REVIEW_REQUIRED, fit,
                "observable needs a ligand; none resolved for this system",
                existing))
            continue
        if existing and _has_fit_provenance(system_dir, obs.key):
            plan.observables.append(ObservablePlan(
                obs.key, obs.tool, AVAILABLE, fit,
                "valid mdfit.xtc-based output already present", existing))
            continue
        if existing:
            # A historical result exists but was computed on raw md.xtc — it is
            # superseded regardless of whether mdfit.xtc has been regenerated yet.
            reason = ("historical output present but computed on raw md.xtc "
                      "(no mdfit.xtc provenance) — recompute on mdfit.xtc")
            if base_blocked:
                reason += f"; first: {base_blocked}"
            plan.observables.append(ObservablePlan(
                obs.key, obs.tool, SUPERSEDED_PENDING, fit, reason, existing))
            continue
        if base_blocked:
            plan.observables.append(ObservablePlan(
                obs.key, obs.tool, REVIEW_REQUIRED, fit, base_blocked, existing))
            continue
        plan.observables.append(ObservablePlan(
            obs.key, obs.tool, PLANNED, fit,
            "mdfit.xtc present; command buildable when execution authorised",
            existing))
    return plan


# names that indicate an apo / cofactor-only system (no small-molecule ligand)
_APO_MARKERS = ("apo", "comp")


def plan_mechanistic(study_root: str | Path, *, reference: Optional[str | Path] = None,
                     gmx: str = "gmx") -> dict:
    """Dry-run mechanistic plan for *study_root* (new systems) plus the optional
    historical *reference* study (A1 / APO / COA / COMP).  Never executes."""
    root = Path(study_root).resolve()
    systems: list[SystemPlan] = []

    if reference:
        ref = Path(reference).resolve()
        for sub in sorted(p for p in ref.iterdir() if p.is_dir()):
            if not (sub / "md.xtc").is_file() and not (sub / FIT_TRAJECTORY_NAME).is_file():
                continue
            has_lig = not any(m in sub.name.lower() for m in ("apo",))
            systems.append(plan_system(
                sub, role="reference", gmx=gmx,
                has_ligand=has_lig and sub.name.lower() != "apo"))

    if root.is_dir():
        manifest = build_manifest(root, inspect_trajectories=False, gmx=gmx)
        for rec in manifest.systems:
            legacy = rec.legacy_study or {}
            sdir = Path(legacy.get("directory") or rec.system_id)
            name = legacy.get("system") or rec.system_id
            fit = _fit_trajectory(rec, sdir if sdir.is_dir() else None)
            raw = next((p for p in rec.production_trajectory_paths
                        if Path(p).is_file()), None)
            has_lig = bool((legacy.get("selections", {}) or {}).get("ligand", {}).get("group_id") is not None) \
                or (sdir / "index.ndx").is_file() and "LIG" in _index_names(sdir / "index.ndx")
            systems.append(plan_system(
                sdir, role="new", gmx=gmx, system_name=name,
                fit_trajectory=fit, raw_trajectory=raw, has_ligand=has_lig))

    counts: dict[str, int] = {}
    for s in systems:
        for o in s.observables:
            counts[o.status] = counts.get(o.status, 0) + 1

    return {
        "profile": PROFILE,
        "protocol_version": PROTOCOL_VERSION,
        "analysis_trajectory": FIT_TRAJECTORY_NAME,
        "raw_md_xtc_is_forbidden_as_analysis_input": True,
        "study_root": str(root),
        "reference_study": str(Path(reference).resolve()) if reference else None,
        "execution": "disabled; dry-run planning only",
        "status_counts": counts,
        "systems": [s.to_dict() for s in systems],
    }


def _index_names(path: Path) -> set[str]:
    try:
        from analysis.campaign.structure.index_groups import parse_index_groups
        return {g.name for g in parse_index_groups(str(path))}
    except Exception:  # noqa: BLE001
        return set()
