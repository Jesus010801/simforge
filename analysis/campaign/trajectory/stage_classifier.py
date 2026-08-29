"""Evidence-based MD-workflow stage classification for a trajectory.

``nvt.xtc`` / ``npt.xtc`` / ``md.xtc`` are *different workflow stages*, not
production segments.  This module classifies each trajectory into a
:class:`TrajectoryStage` from the strongest available combination of:

* filename stem semantics (``nvt`` / ``npt`` / ``md`` / ``prod`` / ``partNNNN`` …)
  — evidence, never absolute truth;
* the paired ``.mdp`` (``mdout.mdp`` or ``<stem>.mdp``): ``gen_vel``,
  pressure/temperature coupling, ``nsteps`` × ``dt`` (total time),
  ``nstxout-compressed`` / ``nstxtcout`` (output spacing);
* PBC-derivative name markers (``nojump`` / ``whole`` / ``center`` / ``fit``)
  → ``DERIVED``;
* frame spacing and duration relative to sibling trajectories in the same
  directory (production is typically the longest with the coarsest spacing).

Filenames alone are never decisive: a paired ``.mdp`` that says
``gen_vel = yes`` + no pressure coupling + 50 ps overrides a stem of ``md``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.models import TrajectoryStage

# stem-prefix -> (stage, weight)
_STEM_RULES: list[tuple[re.Pattern, str, float]] = [
    (re.compile(r"^(nvt|equil[_-]?nvt|eq[_-]?nvt)", re.I), TrajectoryStage.EQUILIBRATION_NVT, 0.55),
    (re.compile(r"^(npt|equil[_-]?npt|eq[_-]?npt)", re.I), TrajectoryStage.EQUILIBRATION_NPT, 0.55),
    (re.compile(r"^(em|min|minim|steep|cg)([_-]|\d|$)", re.I), TrajectoryStage.MINIMIZATION, 0.6),
    (re.compile(r"^(ions?|solv|genion|neutral)", re.I), TrajectoryStage.PREPARATION, 0.6),
    (re.compile(r"^(heat|anneal|therm|warmup|restrain|posre|relax)", re.I),
     TrajectoryStage.EQUILIBRATION_OTHER, 0.5),
    (re.compile(r"^(equil|eq)([_-]|\d|$)", re.I), TrajectoryStage.EQUILIBRATION_OTHER, 0.45),
    (re.compile(r"^(md|prod|production|run|sim|traj|dyn)", re.I), TrajectoryStage.PRODUCTION, 0.45),
]
_DERIVED_MARKERS = (
    "nojump", "nopbc", "no_pbc", "_pbc", "-pbc", "whole", "center", "centered",
    "_fit", "-fit", "fitted", "_mol", "aligned", "wrapped", "skip", "_dt",
)
_PART_RE = re.compile(r"(?:part|seg|chunk|cont|restart|extend)[\-_.]?(\d+)", re.I)


@dataclass
class StageResult:
    stage: str
    confidence: float
    evidence: list[str] = field(default_factory=list)
    part_index: Optional[int] = None
    derived_from: Optional[str] = None
    mdp_dt_out_ps: Optional[float] = None
    mdp_duration_ps: Optional[float] = None


def _find_mdp(traj: Path) -> Optional[Path]:
    d = traj.parent
    stem = traj.stem.lower()
    for cand in (d / f"{traj.stem}.mdp", d / f"{stem}.mdp", d / "mdout.mdp", d / "md.mdp"):
        if cand.is_file():
            # mdout.mdp / md.mdp only count for a production-stem trajectory
            if cand.name in ("mdout.mdp", "md.mdp") and not stem.startswith(("md", "prod", "run")):
                continue
            return cand
    return None


def _parse_mdp(path: Path) -> dict:
    out: dict = {}
    try:
        for raw in path.read_text(errors="replace").splitlines():
            line = raw.split(";", 1)[0].strip()
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip().lower().replace("_", "-")] = v.strip().lower()
    except OSError:
        return {}
    return out


def _mdp_num(mdp: dict, key: str) -> Optional[float]:
    try:
        return float(mdp[key].split()[0])
    except (KeyError, ValueError, IndexError):
        return None


def mdp_timing(mdp: dict) -> tuple[Optional[float], Optional[float]]:
    """(output spacing in ps, total duration in ps) from an mdp dict, or (None, None)."""
    dt = _mdp_num(mdp, "dt")
    nsteps = _mdp_num(mdp, "nsteps")
    nstout = (_mdp_num(mdp, "nstxout-compressed") or _mdp_num(mdp, "nstxtcout")
              or _mdp_num(mdp, "nstxout"))
    duration = nsteps * dt if (nsteps is not None and dt is not None) else None
    spacing = nstout * dt if (nstout is not None and dt is not None) else None
    return spacing, duration


def _mdp_evidence(mdp: dict) -> Optional[tuple[str, float, list[str]]]:
    if not mdp:
        return None
    ev: list[str] = []

    def _num(key: str) -> Optional[float]:
        return _mdp_num(mdp, key)

    integrator = mdp.get("integrator", "")
    if integrator in ("steep", "cg", "l-bfgs"):
        return TrajectoryStage.MINIMIZATION, 0.9, [f"mdp integrator={integrator}"]

    nsteps = _num("nsteps")
    dt = _num("dt")
    total_ps = nsteps * dt if (nsteps is not None and dt is not None) else None
    gen_vel = mdp.get("gen-vel", "no") in ("yes", "true", "1")
    pcoupl = mdp.get("pcoupl", "no")
    pcoupl_on = pcoupl not in ("", "no", "none")

    if total_ps is not None:
        ev.append(f"mdp total time = {total_ps:.0f} ps (nsteps={nsteps:.0f} × dt={dt})")

    # Long run => production
    if total_ps is not None and total_ps >= 2000:
        return TrajectoryStage.PRODUCTION, 0.6, ev + ["long trajectory (≥ 2 ns)"]

    # Fresh velocities => start of equilibration
    if gen_vel:
        ev.append("mdp gen_vel = yes (fresh velocities → equilibration start)")
        if not pcoupl_on:
            return TrajectoryStage.EQUILIBRATION_NVT, 0.7, ev + ["no pressure coupling"]
        return TrajectoryStage.EQUILIBRATION_NPT, 0.6, ev + [f"pcoupl = {pcoupl}"]

    # Short run, no fresh velocities: WEAK equilibration signal only.  A short
    # (<1 ns) run with gen_vel=no could equally be a short production run or an
    # extension, so this must never outweigh a production-suggesting filename
    # stem (weight 0.45).  The directory-relative "much longer sibling exists"
    # check (in classify_trajectory_stage) is the real discriminator.
    if total_ps is not None and total_ps < 1000:
        if pcoupl_on:
            return TrajectoryStage.EQUILIBRATION_NPT, 0.30, ev + [
                f"short run (<1 ns) with pressure coupling ({pcoupl}), no fresh velocities "
                f"— weak NPT-equilibration signal"]
        return TrajectoryStage.EQUILIBRATION_OTHER, 0.25, ev + ["short run (<1 ns), weak signal"]

    return None


def classify_trajectory_stage(
    traj_path: str | Path,
    *,
    dt_ps: Optional[float] = None,
    duration_ps: Optional[float] = None,
    sibling_dts: Optional[list[float]] = None,
    sibling_durations: Optional[list[float]] = None,
) -> StageResult:
    """Classify one trajectory.  ``sibling_*`` are the values for the other
    trajectories in the same directory (used for relative reasoning)."""
    traj = Path(traj_path)
    stem = traj.stem
    low = stem.lower()
    evidence: list[str] = []

    # ── PBC / fit derivative ────────────────────────────────────────────────
    if any(m in low for m in _DERIVED_MARKERS):
        return StageResult(
            TrajectoryStage.DERIVED, 0.85,
            [f"filename '{traj.name}' contains a PBC/fit derivative marker"],
        )

    votes: dict[str, float] = {}
    mdp_dt_out: Optional[float] = None
    mdp_duration: Optional[float] = None

    # ── filename stem ──────────────────────────────────────────────────────
    part_m = _PART_RE.search(stem)
    part_index = int(part_m.group(1)) if part_m else None
    for rx, stage, w in _STEM_RULES:
        if rx.match(low):
            votes[stage] = votes.get(stage, 0.0) + w
            evidence.append(f"filename stem '{stem}' matches {stage} pattern (w={w})")
            break
    if part_index is not None:
        votes[TrajectoryStage.PRODUCTION] = votes.get(TrajectoryStage.PRODUCTION, 0.0) + 0.35
        evidence.append(f"continuation marker part/seg {part_index} → production segment")

    # ── paired mdp ────────────────────────────────────────────────────────
    mdp_path = _find_mdp(traj)
    if mdp_path:
        mdp = _parse_mdp(mdp_path)
        mdp_dt_out, mdp_duration = mdp_timing(mdp)
        mdp_hit = _mdp_evidence(mdp)
        if mdp_hit:
            stage, w, ev = mdp_hit
            votes[stage] = votes.get(stage, 0.0) + w
            evidence.append(f"paired mdp '{mdp_path.name}': " + "; ".join(ev))
    if duration_ps is None:
        duration_ps = mdp_duration
    if dt_ps is None:
        dt_ps = mdp_dt_out

    # ── relative frame spacing / duration ─────────────────────────────────
    if dt_ps is not None and sibling_dts:
        others = [d for d in sibling_dts if d and d > 0]
        if others and dt_ps >= max(others + [dt_ps]) and dt_ps >= 2.0:
            votes[TrajectoryStage.PRODUCTION] = votes.get(TrajectoryStage.PRODUCTION, 0.0) + 0.25
            evidence.append(
                f"coarsest output spacing in directory ({dt_ps} ps vs {sorted(set(others))}) "
                f"→ production-like")
        elif others and dt_ps <= min(others + [dt_ps]) and dt_ps < 1.0:
            votes[TrajectoryStage.EQUILIBRATION_OTHER] = (
                votes.get(TrajectoryStage.EQUILIBRATION_OTHER, 0.0) + 0.15)
            evidence.append(f"finest output spacing in directory ({dt_ps} ps) → equilibration-like")
    if duration_ps is not None:
        if duration_ps >= 5000:
            votes[TrajectoryStage.PRODUCTION] = votes.get(TrajectoryStage.PRODUCTION, 0.0) + 0.3
            evidence.append(f"trajectory duration {duration_ps:.0f} ps (≥ 5 ns) → production-like")
        elif duration_ps <= 250 and sibling_durations and max(sibling_durations) > duration_ps * 3:
            votes[TrajectoryStage.EQUILIBRATION_OTHER] = (
                votes.get(TrajectoryStage.EQUILIBRATION_OTHER, 0.0) + 0.15)
            evidence.append(
                f"short trajectory ({duration_ps:.0f} ps) vs a much longer sibling "
                f"→ equilibration-like")

    if not votes:
        return StageResult(TrajectoryStage.UNKNOWN, 0.0,
                           evidence or ["no stage evidence available"], part_index,
                           mdp_dt_out_ps=mdp_dt_out, mdp_duration_ps=mdp_duration)

    stage = max(votes, key=votes.get)
    # confidence: winning vote squashed to (0,1), penalised if a rival is close
    top = votes[stage]
    rival = max((v for k, v in votes.items() if k != stage), default=0.0)
    conf = max(0.0, min(0.97, top - 0.5 * rival))
    if rival > 0 and top - rival < 0.2:
        evidence.append(f"NOTE: competing stage evidence {dict(sorted(votes.items()))}")
    return StageResult(stage, round(conf, 3), evidence, part_index,
                       mdp_dt_out_ps=mdp_dt_out, mdp_duration_ps=mdp_duration)
