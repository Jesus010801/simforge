"""Lightweight trajectory inspection — runs *before* any observable.

Uses ``gmx check`` (O(1) memory, streams the trajectory once) to record frame
count, timestep, time range, atom count, box presence and precision.  The
reference-structure atom count is used as the topology atom count for the
compatibility check (cheap and reliable; ``gmx dump`` on a multi-GB ``.tpr`` is
avoided per the performance rules).

Nothing here transforms coordinates.  Original files are only read.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from analysis.campaign.gmx import gmx_available, run_gmx
from analysis.campaign.models import (
    CampaignWarning, Severity, TrajectoryInspection, TrajectorySegment,
)
from analysis.campaign.fingerprint import fingerprint_file

_ATOMS_RE = re.compile(r"#\s*Atoms\s+(\d+)")
_PRECISION_RE = re.compile(r"Precision\s+([\d.eE+-]+)")
_FIRST_FRAME_RE = re.compile(r"Reading frame\s+0\s+time\s+([-\d.eE+]+)")
_LAST_FRAME_RE = re.compile(r"Last frame\s+(\d+)\s+time\s+([-\d.eE+]+)")
_ITEM_RE = re.compile(r"^(Coords|Box|Step|Time|Velocities|Forces)\s+(\d+)\s+([\d.eE+-]*)")


def _structure_atom_count(structure_path: Optional[Path]) -> Optional[int]:
    if not structure_path or not Path(structure_path).is_file():
        return None
    p = Path(structure_path)
    try:
        if p.suffix.lower() == ".gro":
            lines = p.read_text(errors="replace").splitlines()
            return int(lines[1].strip()) if len(lines) >= 2 else None
        if p.suffix.lower() == ".pdb":
            n = 0
            for line in p.read_text(errors="replace").splitlines():
                if line[:6].strip() in ("ATOM", "HETATM"):
                    n += 1
            return n or None
    except (ValueError, OSError):
        return None
    return None


def _parse_gmx_check(text: str) -> dict:
    out: dict = {}
    m = _ATOMS_RE.search(text)
    if m:
        out["atom_count"] = int(m.group(1))
    m = _PRECISION_RE.search(text)
    if m:
        out["precision"] = "single"     # gmx only prints Precision for compressed (xtc)
    m = _FIRST_FRAME_RE.search(text)
    if m:
        out["start_time_ps"] = float(m.group(1))
    m = _LAST_FRAME_RE.search(text)
    if m:
        out["last_frame_index"] = int(m.group(1))
        out["end_time_ps"] = float(m.group(2))
    for line in text.splitlines():
        im = _ITEM_RE.match(line.strip())
        if not im:
            continue
        item, nframes, dt = im.group(1), int(im.group(2)), im.group(3)
        if item == "Coords":
            out["n_frames"] = nframes
            if dt:
                out["dt_ps"] = float(dt)
        elif item == "Box":
            out["box_frames"] = nframes
        elif item == "Time" and dt and "dt_ps" not in out:
            out["dt_ps"] = float(dt)
    return out


def inspect_trajectory(
    *,
    trajectory_paths: list[str | Path],
    topology_path: Optional[str | Path] = None,
    structure_path: Optional[str | Path] = None,
    membrane_present: Optional[bool] = None,
    gmx: str = "gmx",
    fingerprint_segments: bool = True,
) -> TrajectoryInspection:
    traj_paths = [str(Path(p).resolve()) for p in trajectory_paths]
    insp = TrajectoryInspection(
        trajectory_paths=traj_paths,
        topology_path=str(Path(topology_path).resolve()) if topology_path else None,
        structure_path=str(Path(structure_path).resolve()) if structure_path else None,
        tool="gmx check",
        membrane_present=membrane_present,
    )

    ref_for_atoms = structure_path
    insp.atom_count_topology = _structure_atom_count(Path(ref_for_atoms) if ref_for_atoms else None)

    if not gmx_available(gmx):
        insp.inspected = False
        insp.warnings.append(CampaignWarning(
            "gmx_unavailable",
            f"'{gmx}' not found; trajectory could not be inspected. Frame count, "
            f"timestep and box presence are unknown.",
            Severity.WARN,
        ))
        # still record segments (paths + fingerprints) so discovery is reproducible
        _record_segments(insp, traj_paths, fingerprint_segments)
        return insp

    per_segment: list[dict] = []
    combined_out: list[str] = []
    for tp in traj_paths:
        res = run_gmx(["check", "-f", tp], gmx=gmx, timeout=1800)
        combined_out.append(f"$ {' '.join(res.argv)}\n{res.stderr}\n{res.stdout}")
        if not res.available:
            insp.inspected = False
            insp.warnings.append(CampaignWarning(
                "gmx_unavailable", f"gmx not runnable: {res.stderr}", Severity.WARN))
            _record_segments(insp, traj_paths, fingerprint_segments)
            return insp
        if res.returncode != 0:
            insp.warnings.append(CampaignWarning(
                "gmx_check_failed",
                f"`gmx check -f {Path(tp).name}` exited {res.returncode}: "
                f"{res.stderr.strip()[-300:]}",
                Severity.WARN,
            ))
        per_segment.append(_parse_gmx_check(res.stderr + "\n" + res.stdout))

    insp.raw_tool_output = "\n\n".join(combined_out)[-20000:]
    insp.inspected = True

    # ── segments ────────────────────────────────────────────────────────────
    _record_segments(insp, traj_paths, fingerprint_segments, per_segment)

    # ── aggregate ───────────────────────────────────────────────────────────
    atom_counts = {d.get("atom_count") for d in per_segment if d.get("atom_count")}
    if len(atom_counts) == 1:
        insp.atom_count_trajectory = atom_counts.pop()
    elif len(atom_counts) > 1:
        insp.warnings.append(CampaignWarning(
            "segment_atom_count_mismatch",
            f"trajectory segments report different atom counts: {sorted(atom_counts)}",
            Severity.ERROR,
        ))

    frame_totals = [d.get("n_frames") for d in per_segment if d.get("n_frames") is not None]
    if frame_totals:
        insp.n_frames = sum(frame_totals)
    dts = {round(d["dt_ps"], 6) for d in per_segment if d.get("dt_ps")}
    if len(dts) == 1:
        insp.dt_ps = dts.pop()
    elif len(dts) > 1:
        insp.warnings.append(CampaignWarning(
            "inconsistent_timestep",
            f"segments have different timesteps: {sorted(dts)} ps", Severity.WARN))

    # per-segment end time: from the "Last frame" line if present, else
    # start + (n_frames - 1) * dt  (gmx check does not always print "Last frame")
    for d in per_segment:
        if d.get("end_time_ps") is None and d.get("n_frames") and d.get("dt_ps") is not None:
            s0 = d.get("start_time_ps", 0.0) or 0.0
            d["end_time_ps"] = s0 + (d["n_frames"] - 1) * d["dt_ps"]
            d["_end_time_derived"] = True

    starts = [d["start_time_ps"] for d in per_segment if d.get("start_time_ps") is not None]
    ends = [d["end_time_ps"] for d in per_segment if d.get("end_time_ps") is not None]
    if starts:
        insp.start_time_ps = min(starts)
    if ends:
        insp.end_time_ps = max(ends)
    if insp.start_time_ps is not None and insp.end_time_ps is not None:
        insp.total_duration_ps = insp.end_time_ps - insp.start_time_ps

    box_ok = [d for d in per_segment if d.get("box_frames") and d.get("n_frames")]
    if box_ok:
        insp.box_present = all(d["box_frames"] == d["n_frames"] for d in box_ok)
        if not insp.box_present:
            insp.warnings.append(CampaignWarning(
                "missing_box",
                "one or more frames have no box vectors — PBC-aware analyses are unsafe",
                Severity.ERROR,
            ))
    precisions = {d.get("precision") for d in per_segment if d.get("precision")}
    if precisions:
        insp.precision = precisions.pop()

    # ── atom-count compatibility ───────────────────────────────────────────
    if insp.atom_count_trajectory and insp.atom_count_topology:
        insp.atom_count_match = insp.atom_count_trajectory == insp.atom_count_topology
        if not insp.atom_count_match:
            insp.warnings.append(CampaignWarning(
                "atom_count_mismatch",
                f"trajectory has {insp.atom_count_trajectory} atoms but the reference "
                f"structure has {insp.atom_count_topology} — topology/trajectory pairing "
                f"is likely wrong",
                Severity.ERROR,
            ))

    # ── segment ordering / overlap ────────────────────────────────────────
    _assess_segments(insp)

    # ── membrane orientation note (best effort, from the reference gro) ────
    if membrane_present and structure_path:
        insp.membrane_orientation_note = _membrane_orientation_note(structure_path)

    return insp


def _record_segments(
    insp: TrajectoryInspection,
    traj_paths: list[str],
    do_fingerprint: bool,
    per_segment: Optional[list[dict]] = None,
) -> None:
    insp.segments = []
    for i, tp in enumerate(traj_paths):
        seg = TrajectorySegment(path=tp)
        if do_fingerprint and Path(tp).is_file():
            seg.fingerprint = fingerprint_file(Path(tp))
        if per_segment and i < len(per_segment):
            d = per_segment[i]
            seg.n_frames = d.get("n_frames")
            seg.start_time_ps = d.get("start_time_ps")
            seg.end_time_ps = d.get("end_time_ps")
            seg.dt_ps = d.get("dt_ps")
        insp.segments.append(seg)


def _assess_segments(insp: TrajectoryInspection) -> None:
    segs = [s for s in insp.segments if s.start_time_ps is not None and s.end_time_ps is not None]
    if len(insp.segments) <= 1:
        insp.segment_ordering = "single"
        return
    if len(segs) != len(insp.segments):
        insp.segment_ordering = "uncertain"
        insp.warnings.append(CampaignWarning(
            "segment_times_unknown",
            "cannot determine time ranges for every segment — ordering is unverified; "
            "segments will not be concatenated",
            Severity.WARN,
        ))
        return
    ordered = sorted(range(len(segs)), key=lambda i: segs[i].start_time_ps)
    overlap = False
    for a, b in zip(ordered, ordered[1:]):
        if segs[a].end_time_ps > segs[b].start_time_ps + 1e-6:
            overlap = True
    insp.segment_overlap = overlap
    for rank, idx in enumerate(ordered):
        insp.segments[idx].order_index = rank
        insp.segments[idx].order_evidence = ["start_time_ps ascending"]
    if overlap:
        insp.segment_ordering = "uncertain"
        insp.warnings.append(CampaignWarning(
            "segment_overlap",
            "trajectory segments overlap in time — concatenation would double-count frames",
            Severity.ERROR,
        ))
    else:
        insp.segment_ordering = "confident"


def _membrane_orientation_note(structure_path: str | Path) -> str:
    """Best-effort: is the bilayer normal along z? (from lipid headgroup spread)."""
    try:
        from analysis.campaign.structure.spatial import analyse_membrane_embedding
        emb = analyse_membrane_embedding(structure_path)
        if emb.membrane_present and emb.membrane_z_min is not None:
            return (f"bilayer located along z in the reference structure "
                    f"({emb.note}); assumed membrane-normal = z")
        if emb.membrane_present:
            return "lipids present but bilayer slab could not be located along z"
    except Exception:  # noqa: BLE001 - best effort only
        pass
    return ""
