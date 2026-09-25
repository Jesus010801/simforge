"""Level 3 — molecule-level periodic image changes in finite multi-molecule components.

A component made of several topological molecules (an oligomeric receptor, a
receptor + peptide complex) can be *split across periodic images* in the stored
coordinates although it is physically intact: one whole molecule sits one
lattice vector away from the others.  The component's aggregate centre then
moves by only a fraction of a box (1 of 4 chains ⇒ ~¼ box) — which the
component-level jump detector cannot see.  This detector looks at every
molecule's own whole-molecule centre.

Observation vs interpretation: a displacement that equals a lattice translation
``n1·a + n2·b + n3·c`` up to a small residual is reported as a *periodic image
change of an intact molecule* — never as dissociation.  Periodically extended
components (membrane, solvent, ions) are not applicable: individual molecules
changing image is expected there and needs domain-specific geometry.
"""
from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

from analysis.campaign.diagnostics.context import DiagnosticContext, ProbeError, box_matrix, box_valid
from analysis.campaign.diagnostics.registry import detector
from analysis.campaign.models import Diagnostic, EvidenceConfidence, SamplingStrategy, Severity

FINITE_ROLES = ("receptor", "peptide", "protein_partner", "ligand", "cofactor",
                "nucleic_acid", "complex")
EXTENDED_ROLES = ("membrane", "water", "ions")
MAX_MOLECULES = 64                    # a finite assembly, not a solvent / membrane
PROBE_SCHEMA = "simforge/diagnostic-probe/molecule-centres/v1"


def lattice_decomposition(d: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(n, lattice_vector, residual)`` with ``d = n @ B + residual`` and integer
    ``n`` the nearest lattice point in fractional coordinates (rows of ``B`` are
    the cell vectors a, b, c — valid for triclinic cells)."""
    f = d @ np.linalg.inv(B)
    n = np.round(f)
    lat = n @ B
    return n.astype(int), lat, d - lat


def image_changes(pos: np.ndarray, boxes, residual_nm: float) -> list[dict]:
    """Consecutive-frame displacements that are a non-zero lattice translation."""
    out = []
    for i in range(len(pos) - 1):
        B = boxes[i + 1] if boxes is not None else None
        if B is None or not box_valid(B):
            continue
        d = pos[i + 1] - pos[i]
        n, lat, res = lattice_decomposition(d, B)
        if n.any() and float(np.linalg.norm(res)) <= residual_nm:
            out.append({"frames": [i, i + 1], "raw_displacement_nm": np.round(d, 4).tolist(),
                        "lattice_translation": n.tolist(),
                        "lattice_vector_nm": np.round(lat, 4).tolist(),
                        "minimum_image_displacement_nm": np.round(res, 4).tolist(),
                        "box_nm": np.round(B, 4).tolist()})
    return out


def split_frames(positions: list[np.ndarray], boxes) -> list[int]:
    """Frames where some molecule's placement relative to molecule 0 is a non-zero
    lattice translation away from its minimum-image placement (assembly split)."""
    n_frames = len(positions[0])
    out = []
    for t in range(n_frames):
        B = boxes[t] if boxes is not None else None
        if B is None or not box_valid(B):
            continue
        ref = positions[0][t]
        if any(lattice_decomposition(p[t] - ref, B)[0].any() for p in positions[1:]):
            out.append(t)
    return out


def _finite_multi(ctx: DiagnosticContext) -> dict:
    """role -> [(start, end, moltype)] for finite roles spanning > 1 molecule."""
    if not ctx.structure_path or not ctx.structure_path.lower().endswith(".tpr"):
        return {}
    from analysis.campaign.structure.index_groups import resolve_index_group
    from analysis.campaign.structure.molecules import molecule_partition
    part = molecule_partition(ctx.structure_path, ctx.gmx)
    if part is None or not ctx.index_path:
        return {}
    out = {}
    for role in FINITE_ROLES:
        name = ctx.groups.get(role)
        grp = resolve_index_group(ctx.index_path, name) if name else None
        if grp is None:
            continue
        sel = part.select(grp.atom_ids)
        if 1 < sel.n_molecules <= MAX_MOLECULES:
            out[role] = [(part.starts[m], part.starts[m + 1] - 1, part.moltypes[m])
                         for m in sel.molecules]
    return out


def _applicable(ctx: DiagnosticContext) -> tuple[bool, str]:
    if not ctx.structure_path or not ctx.structure_path.lower().endswith(".tpr"):
        return False, "needs the .tpr molecule topology"
    ti = ctx.get_time_index()
    if ti.n_frames < 2:
        return False, "fewer than two frames"
    if not _finite_multi(ctx):
        return False, ("no finite component spans more than one topological molecule "
                       "(membrane / solvent / ions are periodically extended: not applicable)")
    return True, ""


def _probe(ctx: DiagnosticContext, mols: list[tuple[int, int, str]]) -> np.ndarray:
    """(n_frames, n_molecules, 3) whole-molecule centres (COM, .tpr masses)."""
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.gmx import gmx_version, run_gmx
    fp = ctx.trajectory_fingerprint()
    key = hashlib.sha256(json.dumps({
        "schema": PROBE_SCHEMA, "trajectory": fp.digest if fp else None,
        "tpr": fingerprint_file(Path(ctx.structure_path)).digest, "molecules": mols,
        "gmx": gmx_version(ctx.gmx)}, sort_keys=True).encode()).hexdigest()[:24]
    cache = Path(ctx.cache_dir) / "molecule_centres" / f"{key}.npy" if ctx.cache_dir else None
    n_frames = ctx.get_time_index().n_frames
    if cache is not None and cache.is_file():
        arr = np.load(cache)
        if arr.shape == (n_frames, len(mols), 3):
            ctx.probe_records.append({"probe": "molecule_centres", "cache_key": key, "from_cache": True})
            return arr
    sel = [f"com of atomnr {a} to {b}" for a, b, _ in mols]
    with tempfile.TemporaryDirectory(prefix="simforge_molc_") as tmp:
        out = Path(tmp) / "centres.xvg"
        res = run_gmx(["trajectory", "-f", str(Path(ctx.trajectory_path).resolve()),
                       "-s", ctx.structure_path, "-rmpbc", "-nopbc", "-ox", str(out),
                       "-select", *sel], gmx=ctx.gmx, timeout=4 * 3600)
        if not res.ok or not out.is_file():
            raise ProbeError(f"gmx trajectory failed (rc={res.returncode}): {res.stderr.strip()[-300:]}")
        rows = [[float(v) for v in l.split()] for l in out.read_text().splitlines()
                if l.strip() and l.lstrip()[0] not in "#@"]
    if len(rows) != n_frames or any(len(r) != 1 + 3 * len(mols) for r in rows):
        raise ProbeError("molecule-centre probe rows/columns do not match the time index")
    arr = np.asarray(rows)[:, 1:].reshape(n_frames, len(mols), 3)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache, arr)
    ctx.probe_records.append({"probe": "molecule_centres", "cache_key": key, "from_cache": False,
                              "command": list(res.argv)})
    return arr


@detector("molecule_periodic_image_change", "1", 3,
          "whole molecules of a finite multi-molecule component changing periodic image "
          "(assembly split across images in the stored coordinates)",
          applicable=_applicable, roles=lambda ctx: sorted(_finite_multi(ctx)),
          sampling=SamplingStrategy.ALL_FRAMES)
def molecule_periodic_image_change(ctx: DiagnosticContext) -> list[Diagnostic]:
    ti = ctx.get_time_index()
    raw = ctx.boxes()
    boxes = [box_matrix(b) for b in raw] if raw is not None and len(raw) == ti.n_frames else None
    p = ctx.params
    out = []
    for role, mols in _finite_multi(ctx).items():
        arr = _probe(ctx, mols)
        per_mol = []
        events = []
        for k, (a, b, mt) in enumerate(mols):
            ev = image_changes(arr[:, k, :], boxes, p.wrap_residual_nm)
            for e in ev:
                e.update({"molecule": k + 1, "moltype": mt, "atoms": [a, b],
                          "times_ps": [ti.times_ps[e["frames"][0]], ti.times_ps[e["frames"][1]]]})
            per_mol.append(len(ev))
            events += ev
        split = split_frames([arr[:, k, :] for k in range(len(mols))], boxes)
        if not events and not split:
            continue
        events.sort(key=lambda e: e["frames"][0])
        first = events[0] if events else None
        out.append(Diagnostic(
            code="molecule_periodic_image_change", severity=Severity.WARN, scope=f"component:{role}",
            message=(f"{role}: {len(mols)} topological molecules; {len(events)} whole-molecule "
                     f"periodic image change(s)"
                     + (f" (first: molecule {first['molecule']} {first['moltype']} at "
                        f"{first['times_ps'][1]} ps, raw Δ {first['raw_displacement_nm']} nm = "
                        f"lattice {first['lattice_translation']} + residual "
                        f"{first['minimum_image_displacement_nm']} nm)" if first else "")
                     + f"; the stored coordinates split the assembly across periodic images in "
                       f"{len(split)} of {ti.n_frames} frames"),
            frames=[split[0], split[-1]] if split else first["frames"],
            times_ps=([ti.times_ps[split[0]], ti.times_ps[split[-1]]] if split
                      else first["times_ps"]),
            measured={"n_molecules": len(mols), "image_changes_per_molecule": per_mol,
                      "n_frames_split": len(split),
                      "events": events[:p.max_events_reported],
                      "events_truncated": len(events) > p.max_events_reported},
            interpretation="periodic_image_change_of_intact_molecules",
            confidence=EvidenceConfidence.STRONG,
            uncertainty=["a molecule's whole-molecule centre moving by a lattice vector is an image "
                         "choice, not dissociation; physical separation would not be a lattice "
                         "translation"],
            candidate_remediations=["cluster_assembly"],
            provenance={"molecules": [{"molecule": k + 1, "moltype": mt, "atoms": [a, b]}
                                      for k, (a, b, mt) in enumerate(mols)],
                        "centre": "centre of mass of each whole molecule (.tpr masses, -rmpbc)",
                        "lattice": "nearest lattice point in fractional (triclinic) coordinates"},
            sampling={"strategy": SamplingStrategy.ALL_FRAMES, "frames_examined": ti.n_frames,
                      "total_frames": ti.n_frames}))
    return out
