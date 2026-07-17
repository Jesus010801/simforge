"""
validators/protein_membrane_interface.py
Phase 9A: Protein–Membrane Interface Evaluator.

Evaluates whether the bilayer actually covers the lipid-facing TM protein
surface.  Does NOT modify coordinates, does NOT move lipids.

Algorithm
---------
1. Parse GRO: TM heavy atoms, all-protein heavy atoms, lipid heavy atoms.
2. Identify TM surface atoms: TM heavy atoms with few protein neighbours
   within surface_neighbor_radius_nm (low local density → exposed surface).
3. For each TM surface atom compute:
     nearest_lipid_dist       (any lipid heavy atom)
     nearest_headgroup_dist   (P / N atoms)
     nearest_tail_dist        (lipid C atoms)
4. Classify coverage per surface atom:
     covered        <= covered_threshold_nm
     weakly_covered covered_threshold_nm < d <= weakly_covered_threshold_nm
     exposed_gap    > weakly_covered_threshold_nm
5. Find contiguous gap clusters (union-find on exposed atoms).
6. Compute per-residue and per-TM-segment metrics.
7. Write protein_membrane_interface_report.json.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    import numpy as np
    _HAS_NP = True
except ImportError:
    _HAS_NP = False


# ── Constants ─────────────────────────────────────────────────────────────────

DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset({
    "DPP", "DPPC", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER",
})

# Solvent / ion resnames to ignore
_SOLVENT_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
    "NA", "CL", "SOD", "MG", "K", "CA",
})

# Headgroup atom names (phosphorus, nitrogen are canonical headgroup markers)
_HEADGROUP_ATOMNAMES: frozenset[str] = frozenset({
    "P", "P1",                             # phosphorus
    "N", "NZ", "NM1", "NM2", "NM3",       # choline / ethanolamine nitrogen
    "O11", "O12", "O13", "O14",            # phosphate oxygens
    "O21", "O22",                          # glycerol-ester oxygens
})


# ── GRO parser ────────────────────────────────────────────────────────────────

@dataclass
class _Atom:
    resid:    int
    resname:  str
    atomname: str
    x: float
    y: float
    z: float


def _parse_gro(gro_path: Path) -> tuple[list[_Atom], tuple[float, float, float]]:
    lines = gro_path.read_text().splitlines()
    n = int(lines[1].strip())
    atoms: list[_Atom] = []
    for ln in lines[2 : 2 + n]:
        if len(ln) < 44:
            continue
        try:
            resid    = int(ln[0:5])
            resname  = ln[5:10].strip()
            atomname = ln[10:15].strip()
            # Skip hydrogen by first character or single-H names
            if atomname.upper().startswith("H"):
                continue
            x = float(ln[20:28])
            y = float(ln[28:36])
            z = float(ln[36:44])
            atoms.append(_Atom(resid, resname, atomname, x, y, z))
        except ValueError:
            pass
    box_parts = lines[2 + n].split()
    box = (float(box_parts[0]), float(box_parts[1]), float(box_parts[2]))
    return atoms, box


# ── Distance helpers ──────────────────────────────────────────────────────────

def _coords(atoms: list[_Atom]):
    """Return (N, 3) float array."""
    if not _HAS_NP:
        raise RuntimeError("numpy is required for protein_membrane_interface")
    return np.array([[a.x, a.y, a.z] for a in atoms], dtype=np.float32)


def _batch_min_distances(
    query:    "np.ndarray",   # (Q, 3)
    database: "np.ndarray",   # (D, 3)
    batch:    int = 256,
) -> "np.ndarray":           # (Q,) minimum distance from each query to any DB point
    """Compute min distance from each query point to any database point."""
    import numpy as np
    Q = len(query)
    min_dists = np.full(Q, np.inf, dtype=np.float32)
    for i in range(0, Q, batch):
        q_batch = query[i : i + batch]               # (B, 3)
        diff = q_batch[:, None, :] - database[None, :, :]  # (B, D, 3)
        dists = np.sqrt((diff ** 2).sum(axis=2))     # (B, D)
        min_dists[i : i + batch] = dists.min(axis=1)
    return min_dists


def _count_neighbors(
    query:    "np.ndarray",   # (Q, 3)
    database: "np.ndarray",   # (D, 3)
    radius:   float,
    batch:    int = 256,
) -> "np.ndarray":           # (Q,) count
    import numpy as np
    Q = len(query)
    counts = np.zeros(Q, dtype=np.int32)
    for i in range(0, Q, batch):
        q_batch = query[i : i + batch]
        diff = q_batch[:, None, :] - database[None, :, :]
        dists = np.sqrt((diff ** 2).sum(axis=2))
        counts[i : i + batch] = (dists < radius).sum(axis=1)
    return counts


# ── TM segment detection ──────────────────────────────────────────────────────

def _find_tm_segments(
    tm_residues: set[int],
    gap_threshold: int = 5,
) -> list[tuple[int, int]]:
    """Split TM residue set into contiguous segments."""
    if not tm_residues:
        return []
    sorted_res = sorted(tm_residues)
    segments: list[tuple[int, int]] = []
    start = prev = sorted_res[0]
    for r in sorted_res[1:]:
        if r - prev > gap_threshold:
            segments.append((start, prev))
            start = r
        prev = r
    segments.append((start, sorted_res[-1]))
    return segments


# ── Union-Find for gap clusters ───────────────────────────────────────────────

def _find_gap_clusters(
    exposed_indices: list[int],
    coords: "np.ndarray",        # coords of ALL surface atoms (indexed by exposed_indices)
    cluster_distance: float,
    batch: int = 256,
) -> list[list[int]]:
    """Return list of clusters; each cluster is a list of surface-atom indices."""
    import numpy as np

    if len(exposed_indices) == 0:
        return []

    ex_coords = coords[exposed_indices]   # (E, 3)
    E = len(exposed_indices)
    parent = list(range(E))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    # Connect pairs closer than cluster_distance
    for i in range(0, E, batch):
        batch_c = ex_coords[i : i + batch]   # (B, 3)
        diff = batch_c[:, None, :] - ex_coords[None, :, :]
        dists = np.sqrt((diff ** 2).sum(axis=2))  # (B, E)
        rows, cols = np.where((dists < cluster_distance) & (dists > 0))
        for r, c in zip(rows, cols):
            _union(i + int(r), int(c))

    # Collect components
    from collections import defaultdict
    comp: dict[int, list[int]] = defaultdict(list)
    for idx in range(E):
        comp[_find(idx)].append(exposed_indices[idx])
    return list(comp.values())


# ── Main evaluator ────────────────────────────────────────────────────────────

def evaluate_protein_membrane_interface(
    gro_path:                       str | Path,
    tm_residues:                    set[int],
    lipid_resnames:                 Optional[set[str]] = None,
    ec_residues:                    Optional[set[int]] = None,
    ic_residues:                    Optional[set[int]] = None,
    tm_segments:                    Optional[list[tuple[int, int]]] = None,
    covered_threshold_nm:           float = 0.45,
    weakly_covered_threshold_nm:    float = 0.70,
    surface_neighbor_radius_nm:     float = 0.50,
    surface_neighbor_threshold:     int   = 18,
    cluster_distance_nm:            float = 0.50,
    min_fraction_covered:           float = 0.75,
    max_fraction_exposed_gap:       float = 0.15,
    max_p90_distance_nm:            float = 0.70,
    output_dir:                     Optional[str | Path] = None,
) -> dict:
    """
    Evaluate the protein–membrane interface quality for a TM protein in a
    lipid bilayer.  Writes protein_membrane_interface_report.json.

    Returns the report dict.
    """
    import numpy as np

    gro_path = Path(gro_path)
    if not gro_path.exists():
        raise FileNotFoundError(f"GRO not found: {gro_path}")

    if output_dir is None:
        output_dir = gro_path.parent
    output_dir = Path(output_dir)

    lipid_set   = frozenset(lipid_resnames or DEFAULT_LIPID_RESNAMES)
    ec_set      = ec_residues or set()
    ic_set      = ic_residues or set()
    warnings: list[str] = []

    # ── Parse atoms ───────────────────────────────────────────────────────────
    all_atoms, box = _parse_gro(gro_path)

    tm_atoms    = [a for a in all_atoms if a.resid in tm_residues
                   and a.resname not in lipid_set and a.resname not in _SOLVENT_RESNAMES]
    prot_atoms  = [a for a in all_atoms
                   if a.resname not in lipid_set and a.resname not in _SOLVENT_RESNAMES]
    lipid_atoms = [a for a in all_atoms if a.resname in lipid_set]
    hg_atoms    = [a for a in lipid_atoms if a.atomname in _HEADGROUP_ATOMNAMES]
    tail_atoms  = [a for a in lipid_atoms
                   if a.atomname.startswith("C") and a.atomname not in _HEADGROUP_ATOMNAMES]

    if not tm_atoms:
        warnings.append("No TM heavy atoms found — interface evaluation skipped")
        report = _empty_report(str(gro_path), tm_residues, warnings)
        _write_report(report, output_dir)
        return report

    if not lipid_atoms:
        warnings.append("No lipid heavy atoms found — interface evaluation skipped")
        report = _empty_report(str(gro_path), tm_residues, warnings)
        _write_report(report, output_dir)
        return report

    tm_coords   = _coords(tm_atoms)     # (T, 3)
    prot_coords = _coords(prot_atoms)   # (P, 3)
    lip_coords  = _coords(lipid_atoms)  # (L, 3)

    # ── Identify TM surface atoms ─────────────────────────────────────────────
    # Count protein-heavy-atom neighbours within surface_neighbor_radius_nm
    neighbor_counts = _count_neighbors(
        tm_coords, prot_coords,
        radius=surface_neighbor_radius_nm,
    )
    # Subtract self (each TM atom counts itself against prot_atoms if prot contains TM)
    # prot_atoms includes TM atoms, so each TM atom is counted once as its own neighbour
    neighbor_counts = np.maximum(0, neighbor_counts - 1)

    surface_mask  = neighbor_counts < surface_neighbor_threshold  # (T,) bool
    surface_idxs  = np.where(surface_mask)[0]

    if len(surface_idxs) == 0:
        warnings.append(
            f"No TM surface atoms found (threshold={surface_neighbor_threshold}); "
            "all TM atoms appear buried — check threshold or TM annotation"
        )
        report = _empty_report(str(gro_path), tm_residues, warnings)
        _write_report(report, output_dir)
        return report

    surface_tm_atoms = [tm_atoms[i] for i in surface_idxs]
    surf_coords      = tm_coords[surface_idxs]   # (S, 3)

    # ── Lipid distance computation ────────────────────────────────────────────
    nearest_lip  = _batch_min_distances(surf_coords, lip_coords)

    if hg_atoms:
        hg_coords   = _coords(hg_atoms)
        nearest_hg  = _batch_min_distances(surf_coords, hg_coords)
    else:
        nearest_hg  = nearest_lip.copy()
        warnings.append("No headgroup (P/N) atoms found; using all-lipid for headgroup distance")

    if tail_atoms:
        tail_coords = _coords(tail_atoms)
        nearest_tail = _batch_min_distances(surf_coords, tail_coords)
    else:
        nearest_tail = nearest_lip.copy()
        warnings.append("No tail (C-chain) atoms found; using all-lipid for tail distance")

    # ── Coverage classification ───────────────────────────────────────────────
    covered_mask  = nearest_lip <= covered_threshold_nm
    weakly_mask   = (nearest_lip > covered_threshold_nm) & (nearest_lip <= weakly_covered_threshold_nm)
    exposed_mask  = nearest_lip > weakly_covered_threshold_nm

    n_surf        = len(surface_idxs)
    n_covered     = int(covered_mask.sum())
    n_weakly      = int(weakly_mask.sum())
    n_exposed     = int(exposed_mask.sum())

    frac_covered      = n_covered  / n_surf if n_surf else 0.0
    frac_weakly       = n_weakly   / n_surf if n_surf else 0.0
    frac_exposed      = n_exposed  / n_surf if n_surf else 0.0
    mean_dist         = float(nearest_lip.mean())
    p90_dist          = float(np.percentile(nearest_lip, 90))
    max_dist          = float(nearest_lip.max())

    # ── Gap clusters ──────────────────────────────────────────────────────────
    exposed_surf_idxs = surface_idxs[exposed_mask].tolist()
    gap_clusters_raw  = _find_gap_clusters(
        exposed_surf_idxs, tm_coords, cluster_distance_nm
    )
    gap_clusters_raw.sort(key=len, reverse=True)
    largest_gap = len(gap_clusters_raw[0]) if gap_clusters_raw else 0

    gap_clusters_report = []
    for cl in gap_clusters_raw:
        cl_atoms = [tm_atoms[i] for i in cl]
        cl_resids = sorted({a.resid for a in cl_atoms})
        cl_coords_np = tm_coords[cl]
        gap_clusters_report.append({
            "size":       len(cl),
            "residues":   cl_resids,
            "centroid_nm": [
                round(float(cl_coords_np[:, 0].mean()), 3),
                round(float(cl_coords_np[:, 1].mean()), 3),
                round(float(cl_coords_np[:, 2].mean()), 3),
            ],
            "max_dist_nm": round(float(nearest_lip[exposed_mask][
                [i for i, idx in enumerate(surface_idxs[exposed_mask]) if idx in cl]
            ].max() if len(cl) > 0 else 0), 3),
        })

    # ── Per-residue metrics ────────────────────────────────────────────────────
    from collections import defaultdict
    resid_surf_count: dict[int, int]   = defaultdict(int)
    resid_cov_count:  dict[int, int]   = defaultdict(int)
    resid_dists:      dict[int, list]  = defaultdict(list)

    for i, (atom, d, cov) in enumerate(zip(
        surface_tm_atoms, nearest_lip, covered_mask
    )):
        resid_surf_count[atom.resid] += 1
        resid_cov_count[atom.resid]  += int(cov)
        resid_dists[atom.resid].append(float(d))

    coverage_by_residue = {
        r: round(resid_cov_count[r] / resid_surf_count[r], 3)
        for r in resid_surf_count
    }
    mean_dist_by_residue = {
        r: round(float(np.mean(resid_dists[r])), 3)
        for r in resid_dists
    }

    # Worst residues (lowest coverage fraction among those with >= 2 surface atoms)
    worst_residues = sorted(
        [r for r in coverage_by_residue if resid_surf_count[r] >= 2],
        key=lambda r: coverage_by_residue[r],
    )[:10]

    # ── Per-TM-segment metrics ────────────────────────────────────────────────
    if tm_segments is None:
        tm_segments = _find_tm_segments(tm_residues)

    seg_metrics = []
    for s_idx, (seg_start, seg_end) in enumerate(tm_segments):
        seg_res = set(range(seg_start, seg_end + 1)) & tm_residues
        seg_surf = sum(resid_surf_count[r] for r in seg_res if r in resid_surf_count)
        seg_cov  = sum(resid_cov_count[r]  for r in seg_res if r in resid_cov_count)
        all_d    = [d for r in seg_res for d in resid_dists.get(r, [])]
        seg_metrics.append({
            "segment_index":    s_idx,
            "residue_range":    [seg_start, seg_end],
            "n_surface_atoms":  seg_surf,
            "coverage_fraction": round(seg_cov / seg_surf, 3) if seg_surf else None,
            "mean_distance_nm": round(float(np.mean(all_d)), 3) if all_d else None,
        })

    worst_seg = None
    scored_segs = [s for s in seg_metrics if s["coverage_fraction"] is not None]
    if scored_segs:
        worst_seg = min(scored_segs, key=lambda s: s["coverage_fraction"])

    # ── Quality gates ─────────────────────────────────────────────────────────
    quality_passed = True
    if frac_covered < min_fraction_covered:
        warnings.append(
            f"TM surface under-covered: fraction_covered={frac_covered:.3f} "
            f"< min_fraction_covered={min_fraction_covered}"
        )
        quality_passed = False
    if frac_exposed > max_fraction_exposed_gap:
        warnings.append(
            f"TM surface gap too large: fraction_exposed_gap={frac_exposed:.3f} "
            f"> max_fraction_exposed_gap={max_fraction_exposed_gap}"
        )
        quality_passed = False
    if p90_dist > max_p90_distance_nm:
        warnings.append(
            f"p90_nearest_lipid_distance={p90_dist:.3f} nm "
            f"> max_p90_distance_nm={max_p90_distance_nm}"
        )
        quality_passed = False

    # ── TM residue ranges for report ──────────────────────────────────────────
    tm_seg_list = [[s, e] for s, e in _find_tm_segments(tm_residues)]

    report = {
        "source_file":                   str(gro_path),
        "tm_residue_ranges":             tm_seg_list,
        "n_tm_surface_atoms":            n_surf,
        "fraction_covered":              round(frac_covered, 4),
        "fraction_weakly_covered":       round(frac_weakly, 4),
        "fraction_exposed_gap":          round(frac_exposed, 4),
        "mean_nearest_lipid_distance":   round(mean_dist, 4),
        "p90_nearest_lipid_distance":    round(p90_dist, 4),
        "max_nearest_lipid_distance":    round(max_dist, 4),
        "mean_nearest_headgroup_distance": round(float(nearest_hg.mean()), 4),
        "mean_nearest_tail_distance":    round(float(nearest_tail.mean()), 4),
        "n_covered":                     n_covered,
        "n_weakly_covered":              n_weakly,
        "n_exposed_gap":                 n_exposed,
        "largest_contiguous_gap_cluster_size": largest_gap,
        "gap_clusters":                  gap_clusters_report,
        "coverage_fraction_by_residue":  {str(k): v for k, v in coverage_by_residue.items()},
        "mean_distance_by_residue":      {str(k): v for k, v in mean_dist_by_residue.items()},
        "coverage_fraction_by_tm_segment": seg_metrics,
        "worst_residues":                worst_residues,
        "worst_tm_segment":              worst_seg,
        "worst_tm_segments":             sorted(scored_segs, key=lambda s: s["coverage_fraction"] or 1.0)[:3],
        "quality_passed":                quality_passed,
        "quality_thresholds": {
            "min_fraction_covered":             min_fraction_covered,
            "max_fraction_exposed_gap":         max_fraction_exposed_gap,
            "max_p90_nearest_lipid_distance_nm": max_p90_distance_nm,
        },
        "warnings": warnings,
    }

    _write_report(report, output_dir)
    return report


# ── Helpers ───────────────────────────────────────────────────────────────────

def _empty_report(source_file: str, tm_residues: set[int], warnings: list[str]) -> dict:
    tm_seg_list = [[s, e] for s, e in _find_tm_segments(tm_residues)]
    return {
        "source_file":                     source_file,
        "tm_residue_ranges":               tm_seg_list,
        "n_tm_surface_atoms":              0,
        "fraction_covered":                None,
        "fraction_weakly_covered":         None,
        "fraction_exposed_gap":            None,
        "mean_nearest_lipid_distance":     None,
        "p90_nearest_lipid_distance":      None,
        "max_nearest_lipid_distance":      None,
        "mean_nearest_headgroup_distance": None,
        "mean_nearest_tail_distance":      None,
        "n_covered":                       0,
        "n_weakly_covered":                0,
        "n_exposed_gap":                   0,
        "largest_contiguous_gap_cluster_size": 0,
        "gap_clusters":                    [],
        "coverage_fraction_by_residue":    {},
        "mean_distance_by_residue":        {},
        "coverage_fraction_by_tm_segment": [],
        "worst_residues":                  [],
        "worst_tm_segment":                None,
        "worst_tm_segments":               [],
        "quality_passed":                  False,
        "quality_thresholds":              {},
        "warnings":                        warnings,
    }


def _write_report(report: dict, output_dir: Path) -> None:
    (output_dir / "protein_membrane_interface_report.json").write_text(
        json.dumps(report, indent=2)
    )
