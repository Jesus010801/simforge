"""
validators/embedding_quality.py

Local TM-aware membrane embedding quality diagnostics.

Metrics computed over the local TM analysis region (NOT the full simulation box):
  - Local lipid count and APL within the expanded TM footprint or radial shell
  - Void fraction: deficit of lipids relative to the reference APL
  - Robust bilayer core boundaries using median and percentile headgroup Z
  - Bilayer thickness (core = inner headgroup boundaries; outer = extremes)

All metrics are diagnostic. This module never blocks or remediates.

Analysis region modes
---------------------
"tm_footprint" (default):
    Circular region of radius (tm_footprint_radius + footprint_margin_nm).
    Includes lipids directly adjacent to and surrounding the TM bundle.

"radial_shell":
    Annular ring from tm_footprint_radius to (tm_footprint_radius + shell_width_nm).
    Excludes the protein core — measures only the immediately adjacent lipid belt.

Void fraction definition
------------------------
    n_expected  = local_region_area / reference_apl_nm2
    void_fraction = max(0, 1 - n_found / n_expected)

A void_fraction of 0 means the region is packed at the reference density.
A positive value means fewer lipids than expected (gap / cavity effect).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

REGION_TM_FOOTPRINT = "tm_footprint"
REGION_RADIAL_SHELL = "radial_shell"

DEFAULT_FOOTPRINT_MARGIN_NM = 0.5
DEFAULT_SHELL_WIDTH_NM      = 1.0


class EmbeddingQualityResult(BaseModel):
    """
    Post-embedding quality metrics within the local TM analysis region.
    All values are diagnostic; none block the workflow.
    """
    source_file:                str = ""
    analysis_region:            str = REGION_TM_FOOTPRINT

    # TM bundle geometry (nm)
    tm_center_x:                Optional[float] = None
    tm_center_y:                Optional[float] = None
    tm_footprint_radius_nm:     Optional[float] = None   # max CA-to-center XY distance
    local_region_area_nm2:      Optional[float] = None   # XY area of analysis region
    footprint_margin_nm:        float = DEFAULT_FOOTPRINT_MARGIN_NM
    shell_width_nm:             float = DEFAULT_SHELL_WIDTH_NM

    # Local lipid density and void metrics
    n_lipids_total:             int = 0
    n_lipids_in_region:         int = 0
    apl_local_nm2:              Optional[float] = None   # local_region_area / n_lipids_in_region
    reference_apl_nm2:          Optional[float] = None
    void_fraction_local:        Optional[float] = None   # 1 − found/expected; 0 = perfect packing
    void_area_local_nm2:        Optional[float] = None   # void_fraction × local_region_area

    # Bilayer Z boundaries — robust statistics (nm)
    midplane_z_median:          Optional[float] = None   # median of tail_atom Z
    midplane_z_mean:            Optional[float] = None
    z_core_top_median:          Optional[float] = None   # median of top-leaflet headgroup Z
    z_core_bot_median:          Optional[float] = None   # median of bot-leaflet headgroup Z
    z_core_top_mean:            Optional[float] = None
    z_core_bot_mean:            Optional[float] = None
    z_core_top_p25:             Optional[float] = None   # 25th pct of top HG Z (inner boundary)
    z_core_bot_p75:             Optional[float] = None   # 75th pct of bot HG Z (inner boundary)
    z_outer_top:                Optional[float] = None   # max top-leaflet HG Z (outermost surface)
    z_outer_bot:                Optional[float] = None   # min bot-leaflet HG Z (outermost surface)
    headgroup_z_iqr_top:        Optional[float] = None   # IQR of top-leaflet HG Z (leaflet disorder)
    headgroup_z_iqr_bot:        Optional[float] = None
    bilayer_core_thickness_nm:  Optional[float] = None   # z_core_top_median − z_core_bot_median
    bilayer_outer_thickness_nm: Optional[float] = None   # z_outer_top − z_outer_bot

    # Custom TM-aware backend metrics
    tm_burial_score:                  Optional[float] = None
    soluble_domain_core_penetration:  Optional[int] = None

    # Counts
    n_headgroup_atoms_found:    int = 0
    n_tail_atoms_found:         int = 0
    n_tm_ca_atoms:              int = 0

    message:                    str = ""
    warnings:                   list[str] = Field(default_factory=list)


def diagnose_embedding_quality(
    gro_path:            Path | str,
    lipid_resname:       str,
    headgroup_atom:      str,
    tail_atom:           str,
    tm_residues:         Optional[set[int]] = None,
    footprint_margin_nm: float = DEFAULT_FOOTPRINT_MARGIN_NM,
    shell_width_nm:      float = DEFAULT_SHELL_WIDTH_NM,
    reference_apl_nm2:   Optional[float] = None,
    analysis_region:     str = REGION_TM_FOOTPRINT,
) -> EmbeddingQualityResult:
    """
    Compute local TM-aware membrane embedding quality metrics.

    The analysis region is restricted to the expanded TM footprint or a radial
    shell around the TM bundle — never the whole simulation box.

    Args:
        gro_path:            Path to embedded system .gro (post shrink-loop convergence).
        lipid_resname:       GRO residue name of lipid (e.g. "DPP" for DPPC OPLS-AA).
        headgroup_atom:      Atom name marking leaflet headgroups (e.g. "O33").
        tail_atom:           Atom name marking bilayer midplane (e.g. "C50").
        tm_residues:         Set of protein residue numbers in TM bundle. None = all protein.
        footprint_margin_nm: Expand TM bounding radius by this amount (tm_footprint mode).
        shell_width_nm:      Width of the annular shell (radial_shell mode only).
        reference_apl_nm2:   Expected APL from literature/force field in nm². Required for
                             void_fraction computation; if None, void metrics are omitted.
        analysis_region:     "tm_footprint" (default) or "radial_shell".
    """
    path = Path(gro_path)
    if not path.exists():
        return EmbeddingQualityResult(
            source_file=str(path),
            message=f"GRO file not found: {path}",
        )

    solvent_names = {"SOL", "HOH", "NA", "CL", "K", "MG", "CA"}

    try:
        lines = path.read_text().splitlines()
    except OSError as exc:
        return EmbeddingQualityResult(
            source_file=str(path),
            message=f"GRO read error: {exc}",
        )

    atom_lines = lines[2:-1]

    tm_ca_coords: list[tuple[float, float, float]] = []
    soluble_coords: list[tuple[float, float, float]] = []
    headgroup_z:  list[float] = []
    tail_z:       list[float] = []
    lipid_coms:   dict[int, list[tuple[float, float, float]]] = {}

    for line in atom_lines:
        if len(line) < 44:
            continue
        resname  = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            resnum = int(line[0:5])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        if resname == lipid_resname:
            lipid_coms.setdefault(resnum, []).append((x, y, z))
            if atomname == headgroup_atom:
                headgroup_z.append(z)
            elif atomname == tail_atom:
                tail_z.append(z)
        elif resname not in solvent_names:
            is_tm = (tm_residues is None or resnum in tm_residues)
            if atomname == "CA" and is_tm:
                tm_ca_coords.append((x, y, z))
            elif tm_residues is not None and not is_tm:
                soluble_coords.append((x, y, z))

    n_lipids_total = len(lipid_coms)
    n_tm_ca        = len(tm_ca_coords)

    if not tm_ca_coords:
        return EmbeddingQualityResult(
            source_file=str(path),
            n_lipids_total=n_lipids_total,
            n_headgroup_atoms_found=len(headgroup_z),
            n_tail_atoms_found=len(tail_z),
            message="No TM CA atoms found — quality diagnostics skipped.",
        )

    # TM bundle XY centroid and footprint radius
    tm_cx = sum(c[0] for c in tm_ca_coords) / n_tm_ca
    tm_cy = sum(c[1] for c in tm_ca_coords) / n_tm_ca
    tm_r  = max(math.sqrt((c[0]-tm_cx)**2 + (c[1]-tm_cy)**2) for c in tm_ca_coords)

    # Analysis region geometry
    if analysis_region == REGION_RADIAL_SHELL:
        inner_r     = tm_r
        outer_r     = tm_r + shell_width_nm
        region_area = math.pi * (outer_r**2 - inner_r**2)
    else:  # tm_footprint
        inner_r     = 0.0
        outer_r     = tm_r + footprint_margin_nm
        region_area = math.pi * outer_r**2

    # Count lipid residues whose XY COM falls within the analysis region
    n_in_region = 0
    for coords in lipid_coms.values():
        com_x = sum(c[0] for c in coords) / len(coords)
        com_y = sum(c[1] for c in coords) / len(coords)
        r = math.sqrt((com_x - tm_cx)**2 + (com_y - tm_cy)**2)
        if analysis_region == REGION_RADIAL_SHELL:
            in_region = inner_r <= r <= outer_r
        else:
            in_region = r <= outer_r
        if in_region:
            n_in_region += 1

    apl_local = region_area / n_in_region if n_in_region > 0 else None

    # Void fraction — only meaningful when reference_apl_nm2 is provided
    void_fraction: Optional[float] = None
    void_area:     Optional[float] = None
    if reference_apl_nm2 and reference_apl_nm2 > 0 and region_area > 0:
        n_expected    = region_area / reference_apl_nm2
        void_fraction = round(max(0.0, 1.0 - n_in_region / n_expected), 4)
        void_area     = round(void_fraction * region_area, 4)

    # Robust bilayer boundaries from headgroup and tail Z coordinates
    bilayer = _compute_robust_leaflet_boundaries(headgroup_z, tail_z)

    # Compute TM burial and soluble core penetration
    z_bot = bilayer.get("z_core_bot_median")
    z_top = bilayer.get("z_core_top_median")
    tm_burial = None
    soluble_penetration = None
    if z_bot is not None and z_top is not None:
        if tm_ca_coords:
            n_buried = sum(1 for c in tm_ca_coords if z_bot <= c[2] <= z_top)
            tm_burial = round(n_buried / len(tm_ca_coords), 4)
        else:
            tm_burial = 0.0
        if soluble_coords:
            soluble_penetration = sum(1 for c in soluble_coords if z_bot <= c[2] <= z_top)
        else:
            soluble_penetration = 0

    # Diagnostic warnings (non-blocking)
    warn: list[str] = []
    if n_in_region == 0 and n_lipids_total > 0:
        warn.append(
            f"No lipids found within {analysis_region} r≤{outer_r:.2f} nm "
            "— footprint_margin_nm may be too small."
        )
    if bilayer.get("n_hg", 0) < 10 and bilayer.get("n_hg", 0) > 0:
        warn.append(
            f"Only {bilayer['n_hg']} headgroup atoms found "
            "— bilayer boundary statistics may be unreliable."
        )
    if void_fraction is not None and void_fraction > 0.3:
        warn.append(
            f"void_fraction_local={void_fraction:.3f} > 0.3 — "
            "significant lipid gap near TM bundle; check convergence."
        )

    region_desc = f"r≤{outer_r:.2f} nm" if analysis_region == REGION_TM_FOOTPRINT \
                  else f"{inner_r:.2f}≤r≤{outer_r:.2f} nm"
    msg = (
        f"{n_in_region}/{n_lipids_total} lipids in {analysis_region} ({region_desc}); "
        + (f"void_fraction={void_fraction:.3f}" if void_fraction is not None
           else "reference_apl not provided — void metrics omitted")
    )

    return EmbeddingQualityResult(
        source_file=str(path),
        analysis_region=analysis_region,
        tm_center_x=round(tm_cx, 4),
        tm_center_y=round(tm_cy, 4),
        tm_footprint_radius_nm=round(tm_r, 4),
        local_region_area_nm2=round(region_area, 4),
        footprint_margin_nm=footprint_margin_nm,
        shell_width_nm=shell_width_nm,
        n_lipids_total=n_lipids_total,
        n_lipids_in_region=n_in_region,
        apl_local_nm2=round(apl_local, 4) if apl_local is not None else None,
        reference_apl_nm2=reference_apl_nm2,
        void_fraction_local=void_fraction,
        void_area_local_nm2=void_area,
        midplane_z_median=bilayer.get("midplane_z_median"),
        midplane_z_mean=bilayer.get("midplane_z_mean"),
        z_core_top_median=bilayer.get("z_core_top_median"),
        z_core_bot_median=bilayer.get("z_core_bot_median"),
        z_core_top_mean=bilayer.get("z_core_top_mean"),
        z_core_bot_mean=bilayer.get("z_core_bot_mean"),
        z_core_top_p25=bilayer.get("z_core_top_p25"),
        z_core_bot_p75=bilayer.get("z_core_bot_p75"),
        z_outer_top=bilayer.get("z_outer_top"),
        z_outer_bot=bilayer.get("z_outer_bot"),
        headgroup_z_iqr_top=bilayer.get("headgroup_z_iqr_top"),
        headgroup_z_iqr_bot=bilayer.get("headgroup_z_iqr_bot"),
        bilayer_core_thickness_nm=bilayer.get("bilayer_core_thickness_nm"),
        bilayer_outer_thickness_nm=bilayer.get("bilayer_outer_thickness_nm"),
        tm_burial_score=tm_burial,
        soluble_domain_core_penetration=soluble_penetration,
        n_headgroup_atoms_found=len(headgroup_z),
        n_tail_atoms_found=len(tail_z),
        n_tm_ca_atoms=n_tm_ca,
        message=msg,
        warnings=warn,
    )


def _compute_robust_leaflet_boundaries(
    headgroup_z: list[float],
    tail_z:      list[float],
) -> dict:
    """
    Compute bilayer Z boundaries using median, percentile, and IQR statistics.

    Returns a dict with keys: midplane_z_median, midplane_z_mean,
    z_core_top_median, z_core_bot_median, z_core_top_mean, z_core_bot_mean,
    z_core_top_p25, z_core_bot_p75, z_outer_top, z_outer_bot,
    headgroup_z_iqr_top, headgroup_z_iqr_bot,
    bilayer_core_thickness_nm, bilayer_outer_thickness_nm, n_hg.
    """
    result: dict = {"n_hg": len(headgroup_z)}
    if not headgroup_z:
        return result

    # Midplane from tail atoms (most reliable; falls back to headgroup midpoint)
    if tail_z:
        midplane_median = _median(tail_z)
        midplane_mean   = sum(tail_z) / len(tail_z)
    else:
        midplane_median = (_percentile(headgroup_z, 50.0))
        midplane_mean   = sum(headgroup_z) / len(headgroup_z)

    result["midplane_z_median"] = round(midplane_median, 4)
    result["midplane_z_mean"]   = round(midplane_mean, 4)

    # Split headgroups by midplane into leaflets
    top_hg = [z for z in headgroup_z if z > midplane_median]
    bot_hg = [z for z in headgroup_z if z <= midplane_median]

    if top_hg:
        result["z_core_top_median"]   = round(_median(top_hg), 4)
        result["z_core_top_mean"]     = round(sum(top_hg) / len(top_hg), 4)
        result["z_core_top_p25"]      = round(_percentile(top_hg, 25.0), 4)  # inner boundary
        result["z_outer_top"]         = round(max(top_hg), 4)
        result["headgroup_z_iqr_top"] = round(_iqr(top_hg), 4)

    if bot_hg:
        result["z_core_bot_median"]   = round(_median(bot_hg), 4)
        result["z_core_bot_mean"]     = round(sum(bot_hg) / len(bot_hg), 4)
        result["z_core_bot_p75"]      = round(_percentile(bot_hg, 75.0), 4)  # inner boundary
        result["z_outer_bot"]         = round(min(bot_hg), 4)
        result["headgroup_z_iqr_bot"] = round(_iqr(bot_hg), 4)

    if top_hg and bot_hg:
        result["bilayer_core_thickness_nm"]  = round(
            result["z_core_top_median"] - result["z_core_bot_median"], 4
        )
        result["bilayer_outer_thickness_nm"] = round(
            result["z_outer_top"] - result["z_outer_bot"], 4
        )

    return result


# ── Pure-Python statistics helpers (no numpy) ─────────────────────────────────

def _percentile(values: list[float], p: float) -> float:
    """p-th percentile (0–100) via linear interpolation on a sorted list."""
    if not values:
        raise ValueError("_percentile: empty list")
    sv  = sorted(values)
    n   = len(sv)
    if n == 1:
        return sv[0]
    idx = (p / 100.0) * (n - 1)
    lo  = int(idx)
    hi  = min(lo + 1, n - 1)
    frac = idx - lo
    return sv[lo] * (1.0 - frac) + sv[hi] * frac


def _median(values: list[float]) -> float:
    return _percentile(values, 50.0)


def _iqr(values: list[float]) -> float:
    return _percentile(values, 75.0) - _percentile(values, 25.0)
