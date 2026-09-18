"""
validators/membrane_water/report_text.py

Human-readable explanation of why each water population was preserved or
removed, generated from a classify_membrane_water() report plus the
resolved removal set.
"""
from __future__ import annotations

from pathlib import Path

_EXPLANATIONS = {
    "bulk": "outside the membrane core (above/below the leaflets) — always preserved.",
    "pore": "inside a solvent-accessible channel that is protein-walled and connects both "
            "sides of the membrane — a genuine transmembrane pore/lumen. Always preserved.",
    "vestibule": "inside a solvent-accessible cavity that is protein-walled and connects to "
                 "only one side of the membrane — an entrance vestibule. Always preserved.",
    "isolated_cavity": "inside a fully enclosed internal cavity not connected to bulk on "
                        "either side. Preserved by default (remove_isolated_internal_water=True "
                        "to remove).",
    "membrane_defect": "inside the membrane core but the surrounding wall is predominantly "
                        "lipid, not protein — a lateral lipid-packing defect or non-physical "
                        "bridge through the hydrophobic core, not a real pore. Removed.",
    "steric_clash": "overlaps a protein/lipid atom's van der Waals volume (a real steric clash). Removed.",
    "ambiguous": "classification confidence was below the safety threshold. Preserved by default "
                 "(fail-safe: never removed without confident classification).",
}


def build_human_report(report: dict, removed_resids: "set[int]") -> str:
    lines: list[str] = []
    lines.append("Membrane water cleanup report")
    lines.append("=" * 32)
    lines.append("")
    tm = report["tm_annotation"]
    lines.append(
        f"TM annotation: source={tm['source']}, residues={tm['residue_count']}, "
        f"segments={tm['segments']}, confidence={tm['confidence']:.2f}"
    )
    lines.append(f"Membrane core Z-range (global fallback slab): {report['membrane_core_z_range']} nm")
    lines.append(f"Cleanup mode used: {report['cleanup_mode']} (requested: {report.get('requested_cleanup_mode', report['cleanup_mode'])})")
    lines.append(f"Overall confidence: {report['confidence']:.2f}")
    if report.get("warnings"):
        lines.append("")
        lines.append("Warnings:")
        for w in report["warnings"]:
            lines.append(f"  - {w}")
    lines.append("")
    lines.append(f"Detected {report['n_regions_detected']} membrane-core region(s): "
                 f"{report['n_transmembrane_pores']} pore(s), "
                 f"{report['n_one_sided_vestibules']} vestibule(s), "
                 f"{report['n_isolated_cavities']} isolated cavity(ies).")
    lines.append("")
    vd = report.get("voxel_diagnostics")
    if vd:
        lines.append(
            "Voxel geometry (computed from protein/lipid atoms only, before any "
            "water position is consulted):"
        )
        lines.append(f"  grid: {vd['grid_shape']} @ {vd['grid_spacing_nm']} nm  "
                      f"({vd['total_voxels']} voxels total)")
        lines.append(f"  protein-solid voxels:            {vd['protein_solid_voxels']}")
        lines.append(f"  lipid-solid voxels:               {vd['lipid_solid_voxels']}")
        lines.append(f"  solvent-accessible void voxels:   {vd['solvent_accessible_void_voxels']}")
        lines.append(f"  core-band void voxels analyzed:   {vd['core_band_void_voxels_analyzed']}")
        lines.append(f"  -> pore voxels:                   {vd['pore_voxels']}")
        lines.append(f"  -> vestibule voxels:              {vd['vestibule_voxels']}")
        lines.append(f"  -> isolated-cavity voxels:        {vd['isolated_cavity_voxels']}")
        lines.append(f"  -> membrane-defect voxels:        {vd['membrane_defect_voxels']}")
        lines.append(f"  -> ambiguous voxels:              {vd['ambiguous_voxels']}")
        lines.append("")
    nb = report.get("non_bulk_water_breakdown")
    if nb:
        lines.append("Non-bulk ('membrane-core') water, by region classification:")
        for category, count in nb.items():
            lines.append(f"  {category}: {count}")
        lines.append("")
    lines.append("Water populations:")
    by_cat = report["water_resid_by_category"]
    for category, explanation in _EXPLANATIONS.items():
        resids = by_cat.get(category, [])
        n_removed = sum(1 for r in resids if r in removed_resids)
        n_preserved = len(resids) - n_removed
        lines.append(f"  {category} ({len(resids)} total): {explanation}")
        lines.append(f"      preserved={n_preserved}, removed={n_removed}")
    lines.append("")
    lines.append(f"Total removed: {len(removed_resids)}")
    lines.append(f"Total preserved: {sum(len(v) for v in by_cat.values()) - len(removed_resids)}")
    return "\n".join(lines) + "\n"


def write_human_report(report: dict, removed_resids: "set[int]", path: "Path | str") -> None:
    Path(path).write_text(build_human_report(report, removed_resids))
