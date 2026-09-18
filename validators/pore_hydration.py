"""Production pore-hydration adapter backed only by membrane-water atlas v2."""
from pathlib import Path
import json
from validators.membrane_water_v2.api import classify_membrane_water, run_cleanup

CATEGORIES=("bulk_water","lipid_core_outside_protein","channel_lumen_water",
            "protein_internal_cavity_water","clash_water","ambiguous")


def _adapt(report):
    mapping={"bulk":"bulk_water","protein_pore":"channel_lumen_water","vestibule":"protein_internal_cavity_water",
       "internal_cavity":"protein_internal_cavity_water","buried_pocket":"protein_internal_cavity_water",
       "membrane_defect_candidate":"lipid_core_outside_protein","lipid_facing_void":"lipid_core_outside_protein",
       "steric_clash":"clash_water","ambiguous":"ambiguous"}
    ids={c:[] for c in CATEGORIES}
    for item in report["water_mappings"]:
        ids[mapping.get(item["classification"],"ambiguous")].append(item["molecule_uid"])
    counts={f"n_{c}":len(ids[c]) for c in CATEGORIES}
    return ids,counts


def classify_pore_hydration(gro_path,*,tm_residues=None,lipid_resnames=None,output_dir=None,enabled=True,
        grid_spacing_nm=.15,protein_padding_nm=.14,protein_clash_nm=.20,lipid_clash_nm=.18):
    out=Path(output_dir or Path(gro_path).parent);out.mkdir(parents=True,exist_ok=True)
    if not enabled:
        result={"enabled":False,"water_ids_by_class":{c:[] for c in CATEGORIES},"warnings":[],"coordinate_file_modified":False}
        (out/"pore_hydration_report.json").write_text(json.dumps(result,indent=2));return result
    report=classify_membrane_water(gro_path,tm_residues=tm_residues,lipid_resnames=lipid_resnames,
       grid_spacing_nm=grid_spacing_nm,solvent_probe_radius_nm=protein_padding_nm)
    ids,counts=_adapt(report)
    result={"enabled":True,"channel_like_system_detected":bool(counts["n_channel_lumen_water"]),
       "n_waters_total":sum(counts.values()),**counts,"membrane_slab_z_min":None,"membrane_slab_z_max":None,
       "pore_hydration_score":0.0,"lumen_continuity_score":0.0,"warnings":report["warnings"],
       "coordinate_file_modified":False,"water_ids_by_class":ids,"region_count":len(report["regions"]),
       "region_graph":report["region_graph"],"pathways":report["pathways"]}
    lumen=len(ids["channel_lumen_water"]);internal=lumen+len(ids["protein_internal_cavity_water"])
    result["lumen_continuity_score"]=lumen/internal if internal else 0.0
    (out/"pore_hydration_report.json").write_text(json.dumps(result,indent=2))
    return result


def clean_water_channel_aware(gro_in,gro_out,*,tm_residues=None,output_dir=None,topol_in=None,topol_out=None):
    out=Path(output_dir or Path(gro_out).parent);out.mkdir(parents=True,exist_ok=True)
    result=run_cleanup(gro_in,gro_out,tm_residues=tm_residues,output_dir=out,
        topol_in=topol_in,topol_out=topol_out,cleanup_mode="conservative")
    report=result["classification"];ids,counts=_adapt(report)
    cleaned={"input_water_molecules":sum(counts.values()),
      "n_water_molecules_removed":result["n_water_molecules_removed"],
      "n_water_atoms_removed":result["n_water_atoms_removed"],
      "n_water_oxygens_remaining_in_core":result["n_water_oxygens_remaining_in_core"],
      "n_pore_waters_preserved":counts["n_channel_lumen_water"],
      "n_internal_cavity_waters_preserved":counts["n_protein_internal_cavity_water"],
      "output_gro_path":str(gro_out),"topology_updated":result["topology_updated"],
      "cleanup_passed":result["cleanup_passed"],"pore_aware":True,"warnings":result["warnings"],
      "schema_version":result["schema_version"],
      "confirmed_hard_clashes_remaining":result["confirmed_hard_clashes_remaining"],
      "supported_membrane_defects_remaining":result["supported_membrane_defects_remaining"]}
    (out/"clean_water_report.json").write_text(json.dumps(cleaned,indent=2))
    (out/"water_report.json").write_text(json.dumps({"passed":result["cleanup_passed"],
      "n_water_molecules_removed":result["n_water_molecules_removed"],
      "n_waters_remaining":sum(result["post_cleanup_classification"]["water_counts"].values()),
      "errors":result["validation"]["errors"],"warnings":result["warnings"],"confidence":report["confidence"]},indent=2))
    hyd={"enabled":True,"channel_like_system_detected":bool(counts["n_channel_lumen_water"]),
       **counts,"water_ids_by_class":ids,"coordinate_file_modified":bool(result["n_water_molecules_removed"]),
       "n_water_molecules_removed":result["n_water_molecules_removed"],"output_gro_path":str(gro_out),
       "warnings":result["warnings"],"region_graph":report["region_graph"],"pathways":report["pathways"]}
    (out/"pore_hydration_report.json").write_text(json.dumps(hyd,indent=2))
    return hyd
