"""JSON/text and region point-cloud exports for audit and VMD review."""
import json
from pathlib import Path
import numpy as np


def atlas_report(atlas, mappings, decisions, runtime_seconds=None):
    counts={}
    water_by_region={}
    for m in mappings:
        counts[m.classification]=counts.get(m.classification,0)+1
        key=str(m.region_id) if m.region_id is not None else "bulk"
        water_by_region[key]=water_by_region.get(key,0)+1
    regions=[]
    for r in atlas.regions:
        regions.append({"id":r.id,"classification":r.classification,"attributes":r.attributes,
            "evidence":r.evidence,"uncertainty":r.uncertainty,"voxel_count":r.voxel_count,
            "water_count":water_by_region.get(str(r.id),0)})
    removed={}
    for d in decisions:
        if d.action=="REMOVE": removed[d.reason]=removed.get(d.reason,0)+1
    return {"schema_version":atlas.schema_version,"regions":regions,
        "components":[vars(x) for x in atlas.components],
        "hierarchy":{"aggregates":atlas.hierarchy.aggregates,"scale_levels_nm":atlas.hierarchy.scale_levels_nm},
        "portals":[vars(x) for x in atlas.graph.portals],
        "pathways":atlas.pathways.networks,"water_counts_by_region_class":counts,
        "water_counts_by_region":water_by_region,
        "cleanup":{"policy":decisions[0].policy if decisions else "conservative",
                   "preserved":sum(d.action=="PRESERVE" for d in decisions),
                   "removed":sum(d.action=="REMOVE" for d in decisions),"removed_categories":removed},
        "uncertain_water_count":sum(bool(m.uncertainty) for m in mappings),
        "geometry":{"grid_shape":[int(x) for x in atlas.field.lattice.shape],"grid_steps_nm":atlas.field.lattice.steps.tolist(),
          "nominal_free_voxels":int(atlas.field.nominal_free.sum()),
          "certified_free_voxels":int(atlas.field.certified_free.sum()),
          "possible_free_voxels":int(atlas.field.possible_free.sum()),
          "membrane_axis":atlas.membrane.normal_axis,"supported_patch_fraction":float(atlas.membrane.supported.mean())},
        "warnings":atlas.warnings,"timings":{"geometry_seconds":atlas.timings.get("geometry_seconds"),
          "total_seconds":runtime_seconds}}


def write_outputs(report, atlas, system, mappings, decisions, output_dir, cleaned_path):
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    (out/"membrane_water_v2_report.json").write_text(json.dumps(report,indent=2,default=lambda o: o.item() if hasattr(o,"item") else list(o) if isinstance(o,tuple) else str(o)))
    lines=["SimForge membrane-water region atlas v2",f"Regions: {len(report['regions'])}",
           f"Water mappings: {sum(report['water_counts_by_region_class'].values())}",
           f"Preserved: {report['cleanup']['preserved']}; removed: {report['cleanup']['removed']}",
           f"Grid: {report['geometry']['grid_shape']} at steps {report['geometry']['grid_steps_nm']} nm"]
    lines += [f"region_{r['id']} {r['classification']} voxels={r['voxel_count']} waters={r['water_count']} evidence={r['evidence']} uncertainty={r['uncertainty']}" for r in report["regions"]]
    lines += ["Warnings:"]+["- "+str(x) for x in report["warnings"]]
    (out/"membrane_water_v2_report.txt").write_text("\n".join(lines)+"\n")
    classes={"protein_pore":"region_pore.pdb","vestibule":"region_vestibule.pdb",
      "internal_cavity":"region_cavity.pdb","lipid_facing_void":"region_lipid_void.pdb",
      "membrane_defect_candidate":"region_membrane_defect.pdb","ambiguous":"region_ambiguous.pdb"}
    for cls,name in classes.items():
        records=[]; serial=1
        for r in atlas.regions:
            if r.classification!=cls or r.voxel_indices is None: continue
            pts=atlas.field.lattice.points(r.voxel_indices)
            for x,y,z in pts:
                records.append(f"HETATM{serial%100000:5d}  C   REG A{r.id%10000:4d}    {x*10:8.3f}{y*10:8.3f}{z*10:8.3f}{1.:6.2f}{0.:6.2f}          C \n")
                serial+=1
        (out/name).write_text("".join(records)+"END\n")
    # category-specific molecule identity exports are machine-readable.
    cat={k:[] for k in classes}
    for m in mappings: cat.setdefault(m.classification,[]).append(m.molecule_uid)
    (out/"removed_water_categories.json").write_text(json.dumps(
      {"removed_molecule_uids":[d.molecule_id for d in decisions if d.action=="REMOVE"],
       "removed_by_category":{k:[m.molecule_uid for m,d in zip(mappings,decisions) if d.action=="REMOVE" and m.classification==k] for k in {m.classification for m,d in zip(mappings,decisions) if d.action=="REMOVE"}},
       "mappings_by_category":cat},indent=2))

