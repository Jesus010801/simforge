#!/usr/bin/env python3
"""Analysis-only multiscale atlas report; never writes cleaned coordinates."""
import argparse
import json
from pathlib import Path
from time import perf_counter
import resource

from validators.membrane_water_v2.engine import build_atlas
from validators.membrane_water_v2.stability import analyse_multiscale, stability_dict


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("structure",type=Path)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--spacings",type=float,nargs="+",default=[.18,.20,.22,.24])
    args=parser.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    start=perf_counter(); atlases=[];systems=[];scale_rows=[]
    for spacing in args.spacings:
        t=perf_counter();system,atlas=build_atlas(args.structure,spacing=spacing)
        systems.append(system);atlases.append(atlas)
        scale_rows.append({"spacing_nm":spacing,"region_count":len(atlas.regions),
          "classification_counts":{c:sum(r.classification==c for r in atlas.regions) for c in sorted({r.classification for r in atlas.regions})},
          "portal_count":len(atlas.graph.portals),"nominal_free_voxels":int(atlas.field.nominal_free.sum()),
          "geometry_seconds":perf_counter()-t,"grid_shape":list(atlas.field.lattice.shape)})
    stability,water_counts=analyse_multiscale(atlases,args.spacings,systems)
    data=stability_dict(stability,water_counts)
    data.update({"input":str(args.structure),"scale_summaries":scale_rows,
      "total_runtime_seconds":perf_counter()-start,"peak_rss_kb":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      "cleanup_policy_invoked":False,
      "lateral_membrane_water_assessment":_lateral_assessment(data,scale_rows)})
    jp=args.out_dir/"membrane_water_v2_stability.json"
    jp.write_text(json.dumps(data,indent=2,sort_keys=True,default=_json_default)+"\n")
    (args.out_dir/"membrane_water_v2_stability.txt").write_text(_text(data))
    print(f"JSON: {jp}\nTXT: {args.out_dir/'membrane_water_v2_stability.txt'}")


def _json_default(value):
    if hasattr(value, "item"): return value.item()
    if hasattr(value, "tolist"): return value.tolist()
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")

def _lateral_assessment(data,scale_rows):
    classes={"lipid_facing_void","membrane_defect_candidate","membrane_defect"}
    candidate=[]
    for r in data["regions"]:
        if r["representative_classification"] in classes:
            candidate.append(r)
    if not candidate:
        return {"outcome":"unresolved_geometry","reason":"No persistent lipid-facing or membrane-defect class was independently supported by scale matching; broad lateral water remains ambiguous unless mapped to a candidate below.","candidate_region_ids":[]}
    stable=[r for r in candidate if r["status"]=="stable_region"]
    unstable=[r for r in candidate if r["status"]=="unstable_subdivision"]
    unresolved=[r for r in candidate if r["status"]=="unresolved_boundary"]
    outcome="unresolved_geometry" if unresolved or not stable else "unstable_subdivisions_of_larger_defect" if unstable else "stable_lipid_facing_or_defect_regions"
    return {"outcome":outcome,"candidate_region_ids":[r["id"] for r in candidate],
      "stable_ids":[r["id"] for r in stable],"unstable_ids":[r["id"] for r in unstable],"unresolved_ids":[r["id"] for r in unresolved]}


def _text(data):
    lines=["SimForge membrane-water multiscale stability analysis",f"Input: {data['input']}",
      f"Scales (nm): {', '.join(map(str,data['spacings_nm']))}",
      f"Runtime: {data['total_runtime_seconds']:.1f} s; peak RSS: {data['peak_rss_kb']} kB",
      "Cleanup policy invoked: no", "", "Per-resolution atlas:"]
    for row in data["scale_summaries"]:
        lines.append(f"  {row['spacing_nm']:.2f} nm: {row['region_count']} elementary regions; {row['portal_count']} portals; {row['classification_counts']}; {row['geometry_seconds']:.1f} s")
    lines += ["",f"Persistent tracks: {len(data['regions'])}",
      f"Stable: {sum(r['status']=='stable_region' for r in data['regions'])}; unstable subdivisions: {sum(r['status']=='unstable_subdivision' for r in data['regions'])}; unresolved: {sum(r['status']=='unresolved_boundary' for r in data['regions'])}",
      "", "Regions:"]
    for r in data["regions"]:
        lines.append(f"  {r['id']} {r['status']} class={r['representative_classification']} persistence={r['persistence']:.2f} overlap={r['spatial_overlap']:.2f} class_stability={r['classification_stability']:.2f} lining={r['lining_stability']:.2f} portals={r['portal_stability']:.2f} adjacency={r['adjacency_stability']:.2f} volume={r['volume_stability']:.2f} centroid_drift_nm={r['centroid_drift_nm']:.3f} bottleneck_drift_nm={r['bottleneck_drift_nm']:.3f} waters={r['water_count']}")
        for warning in r["uncertainty"]:lines.append(f"    uncertainty: {warning}")
    lines += ["", "Scale transitions:"]
    lines += [f"  {t}" for t in data["transitions"]] or ["  none"]
    lines += ["", "Merge hierarchy:"]
    lines += [f"  {h['parent']} ({h['status']}): {h['members']}" for h in data["merge_hierarchy"]]
    lines += ["", "Lateral membrane-water assessment:",f"  {data['lateral_membrane_water_assessment']}",""]
    return "\n".join(lines)


if __name__=="__main__":main()
