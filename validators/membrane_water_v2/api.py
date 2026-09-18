"""Stable public API for v2 atlas and cleanup policy."""
from pathlib import Path
from collections import defaultdict
from .engine import build_atlas
from .models import ExclusionModel
from .water_mapping import map_waters
from .policy import decide
from .cleanup import clean_gro


def _classify(gro_path, tm_residues=None, lipid_resnames=None, grid_spacing_nm=.15,
              solvent_probe_radius_nm=.14, phase=(.5,.5,.5), cleanup_mode="conservative"):
    system, atlas = build_atlas(gro_path, spacing=grid_spacing_nm, phase=phase,
        exclusion=ExclusionModel(probe_radius_nm=solvent_probe_radius_nm))
    mappings = map_waters(system, atlas, atlas.exclusion)
    categories=defaultdict(list)
    for m in mappings:
        categories[m.classification].append(m.molecule_uid)
    classes=("bulk","protein_pore","vestibule","internal_cavity","buried_pocket","lipid_facing_void",
             "membrane_defect_candidate","ambiguous","steric_clash")
    counts={c:len(categories[c]) for c in classes}
    # Compatibility lipid report is derived from v2 pore-region labels at
    # actual lipid atom positions; it does not classify pores from TM circles.
    pore_ids={r.id for r in atlas.regions if r.classification=="protein_pore"}
    lipid_molecules={}
    for gi,mid in enumerate(system.geometry.molecule_ids):
        if system.geometry.species[gi]==2:
            lipid_molecules.setdefault(int(mid),[]).append(gi)
    forbidden=[]
    for mid,indices in lipid_molecules.items():
        pts=system.geometry.coordinates[indices]
        ijk=atlas.field.lattice.indices(pts)
        hit=False
        for ix,iy,iz in ijk:
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    for dz in (-1,0,1):
                        label=int(atlas.labels[(ix+dx)%atlas.labels.shape[0],(iy+dy)%atlas.labels.shape[1],max(0,min(int(iz+dz),atlas.labels.shape[2]-1))])
                        if label in pore_ids:
                            hit=True;break
                    if hit:break
                if hit:break
            if hit:break
        if hit:forbidden.append(mid)
    lipid_count=len(lipid_molecules)
    public={"schema_version":atlas.schema_version,"engine":"membrane_water_v2",
        "regions":[{"region_id":r.id,"classification":r.classification,"voxel_count":r.voxel_count,
                    "attributes":r.attributes,"evidence":r.evidence,"uncertainty":r.uncertainty} for r in atlas.regions],
        "region_graph":{"nodes":atlas.graph.nodes,"portals":[vars(p) for p in atlas.graph.portals]},
        "hierarchy":{"components":[vars(c) for c in atlas.hierarchy.components],"aggregates":atlas.hierarchy.aggregates},
        "pathways":atlas.pathways.networks,
        "water_mappings":[{"molecule_uid":m.molecule_uid,"region_id":m.region_id,"classification":m.classification,
                           "hard_clash":m.hard_clash,"membrane_interior":m.membrane_interior,"uncertainty":m.uncertainty} for m in mappings],
        "water_counts":counts,"water_resid_by_category":dict(categories),
        "cleanup_mode":cleanup_mode,"confidence":1.0 if not atlas.warnings else .5,
        "warnings":atlas.warnings}
    # Compatibility aliases use unique water-instance IDs, never wrapped GRO residue numbers.
    n_pore=sum(r.classification=="protein_pore" for r in atlas.regions)
    n_vest=sum(r.classification=="vestibule" for r in atlas.regions)
    n_cavity=sum(r.classification in {"internal_cavity","buried_pocket"} for r in atlas.regions)
    uid_to_resid={w.uid:w.residue_number for w in system.waters}
    aliases={name:[uid_to_resid[x] for x in values] for name,values in categories.items()}
    aliases["pore"]=list(aliases["protein_pore"])
    aliases["isolated_cavity"]=list(aliases["internal_cavity"])+list(aliases["buried_pocket"])
    aliases["membrane_defect"]=list(aliases["membrane_defect_candidate"])
    public["water_counts"]["pore"]=counts["protein_pore"]
    public["water_counts"]["membrane_defect"]=counts["membrane_defect_candidate"]
    core_range=None
    if atlas.membrane.normal_axis==2:
        center=atlas.membrane.center_fraction*atlas.geometry.cell.lengths[2]
        core_range=[float(center+atlas.membrane.lower.min()+.35),float(center+atlas.membrane.upper.max()-.35)]
    public["membrane_core_z_range"]=core_range
    public["voxel_diagnostics"]={"grid_shape":[int(x) for x in atlas.field.lattice.shape],
      "nominal_free":int(atlas.field.nominal_free.sum()),"hard_excluded":int(atlas.field.hard_excluded.sum()),
      "region_voxels":int((atlas.labels>0).sum())}
    public.update({"enabled":True,"grid_spacing_nm":grid_spacing_nm,"protein_padding_nm":solvent_probe_radius_nm,
      "membrane_core_z_range":core_range,"n_regions_detected":len(atlas.regions),"n_transmembrane_pores":n_pore,
      "n_one_sided_vestibules":n_vest,"n_isolated_cavities":n_cavity,"n_bulk_waters":counts["bulk"],
      "n_pore_waters":counts["protein_pore"],"n_membrane_core_waters":sum(m.membrane_interior for m in mappings),
      "n_isolated_internal_waters":counts["internal_cavity"],"n_valid_external_lipids":lipid_count-len(forbidden),
      "n_forbidden_pore_lipids":len(forbidden),"forbidden_pore_lipid_resids":forbidden,"pore_water_resids":aliases["pore"],
      "water_resid_by_category":aliases,"removed_water_resids_candidate":sorted({uid_to_resid[m.molecule_uid] for m in mappings if m.hard_clash or (m.classification=="membrane_defect_candidate" and not m.uncertainty)})})
    return public,system,atlas,mappings


def classify_membrane_water(gro_path, *, tm_residues=None, lipid_resnames=None, grid_spacing_nm=.15,
                            solvent_probe_radius_nm=.14, phase=(.5,.5,.5), cleanup_mode="conservative", **kwargs):
    return _classify(gro_path,tm_residues,lipid_resnames,grid_spacing_nm,solvent_probe_radius_nm,phase,cleanup_mode)[0]


def run_cleanup(gro_in, gro_out, *, tm_residues=None, lipid_resnames=None, lipid="DPPC", forcefield="opls-aa",
                topol_in=None, topol_out=None, output_dir=None, cleanup_mode="conservative",
                remove_membrane_core_water=True, remove_isolated_internal_water=False,
                grid_spacing_nm=.15, solvent_probe_radius_nm=.14, **kwargs):
    from time import perf_counter
    from .report import atlas_report, write_outputs
    t=perf_counter()
    report,system,atlas,mappings=_classify(gro_in,tm_residues,lipid_resnames,grid_spacing_nm,solvent_probe_radius_nm,
                                           cleanup_mode=cleanup_mode)
    decisions=[decide(m,cleanup_mode) for m in mappings]
    clean_gro(gro_in,gro_out,mappings,cleanup_mode)
    topology_updated=False;topology_report=None;ion_sync_reports={}
    if topol_in and topol_out:
        from shutil import copy2
        from validators.topology_sync import sync_topology_molecule_count
        copy2(topol_in,topol_out)
        from .identity import load_system
        retained_water_count=len(load_system(gro_out,lipid_resnames).waters)
        topology_report=sync_topology_molecule_count(gro_out,topol_out,"SOL",
            report_path=(Path(output_dir or Path(gro_out).parent)/"topology_sync_water_report.json"),
            molecule_count_override=retained_water_count)
        topology_updated=bool(topology_report.get("synchronized"))
        # Legacy system topologies can carry stale monatomic-ion counts (for
        # example, an ionized topology paired with a pre-ionization GRO). Keep
        # the output topology coordinate-consistent for recognized ions too.
        from .identity import IONS
        for ion_name in sorted(IONS):
            ion_sync_reports[ion_name]=sync_topology_molecule_count(
                gro_out,topol_out,ion_name,
                report_path=Path(output_dir or Path(gro_out).parent)/f"topology_sync_{ion_name}_report.json",
                atoms_per_molecule=1)
            if ion_sync_reports[ion_name].get("old_count") is not None:
                topology_updated=topology_updated and bool(ion_sync_reports[ion_name].get("synchronized"))
    report_out=atlas_report(atlas,mappings,decisions,perf_counter()-t)
    if output_dir:
        write_outputs(report_out,atlas,system,mappings,decisions,output_dir,gro_out)
    post_counts={c:sum(1 for m,d in zip(mappings,decisions) if d.action=="PRESERVE" and m.classification==c)
                 for c in report["water_counts"]}
    errors=[]
    from utils.gro_parser import parse_gro
    before=parse_gro(gro_in);after=parse_gro(gro_out)
    expected=len(before.atoms)-sum(len(w.atom_indices) for w in system.waters if any(d.molecule_id==w.uid and d.action=="REMOVE" for d in decisions))
    if len(after.atoms)!=expected: errors.append("atom count mismatch after cleanup")
    validation={"passed":not errors,"errors":errors,"warnings":[]}
    remaining_core=sum(1 for m,d in zip(mappings,decisions) if d.action=="PRESERVE" and m.membrane_interior)
    return {"classification":report,"post_cleanup_classification":{"water_counts":post_counts},
       "removed_water_resids":[d.molecule_id for d in decisions if d.action=="REMOVE"],
       "n_water_molecules_removed":sum(d.action=="REMOVE" for d in decisions),
       "n_water_atoms_removed":sum(len(w.atom_indices) for w in system.waters if any(d.molecule_id==w.uid and d.action=="REMOVE" for d in decisions)),
       "n_water_oxygens_remaining_in_core":remaining_core,"output_gro_path":str(gro_out),
       "topology_updated":topology_updated,"topology_sync_report":topology_report,"topology_ion_sync_reports":ion_sync_reports,
       "validation":validation,"cleanup_passed":validation["passed"],"cleanup_mode":cleanup_mode,
       "warnings":report["warnings"],"schema_version":atlas.schema_version,
       "confirmed_hard_clashes_remaining":sum(m.hard_clash and d.action=="PRESERVE" for m,d in zip(mappings,decisions)),
       "supported_membrane_defects_remaining":sum(m.classification=="membrane_defect_candidate" and d.action=="PRESERVE" and not m.uncertainty for m,d in zip(mappings,decisions)),
       "atlas_report":report_out}

