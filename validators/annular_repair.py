"""Conservative post-shrink TM/slab annular lipid relocation."""
from __future__ import annotations
import json, math
from pathlib import Path
from collections import defaultdict
from .protein_membrane_interface import evaluate_protein_membrane_interface
from .membrane_validators import detect_trapped_lipids


def _parse(path):
    lines=Path(path).read_text().splitlines(); atoms=[]
    for i,l in enumerate(lines[2:-1],2):
        if len(l)<44: continue
        try: atoms.append((i,l,int(l[:5]),l[5:10].strip(),l[10:15].strip(),float(l[20:28]),float(l[28:36]),float(l[36:44])))
        except ValueError: pass
    return lines,atoms

def repair_annular_packing(gro_path, tm_residues, output_path=None, report_path=None, enabled=False, allow_insertion=False):
    gro=Path(gro_path); out=Path(output_path or gro); rep=Path(report_path or gro.parent/'annular_repair_report.json')
    before= evaluate_protein_membrane_interface(gro, set(tm_residues), output_dir=gro.parent)
    diag_before=detect_trapped_lipids(gro, tm_residues=set(tm_residues))
    report={'enabled':bool(enabled),'mode':'insert_allowed' if allow_insertion else 'relocate_only','n_exposed_clusters_before':len(before.get('gap_clusters',[])),'n_exposed_clusters_after':len(before.get('gap_clusters',[])), 'fraction_covered_before':before.get('fraction_covered'),'fraction_covered_after':before.get('fraction_covered'),'exposed_gap_fraction_before':before.get('fraction_exposed_gap'),'exposed_gap_fraction_after':before.get('fraction_exposed_gap'),'p90_nearest_lipid_distance_before':before.get('p90_nearest_lipid_distance'),'p90_nearest_lipid_distance_after':before.get('p90_nearest_lipid_distance'),'n_relocated_lipids':0,'n_inserted_lipids':0,'n_rejected_candidates':0,'rejection_reasons':{},'true_cavity_trapped_before':diag_before.n_true_cavity_trapped,'true_cavity_trapped_after':diag_before.n_true_cavity_trapped,'coordinate_file_modified':False,'topology_modified':False}
    if not enabled or not before.get('gap_clusters'):
        rep.write_text(json.dumps(report,indent=2)); return report
    lines,atoms=_parse(gro); lip=defaultdict(list); prot=[]
    for idx,l,res, rn, an,x,y,z in atoms:
        if rn in {'DPP','DPPC'}: lip[res].append((idx,x,y,z))
        elif rn not in {'SOL','NA','CL','K','MG','CA'}: prot.append((x,y,z))
    # choose one outer lipid per largest cluster; move only if COM is far from receptor
    used=set(); moved=0
    for cl in sorted(before['gap_clusters'], key=lambda c:c.get('size',0), reverse=True)[:3]:
        cx,cy,cz=cl['centroid_nm']; best=None
        for res,aa in lip.items():
            if res in used: continue
            ox=sum(a[1] for a in aa)/len(aa); oy=sum(a[2] for a in aa)/len(aa); oz=sum(a[3] for a in aa)/len(aa)
            d=math.hypot(ox-cx,oy-cy)
            if d<2.0: continue
            score=d
            if best is None or score>best[0]: best=(score,res,ox,oy,oz)
        if best is None: report['n_rejected_candidates']+=1; report['rejection_reasons']['no_outer_annular_candidate']=report['rejection_reasons'].get('no_outer_annular_candidate',0)+1; continue
        _,res,ox,oy,oz=best; dx,dy=cx-ox,cy-oy
        # reject if translated lipid atom approaches any protein atom below 0.22 nm
        aa=lip[res]; shifted=[(x+dx,y+dy,z) for _,x,y,z in aa]
        if any(math.dist(q,p)<0.22 for q in shifted for p in prot):
            report['n_rejected_candidates']+=1; report['rejection_reasons']['tm_or_protein_hard_clash']=report['rejection_reasons'].get('tm_or_protein_hard_clash',0)+1; continue
        for (idx,_,_,_),(_,nx,ny,nz) in zip(aa,shifted):
            old=lines[idx]; lines[idx]=old[:20]+f'{nx:8.3f}{ny:8.3f}{nz:8.3f}'+old[44:]
        used.add(res); moved+=1
    if moved:
        Path(out).write_text('\n'.join(lines)+'\n'); report['coordinate_file_modified']=True; report['n_relocated_lipids']=moved
        after=evaluate_protein_membrane_interface(out,set(tm_residues),output_dir=gro.parent); da=detect_trapped_lipids(out,tm_residues=set(tm_residues))
        report.update(n_exposed_clusters_after=len(after.get('gap_clusters',[])),fraction_covered_after=after.get('fraction_covered'),exposed_gap_fraction_after=after.get('fraction_exposed_gap'),p90_nearest_lipid_distance_after=after.get('p90_nearest_lipid_distance'),true_cavity_trapped_after=da.n_true_cavity_trapped)
    rep.write_text(json.dumps(report,indent=2)); return report
