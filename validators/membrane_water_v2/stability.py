"""Multiscale geometric stability analysis for solvent-space regions.

This module is deliberately observational: it creates scale-matched region
tracks and uncertainty metadata, but does not make cleanup decisions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import defaultdict
from time import perf_counter
import numpy as np
from scipy.spatial import cKDTree

SCHEMA = "membrane-water-stability/1.0"


@dataclass
class RegionScaleMetric:
    spacing_nm: float
    source_region_id: int
    overlap_score: float
    classification: str
    volume_nm3: float
    centroid_fractional: tuple[float, float, float]
    bottleneck_nm: float
    portal_signature: tuple
    adjacency_signature: tuple
    lining_signature: tuple


@dataclass
class PersistentRegion:
    id: str
    members: list[RegionScaleMetric]
    persistence: float
    spatial_overlap: float
    classification_stability: float
    lining_stability: float
    portal_stability: float
    adjacency_stability: float
    volume_stability: float
    centroid_drift_nm: float
    bottleneck_drift_nm: float
    status: str
    representative_classification: str
    destructive_eligible: bool = False
    uncertainty: list[str] = field(default_factory=list)


@dataclass
class MultiScaleStability:
    spacings_nm: list[float]
    regions: list[PersistentRegion]
    merge_hierarchy: list[dict]
    transitions: list[dict]
    schema_version: str = SCHEMA
    timings: dict = field(default_factory=dict)


def _fractional_centroid(atlas, points):
    f = atlas.geometry.cell.fractional(points) % 1.0
    # Circular means remain continuous for clouds crossing the periodic seam.
    z = np.exp(2j*np.pi*f).mean(axis=0)
    return (np.angle(z) / (2*np.pi)) % 1.0


def _portal_signature(atlas, rid):
    result=[]
    for p in atlas.graph.portals:
        if rid in p.regions:
            other = p.regions[1] if p.regions[0] == rid else p.regions[0]
            result.append(("bulk" if other in (0,-1) else "region", round(float(p.area_nm2),2), bool(p.certified)))
    return tuple(sorted(result))


def _adjacency_signature(atlas, rid):
    return tuple(sorted(("bulk" if (p.regions[1] if p.regions[0]==rid else p.regions[0]) in (0,-1)
                         else atlas.regions[(p.regions[1] if p.regions[0]==rid else p.regions[0])-1].classification)
                        for p in atlas.graph.portals if rid in p.regions))


def _metrics(atlas, spacing):
    result=[]
    cell=atlas.geometry.cell
    voxel_volume=abs(np.linalg.det(cell.matrix))/np.prod(atlas.field.lattice.shape)
    normal=atlas.membrane.normal_axis
    for r in atlas.regions:
        indices=np.asarray(r.voxel_indices,dtype=np.int64)
        points=atlas.field.lattice.points(indices)
        frac=_fractional_centroid(atlas,points)
        cl=atlas.field.clearance.ravel()[indices]
        walls=r.attributes.get("wall_species_counts",{})
        denom=sum(walls.values()) or 1
        # Preserve the species-specific evidence; stability compares it as a
        # normalized composition rather than allowing sample count to dominate.
        lining=tuple(sorted((str(k),round(v/denom,3)) for k,v in walls.items()))
        result.append(RegionScaleMetric(float(spacing),r.id,1.0,r.classification,
            float(len(indices)*voxel_volume),tuple(map(float,frac)),
            float(np.quantile(cl,.05)) if len(cl) else 0.0,
            _portal_signature(atlas,r.id),_adjacency_signature(atlas,r.id),lining))
    return result


def _cloud(atlas, region):
    points=atlas.field.lattice.points(region.voxel_indices)
    return atlas.geometry.cell.fractional(points)%1.0


def _pair_overlap(aa,ta,bb,tb,tolerance_nm):
    """Periodic bidirectional coverage, independent of grid/component IDs."""
    da=ta.query(bb,k=1)[0];db=tb.query(aa,k=1)[0]
    return float(np.mean(da<=tolerance_nm)),float(np.mean(db<=tolerance_nm))


def _fractional_radius(atlas,region):
    points=atlas.field.lattice.points(region.voxel_indices)
    center=_fractional_centroid(atlas,points)
    delta=(atlas.geometry.cell.fractional(points)-center+.5)%1-.5
    return float(np.linalg.norm(delta*atlas.geometry.cell.lengths,axis=1).max(initial=0.0))

def _jaccard_signature(a,b):
    aa,bb=set(a),set(b)
    return 1.0 if not aa and not bb else len(aa&bb)/len(aa|bb) if aa|bb else 1.0


def analyse_multiscale(atlases, spacings_nm, systems=None, mapping_factory=None):
    """Match independently discovered atlas regions across resolutions.

    ``atlases`` must all be generated from the same molecular coordinates and
    cell. ``systems`` is optional and is used only after geometry matching to
    count waters in the canonical (closest-to-0.20 nm) atlas.
    """
    start=perf_counter()
    if len(atlases)<2 or len(atlases)!=len(spacings_nm):
        raise ValueError("Provide matching atlases and at least two resolutions")
    order=np.argsort(spacings_nm);atlases=[atlases[i] for i in order];spacings_nm=[float(spacings_nm[i]) for i in order]
    if systems is not None: systems=[systems[i] for i in order]
    nscale=len(atlases);records=[];by_node={}
    for si,(atlas,spacing) in enumerate(zip(atlases,spacings_nm)):
        mm=_metrics(atlas,spacing)
        for r,m in zip(atlas.regions,mm):
            node=(si,r.id);by_node[node]=(atlas,r,m);records.append(node)
    edges=[];transitions=[];tol=max(spacings_nm)*np.sqrt(3)/2
    lengths=atlases[0].geometry.cell.lengths
    if any(not np.allclose(lengths,a.geometry.cell.lengths,atol=1e-6) for a in atlases[1:]):
        raise ValueError("Multiscale atlases must share one periodic cell")
    spatial={}
    for node in records:
        atlas,r,_=by_node[node];cloud=_cloud(atlas,r)*lengths
        spatial[node]=(cloud,cKDTree(cloud,boxsize=lengths),_fractional_radius(atlas,r))
    # Adjacent scales detect ordinary persistence; all scale pairs also catch
    # temporary subdivisions that disappear before the next sampled level.
    for si in range(nscale-1):
        left=[n for n in records if n[0]==si];right=[n for n in records if n[0]==si+1]
        scored=[]
        for na in left:
            ca0=np.asarray(by_node[na][2].centroid_fractional)*lengths
            cloud_a,tree_a,radius_a=spatial[na]
            for nb in right:
                cb0=np.asarray(by_node[nb][2].centroid_fractional)*lengths
                delta=(ca0-cb0+lengths/2)%lengths-lengths/2
                cloud_b,tree_b,radius_b=spatial[nb]
                if np.linalg.norm(delta)>radius_a+radius_b+tol: continue
                ca,cb=_pair_overlap(cloud_a,tree_a,cloud_b,tree_b,tol)
                # Strong one-sided containment links fine subdivisions to a coarser parent.
                score=max(ca,cb)
                if score>=.65:
                    scored.append((na,nb,score,ca,cb))
        # A split or merge remains represented as a multi-parent hierarchy.
        for na,nb,score,ca,cb in scored: edges.append((na,nb,score))
        children=defaultdict(list);parents=defaultdict(list)
        for na,nb,*_ in scored:children[na].append(nb);parents[nb].append(na)
        for n,vals in children.items():
            if len(vals)>1: transitions.append({"kind":"split","scale_nm":spacings_nm[si],"source":n[1],"source_node":n,"targets":[x[1] for x in vals],"target_nodes":vals})
        for n,vals in parents.items():
            if len(vals)>1: transitions.append({"kind":"merge","scale_nm":spacings_nm[si+1],"sources":[x[1] for x in vals],"source_nodes":vals,"target":n[1],"target_node":n})
    # Also link nonadjacent resolutions. This prevents a one-scale bridge from
    # collapsing two structures that are separate at the finer and coarser
    # observations. Such a bridge remains an explicit unresolved node.
    for si in range(nscale-2):
        for sj in range(si+2,nscale):
            left=[n for n in records if n[0]==si];right=[n for n in records if n[0]==sj]
            for na in left:
                c0=np.asarray(by_node[na][2].centroid_fractional)*lengths
                ca,ta,ra=spatial[na]
                for nb in right:
                    d0=np.asarray(by_node[nb][2].centroid_fractional)*lengths
                    cb,tb,rb=spatial[nb]
                    delta=(c0-d0+lengths/2)%lengths-lengths/2
                    if np.linalg.norm(delta)>ra+rb+tol: continue
                    x,y=_pair_overlap(ca,ta,cb,tb,tol)
                    score=max(x,y)
                    if score>=.65: edges.append((na,nb,score))
    # Tracks are connected components in the cross-resolution correspondence
    # graph. Split/merge nodes are retained as hierarchy events below.
    parent={n:n for n in records}
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]];x=parent[x]
        return x
    split_sources={t["source_node"] for t in transitions if t["kind"]=="split"}
    for a,b,_ in edges:
        if (a in split_sources or b in split_sources) and abs(a[0]-b[0])==1: continue
        x,y=find(a),find(b)
        if x!=y: parent[y]=x
    groups=defaultdict(list)
    for n in records:groups[find(n)].append(n)
    output=[];hierarchy=[]
    for ix,nodes in enumerate(sorted(groups.values(),key=lambda g:min(g))):
        members=[by_node[n][2] for n in sorted(nodes)]
        by_scale=defaultdict(list)
        for m in members: by_scale[m.spacing_nm].append(m)
        # A split component set at one resolution is one observation in the
        # persistence series. Sum volume and compute a periodic volume-weighted
        # centroid before comparing scales.
        observations=[]
        for spacing,parts in sorted(by_scale.items()):
            weights=np.array([max(m.volume_nm3,1e-12) for m in parts]); coords0=np.array([m.centroid_fractional for m in parts])
            z=np.sum(weights[:,None]*np.exp(2j*np.pi*coords0),axis=0)
            center=(np.angle(z)/(2*np.pi))%1
            classes0=[m.classification for m in parts]
            observations.append(RegionScaleMetric(spacing,-1,1.0,max(set(classes0),key=classes0.count),float(weights.sum()),tuple(center),
                float(min(m.bottleneck_nm for m in parts)),tuple(sorted(x for m in parts for x in m.portal_signature)),
                tuple(sorted(x for m in parts for x in m.adjacency_signature)),tuple(sorted(x for m in parts for x in m.lining_signature))))
        series=observations
        present=len(series)/nscale
        pair_scores=[]
        for a,b,s in edges:
            if a in nodes and b in nodes:pair_scores.append(s)
        overlap=float(np.mean(pair_scores)) if pair_scores else 0.0
        classes=[m.classification for m in series]
        cls_stab=max(classes.count(c) for c in set(classes))/len(classes)
        lining_stab=np.mean([_jaccard_signature(series[i].lining_signature,series[i+1].lining_signature)
                             for i in range(len(series)-1)]) if len(series)>1 else 0.0
        portal_stab=np.mean([_jaccard_signature(series[i].portal_signature,series[i+1].portal_signature)
                             for i in range(len(series)-1)]) if len(series)>1 else 0.0
        adjacency_stab=np.mean([_jaccard_signature(series[i].adjacency_signature,series[i+1].adjacency_signature)
                                for i in range(len(series)-1)]) if len(series)>1 else 0.0
        vols=np.array([m.volume_nm3 for m in series]);vol_stab=float(np.exp(-np.std(np.log(np.maximum(vols,1e-12))))) if len(vols)>1 else 0.0
        coords=np.array([m.centroid_fractional for m in series]);diff=[]
        lengths=atlases[0].geometry.cell.lengths
        for a,b in zip(coords[:-1],coords[1:]):
            df=(b-a+.5)%1-.5;diff.append(np.linalg.norm(df*lengths))
        drift=float(max(diff,default=0.0))
        bottleneck=float(np.ptp([m.bottleneck_nm for m in series])) if len(series)>1 else 0.0
        # Eligibility is intentionally stringent and informational only.
        unresolved_track=any(t.get("source_node") in nodes for t in transitions if t["kind"]=="split")
        subdivided=any(len(parts)>1 for parts in by_scale.values())
        stable=(not unresolved_track and not subdivided and present>=.75 and overlap>=.55 and cls_stab==1.0 and lining_stab>=.65
                and portal_stab>=.5 and adjacency_stab>=.5 and vol_stab>=.55
                and drift<=2*tol and bottleneck<=2*max(spacings_nm))
        status="unresolved_boundary" if unresolved_track else "stable_region" if stable else "unstable_subdivision" if present>=.5 or subdivided else "unresolved_boundary"
        warnings=[]
        if len(set(classes))>1:warnings.append("classification changes across resolutions")
        if len({m.lining_signature for m in series})>1:warnings.append("lining composition changes across resolutions")
        if present<.75:warnings.append("region absent at multiple resolutions")
        if unresolved_track:warnings.append("transient split topology; boundary unresolved")
        if subdivided:warnings.append("subregions merge into a persistent larger region")
        region=PersistentRegion(f"persistent_{ix+1:04d}",members,present,overlap,cls_stab,float(lining_stab),
             float(portal_stab),float(adjacency_stab),vol_stab,drift,bottleneck,status,
             max(set(classes),key=classes.count),False,warnings)
        output.append(region)
        hierarchy.append({"parent":region.id,"status":status,"members":[{"spacing_nm":m.spacing_nm,"region_id":m.source_region_id} for m in members]})
    # A transient coarse-scale merge of multiple otherwise persistent tracks is
    # an unresolved boundary; keep those tracks separate and flag the event.
    for t in transitions:
        if t["kind"]=="split": t["status"]="unresolved_boundary"
        else: t["status"]="unstable_subdivision"
    water_counts={}
    if systems is not None:
        from .water_mapping import map_waters
        ci=int(np.argmin(np.abs(np.asarray(spacings_nm)-.20)))
        amap=atlases[ci]; mappings=map_waters(systems[ci],amap,amap.exclusion)
        node_to_track={n:find(n) for n in records};root_to_id={find(n):f"persistent_{i+1:04d}" for i,g in enumerate(sorted(groups.values(),key=lambda g:min(g))) for n in g}
        region_to_track={rid:root_to_id[node_to_track[(ci,rid)]] for rid in [r.id for r in amap.regions]}
        for mapping in mappings:
            key=region_to_track.get(mapping.region_id,"unresolved")
            water_counts[key]=water_counts.get(key,0)+1
    stability=MultiScaleStability(spacings_nm,output,hierarchy,transitions,timings={"analysis_seconds":perf_counter()-start})
    return stability,water_counts


def stability_dict(stability,water_counts=None):
    water_counts=water_counts or {}
    return {"schema_version":stability.schema_version,"spacings_nm":stability.spacings_nm,
      "regions":[{**vars(r),"members":[vars(m) for m in r.members],"water_count":int(water_counts.get(r.id,0))} for r in stability.regions],
      "merge_hierarchy":stability.merge_hierarchy,"transitions":[{k:v for k,v in t.items() if not k.endswith("_node") and not k.endswith("_nodes")} for t in stability.transitions],"timings":stability.timings,
      "matching_parameters":{"spatial_tolerance_nm":max(stability.spacings_nm)*float(np.sqrt(3)/2),"minimum_containment_overlap":0.65,"persistence_required_for_stable":0.75,"automatic_destructive_eligibility":False}}
