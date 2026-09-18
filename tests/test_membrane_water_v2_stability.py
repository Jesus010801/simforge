"""Multiscale region matching invariants, independent of cleanup policy."""
from types import SimpleNamespace
import numpy as np

from validators.membrane_water_v2.cell import PeriodicCell, Lattice
from validators.membrane_water_v2.models import Region, RegionGraph, Portal
from validators.membrane_water_v2.stability import analyse_multiscale


def atlas(regions, spacing, lengths=(4.,4.,4.)):
    cell=PeriodicCell(np.diag(lengths))
    lattice=Lattice.create(cell,spacing)
    shape=lattice.shape
    labels=np.zeros(shape,np.int32)
    items=[]
    for rid,vox in enumerate(regions,1):
        vox=np.asarray(vox,dtype=int)
        labels.ravel()[vox]=rid
        items.append(Region(rid,rid,len(vox),"internal_cavity",{"wall_species_counts":{"1":20},
            "upper_access":False,"lower_access":False},voxel_indices=vox))
    clear=np.full(shape,.3)
    return SimpleNamespace(geometry=SimpleNamespace(cell=cell),field=SimpleNamespace(lattice=lattice,clearance=clear),
        regions=items,graph=RegionGraph([r.id for r in items],[]),membrane=SimpleNamespace(normal_axis=2),labels=labels,exclusion=None)


def box(lattice,lo,hi):
    x,y,z=np.indices(lattice.shape)
    grids=(x,y,z);m=np.ones(lattice.shape,bool)
    for q,(a,b) in enumerate(zip(lo,hi)):
        f=(grids[q]+.5)/lattice.shape[q]
        m &= (f>=a)&(f<b)
    return np.flatnonzero(m)


def test_fine_fragmented_cavity_consolidates_under_persistent_parent():
    atl=[];sp=[.18,.20,.22,.24]
    for i,s in enumerate(sp):
        lat=Lattice.create(PeriodicCell(np.diag([4.,4.,4.])),s)
        whole=box(lat,(.4,.4,.4),(.6,.6,.6))
        if i==0:
            # Fine segmentation reports many micro-regions within one physical cavity.
            chunks=[box(lat,(.4+dx*.1,.4+dy*.1,.4+dz*.1),(.4+(dx+1)*.1,.4+(dy+1)*.1,.4+(dz+1)*.1)) for dx in range(2) for dy in range(2) for dz in range(2)]
            atl.append(atlas(chunks,s))
        else:atl.append(atlas([whole],s))
    result,_=analyse_multiscale(atl,sp)
    assert len(result.regions)==1
    assert result.regions[0].status=="unstable_subdivision"
    assert result.regions[0].persistence==1.0
    assert len(result.merge_hierarchy[0]["members"])>4


def test_two_distinct_channels_keep_separate_tracks_across_scales():
    atl=[];sp=[.18,.20,.22,.24]
    for s in sp:
        lat=Lattice.create(PeriodicCell(np.diag([8.,4.,4.])),s)
        atl.append(atlas([box(lat,(.1,.2,.1),(.2,.8,.9)),box(lat,(.6,.2,.1),(.7,.8,.9))],s,(8.,4.,4.)))
    result,_=analyse_multiscale(atl,sp)
    assert len(result.regions)==2
    assert all(r.persistence==1.0 for r in result.regions)
    assert result.regions[0].id!=result.regions[1].id


def test_one_scale_connection_is_reported_as_unresolved_topology():
    sp=[.18,.20,.22,.24];atl=[]
    for i,s in enumerate(sp):
        lat=Lattice.create(PeriodicCell(np.diag([8.,4.,4.])),s)
        left=box(lat,(.1,.2,.1),(.2,.8,.9));right=box(lat,(.6,.2,.1),(.7,.8,.9))
        if i==1: atl.append(atlas([np.concatenate([left,right])],s,(8.,4.,4.)))
        else: atl.append(atlas([left,right],s,(8.,4.,4.)))
    result,_=analyse_multiscale(atl,sp)
    assert any(t.get("status")=="unresolved_boundary" for t in result.transitions)
    assert sum(r.status=="stable_region" for r in result.regions)==2
    assert sum(r.status=="unresolved_boundary" for r in result.regions)==1
    assert any("transient split topology" in w for r in result.regions for w in r.uncertainty)
