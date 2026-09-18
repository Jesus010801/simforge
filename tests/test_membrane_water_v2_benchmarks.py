"""Topology-first scientific fixtures A-P for the water-independent atlas."""
from types import SimpleNamespace
import numpy as np
import pytest
from validators.membrane_water_v2.cell import PeriodicCell,Lattice
from validators.membrane_water_v2.models import MembraneAtlas
from validators.membrane_water_v2.regions import discover_regions,local_classify
from validators.membrane_water_v2.region_graph import build_graph
from validators.membrane_water_v2.pathways import build_pathways


def fixture(tubes=(), defect=(), cavities=(), curvature=False):
    cell=PeriodicCell(np.diag([12.0,3.2,3.0])); lattice=Lattice.create(cell,.1)
    nx,ny,nz=lattice.shape
    lower=np.full((nx,ny),-.9); upper=np.full((nx,ny),.9)
    if curvature:
        xx=np.arange(nx)[:,None]*.1
        lower += .12*np.sin(2*np.pi*xx/12.0)
        upper += .12*np.sin(2*np.pi*xx/12.0)
    support=np.ones((nx,ny),bool)
    membrane=MembraneAtlas(2,(0,1),.5,lower,upper,support,np.array([0,0,1.]))
    free=np.zeros(lattice.shape,bool); species=np.zeros(lattice.shape,np.uint8)
    for x,y,cls,shift in list(tubes)+list(defect):
        for z in range(0,nz):
            xx=(x+int((z-15)*shift))%nx
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    free[(xx+dx)%nx,(y+dy)%ny,z]=True
                    species[(xx+dx)%nx,(y+dy)%ny,z]=cls
    for x,y,z in cavities:
        for dx in (-2,-1,0,1,2):
            for dy in (-2,-1,0,1,2):
                for dz in (-2,-1,0,1,2):
                    if dx*dx+dy*dy+dz*dz<=4:
                        free[(x+dx)%nx,(y+dy)%ny,z+dz]=True
                        species[(x+dx)%nx,(y+dy)%ny,z+dz]=1
    field=SimpleNamespace(lattice=lattice,nominal_free=free,clearance=free.astype(float)*.3,
                          nearest_species=species,uncertain=np.zeros_like(free))
    return field,membrane


def analyse(*args,**kwargs):
    field,mem=fixture(*args,**kwargs)
    labels,regions,components,sides=discover_regions(field,mem)
    local_classify(regions,mem,min_wall_samples=1)
    graph=build_graph(labels,regions,field,sides)
    return field,labels,regions,components,graph


@pytest.mark.parametrize("case",list("ABCDEFGHIJKLMNOP"))
def test_synthetic_topology_matrix(case):
    if case in "AB":
        _,labels,regions,_,_=analyse()
        assert len(regions)==0
    elif case in "CDLMO":
        x=1 if case=="M" else 60
        shift=1 if case=="L" else 0
        _,labels,regions,_,_=analyse(tubes=[(x,16,1,shift)])
        assert len(regions)==1
        assert regions[0].attributes["upper_access"] and regions[0].attributes["lower_access"]
        assert regions[0].classification=="protein_pore"
    elif case=="E":
        _,_,regions,_,_=analyse(tubes=[(20,8,1,0),(70,23,1,0)])
        assert len(regions)==2
        assert all(r.classification=="protein_pore" for r in regions)
    elif case=="F":
        _,_,regions,_,_=analyse(tubes=[(20,9,1,0)],defect=[(70,23,2,0)])
        assert len(regions)==2
    elif case=="G":
        _,_,regions,_,_=analyse(tubes=[(20,7,1,0),(60,16,1,0)],defect=[(100,25,2,0)])
        assert len(regions)==3
    elif case=="H":
        field,mem=fixture()
        for z in range(14,30):
            field.nominal_free[60,16,z]=True;field.nearest_species[60,16,z]=1
        labels,regions,_,_=discover_regions(field,mem);local_classify(regions,mem,min_wall_samples=1)
        assert len(regions)==1 and regions[0].attributes["upper_access"] != regions[0].attributes["lower_access"]
        assert regions[0].classification=="vestibule"
    elif case=="I":
        _,_,regions,_,_=analyse(cavities=[(60,16,15)])
        assert len(regions)==1 and regions[0].classification=="internal_cavity"
    elif case=="J":
        field,labels,regions,_,graph=analyse(tubes=[(30,12,1,0),(80,20,1,0)])
        paths=build_pathways(regions,graph).networks
        assert len(regions)==2 and len(paths)==2
        upper=[p for p in graph.portals if 0 in p.regions]
        assert len(upper)==2 and upper[0].regions[1]==upper[1].regions[1]
    elif case=="K":
        _,_,regions,_,_=analyse(tubes=[(60,16,1,0)],curvature=True)
        assert len(regions)==1 and regions[0].attributes["upper_access"] and regions[0].attributes["lower_access"]
    elif case=="N":
        _,labels,regions,_,_=analyse(tubes=[(0,16,1,0)])
        assert len(regions)==1
    elif case=="P":
        _,_,regions,_,_=analyse(cavities=[(60,16,12),(60,16,18)])
        assert len(regions)==2
        assert not any(r.attributes["upper_access"] and r.attributes["lower_access"] for r in regions)


def test_two_channels_five_nm_apart_are_independent():
    _,_,regions,_,_=analyse(tubes=[(20,8,1,0),(70,8,1,0)])
    assert len(regions)==2 and len({r.id for r in regions})==2


def test_periodic_translation_preserves_component_topology():
    field,mem=fixture(tubes=[(0,12,1,0),(90,20,1,0)])
    first=discover_regions(field,mem)[1]
    shifted=SimpleNamespace(**vars(field));shifted.nominal_free=np.roll(field.nominal_free,7,axis=0)
    shifted.clearance=np.roll(field.clearance,7,axis=0);shifted.nearest_species=np.roll(field.nearest_species,7,axis=0)
    second=discover_regions(shifted,mem)[1]
    assert sorted(r.voxel_count for r in first)==sorted(r.voxel_count for r in second)

