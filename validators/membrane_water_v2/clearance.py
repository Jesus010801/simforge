"""Exact periodic distance to unions of unequal atomic spheres."""
import numpy as np
from scipy.spatial import cKDTree
from .models import FreeSpaceField
from .cell import Lattice

class ClearanceIndex:
    def __init__(self,geometry):
        self.geometry=geometry;self.groups=[];geometry.cell.require_orthogonal()
        coords=geometry.cell.local(geometry.coordinates)
        for radius in np.unique(geometry.radii):
            indices=np.flatnonzero(geometry.radii==radius)
            self.groups.append((float(radius),indices,cKDTree(coords[indices],boxsize=geometry.cell.lengths)))

    def query(self,points,model):
        local=self.geometry.cell.local(points)
        clearance=np.full(len(points),np.inf);hard=np.full(len(points),np.inf)
        nearest=np.full(len(points),-1,dtype=np.int32)
        for radius,indices,tree in self.groups:
            d,j=tree.query(local,workers=1);candidate=d-radius
            better=candidate<clearance
            nearest[better]=indices[j[better]]
            clearance=np.minimum(clearance,candidate)
            hard=np.minimum(hard,d-model.hard_contact_scale*(radius+model.oxygen_radius_nm))
        return clearance,hard,nearest


def build_field(geometry,model,spacing=.15,phase=(.5,.5,.5),max_voxels=30_000_000):
    lattice=Lattice.create(geometry.cell,spacing,phase);size=int(np.prod(lattice.shape))
    if size>max_voxels:
        raise ValueError(f"Requested {size} cells exceeds resource limit; preserve water, refine explicitly")
    index=ClearanceIndex(geometry);clearance=np.empty(size,np.float32);hard=np.empty(size,bool)
    species=np.zeros(size,np.uint8)
    for start in range(0,size,100_000):
        stop=min(size,start+100_000)
        c,hm,nearest=index.query(lattice.points(np.arange(start,stop)),model)
        clearance[start:stop]=c;hard[start:stop]=hm < -lattice.error_bound
        valid=nearest>=0;species[start:stop][valid]=geometry.species[nearest[valid]]
    c=clearance.reshape(lattice.shape);error=lattice.error_bound+model.radius_uncertainty_nm
    return FreeSpaceField(lattice,c,species.reshape(lattice.shape),hard.reshape(lattice.shape),
        c>model.probe_radius_nm,c>model.probe_radius_nm+error,c>model.probe_radius_nm-error,
        np.abs(c-model.probe_radius_nm)<=error)

