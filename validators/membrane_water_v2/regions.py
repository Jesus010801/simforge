"""Exhaustive membrane-interior components and reservoir reachability."""
import numpy as np
from scipy import ndimage
from .models import Region, VoidComponent
from .void_components import components as periodic_components

_CONN=ndimage.generate_binary_structure(3,1)


def discover_regions(field,membrane):
    """Find every core component and test real full-void paths to each reservoir."""
    upper_seed,lower_seed,core,supported,upper_face,lower_face=_membrane_fields(membrane,field.lattice)
    # Reservoir accessibility is flood connectivity in the complete free-space
    # geometry, with periodicity along the membrane but no wrap across its normal.
    periodic=tuple(i!=membrane.normal_axis for i in range(3))
    full_labels,nfull=periodic_components(field.nominal_free,periodic)
    upper_ids=set(np.unique(full_labels[upper_seed]).tolist())-{0}
    lower_ids=set(np.unique(full_labels[lower_seed]).tolist())-{0}
    upper_reached=np.isin(full_labels,list(upper_ids)) if upper_ids else np.zeros_like(field.nominal_free)
    lower_reached=np.isin(full_labels,list(lower_ids)) if lower_ids else np.zeros_like(field.nominal_free)

    candidate=field.nominal_free & core & supported
    labels,n=ndimage.label(candidate,structure=_CONN)
    parent=np.arange(n+1)
    def root(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for axis in range(3):
        a,b=np.take(labels,0,axis),np.take(labels,-1,axis)
        for u,v in np.unique(np.stack((a.ravel(),b.ravel()),axis=1),axis=0):
            if u and v:
                ru,rv=root(int(u)),root(int(v));parent[max(ru,rv)]=min(ru,rv)
    roots=sorted({root(int(i)) for i in np.unique(labels) if i})
    remap={r:i+1 for i,r in enumerate(roots)}
    lut=np.zeros(n+1,dtype=np.int32)
    for old in range(1,n+1):
        lut[old]=remap[root(old)]
    # One dense lookup remaps all voxels in O(grid size), independent of the
    # number of disconnected components.
    labels=lut[labels];regions=[];component_records=[]
    # Group labeled voxel indices once. A full-volume labels==cid scan for every
    # component is quadratic on fragmented large-system grids.
    flat_labels=labels.ravel();flat_voxels=np.flatnonzero(flat_labels)
    order=np.argsort(flat_labels[flat_voxels],kind="stable")
    sorted_voxels=flat_voxels[order];sorted_labels=flat_labels[sorted_voxels]
    starts=np.searchsorted(sorted_labels,np.arange(1,len(roots)+1),side="left")
    stops=np.searchsorted(sorted_labels,np.arange(1,len(roots)+1),side="right")
    upper_flat=upper_reached.ravel();lower_flat=lower_reached.ravel()
    clearance_flat=field.clearance.ravel();species_flat=field.nearest_species.ravel()
    cutoff=field.lattice.steps.max()*2+.14
    for cid,(begin,end) in enumerate(zip(starts,stops),1):
        voxels=sorted_voxels[begin:end]
        if not len(voxels):continue
        up=bool(upper_flat[voxels].any());down=bool(lower_flat[voxels].any())
        clear=clearance_flat[voxels]
        shell_species=species_flat[voxels[clear<=cutoff]]
        nearest=shell_species if len(shell_species) else species_flat[voxels]
        regions.append(Region(cid,cid,int(len(voxels)),attributes={
          "membrane_interior":True,"upper_access":up,"lower_access":down,
          "clearance_mean_nm":float(clear.mean()),"clearance_max_nm":float(clear.max()),
          "wall_species_counts":{str(int(sp)):int((nearest==sp).sum()) for sp in np.unique(nearest)},
        },evidence={"connectivity":"nominal full-void reservoir reachability",
                    "core_supported_fraction":1.0},voxel_indices=voxels))
        component_records.append(VoidComponent(cid,int(len(voxels)),[cid],up,down))
    return labels,regions,component_records,(upper_reached,lower_reached,core,supported,upper_face,lower_face)


def _membrane_fields(membrane,lattice):
    axis=membrane.normal_axis
    q=(np.arange(lattice.shape[axis])+lattice.phase[axis])/lattice.shape[axis]
    height=((q-membrane.center_fraction+.5)%1-.5)*lattice.cell.lengths[axis]
    shape=[1,1,1];shape[axis]=len(height);h=height.reshape(shape)
    lo,hi=np.expand_dims(membrane.lower,axis),np.expand_dims(membrane.upper,axis)
    supported=np.broadcast_to(np.expand_dims(membrane.supported,axis),lattice.shape)
    upper=np.broadcast_to(h>hi+.25,lattice.shape)
    lower=np.broadcast_to(h<lo-.25,lattice.shape)
    core=np.broadcast_to((h>lo+.35)&(h<hi-.35),lattice.shape)
    upper_face=np.broadcast_to((h>hi-.65)&(h<hi-.35),lattice.shape)
    lower_face=np.broadcast_to((h<lo+.65)&(h>lo+.35),lattice.shape)
    return upper,lower,core,supported,upper_face,lower_face


def local_classify(regions,membrane,min_wall_samples=8):
    for r in regions:
        up,down=r.attributes["upper_access"],r.attributes["lower_access"]
        walls=r.attributes["wall_species_counts"];protein,lipid=walls.get("1",0),walls.get("2",0)
        denom=protein+lipid;pfrac=protein/denom if denom else 0.0
        r.evidence.update({"protein_wall_fraction":pfrac,"wall_sample_count":denom})
        r.attributes["path_class"]="conduit" if up and down else "vestibule" if up or down else "enclosed"
        r.attributes["boundary_class"]="protein" if pfrac>=.65 else "lipid" if pfrac<=.25 else "mixed"
        if denom<min_wall_samples:
            r.classification="ambiguous";r.uncertainty.append("insufficient molecular boundary samples")
        elif up and down and pfrac>=.65:r.classification="protein_pore"
        elif (up or down) and pfrac>=.55:r.classification="vestibule"
        elif up and down and pfrac<=.25:r.classification="membrane_defect_candidate"
        elif not up and not down and pfrac>=.55:r.classification="internal_cavity"
        elif not up and not down and pfrac<=.25:r.classification="lipid_facing_void"
        else:r.classification="ambiguous";r.uncertainty.append("mixed or conflicting boundary evidence")
        r.evidence["classification_basis"]=["local clearance","core location","full-void reservoir connectivity","local wall composition"]
    return regions

