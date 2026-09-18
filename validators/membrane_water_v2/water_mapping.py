"""Water mapping runs strictly after geometry discovery and classification."""
import numpy as np
from .models import WaterMapping
from .clearance import ClearanceIndex


def map_waters(system,atlas,exclusion):
    geometry=system.geometry;index=ClearanceIndex(geometry)
    out=[];labels=atlas.labels;field=atlas.field;by_id={r.id:r for r in atlas.regions}
    mem=atlas.membrane;axis=mem.normal_axis
    for water in system.waters:
        atom=system.gro.atoms[water.oxygen_index]
        point=np.array([[atom.x,atom.y,atom.z]])
        _,hardmargin,_=index.query(point,exclusion)
        ijk=tuple(field.lattice.indices(point)[0]);rid=int(labels[ijk]);region=by_id.get(rid)
        fractional=geometry.cell.fractional(point)[0]%1
        along=((fractional[axis]-mem.center_fraction+.5)%1-.5)*geometry.cell.lengths[axis]
        tang=tuple(int(ijk[a]) for a in mem.tangent_axes)
        supported=bool(mem.supported[tang])
        local=tuple(fractional[a]*geometry.cell.lengths[a] for a in mem.tangent_axes)
        # Bilinear-free nearest patch query; atlas stores support and per-patch surfaces.
        lower=float(mem.lower[tang]);upper=float(mem.upper[tang])
        in_core=bool(supported and lower+.35<along<upper-.35)
        cls=region.classification if region else ("ambiguous" if in_core else "bulk")
        clash=bool(hardmargin[0]<0)
        if clash: cls="steric_clash"
        uncertainty=[] if water.valid else ["malformed water identity; preserve"]
        if not supported: uncertainty.append("local membrane surface unsupported")
        if field.uncertain is not None and field.uncertain[ijk]:
            uncertainty.append("water lies in uncertain geometric boundary")
        out.append(WaterMapping(water.uid,rid if region else None,cls,clash,in_core,uncertainty,
                                molecule_uid=water.uid,atom_indices=water.atom_indices))
    return out

