"""Portals between classified regions and side reservoirs."""
from .models import Portal,RegionGraph
import numpy as np


def build_graph(labels,regions,field,side_masks,normal_axis=2):
    upper_reached,lower_reached=side_masks[:2]
    upper_face=side_masks[4] if len(side_masks)>4 else upper_reached
    lower_face=side_masks[5] if len(side_masks)>5 else lower_reached
    axis=field.lattice.cell
    steps=field.lattice.steps
    portals=[];nodes=[0,-1]+[r.id for r in regions]
    normal=normal_axis
    # The atlas carries per-region reservoir connectivity from full-space flood.
    for r in regions:
        mask=labels==r.id
        for side_id,opened,face in ((0,r.attributes["upper_access"],upper_face),
                                    (-1,r.attributes["lower_access"],lower_face)):
            if not opened:continue
            count=int((mask&face).sum())
            # Report a one-face portal at minimum when a path reaches the
            # reservoir but the discrete surface strip is thinner than a cell.
            count=max(1,count)
            tangent_area=float(np.prod(np.delete(steps,normal)))
            portals.append(Portal(len(portals)+1,(r.id,side_id),count,
                count*tangent_area,float(field.clearance[mask].max()),True,[(0,0,0)]))
    # Elementary regions touching in free space share explicit portals. This is
    # normally used after watershed subdivision, not component discovery.
    for i,a in enumerate(regions):
        ma=labels==a.id
        for b in regions[i+1:]:
            mb=labels==b.id
            contact=sum(int((ma&np.roll(mb,step,axis)).sum()) for axis in range(3) for step in (-1,1))
            if contact:
                portals.append(Portal(len(portals)+1,(a.id,b.id),contact,
                    contact*float(np.prod(np.delete(steps,normal))),
                    float(min(field.clearance[ma].max(),field.clearance[mb].max())),True,[(0,0,0)]))
    return RegionGraph(nodes,portals)

