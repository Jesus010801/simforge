"""Periodic leaflet patches in an intrinsic cell frame.

This backend represents single-valued, gently curved leaflets along a detected
cell axis. Unsupported/multiply-valued patches are marked uncertain. No fixed
Cartesian Z or synthetic fallback slab authorizes removal.
"""
import numpy as np
from scipy.spatial import cKDTree
from .identity import LIPID
from .models import MembraneAtlas


def build_membrane(geometry, lattice):
    lipid = np.flatnonzero(geometry.species == LIPID)
    # One reference site per lipid instance; avoid mixing P/N/O heights.
    # One pass over lipid atoms avoids rescanning the complete lipid array for
    # each molecule (quadratic for large bilayers). Preserve explicit site
    # preference while grouping by the independent molecule-instance ID.
    preference = {name: rank for rank, name in enumerate(("P", "P1", "O33", "N", "O11"))}
    chosen = {}
    for i in lipid:
        name = geometry.atom_names[i]
        rank = preference.get(name)
        if rank is None:
            continue
        mid = int(geometry.molecule_ids[i])
        previous = chosen.get(mid)
        if previous is None or rank < previous[0]:
            chosen[mid] = (rank, int(i))
    heads = [x[1] for x in chosen.values()]
    source = geometry.cell.fractional(geometry.coordinates[heads]) % 1 if heads else np.empty((0,3))
    concentration = np.abs(np.exp(2j*np.pi*source).mean(axis=0)) if len(source) else np.array([0,0,1])
    axis = int(np.argmax(concentration))
    tangents = tuple(i for i in range(3) if i != axis)
    shape = tuple(lattice.shape[i] for i in tangents)
    center = float(np.angle(np.exp(2j*np.pi*source[:,axis]).mean()) / (2*np.pi)) % 1 if len(source) else .5
    normal = geometry.cell.basis[:,axis]
    if len(source) < 6:
        return MembraneAtlas(axis, tangents, center, np.zeros(shape), np.zeros(shape), np.zeros(shape,bool), normal,
                             ["Insufficient leaflet references; membrane interpretation unresolved"])
    height = ((source[:,axis] - center + .5) % 1 - .5) * geometry.cell.lengths[axis]
    # Split relative to circular bilayer center, not Cartesian Z.
    surfaces, supports = [], []
    grids = np.meshgrid(*[(np.arange(lattice.shape[i])+lattice.phase[i])*lattice.steps[i] for i in tangents], indexing="ij")
    query = np.stack(grids, axis=-1).reshape(-1,2) % geometry.cell.lengths[list(tangents)]
    for side in (-1, 1):
        selected = height * side > 0
        if selected.sum() < 3:
            surfaces.append(np.zeros(shape)); supports.append(np.zeros(shape,bool)); continue
        tree = cKDTree(source[selected][:,tangents]*geometry.cell.lengths[list(tangents)], boxsize=geometry.cell.lengths[list(tangents)])
        distance, near = tree.query(query, k=min(6, int(selected.sum())))
        weights = 1 / np.maximum(distance, .15)**2
        values = height[selected][near]
        avg = (values*weights).sum(axis=1) / weights.sum(axis=1)
        variance = ((values-avg[:,None])**2*weights).sum(axis=1) / weights.sum(axis=1)
        surfaces.append(avg.reshape(shape))
        supports.append(((distance[:,0] < 5.0) & (variance < .5**2)).reshape(shape))
    lower, upper = surfaces
    support = supports[0] & supports[1] & (upper-lower > .8) & (upper-lower < 7.)
    for surf in surfaces:
        slopes = [(np.roll(surf,-1,a)-np.roll(surf,1,a))/(2*lattice.steps[tangents[a]]) for a in range(2)]
        support &= np.hypot(*slopes) < 1.0
    warnings = []
    if not support.all():
        warnings.append(f"{1-support.mean():.1%} of leaflet patches unresolved; no defect removal there")
    return MembraneAtlas(axis, tangents, center, lower, upper, support, normal, warnings)


def membrane_fields(membrane, lattice, interface_margin=.25):
    axis = membrane.normal_axis
    q = (np.arange(lattice.shape[axis])+lattice.phase[axis])/lattice.shape[axis]
    height = ((q-membrane.center_fraction+.5)%1-.5)*lattice.cell.lengths[axis]
    view = [1,1,1]; view[axis] = len(height)
    h = height.reshape(view)
    low = np.expand_dims(membrane.lower, axis)
    high = np.expand_dims(membrane.upper, axis)
    supported = np.broadcast_to(np.expand_dims(membrane.supported,axis),lattice.shape)
    above = np.broadcast_to(h > high + interface_margin, lattice.shape)
    below = np.broadcast_to(h < low - interface_margin, lattice.shape)
    # A headgroup-to-headgroup band is not all hydrocarbon: inset it.
    core = np.broadcast_to((h > low+.35)&(h < high-.35),lattice.shape)
    return above, below, core, supported
