"""Single v2 scientific source of truth."""
from time import perf_counter
from .identity import load_system
from .models import ExclusionModel, RegionHierarchy, Atlas
from .clearance import build_field
from .membrane_atlas import build_membrane
from .regions import discover_regions, local_classify
from .region_graph import build_graph
from .pathways import build_pathways


def build_atlas(path, spacing=.15, phase=(.5, .5, .5), exclusion=None, max_voxels=30_000_000):
    started = perf_counter(); system = load_system(path); geometry = system.geometry
    exclusion = exclusion or ExclusionModel()
    field = build_field(geometry, exclusion, spacing, phase, max_voxels)
    membrane = build_membrane(geometry, field.lattice)
    labels, regions, components, sides = discover_regions(field, membrane)
    local_classify(regions, membrane)
    graph = build_graph(labels, regions, field, sides, membrane.normal_axis)
    pathways = build_pathways(regions, graph)
    hierarchy = RegionHierarchy(components, [{"id": r.id, "members": [r.id], "kind": "elementary_region"} for r in regions],
                                [exclusion.probe_radius_nm], [])
    warnings = geometry.warnings + membrane.warnings
    atlas = Atlas(geometry, exclusion, field, membrane, labels, regions, graph, hierarchy, pathways, warnings,
                  timings={"geometry_seconds": perf_counter()-started}, components=components)
    return system, atlas

