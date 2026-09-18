"""Portal graph network summaries; no single best corridor defines discovery."""
from collections import defaultdict
from .models import PathwayNetwork


def build_pathways(regions, graph):
    adjacency = defaultdict(list)
    for portal in graph.portals:
        a, b = portal.regions
        if a in (0, -1) or b in (0, -1):
            continue
        adjacency[a].append((b, portal.id)); adjacency[b].append((a, portal.id))
    networks, seen = [], set()
    for region in regions:
        if region.id in seen:
            continue
        stack, members = [region.id], set()
        while stack:
            u = stack.pop()
            if u in members:
                continue
            members.add(u); seen.add(u)
            stack.extend(v for v, _ in adjacency[u] if v not in members)
        upper = any((p.regions[0] in members and p.regions[1] == 0) or
                    (p.regions[1] in members and p.regions[0] == 0) for p in graph.portals)
        lower = any((p.regions[0] in members and p.regions[1] == -1) or
                    (p.regions[1] in members and p.regions[0] == -1) for p in graph.portals)
        networks.append({"region_ids": sorted(members), "has_upper_portal": upper, "has_lower_portal": lower,
                         "shared_junctions": sorted(x for x in members if len(adjacency[x]) > 2),
                         "pathways": "represented by the full region and portal graph"})
    return PathwayNetwork(networks)

