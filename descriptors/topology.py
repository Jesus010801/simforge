# descriptors/topology.py
"""
Descriptor de topología molecular.

Es la base del descriptor engine: construye el grafo molecular y detecta
anillos. Todos los demás descriptors dependen de este.

Interfaz pública:
    compute_topology(atoms: list[dict], bonds: list[dict]) -> TopologyDescriptor

El grafo interno usa índices de átomo (1-based, como vienen de los parsers).
No depende de coordenadas — solo de la tabla de átomos y enlaces.
"""

from __future__ import annotations
from collections import defaultdict
from pydantic import BaseModel


# ─── Contrato de salida ───────────────────────────────────────────────────────

class RingInfo(BaseModel):
    atom_indices: list[int]    # índices (1-based) de los átomos del anillo
    size:         int          # 5, 6, 7...
    is_aromatic:  bool = False # se rellena por aromaticity.py
    is_fused:     bool = False # comparte ≥2 átomos con otro anillo


class TopologyDescriptor(BaseModel):
    # Conteo de átomos pesados (no-H)
    n_heavy_atoms:     int          = 0
    n_hydrogens:       int          = 0
    n_bonds:           int          = 0

    # Elementos presentes
    elements:          list[str]    = []   # únicos, ordenados
    element_counts:    dict         = {}   # {"C": 12, "O": 3, ...}

    # Vecinos por átomo (grafo)
    adjacency:         dict         = {}   # {atom_idx: [neighbor_idx, ...]}

    # Anillos detectados (SSSR — smallest set of smallest rings)
    rings:             list[RingInfo] = []
    n_rings:           int            = 0
    ring_atom_indices: set            = set()   # set flat de todos los átomos en anillo

    # Clasificación de átomos
    terminal_atoms:    list[int]    = []   # átomos con un solo vecino pesado
    branching_atoms:   list[int]    = []   # átomos con ≥3 vecinos pesados

    class Config:
        arbitrary_types_allowed = True


# ─── Construcción del grafo ───────────────────────────────────────────────────

def _build_adjacency(
    atoms: list[dict],
    bonds: list[dict],
) -> dict[int, list[int]]:
    """
    Construye grafo de adyacencia.
    Solo átomos no-H en los vecinos (para ring detection y flexibilidad).
    Los H sí están como nodos pero no como vecinos en este grafo.
    """
    h_indices: set[int] = {a["index"] for a in atoms if a["element"] == "H"}

    adj: dict[int, list[int]] = defaultdict(list)
    for a in atoms:
        adj[a["index"]]  # asegurar que todos los nodos existen

    for b in bonds:
        i, j = b["atom1"], b["atom2"]
        # Grafo pesado: solo enlaces entre átomos no-H
        if i not in h_indices and j not in h_indices:
            adj[i].append(j)
            adj[j].append(i)

    return dict(adj)


# ─── Ring detection (SSSR simplificado) ──────────────────────────────────────

def _find_rings(
    adj:       dict[int, list[int]],
    max_size:  int = 8,
) -> list[list[int]]:
    """
    Encuentra todos los anillos simples de tamaño ≤ max_size.

    Algoritmo: DFS con backtracking desde cada átomo.
    Retorna lista de ciclos únicos (normalizados como frozenset para dedup).
    Solo considera átomos en el grafo pesado.

    No es un SSSR completo (para eso necesitaríamos álgebra lineal sobre
    el espacio de ciclos), pero es suficiente para anillos de 5-6 en
    moléculas pequeñas típicas de ligandos.
    """
    found:   list[frozenset] = []
    visited: set[frozenset]  = set()

    heavy_atoms = [idx for idx, neighbors in adj.items() if neighbors]

    def dfs(start: int, current: int, path: list[int], depth: int):
        if depth > max_size:
            return
        for neighbor in adj.get(current, []):
            if neighbor == start and depth >= 3:
                # Ciclo cerrado
                key = frozenset(path)
                if key not in visited:
                    visited.add(key)
                    found.append(key)
            elif neighbor not in path:
                dfs(start, neighbor, path + [neighbor], depth + 1)

    for atom in heavy_atoms:
        dfs(atom, atom, [atom], 1)

    # Convertir de frozenset a list (con orden consistente)
    result = []
    for ring_set in found:
        result.append(sorted(ring_set))

    # Deduplicar y ordenar por tamaño
    seen = set()
    unique = []
    for ring in result:
        key = frozenset(ring)
        if key not in seen:
            seen.add(key)
            unique.append(ring)

    unique.sort(key=len)
    return unique


def _detect_fused_rings(rings: list[list[int]]) -> list[bool]:
    """
    Marca cada anillo como fusionado si comparte ≥2 átomos con otro anillo.
    (Compartir 1 átomo = espirocíclico, no fusionado.)
    """
    fused = [False] * len(rings)
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            shared = set(rings[i]) & set(rings[j])
            if len(shared) >= 2:
                fused[i] = True
                fused[j] = True
    return fused


# ─── Interfaz pública ────────────────────────────────────────────────────────

def compute_topology(
    atoms: list[dict],
    bonds: list[dict],
) -> TopologyDescriptor:
    """
    Calcula el descriptor de topología molecular.

    atoms: lista de dicts con keys: index, element, x, y, z, charge
    bonds: lista de dicts con keys: atom1, atom2, bond_type
    """
    # ─── Conteos básicos ─────────────────────────────────────────────────────
    elements_raw  = [a["element"] for a in atoms]
    element_count: dict[str, int] = defaultdict(int)
    for el in elements_raw:
        element_count[el] += 1

    n_hydrogens  = element_count.get("H", 0)
    n_heavy      = len(atoms) - n_hydrogens
    elements_uniq = sorted(set(elements_raw))

    # ─── Grafo pesado ────────────────────────────────────────────────────────
    adj = _build_adjacency(atoms, bonds)

    # ─── Anillos ─────────────────────────────────────────────────────────────
    raw_rings = _find_rings(adj, max_size=8)
    fused_flags = _detect_fused_rings(raw_rings)

    ring_objects: list[RingInfo] = []
    ring_atom_indices: set[int]  = set()

    for ring_atoms, is_fused in zip(raw_rings, fused_flags):
        ring_objects.append(RingInfo(
            atom_indices = ring_atoms,
            size         = len(ring_atoms),
            is_fused     = is_fused,
        ))
        for idx in ring_atoms:
            ring_atom_indices.add(idx)

    # ─── Clasificación de átomos ──────────────────────────────────────────────
    terminal_atoms  = []
    branching_atoms = []

    for atom in atoms:
        if atom["element"] == "H":
            continue
        idx      = atom["index"]
        n_heavy_neighbors = len([
            nb for nb in adj.get(idx, [])
        ])   # adj ya es solo pesados
        if n_heavy_neighbors == 1:
            terminal_atoms.append(idx)
        elif n_heavy_neighbors >= 3:
            branching_atoms.append(idx)

    return TopologyDescriptor(
        n_heavy_atoms    = n_heavy,
        n_hydrogens      = n_hydrogens,
        n_bonds          = len(bonds),
        elements         = elements_uniq,
        element_counts   = dict(element_count),
        adjacency        = {k: list(v) for k, v in adj.items()},
        rings            = ring_objects,
        n_rings          = len(ring_objects),
        ring_atom_indices = ring_atom_indices,
        terminal_atoms   = terminal_atoms,
        branching_atoms  = branching_atoms,
    )
