# descriptors/aromaticity.py
"""
Descriptor de aromaticidad molecular.

Clasifica cada anillo detectado por topology.py según su tipo aromático,
y enriquece el TopologyDescriptor con esa información.

Jerarquía de clasificación (de más a menos confiable sin RDKit):
    1. Tipo 4 explícito en bonds (SDF con aromaticidad anotada)
    2. Geometría: planitud del anillo (RMS desviación < umbral)
    3. Composición: solo C/N/O/S + tamaño 5-6 + todos los átomos en anillos fusionados

RDKit futuro: reemplazar _classify_ring() → percepción SSSR exacta + Hückel.

Interfaz pública:
    compute_aromaticity(
        atoms:    list[dict],
        bonds:    list[dict],
        topology: TopologyDescriptor,
    ) -> AromaticityDescriptor
"""

from __future__ import annotations
import math
from pydantic import BaseModel
from descriptors.topology import TopologyDescriptor, RingInfo


# ─── Contrato de salida ───────────────────────────────────────────────────────

class RingAromaticityInfo(BaseModel):
    atom_indices:    list[int]
    size:            int
    is_aromatic:     bool
    is_fused:        bool
    classification:  str    # "aromatic" | "antiaromatic" | "nonaromatic" | "unknown"
    evidence:        str    # "bond_type_4" | "planar_geometry" | "composition_heuristic" | "none"
    planarity_rms:   float = -1.0   # Å, -1 si no calculado


class AromaticityDescriptor(BaseModel):
    # Por anillo
    rings:               list[RingAromaticityInfo] = []

    # Resumen
    n_aromatic_rings:    int        = 0
    n_fused_aromatic:    int        = 0   # anillos aromáticos fusionados entre sí
    aromatic_atom_count: int        = 0
    aromatic_atom_set:   list[int]  = []  # índices (1-based)

    # Sistemas conjugados fusionados (ej: naftaleno = 1 sistema de 2 anillos)
    fused_aromatic_systems: list[list[int]] = []   # lista de grupos de ring_idx


# ─── Planitud ────────────────────────────────────────────────────────────────

def _planarity_rms(atom_indices: list[int], atoms: list[dict]) -> float:
    """
    Calcula la desviación RMS de planitud para un conjunto de átomos.
    Retorna -1.0 si hay menos de 3 átomos o error numérico.
    """
    coords = []
    atom_map = {a["index"]: a for a in atoms}
    for idx in atom_indices:
        a = atom_map.get(idx)
        if a:
            coords.append((a["x"], a["y"], a["z"]))

    if len(coords) < 3:
        return -1.0

    n  = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n

    # Normal del plano por producto cruzado de los dos primeros vectores
    v1 = (coords[1][0]-coords[0][0], coords[1][1]-coords[0][1], coords[1][2]-coords[0][2])
    v2 = (coords[2][0]-coords[0][0], coords[2][1]-coords[0][1], coords[2][2]-coords[0][2])

    normal = (
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0],
    )
    norm_len = math.sqrt(sum(x**2 for x in normal))
    if norm_len < 1e-8:
        return 0.0   # puntos colineales → plano degenerado

    normal = tuple(x / norm_len for x in normal)

    deviations = []
    for c in coords:
        d = (c[0]-cx, c[1]-cy, c[2]-cz)
        dist = abs(sum(d[i]*normal[i] for i in range(3)))
        deviations.append(dist)

    return math.sqrt(sum(d**2 for d in deviations) / len(deviations))


# ─── Clasificación de aromaticidad por anillo ─────────────────────────────────

_AROMATIC_ELEMENTS = {"C", "N", "O", "S"}

# Umbrales de planitud por tamaño de anillo
# Los anillos de 6 son más rígidos; los de 5 (furano, pirrolo) algo más flexibles
_PLANARITY_THRESHOLD = {
    5: 0.10,   # Å RMS
    6: 0.08,
    7: 0.12,   # tropilium y similares
}
_PLANARITY_DEFAULT = 0.10


def _has_aromatic_bond_type(
    atom_indices: list[int],
    bonds:        list[dict],
) -> bool:
    """Retorna True si hay al menos un enlace tipo 4 dentro del anillo."""
    ring_set = set(atom_indices)
    for b in bonds:
        if b["bond_type"] == 4 and b["atom1"] in ring_set and b["atom2"] in ring_set:
            return True
    return False


def _composition_is_aromatic_candidate(
    atom_indices: list[int],
    atoms:        list[dict],
) -> bool:
    """
    Heurística de composición: ¿podría ser aromático?
    Condiciones: tamaño 5-6, solo elementos {C,N,O,S}, sin H en el anillo.
    """
    if len(atom_indices) not in (5, 6):
        return False
    atom_map = {a["index"]: a for a in atoms}
    for idx in atom_indices:
        a = atom_map.get(idx)
        if not a:
            return False
        if a["element"] not in _AROMATIC_ELEMENTS:
            return False
    return True


def _classify_ring(
    ring:   RingInfo,
    atoms:  list[dict],
    bonds:  list[dict],
) -> RingAromaticityInfo:
    """
    Clasifica un anillo como aromático o no.

    Prioridad de evidencia:
      1. Bond type 4 explícito → más confiable
      2. Planitud geométrica + composición → segundo nivel
      3. Solo composición (sin coordenadas suficientes) → heurística débil
    """
    idx_list = ring.atom_indices

    # ── Evidencia 1: bond_type == 4 ──────────────────────────────────────────
    if _has_aromatic_bond_type(idx_list, bonds):
        rms = _planarity_rms(idx_list, atoms)
        return RingAromaticityInfo(
            atom_indices   = idx_list,
            size           = ring.size,
            is_aromatic    = True,
            is_fused       = ring.is_fused,
            classification = "aromatic",
            evidence       = "bond_type_4",
            planarity_rms  = rms,
        )

    # ── Evidencia 2: planitud geométrica ─────────────────────────────────────
    rms = _planarity_rms(idx_list, atoms)
    threshold = _PLANARITY_THRESHOLD.get(ring.size, _PLANARITY_DEFAULT)

    if rms >= 0 and rms <= threshold:
        if _composition_is_aromatic_candidate(idx_list, atoms):
            return RingAromaticityInfo(
                atom_indices   = idx_list,
                size           = ring.size,
                is_aromatic    = True,
                is_fused       = ring.is_fused,
                classification = "aromatic",
                evidence       = "planar_geometry",
                planarity_rms  = rms,
            )

    # ── Evidencia 3: composición sola (sin coordenadas o no planar) ──────────
    if _composition_is_aromatic_candidate(idx_list, atoms) and ring.is_fused:
        # En sistemas fusionados la planitud global puede superar el umbral local
        # pero el anillo sigue siendo aromático (ej: coroneno, acridina)
        return RingAromaticityInfo(
            atom_indices   = idx_list,
            size           = ring.size,
            is_aromatic    = True,
            is_fused       = ring.is_fused,
            classification = "aromatic",
            evidence       = "composition_heuristic",
            planarity_rms  = rms,
        )

    return RingAromaticityInfo(
        atom_indices   = idx_list,
        size           = ring.size,
        is_aromatic    = False,
        is_fused       = ring.is_fused,
        classification = "nonaromatic",
        evidence       = "none",
        planarity_rms  = rms,
    )


# ─── Sistemas fusionados aromáticos ──────────────────────────────────────────

def _find_fused_aromatic_systems(
    aromatic_rings: list[tuple[int, RingAromaticityInfo]],
) -> list[list[int]]:
    """
    Agrupa anillos aromáticos que comparten al menos 2 átomos en sistemas fusionados.
    Retorna lista de grupos, cada grupo es una lista de índices de ring en aromatic_rings.

    Ejemplo: naftaleno → [[0, 1]], antraceno → [[0, 1, 2]]
    """
    if not aromatic_rings:
        return []

    n = len(aromatic_rings)
    # Union-Find simple
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for i in range(n):
        for j in range(i + 1, n):
            set_i = set(aromatic_rings[i][1].atom_indices)
            set_j = set(aromatic_rings[j][1].atom_indices)
            if len(set_i & set_j) >= 2:
                union(i, j)

    groups: dict[int, list[int]] = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(i)

    return list(groups.values())


# ─── Interfaz pública ────────────────────────────────────────────────────────

def compute_aromaticity(
    atoms:    list[dict],
    bonds:    list[dict],
    topology: TopologyDescriptor,
) -> AromaticityDescriptor:
    """
    Clasifica la aromaticidad de cada anillo detectado por topology.

    Enriquece los RingInfo del TopologyDescriptor con is_aromatic=True
    donde corresponda (mutación in-place de los RingInfo).
    """
    ring_classifications: list[RingAromaticityInfo] = []
    aromatic_atom_set: set[int] = set()

    for ring in topology.rings:
        info = _classify_ring(ring, atoms, bonds)
        ring_classifications.append(info)
        # Propagar is_aromatic al RingInfo del topology
        ring.is_aromatic = info.is_aromatic
        if info.is_aromatic:
            for idx in info.atom_indices:
                aromatic_atom_set.add(idx)

    # Sistemas fusionados aromáticos
    aromatic_indexed = [
        (i, rc) for i, rc in enumerate(ring_classifications) if rc.is_aromatic
    ]
    fused_systems = _find_fused_aromatic_systems(aromatic_indexed)

    n_aromatic = sum(1 for rc in ring_classifications if rc.is_aromatic)
    n_fused_aromatic = sum(
        len(g) for g in fused_systems if len(g) > 1
    )

    return AromaticityDescriptor(
        rings                = ring_classifications,
        n_aromatic_rings     = n_aromatic,
        n_fused_aromatic     = n_fused_aromatic,
        aromatic_atom_count  = len(aromatic_atom_set),
        aromatic_atom_set    = sorted(aromatic_atom_set),
        fused_aromatic_systems = fused_systems,
    )
