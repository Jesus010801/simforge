# descriptors/geometry.py
"""
Descriptor de geometría molecular 3D.

Trabaja con coordenadas explícitas de los átomos para inferir propiedades
espaciales que no son accesibles desde conectividad sola.

Depende de:
    topology.py    — grafo molecular, anillos
    aromaticity.py — qué anillos son aromáticos (para planitud por sistema)

Interfaz pública:
    compute_geometry(
        atoms:       list[dict],
        bonds:       list[dict],
        topology:    TopologyDescriptor,
        aromaticity: AromaticityDescriptor,
    ) -> GeometryDescriptor
"""

from __future__ import annotations
import math
from pydantic import BaseModel
from descriptors.topology    import TopologyDescriptor
from descriptors.aromaticity import AromaticityDescriptor


# ─── Contrato de salida ───────────────────────────────────────────────────────

class BoundingBox(BaseModel):
    x_min: float; x_max: float
    y_min: float; y_max: float
    z_min: float; z_max: float
    x_len: float; y_len: float; z_len: float


class ChiralCenter(BaseModel):
    atom_index:  int
    element:     str
    n_neighbors: int


class GeometryDescriptor(BaseModel):
    # Extensión espacial
    max_distance:      float = 0.0   # Å entre los dos átomos pesados más lejanos
    bounding_box:      BoundingBox | None = None
    radius_of_gyration: float = 0.0  # Å desde el centroide

    # Shape
    shape_class:       str = "unknown"   # flat / elongated / globular
    anisotropy:        float = 0.0       # 0=esférico, 1=totalmente lineal

    # Planaridad global
    global_planarity_rms: float = 0.0   # Å RMS de todos los átomos pesados al plano
    is_globally_planar:   bool  = False  # RMS < 0.5Å

    # Planaridad por sistema aromático fusionado
    aromatic_system_planarity: list[float] = []   # RMS por sistema

    # Centros quirales (heurística: C sp3 con 4 vecinos pesados distintos)
    chiral_centers:    list[ChiralCenter] = []
    n_chiral_centers:  int = 0

    # Compactación
    n_heavy_atoms:     int   = 0
    volume_estimate:   float = 0.0   # Å³ aproximado (esfera de radio Rg)


# ─── Utilidades geométricas ───────────────────────────────────────────────────

def _centroid(coords: list[tuple]) -> tuple:
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


def _distance_3d(a: tuple, b: tuple) -> float:
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))


def _planarity_rms(coords: list[tuple]) -> float:
    """RMS de desviación al plano de mínimos cuadrados (producto cruzado)."""
    if len(coords) < 3:
        return 0.0
    cx, cy, cz = _centroid(coords)
    v1 = (coords[1][0]-coords[0][0], coords[1][1]-coords[0][1], coords[1][2]-coords[0][2])
    v2 = (coords[2][0]-coords[0][0], coords[2][1]-coords[0][1], coords[2][2]-coords[0][2])
    normal = (
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0],
    )
    norm_len = math.sqrt(sum(x**2 for x in normal))
    if norm_len < 1e-8:
        return 0.0
    normal = tuple(x / norm_len for x in normal)
    deviations = [
        abs((c[0]-cx)*normal[0] + (c[1]-cy)*normal[1] + (c[2]-cz)*normal[2])
        for c in coords
    ]
    return math.sqrt(sum(d**2 for d in deviations) / len(deviations))


# ─── Extensión espacial ───────────────────────────────────────────────────────

def _compute_max_distance(heavy_coords: list[tuple]) -> float:
    """Distancia máxima entre cualquier par de átomos pesados."""
    max_d = 0.0
    n = len(heavy_coords)
    for i in range(n):
        for j in range(i+1, n):
            d = _distance_3d(heavy_coords[i], heavy_coords[j])
            if d > max_d:
                max_d = d
    return round(max_d, 3)


def _compute_bounding_box(heavy_coords: list[tuple]) -> BoundingBox:
    xs = [c[0] for c in heavy_coords]
    ys = [c[1] for c in heavy_coords]
    zs = [c[2] for c in heavy_coords]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)
    return BoundingBox(
        x_min=round(x_min,3), x_max=round(x_max,3),
        y_min=round(y_min,3), y_max=round(y_max,3),
        z_min=round(z_min,3), z_max=round(z_max,3),
        x_len=round(x_max-x_min,3),
        y_len=round(y_max-y_min,3),
        z_len=round(z_max-z_min,3),
    )


def _compute_radius_of_gyration(heavy_coords: list[tuple]) -> float:
    """Rg = sqrt( mean( |r_i - centroid|^2 ) )"""
    if not heavy_coords:
        return 0.0
    c = _centroid(heavy_coords)
    return round(math.sqrt(
        sum(_distance_3d(coord, c)**2 for coord in heavy_coords) / len(heavy_coords)
    ), 3)


# ─── Shape ───────────────────────────────────────────────────────────────────

def _compute_shape(bb: BoundingBox) -> tuple[str, float]:
    """
    Clasifica la forma molecular usando el bounding box.

    Anisotropía = 1 - (dim_min / dim_max)
        0   → esfera perfecta (globular)
        1   → vara infinita (elongated)

    Shape:
        flat      → una dimensión << otras dos (como xantona)
        elongated → una dimensión >> otras dos (como cadena larga)
        globular  → las tres dimensiones similares
    """
    dims = sorted([bb.x_len, bb.y_len, bb.z_len])
    d_min, d_mid, d_max = dims

    anisotropy = round(1.0 - (d_min / d_max) if d_max > 0 else 0.0, 3)

    ratio_min_mid = d_min / d_mid if d_mid > 0 else 1.0
    ratio_mid_max = d_mid / d_max if d_max > 0 else 1.0

    if ratio_min_mid < 0.35:
        shape = "flat"        # una dimensión muy comprimida
    elif ratio_mid_max < 0.45:
        shape = "elongated"   # una dimensión muy extendida
    else:
        shape = "globular"

    return shape, anisotropy


# ─── Centros quirales ─────────────────────────────────────────────────────────

def _detect_chiral_centers(
    atoms:    list[dict],
    topology: TopologyDescriptor,
) -> list[ChiralCenter]:
    """
    Heurística: un carbono sp3 con exactamente 4 vecinos pesados distintos
    es candidato a centro quiral.

    Limitación: sin RDKit no podemos verificar si los 4 sustituyentes son
    realmente distintos (eso requiere CIP ranking). Esto detecta candidatos.
    """
    atom_map = {a["index"]: a for a in atoms}
    adj      = topology.adjacency
    centers  = []

    for atom in atoms:
        if atom["element"] != "C":
            continue
        idx = atom["index"]
        heavy_neighbors = adj.get(idx, [])
        if len(heavy_neighbors) != 4:
            continue
        # Verificar que los 4 vecinos tienen elementos distintos entre sí
        # o al menos que hay 4 vecinos (condición necesaria, no suficiente)
        neighbor_elements = [atom_map[nb]["element"] for nb in heavy_neighbors if nb in atom_map]
        # Si los 4 vecinos son todos iguales → no quiral (ej: CH4)
        if len(set(neighbor_elements)) < 2:
            continue
        centers.append(ChiralCenter(
            atom_index  = idx,
            element     = "C",
            n_neighbors = len(heavy_neighbors),
        ))

    return centers


# ─── Planaridad por sistema aromático ─────────────────────────────────────────

def _aromatic_system_planarity(
    atoms:       list[dict],
    aromaticity: AromaticityDescriptor,
) -> list[float]:
    """
    Calcula la planitud RMS para cada sistema aromático fusionado.
    Un sistema fusionado = grupo de anillos aromáticos que comparten átomos.
    """
    atom_map = {a["index"]: a for a in atoms}
    rms_values = []

    for system_ring_indices in aromaticity.fused_aromatic_systems:
        # Recolectar todos los átomos del sistema fusionado
        system_atom_set: set[int] = set()
        for ring_idx in system_ring_indices:
            if ring_idx < len(aromaticity.rings):
                for ai in aromaticity.rings[ring_idx].atom_indices:
                    system_atom_set.add(ai)

        coords = []
        for idx in system_atom_set:
            a = atom_map.get(idx)
            if a:
                coords.append((a["x"], a["y"], a["z"]))

        if len(coords) >= 3:
            rms_values.append(round(_planarity_rms(coords), 4))

    # Anillos aromáticos solitarios (no fusionados con otros)
    for i, ring_info in enumerate(aromaticity.rings):
        if not ring_info.is_aromatic:
            continue
        # Ver si este anillo ya está en un sistema fusionado
        in_fused = any(i in sys for sys in aromaticity.fused_aromatic_systems)
        if not in_fused:
            coords = []
            for idx in ring_info.atom_indices:
                a = atom_map.get(idx)
                if a:
                    coords.append((a["x"], a["y"], a["z"]))
            if len(coords) >= 3:
                rms_values.append(round(_planarity_rms(coords), 4))

    return rms_values


# ─── Interfaz pública ────────────────────────────────────────────────────────

def compute_geometry(
    atoms:       list[dict],
    bonds:       list[dict],
    topology:    TopologyDescriptor,
    aromaticity: AromaticityDescriptor,
) -> GeometryDescriptor:
    """
    Calcula el descriptor de geometría molecular 3D.
    """
    atom_map = {a["index"]: a for a in atoms}

    # Solo átomos pesados para la mayoría de cálculos
    heavy_atoms = [a for a in atoms if a["element"] != "H"]
    heavy_coords = [(a["x"], a["y"], a["z"]) for a in heavy_atoms]

    if not heavy_coords:
        return GeometryDescriptor()

    # ── Extensión ─────────────────────────────────────────────────────────────
    max_dist = _compute_max_distance(heavy_coords)
    bbox     = _compute_bounding_box(heavy_coords)
    rg       = _compute_radius_of_gyration(heavy_coords)

    # ── Shape ─────────────────────────────────────────────────────────────────
    shape, anisotropy = _compute_shape(bbox)

    # ── Planaridad global ─────────────────────────────────────────────────────
    global_rms      = round(_planarity_rms(heavy_coords), 4)
    is_planar       = global_rms < 0.5   # Å — umbral para molécula "plana"

    # ── Planaridad por sistema aromático ─────────────────────────────────────
    arom_planarity  = _aromatic_system_planarity(atoms, aromaticity)

    # ── Centros quirales ──────────────────────────────────────────────────────
    chirals         = _detect_chiral_centers(atoms, topology)

    # ── Volumen estimado (esfera de radio Rg) ────────────────────────────────
    volume = round((4/3) * math.pi * rg**3, 2)

    return GeometryDescriptor(
        max_distance             = max_dist,
        bounding_box             = bbox,
        radius_of_gyration       = rg,
        shape_class              = shape,
        anisotropy               = anisotropy,
        global_planarity_rms     = global_rms,
        is_globally_planar       = is_planar,
        aromatic_system_planarity = arom_planarity,
        chiral_centers           = chirals,
        n_chiral_centers         = len(chirals),
        n_heavy_atoms            = len(heavy_atoms),
        volume_estimate          = volume,
    )
