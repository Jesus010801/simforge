# descriptors/flexibility.py
"""
Descriptor de flexibilidad conformacional.

Mejora sustancial sobre el conteo crudo de enlaces rotables en _analyze_chemistry():
    - Excluye enlaces dentro de anillos (no son rotables)
    - Excluye enlaces en sistemas aromáticos
    - Excluye enlaces terminales (=O, -OH en ácidos, -NH2) que tienen
      restricciones electrónicas o son demasiado cortos para afectar
      la búsqueda conformacional
    - Penaliza por número de anillos fusionados (rigidez de scaffold)
    - Calcula flexibilidad efectiva, no solo conteo crudo

La distinción rigidez/flexibilidad es cualitativa pero científicamente informada:
    rigid    : ≤ 3 rot. bonds  (exploración conformacional trivial)
    moderate : 4-7             (muestreo estándar suficiente)
    flexible : 8-14            (requiere muestreo aumentado o docking)
    very_flexible: ≥ 15        (requiere estrategia especial: REST2, metadinámica)

Interfaz pública:
    compute_flexibility(
        atoms:        list[dict],
        bonds:        list[dict],
        topology:     TopologyDescriptor,
        aromaticity:  AromaticityDescriptor,
    ) -> FlexibilityDescriptor
"""

from __future__ import annotations
from pydantic import BaseModel
from descriptors.topology     import TopologyDescriptor
from descriptors.aromaticity  import AromaticityDescriptor


# ─── Contrato de salida ───────────────────────────────────────────────────────

class RotatableBondInfo(BaseModel):
    atom1:       int
    atom2:       int
    bond_type:   int
    reason:      str   # por qué se considera rotable
    is_excluded: bool = False
    exclusion_reason: str = ""


class FlexibilityDescriptor(BaseModel):
    # Conteo
    n_rotatable_bonds:     int   = 0
    n_excluded_bonds:      int   = 0

    # Detalle
    rotatable_bonds:       list[RotatableBondInfo] = []
    excluded_bonds:        list[RotatableBondInfo] = []

    # Métricas de rigidez del scaffold
    n_rings:               int   = 0
    n_fused_aromatic_systems: int = 0
    scaffold_rigidity:     str   = "unknown"   # "rigid_core" | "flexible_core" | "mixed"

    # Clasificación final
    flexibility_class:     str   = "unknown"   # rigid / moderate / flexible / very_flexible
    flexibility_score:     float = 0.0         # métrica continua para futura comparación

    # Recomendación de muestreo
    sampling_recommendation: str = ""


# ─── Criterios de exclusión ───────────────────────────────────────────────────

def _is_in_ring(atom1: int, atom2: int, ring_atom_set: set[int]) -> bool:
    """Ambos átomos en anillos Y comparten un anillo (no basta estar ambos en anillos)."""
    return atom1 in ring_atom_set and atom2 in ring_atom_set


def _both_in_same_ring(
    atom1:    int,
    atom2:    int,
    topology: TopologyDescriptor,
) -> bool:
    for ring in topology.rings:
        ring_set = set(ring.atom_indices)
        if atom1 in ring_set and atom2 in ring_set:
            return True
    return False


def _is_terminal_bond(
    atom1:    int,
    atom2:    int,
    topology: TopologyDescriptor,
) -> tuple[bool, str]:
    """
    Un enlace es 'terminal' (no rotable conformacionalmente relevante) si:
      - Uno de los átomos tiene valencia 1 en el grafo pesado (=O, =S, -F)
      - Es un enlace C=O de carbonilo (no rota)
      - Es un enlace C-O de éster/ácido si el O es terminal (solo 1 vecino pesado)

    Retorna (is_terminal, reason).
    """
    adj = topology.adjacency
    n1 = len(adj.get(atom1, []))
    n2 = len(adj.get(atom2, []))

    if n1 <= 1 or n2 <= 1:
        return True, "terminal_atom"

    return False, ""


# ─── Clasificación de enlaces ─────────────────────────────────────────────────

def _classify_bond(
    bond:        dict,
    atoms:       list[dict],
    topology:    TopologyDescriptor,
    aromatic_set: set[int],
) -> RotatableBondInfo:
    """
    Clasifica un enlace individual como rotable o excluido.
    """
    a1_idx = bond["atom1"]
    a2_idx = bond["atom2"]
    btype  = bond["bond_type"]

    atom_map = {a["index"]: a for a in atoms}
    a1 = atom_map.get(a1_idx)
    a2 = atom_map.get(a2_idx)

    # ── Excluir: involucra hidrógeno ──────────────────────────────────────────
    if not a1 or not a2:
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason="atom_not_found"
        )
    if a1["element"] == "H" or a2["element"] == "H":
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason="hydrogen_bond"
        )

    # ── Excluir: enlace aromático (tipo 4) ────────────────────────────────────
    if btype == 4:
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason="aromatic_bond"
        )

    # ── Excluir: enlace doble o triple (sin rotación libre) ──────────────────
    if btype in (2, 3):
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason="multiple_bond"
        )

    # ── Excluir: ambos átomos en el mismo anillo ──────────────────────────────
    if _both_in_same_ring(a1_idx, a2_idx, topology):
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason="ring_bond"
        )

    # ── Excluir: enlace terminal ──────────────────────────────────────────────
    is_term, term_reason = _is_terminal_bond(a1_idx, a2_idx, topology)
    if is_term:
        return RotatableBondInfo(
            atom1=a1_idx, atom2=a2_idx, bond_type=btype,
            reason="", is_excluded=True, exclusion_reason=term_reason
        )

    # ── Rotable ───────────────────────────────────────────────────────────────
    reason = f"{a1['element']}-{a2['element']} single bond"
    if a1_idx in topology.ring_atom_indices or a2_idx in topology.ring_atom_indices:
        reason += " (exocyclic)"

    return RotatableBondInfo(
        atom1=a1_idx, atom2=a2_idx, bond_type=btype,
        reason=reason, is_excluded=False
    )


# ─── Clasificación de scaffold ────────────────────────────────────────────────

def _classify_scaffold(
    topology:    TopologyDescriptor,
    aromaticity: AromaticityDescriptor,
) -> str:
    if topology.n_rings == 0:
        return "acyclic"
    if aromaticity.n_fused_aromatic > 0:
        if topology.n_rings >= 3:
            return "rigid_core"    # sistema fusionado grande
        return "mixed"
    if topology.n_rings >= 2:
        return "mixed"
    return "flexible_core"


# ─── Interfaz pública ────────────────────────────────────────────────────────

def compute_flexibility(
    atoms:       list[dict],
    bonds:       list[dict],
    topology:    TopologyDescriptor,
    aromaticity: AromaticityDescriptor,
) -> FlexibilityDescriptor:
    """
    Calcula el descriptor de flexibilidad conformacional.
    """
    aromatic_set = set(aromaticity.aromatic_atom_set)

    rotatable: list[RotatableBondInfo] = []
    excluded:  list[RotatableBondInfo] = []

    for bond in bonds:
        info = _classify_bond(bond, atoms, topology, aromatic_set)
        if info.is_excluded:
            excluded.append(info)
        else:
            rotatable.append(info)

    n_rot = len(rotatable)

    # ── Clasificación de flexibilidad ─────────────────────────────────────────
    if n_rot <= 3:
        flex_class = "rigid"
    elif n_rot <= 7:
        flex_class = "moderate"
    elif n_rot <= 14:
        flex_class = "flexible"
    else:
        flex_class = "very_flexible"

    # ── Score continuo (para comparación futura entre ligandos) ──────────────
    # Normalizado al número de átomos pesados para ser tamaño-independiente
    n_heavy = topology.n_heavy_atoms or 1
    flex_score = round(n_rot / n_heavy, 3)

    # ── Scaffold ──────────────────────────────────────────────────────────────
    scaffold = _classify_scaffold(topology, aromaticity)

    # ── Recomendación de muestreo ─────────────────────────────────────────────
    if flex_class == "rigid":
        sampling = (
            "Muestreo estándar suficiente. "
            "El ligando es rígido — la pose de docking es representativa."
        )
    elif flex_class == "moderate":
        sampling = (
            "Muestreo MD estándar suficiente. "
            "Recomendable verificar 2-3 poses de docking antes de producción."
        )
    elif flex_class == "flexible":
        sampling = (
            f"Ligando flexible ({n_rot} rot. bonds): considerar docking exhaustivo "
            "antes de MD. Equilibración extendida recomendada (≥10 ns)."
        )
    else:
        sampling = (
            f"Ligando muy flexible ({n_rot} rot. bonds): requiere estrategia especial. "
            "Opciones: REST2, metadinámica, o docking con RMSD clustering previo. "
            "MD estándar puede no explorar suficiente espacio conformacional."
        )

    return FlexibilityDescriptor(
        n_rotatable_bonds         = n_rot,
        n_excluded_bonds          = len(excluded),
        rotatable_bonds           = rotatable,
        excluded_bonds            = excluded,
        n_rings                   = topology.n_rings,
        n_fused_aromatic_systems  = aromaticity.n_fused_aromatic,
        scaffold_rigidity         = scaffold,
        flexibility_class         = flex_class,
        flexibility_score         = flex_score,
        sampling_recommendation   = sampling,
    )
