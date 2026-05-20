# descriptors/polarity.py
"""
Descriptor de polaridad molecular.

Infiere propiedades electrónicas y de solvatación desde elementos,
conectividad y carga formal — sin RDKit.

Depende de:
    topology.py — grafo molecular, element_counts

Interfaz pública:
    compute_polarity(
        atoms:    list[dict],
        bonds:    list[dict],
        topology: TopologyDescriptor,
    ) -> PolarityDescriptor
"""

from __future__ import annotations
import math
from pydantic import BaseModel
from descriptors.topology import TopologyDescriptor


# ─── Contrato de salida ───────────────────────────────────────────────────────

class FunctionalGroup(BaseModel):
    name:         str         # "hydroxyl", "carbonyl", "phosphate", etc.
    atom_indices: list[int]   # átomos que lo componen
    smarts_like:  str         # descripción legible del patrón


class PolarityDescriptor(BaseModel):
    # Carga
    net_charge:         int   = 0
    n_charged_atoms:    int   = 0

    # H-bond (reglas Lipinski extendidas)
    hbd_count:          int   = 0    # donors: NH, OH
    hba_count:          int   = 0    # acceptors: N, O (no en anillo aromático contados igual)
    lipinski_compliant: bool  = True  # MW<500, logP<5, HBD<5, HBA<10

    # Grupos funcionales detectados
    functional_groups:  list[FunctionalGroup] = []
    group_names:        list[str]             = []   # lista plana para acceso rápido

    # logP estimado (método fragmentario de Wildman-Crippen simplificado)
    logp_estimate:      float = 0.0
    logp_class:         str   = "unknown"   # hydrophilic / moderate / lipophilic / very_lipophilic

    # Polaridad global
    polarity_class:     str   = "unknown"   # nonpolar / polar / charged
    estimated_solubility: str = "unknown"   # poor / moderate / good


# ─── Electronegatividades de Pauling ─────────────────────────────────────────
# Usadas para estimar contribución a polaridad

_ELECTRONEGATIVITY: dict[str, float] = {
    "H":  2.20, "C":  2.55, "N":  3.04, "O":  3.44,
    "S":  2.58, "P":  2.19, "F":  3.98, "Cl": 3.16,
    "Br": 2.96, "I":  2.66, "Se": 2.55, "B":  2.04,
}
_DEFAULT_EN = 2.0


# ─── logP fragmentario (Wildman-Crippen simplificado) ────────────────────────
# Contribuciones por tipo de átomo (muy simplificadas, sin RDKit).
# Valores aproximados de la literatura para estimación rápida.

_LOGP_CONTRIBUTIONS: dict[str, float] = {
    "C_aromatic":    0.13,
    "C_aliphatic":   0.20,
    "N_aromatic":   -0.96,
    "N_aliphatic":  -1.03,
    "O_carbonyl":   -0.45,
    "O_hydroxyl":   -0.67,
    "O_ether":      -0.44,
    "S_thioether":   0.48,
    "S_thiol":       0.52,
    "P":            -1.50,
    "F":             0.14,
    "Cl":            0.60,
    "Br":            0.88,
    "I":             1.35,
    "H_on_C":        0.12,
    "H_on_N":       -0.20,
    "H_on_O":       -0.20,
}


# ─── Detección de grupos funcionales ─────────────────────────────────────────

def _detect_functional_groups(
    atoms:    list[dict],
    bonds:    list[dict],
    topology: TopologyDescriptor,
) -> list[FunctionalGroup]:
    """
    Detecta grupos funcionales comunes por patrones de conectividad.
    Cada patrón busca un átomo central con vecinos de elementos específicos
    y tipos de enlace específicos.
    """
    atom_map  = {a["index"]: a for a in atoms}
    adj       = topology.adjacency
    ring_set  = topology.ring_atom_indices

    # Construir índice de bonds por átomo
    bonds_by_atom: dict[int, list[dict]] = {}
    for b in bonds:
        bonds_by_atom.setdefault(b["atom1"], []).append(b)
        bonds_by_atom.setdefault(b["atom2"], []).append(b)

    def get_bonds(idx: int) -> list[dict]:
        return bonds_by_atom.get(idx, [])

    def bond_type_to(idx: int, neighbor: int) -> int:
        for b in get_bonds(idx):
            other = b["atom2"] if b["atom1"] == idx else b["atom1"]
            if other == neighbor:
                return b["bond_type"]
        return 0

    def neighbors_of_element(idx: int, element: str) -> list[int]:
        return [nb for nb in adj.get(idx, []) if atom_map.get(nb, {}).get("element") == element]

    groups: list[FunctionalGroup] = []

    for atom in atoms:
        idx = atom["index"]
        el  = atom["element"]

        # ── Hydroxyl: O con 1 H y 1 C (no C=O) ─────────────────────────────
        if el == "O":
            h_neighbors = [nb for nb in bonds_by_atom.get(idx, [])
                           if atom_map.get(
                               nb["atom2"] if nb["atom1"]==idx else nb["atom1"],
                               {}).get("element") == "H"]
            c_neighbors = neighbors_of_element(idx, "C")
            if len(h_neighbors) >= 1 and len(c_neighbors) >= 1:
                # Verificar que no es C=O (carbonilo)
                c_idx = c_neighbors[0]
                if bond_type_to(idx, c_idx) == 1:
                    h_idx = (h_neighbors[0]["atom2"]
                             if h_neighbors[0]["atom1"] == idx
                             else h_neighbors[0]["atom1"])
                    groups.append(FunctionalGroup(
                        name="hydroxyl",
                        atom_indices=[idx, h_idx],
                        smarts_like="[OH]",
                    ))

        # ── Carbonyl: C=O ────────────────────────────────────────────────────
        if el == "C":
            for nb in adj.get(idx, []):
                if atom_map.get(nb, {}).get("element") == "O":
                    if bond_type_to(idx, nb) == 2:
                        groups.append(FunctionalGroup(
                            name="carbonyl",
                            atom_indices=[idx, nb],
                            smarts_like="[C]=O",
                        ))
                        break

        # ── Amine: N con H(s) ────────────────────────────────────────────────
        if el == "N":
            all_bonds = bonds_by_atom.get(idx, [])
            h_nbs = [b for b in all_bonds
                     if atom_map.get(
                         b["atom2"] if b["atom1"]==idx else b["atom1"],
                         {}).get("element") == "H"]
            if h_nbs:
                h_indices = [
                    b["atom2"] if b["atom1"]==idx else b["atom1"]
                    for b in h_nbs
                ]
                name = "primary_amine" if len(h_nbs) >= 2 else "secondary_amine"
                groups.append(FunctionalGroup(
                    name=name,
                    atom_indices=[idx] + h_indices,
                    smarts_like="[NH2]" if len(h_nbs) >= 2 else "[NH]",
                ))

        # ── Phosphate: P con ≥2 O ────────────────────────────────────────────
        if el == "P":
            o_neighbors = neighbors_of_element(idx, "O")
            if len(o_neighbors) >= 2:
                groups.append(FunctionalGroup(
                    name="phosphate",
                    atom_indices=[idx] + o_neighbors,
                    smarts_like="[PO4]",
                ))

        # ── Thiol: S con H ───────────────────────────────────────────────────
        if el == "S":
            all_bonds = bonds_by_atom.get(idx, [])
            h_nbs = [b for b in all_bonds
                     if atom_map.get(
                         b["atom2"] if b["atom1"]==idx else b["atom1"],
                         {}).get("element") == "H"]
            if h_nbs:
                groups.append(FunctionalGroup(
                    name="thiol",
                    atom_indices=[idx],
                    smarts_like="[SH]",
                ))

        # ── Carboxyl: C(=O)OH ────────────────────────────────────────────────
        if el == "C":
            o_neighbors = neighbors_of_element(idx, "O")
            if len(o_neighbors) >= 2:
                double_o = [nb for nb in o_neighbors if bond_type_to(idx, nb) == 2]
                single_o = [nb for nb in o_neighbors if bond_type_to(idx, nb) == 1]
                if double_o and single_o:
                    groups.append(FunctionalGroup(
                        name="carboxyl",
                        atom_indices=[idx] + o_neighbors,
                        smarts_like="[C](=O)O",
                    ))

    return groups


# ─── H-bond donors y acceptors ───────────────────────────────────────────────

def _count_hbond(
    atoms:    list[dict],
    bonds:    list[dict],
    groups:   list[FunctionalGroup],
) -> tuple[int, int]:
    """
    Cuenta H-bond donors (HBD) y acceptors (HBA).

    HBD: NH o OH (N o O con H unido)
    HBA: N u O con par solitario (todos los N y O no cargados positivamente)

    Reglas de Lipinski simplificadas.
    """
    atom_map = {a["index"]: a for a in atoms}

    hbd = 0
    hba = 0

    # Donors: grupos con OH o NH
    for g in groups:
        if g.name in ("hydroxyl", "primary_amine", "secondary_amine", "thiol"):
            hbd += 1

    # Acceptors: contar todos los N y O
    for atom in atoms:
        if atom["element"] in ("N", "O"):
            hba += 1

    return hbd, hba


# ─── logP estimado ────────────────────────────────────────────────────────────

def _estimate_logp(
    atoms:       list[dict],
    bonds:       list[dict],
    topology:    TopologyDescriptor,
    groups:      list[FunctionalGroup],
) -> float:
    """
    Estimación fragmentaria simplificada de logP.
    No es Crippen exacto — es una aproximación de orden de magnitud
    útil para clasificar hidrofílico/lipofílico.
    """
    logp = 0.0
    group_names = {g.name for g in groups}
    ring_set    = topology.ring_atom_indices

    adj = topology.adjacency

    for atom in atoms:
        el  = atom["element"]
        idx = atom["index"]

        if el == "C":
            contrib = (_LOGP_CONTRIBUTIONS["C_aromatic"]
                       if idx in ring_set
                       else _LOGP_CONTRIBUTIONS["C_aliphatic"])
        elif el == "N":
            contrib = (_LOGP_CONTRIBUTIONS["N_aromatic"]
                       if idx in ring_set
                       else _LOGP_CONTRIBUTIONS["N_aliphatic"])
        elif el == "O":
            # Distinguir carbonilo vs hidroxilo vs éter
            if "carbonyl" in group_names:
                contrib = _LOGP_CONTRIBUTIONS["O_carbonyl"]
            elif "hydroxyl" in group_names:
                contrib = _LOGP_CONTRIBUTIONS["O_hydroxyl"]
            else:
                contrib = _LOGP_CONTRIBUTIONS["O_ether"]
        elif el == "S":
            contrib = (_LOGP_CONTRIBUTIONS["S_thiol"]
                       if "thiol" in group_names
                       else _LOGP_CONTRIBUTIONS["S_thioether"])
        elif el == "P":
            contrib = _LOGP_CONTRIBUTIONS["P"]
        elif el in ("F", "Cl", "Br", "I"):
            contrib = _LOGP_CONTRIBUTIONS.get(el, 0.0)
        elif el == "H":
            # Buscar átomo al que está unido para determinar contribución
            contrib = _LOGP_CONTRIBUTIONS["H_on_C"]   # default
            for b in bonds:
                other = None
                if b["atom1"] == idx:
                    other = b["atom2"]
                elif b["atom2"] == idx:
                    other = b["atom1"]
                if other is not None:
                    other_el = next(
                        (a["element"] for a in atoms if a["index"] == other), "C"
                    )
                    if other_el == "N":
                        contrib = _LOGP_CONTRIBUTIONS["H_on_N"]
                    elif other_el == "O":
                        contrib = _LOGP_CONTRIBUTIONS["H_on_O"]
                    break
        else:
            contrib = 0.0

        logp += contrib

    return round(logp, 2)


# ─── Clasificaciones finales ──────────────────────────────────────────────────

def _logp_class(logp: float) -> str:
    if logp < 0:
        return "hydrophilic"
    elif logp < 2:
        return "moderate"
    elif logp < 5:
        return "lipophilic"
    else:
        return "very_lipophilic"


def _polarity_class(net_charge: int, hba: int, hbd: int) -> str:
    if net_charge != 0:
        return "charged"
    elif hba > 0 or hbd > 0:
        return "polar"
    else:
        return "nonpolar"


def _solubility_estimate(logp: float, hbd: int, hba: int) -> str:
    """
    Estimación cualitativa de solubilidad en agua.
    Basada en reglas de Lipinski + logP.
    """
    if logp < 1 and (hbd + hba) >= 3:
        return "good"
    elif logp < 3 and (hbd + hba) >= 1:
        return "moderate"
    else:
        return "poor"


def _lipinski_check(
    n_heavy: int,
    logp:    float,
    hbd:     int,
    hba:     int,
) -> bool:
    """
    Regla de Lipinski para drug-likeness (Ro5).
    MW aproximado: n_heavy_atoms × 6.7 (media empírica para C,N,O,H mix)
    """
    mw_approx = n_heavy * 6.7
    return (mw_approx <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)


# ─── Interfaz pública ────────────────────────────────────────────────────────

def compute_polarity(
    atoms:    list[dict],
    bonds:    list[dict],
    topology: TopologyDescriptor,
) -> PolarityDescriptor:
    """
    Calcula el descriptor de polaridad molecular.
    """
    # ── Carga ─────────────────────────────────────────────────────────────────
    net_charge     = sum(a.get("charge", 0) for a in atoms)
    n_charged      = sum(1 for a in atoms if a.get("charge", 0) != 0)

    # ── Grupos funcionales ────────────────────────────────────────────────────
    groups     = _detect_functional_groups(atoms, bonds, topology)
    group_names = list({g.name for g in groups})

    # ── H-bond ────────────────────────────────────────────────────────────────
    hbd, hba = _count_hbond(atoms, bonds, groups)

    # ── logP ──────────────────────────────────────────────────────────────────
    logp      = _estimate_logp(atoms, bonds, topology, groups)
    lclass    = _logp_class(logp)

    # ── Clasificaciones ───────────────────────────────────────────────────────
    pclass    = _polarity_class(net_charge, hba, hbd)
    solubility = _solubility_estimate(logp, hbd, hba)
    lipinski  = _lipinski_check(topology.n_heavy_atoms, logp, hbd, hba)

    return PolarityDescriptor(
        net_charge           = net_charge,
        n_charged_atoms      = n_charged,
        hbd_count            = hbd,
        hba_count            = hba,
        lipinski_compliant   = lipinski,
        functional_groups    = groups,
        group_names          = group_names,
        logp_estimate        = logp,
        logp_class           = lclass,
        polarity_class       = pclass,
        estimated_solubility = solubility,
    )
