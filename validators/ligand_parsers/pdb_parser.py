# validators/ligand_parsers/pdb_parser.py
"""
Parser PDB para ligandos pequeños (sin sección CONECT).

Los PDB de ligandos producidos por herramientas como Avogadro, ChemDraw
o extraídos de estructuras cristalográficas frecuentemente no incluyen
la sección CONECT, por lo que la conectividad debe inferirse por geometría.

Estrategia de inferencia:
    - Enlace: distancia entre átomos < umbral por par de elementos
    - Aromaticidad: anillo de 5-6 átomos plano (desviación < 0.1Å del plano)
      compuesto por elementos típicos (C, N, O, S)
    - Tipo de enlace: heurística por distancia (doble < 1.35Å C-C, etc.)
      — no es química exacta; RDKit lo resolverá correctamente después

Interfaz pública:
    parse_pdb_ligand(path: Path) -> dict

Retorna el mismo contrato que sdf_parser.parse_sdf():
    {
        "mol_name": str,
        "n_atoms":  int,
        "n_bonds":  int,
        "atoms":    list[dict],   # index, element, x, y, z, charge=0
        "bonds":    list[dict],   # atom1, atom2, bond_type (1/2/4)
        "error":    str | None,
    }
"""

from __future__ import annotations
from pathlib import Path
import math


# ─── Umbrales de distancia de enlace (Å) ─────────────────────────────────────
# Conservadores: si dos átomos están más cerca que el umbral → enlace.
# Valores de radio de covalencia extendidos ~10% para tolerar geometrías
# cristalográficas imperfectas.

_COVALENT_RADII: dict[str, float] = {
    "H":  0.31,
    "C":  0.76,
    "N":  0.71,
    "O":  0.66,
    "S":  1.05,
    "P":  1.07,
    "F":  0.57,
    "Cl": 1.02,
    "Br": 1.20,
    "I":  1.39,
}

_DEFAULT_RADIUS = 0.90   # fallback para elementos no listados
_BOND_TOLERANCE = 0.40   # Å de tolerancia extra sobre la suma de radios


def _bond_threshold(el1: str, el2: str) -> float:
    r1 = _COVALENT_RADII.get(el1, _DEFAULT_RADIUS)
    r2 = _COVALENT_RADII.get(el2, _DEFAULT_RADIUS)
    return r1 + r2 + _BOND_TOLERANCE


def _distance(a1: dict, a2: dict) -> float:
    return math.sqrt(
        (a1["x"] - a2["x"]) ** 2 +
        (a1["y"] - a2["y"]) ** 2 +
        (a1["z"] - a2["z"]) ** 2
    )


# ─── Estimación de tipo de enlace por distancia ───────────────────────────────
# Heurística simplificada. Los umbrales son para pares C-C / C-N / C-O.
# Para el resto se usa solo enlace simple.

_DOUBLE_BOND_THRESHOLDS: dict[frozenset, float] = {
    frozenset({"C", "C"}): 1.42,   # C=C aromático / doble
    frozenset({"C", "N"}): 1.38,
    frozenset({"C", "O"}): 1.35,
    frozenset({"C", "S"}): 1.65,
    frozenset({"N", "N"}): 1.28,
    frozenset({"P", "O"}): 1.55,
}


def _estimate_bond_type(el1: str, el2: str, dist: float) -> int:
    """
    Retorna 1 (simple) o 2 (doble).
    La aromaticidad (tipo 4) se asigna después en _detect_aromatic_rings().
    """
    key = frozenset({el1, el2})
    threshold = _DOUBLE_BOND_THRESHOLDS.get(key)
    if threshold is not None and dist <= threshold:
        return 2
    return 1


# ─── Detección de anillos aromáticos ─────────────────────────────────────────

_AROMATIC_ELEMENTS = {"C", "N", "O", "S"}


def _detect_aromatic_rings(atoms: list[dict], bonds: list[dict]) -> set[int]:
    """
    Detecta átomos en anillos aromáticos por:
      1. Estar en un anillo de 5 o 6 miembros cerrado.
      2. El anillo es plano (desviación RMS del plano < 0.15Å).
      3. Solo elementos de _AROMATIC_ELEMENTS.

    Retorna set de índices (1-based) de átomos aromáticos.
    """
    from collections import defaultdict

    # Construir grafo de adyacencia (solo átomos no-H)
    adj: dict[int, list[int]] = defaultdict(list)
    atom_by_index = {a["index"]: a for a in atoms}

    for b in bonds:
        i, j = b["atom1"], b["atom2"]
        ai = atom_by_index.get(i)
        aj = atom_by_index.get(j)
        if ai and aj and ai["element"] != "H" and aj["element"] != "H":
            adj[i].append(j)
            adj[j].append(i)

    aromatic_atoms: set[int] = set()

    # DFS para encontrar ciclos de tamaño 5-6
    def find_cycles(start: int) -> list[list[int]]:
        cycles = []
        stack  = [(start, [start])]
        while stack:
            node, path = stack.pop()
            for neighbor in adj[node]:
                if neighbor == start and len(path) in (5, 6):
                    cycles.append(path[:])
                elif neighbor not in path and len(path) < 6:
                    stack.append((neighbor, path + [neighbor]))
        return cycles

    visited_rings: set[frozenset] = set()

    for atom in atoms:
        if atom["element"] == "H":
            continue
        for cycle in find_cycles(atom["index"]):
            key = frozenset(cycle)
            if key in visited_rings:
                continue
            visited_rings.add(key)

            # Solo elementos aromáticos
            elements_in_ring = {atom_by_index[i]["element"] for i in cycle}
            if not elements_in_ring.issubset(_AROMATIC_ELEMENTS):
                continue

            # Verificar planaridad
            coords = [(atom_by_index[i]["x"],
                       atom_by_index[i]["y"],
                       atom_by_index[i]["z"]) for i in cycle]
            if _is_planar(coords, threshold=0.15):
                for i in cycle:
                    aromatic_atoms.add(i)

    return aromatic_atoms


def _is_planar(coords: list[tuple], threshold: float = 0.15) -> bool:
    """
    Verifica si un conjunto de puntos 3D es aproximadamente plano.
    Calcula el plano de mínimos cuadrados y mide la desviación RMS.
    """
    if len(coords) < 3:
        return True

    # Centroide
    n   = len(coords)
    cx  = sum(c[0] for c in coords) / n
    cy  = sum(c[1] for c in coords) / n
    cz  = sum(c[2] for c in coords) / n

    # Matriz de covarianza 3x3 (manual, sin numpy)
    cov = [[0.0] * 3 for _ in range(3)]
    for c in coords:
        d = (c[0] - cx, c[1] - cy, c[2] - cz)
        for i in range(3):
            for j in range(3):
                cov[i][j] += d[i] * d[j]

    # Para moléculas pequeñas (≤6 átomos) usamos el producto cruzado de
    # dos vectores del anillo como normal del plano — suficiente para
    # la heurística sin necesidad de eigenvalores completos.
    v1 = (
        coords[1][0] - coords[0][0],
        coords[1][1] - coords[0][1],
        coords[1][2] - coords[0][2],
    )
    v2 = (
        coords[2][0] - coords[0][0],
        coords[2][1] - coords[0][1],
        coords[2][2] - coords[0][2],
    )
    # Normal = v1 × v2
    normal = (
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    )
    norm_len = math.sqrt(sum(x**2 for x in normal))
    if norm_len < 1e-6:
        return True   # puntos colineales → plano degenerado, aceptar

    normal = tuple(x / norm_len for x in normal)

    # Distancia de cada punto al plano
    deviations = []
    for c in coords:
        d = (c[0] - cx, c[1] - cy, c[2] - cz)
        dist = abs(sum(d[i] * normal[i] for i in range(3)))
        deviations.append(dist)

    rms = math.sqrt(sum(d**2 for d in deviations) / len(deviations))
    return rms < threshold


# ─── Extracción de elemento desde nombre de átomo PDB ────────────────────────
#
# Regla del formato PDB (wwPDB v3):
#   line[12] == ' ' → elemento de 1 letra, nombre empieza en col 13 (idx 13)
#                     → CA, CB, CD, CG, HA, HB, OC, SC, PF, NA...
#                        primer carácter del nombre = elemento garantizado
#   line[12] != ' ' → nombre ocupa col 12 (puede ser metal de 2 letras o H largo)
#                     → FE, ZN, MG, CA(ion), HD11, HG21...
#                        intentar metal bio-relevante, sino H o primera letra
#
# Esto resuelve la ambigüedad entre CA (alpha-carbon, elemento C)
# y Ca (calcio, elemento Ca), o CD (delta-carbon) y Cd (cadmio).

_BIO_METALS: set[str] = {
    "Fe", "Zn", "Mg", "Ca", "Mn", "Cu", "Co", "Ni", "Mo",
    "Na", "Cl", "Br", "Se", "Al", "Si", "Rb", "Sr", "Ba",
}


def _element_from_atom_name(atom_name: str, col12_char: str = " ") -> str:
    """
    Infiere el elemento desde el nombre de átomo PDB usando la posición col12.

    atom_name  : contenido de line[12:16].strip()
    col12_char : line[12] — el carácter crudo en columna 12 (sin strip)

    Si col12_char == ' ', el elemento es de 1 letra (convención PDB para C,N,O,H,S,P).
    Si col12_char != ' ', puede ser metal de 2 letras o H con nombre largo.
    """
    name = atom_name.strip()
    if not name:
        return "X"

    # Caso 1: col 12 en blanco → elemento de 1 letra garantizado por formato PDB
    if col12_char == " ":
        return name[0].upper()

    # Caso 2: col 12 ocupada → intentar metal bio-relevante
    if len(name) >= 2:
        candidate = name[0].upper() + name[1].lower()
        if candidate in _BIO_METALS:
            return candidate

    # Col 12 ocupada pero no es metal conocido (ej: HD11, HG21) → H
    if name[0].upper() == "H":
        return "H"

    # Fallback
    return name[0].upper()


# ─── Parser principal ─────────────────────────────────────────────────────────

def parse_pdb_ligand(path: Path) -> dict:
    """
    Parsea un PDB de ligando pequeño sin sección CONECT.

    Lee registros HETATM (y ATOM como fallback si no hay HETATM).
    Si hay múltiples residuos en el archivo, usa solo el primero
    (los PDB de ligandos exportados por herramientas a veces incluyen
    agua u otros registros extra).

    Infiere conectividad por distancia.
    Detecta aromaticidad por geometría.
    """
    atoms: list[dict] = []
    mol_name = path.stem
    first_resname: str | None = None   # residuo objetivo (el primero no-HOH)
    first_chain:   str | None = None
    first_resseq:  int | None = None

    raw_lines: list[str] = []

    with open(path, "r") as f:
        raw_lines = f.readlines()

    # ─── Primera pasada: detectar residuo objetivo ────────────────────────────
    # Preferir HETATM; si no hay, usar ATOM.
    has_hetatm = any(l[:6].strip() == "HETATM" for l in raw_lines)
    target_record = "HETATM" if has_hetatm else "ATOM"

    for line in raw_lines:
        record = line[:6].strip()
        if record != target_record:
            continue
        resname = line[17:20].strip()
        if resname in ("HOH", "WAT", "TIP"):   # ignorar agua
            continue
        chain  = line[21] if len(line) > 21 else " "
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            resseq = 0
        first_resname = resname
        first_chain   = chain
        first_resseq  = resseq
        mol_name      = resname
        break

    if first_resname is None:
        return {"error": "No se encontraron registros HETATM/ATOM válidos en el PDB"}

    # ─── Segunda pasada: leer solo átomos del residuo objetivo ───────────────
    for line in raw_lines:
        record = line[:6].strip()
        if record != target_record:
            continue

        resname = line[17:20].strip()
        chain   = line[21] if len(line) > 21 else " "
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            resseq = 0

        # Filtrar: solo el residuo identificado en la primera pasada
        if resname != first_resname or chain != first_chain or resseq != first_resseq:
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            # Columnas 76-78: elemento explícito (estándar PDB v3)
            element = line[76:78].strip() if len(line) > 76 else ""

            if not element:
                # Fallback: inferir desde nombre del átomo (cols 12-16)
                # Se pasa line[12] para distinguir átomos de 1 letra (col12=' ')
                # de posibles metales de 2 letras (col12!=' ').
                col12_char = line[12] if len(line) > 12 else " "
                atom_name  = line[12:16].strip()
                element    = _element_from_atom_name(atom_name, col12_char)

            # Normalizar capitalización del elemento
            if len(element) == 2:
                element = element[0].upper() + element[1].lower()
            else:
                element = element.upper()

            atoms.append({
                "index":   len(atoms) + 1,
                "element": element,
                "x": x, "y": y, "z": z,
                "charge":  0,
            })
        except (ValueError, IndexError) as e:
            return {"error": f"Error parseando línea ATOM/HETATM: {e}\n  → {line.rstrip()}"}

    if not atoms:
        return {"error": "No se encontraron átomos HETATM/ATOM en el PDB"}

    # ─── Inferir conectividad por distancia ────────────────────────────────────────────
    # Límite de valencia: guardia contra falsos positivos.
    _MAX_VALENCE: dict[str, int] = {
        "H": 1, "C": 4, "N": 4, "O": 2, "S": 6, "P": 5,
        "F": 1, "Cl": 1, "Br": 1, "I": 1,
    }
    _DEFAULT_VALENCE = 6

    bonds: list[dict] = []
    valence_count: dict[int, int] = {a['index']: 0 for a in atoms}
    n = len(atoms)

    # Ordenar pares por distancia para priorizar enlaces más cortos.
    candidates: list[tuple] = []
    for i in range(n):
        for j in range(i + 1, n):
            a1 = atoms[i]
            a2 = atoms[j]
            dist = _distance(a1, a2)
            threshold = _bond_threshold(a1["element"], a2["element"])
            if dist <= threshold:
                candidates.append((dist, a1, a2))

    candidates.sort(key=lambda x: x[0])

    for dist, a1, a2 in candidates:
        max1 = _MAX_VALENCE.get(a1["element"], _DEFAULT_VALENCE)
        max2 = _MAX_VALENCE.get(a2["element"], _DEFAULT_VALENCE)
        if valence_count[a1["index"]] >= max1:
            continue
        if valence_count[a2["index"]] >= max2:
            continue
        bond_type = _estimate_bond_type(a1["element"], a2["element"], dist)
        bonds.append({
            "atom1":     a1["index"],
            "atom2":     a2["index"],
            "bond_type": bond_type,
        })
        valence_count[a1["index"]] += 1
        valence_count[a2["index"]] += 1

    # ─── Detectar aromaticidad y reclasificar enlaces ─────────────────────────
    aromatic_atom_indices = _detect_aromatic_rings(atoms, bonds)
    if aromatic_atom_indices:
        for b in bonds:
            if b["atom1"] in aromatic_atom_indices and b["atom2"] in aromatic_atom_indices:
                b["bond_type"] = 4   # aromático

    return {
        "mol_name": mol_name,
        "n_atoms":  len(atoms),
        "n_bonds":  len(bonds),
        "atoms":    atoms,
        "bonds":    bonds,
        "error":    None,
    }