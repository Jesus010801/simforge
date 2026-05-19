# validators/protein_validator.py
"""
Validator estructural de archivos PDB.

Implementación actual: parser propio mínimo (línea por línea).
Implementación futura: MDAnalysis (sin cambiar la interfaz pública).

Interfaz pública garantizada:
    validate_protein(path: str | Path) -> ProteinValidationResult

Agregar MDAnalysis = reemplazar _parse_pdb() únicamente.
"""

from __future__ import annotations
from pathlib import Path
from pydantic import BaseModel
from typing import Optional
from pydantic import Field
from core.models import Warning, Risk, Recommendation, Severity


# ─── Resultado de validación ─────────────────────────────────────────────────

class ResidueInfo(BaseModel):
    chain:   str
    number:  int
    name:    str


class ProteinValidationResult(BaseModel):
    
    source_file:      str

    # Estructura detectada
    chains:           list[str]         = []
    total_residues:   int               = 0
    has_hydrogens:    bool              = False
    has_hetatm:       list[str]         = []   # nombres de HETATM encontrados

    # Problemas detectados
    missing_residues: list[ResidueInfo] = []
    exposed_termini:  list[str]         = []   # ej: ["A_N", "A_C", "B_N"]

    # ─── Oligomerización ─────────────────

    oligomeric_state: Optional[str] = None

    chain_sizes: dict = Field(
        default_factory=dict
    )

    likely_oligomer: bool = False

    # Inferencia
    warnings:         list[Warning]         = []
    risks:            list[Risk]            = []
    recommendations:  list[Recommendation] = []
   

# ─── Parser PDB mínimo (implementación actual) ───────────────────────────────

def _parse_pdb(path: Path) -> dict:
    chains         = set()
    residues       = {}   # (chain, resnum) → resname
    hetatm_names   = set()
    missing        = []
    has_hydrogens  = False

    with open(path, "r") as f:
        for line in f:

            record = line[:6].strip()

            # ─── ATOM ────────────────────────────────────────────────────────
            if record == "ATOM":
                chain  = line[21].strip()
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                atom_name = line[12:16].strip()

                chains.add(chain)
                residues[(chain, resnum)] = resname

                if atom_name.startswith("H"):
                    has_hydrogens = True

            # ─── HETATM ──────────────────────────────────────────────────────
            elif record == "HETATM":
                resname = line[17:20].strip()
                if resname != "HOH":   # agua es esperada, ignorar
                    hetatm_names.add(resname)

            # ─── REMARK 465 (residuos faltantes) ─────────────────────────────
            elif record == "REMARK":
                remark_num = line[6:10].strip()
                if remark_num == "465":
                    content = line[10:].strip()
                    # Formato típico: "  M RES C SSSS"
                    parts = content.split()
                    if len(parts) >= 3:
                        try:
                            resname = parts[-3] if len(parts) >= 3 else ""
                            chain   = parts[-2] if len(parts) >= 2 else ""
                            resnum  = int(parts[-1])
                            missing.append({
                                "chain":  chain,
                                "number": resnum,
                                "name":   resname,
                            })
                        except (ValueError, IndexError):
                            pass

    return {
        "chains":        sorted(chains),
        "residues":      residues,
        "hetatm_names":  list(hetatm_names),
        "missing":       missing,
        "has_hydrogens": has_hydrogens,
    }


# ─── Detección de terminales expuestos ───────────────────────────────────────

def _detect_exposed_termini(residues: dict) -> list[str]:
    """
    Detecta terminales N y C por cadena.
    Retorna lista de strings: ["A_N", "A_C", "B_N", "B_C"]
    """
    from collections import defaultdict

    by_chain = defaultdict(list)
    for (chain, resnum), _ in residues.items():
        by_chain[chain].append(resnum)

    termini = []
    for chain, resnums in by_chain.items():
        resnums_sorted = sorted(resnums)
        termini.append(f"{chain}_N")   # terminal N = primer residuo
        termini.append(f"{chain}_C")   # terminal C = último residuo

    return termini

# ─── Detección de oligomerización ────────────────────────────────────────────

def _detect_oligomeric_state(residues: dict, chains: list[str]) -> dict:
    """
    Detecta posible oligomerización comparando tamaño de cadenas.
    Usa tolerancia del 5% para cubrir asimetrías cristalográficas reales.
    """
    from collections import defaultdict

    by_chain = defaultdict(list)
    for (chain, resnum), resname in residues.items():
        by_chain[chain].append(resnum)

    chain_sizes = {chain: len(resnums) for chain, resnums in by_chain.items()}

    if len(chain_sizes) < 2:
        return {
            "chain_sizes":       chain_sizes,
            "oligomeric_groups": {},
            "likely_oligomer":   False,
        }

    # Agrupar cadenas por tamaño similar (tolerancia 5%)
    chain_list   = list(chain_sizes.items())
    visited      = set()
    similar_groups = []

    for i, (chain_i, size_i) in enumerate(chain_list):
        if chain_i in visited:
            continue
        group = [chain_i]
        for j, (chain_j, size_j) in enumerate(chain_list):
            if i == j or chain_j in visited:
                continue
            tolerance = 0.05 * max(size_i, size_j)
            if abs(size_i - size_j) <= tolerance:
                group.append(chain_j)
        if len(group) > 1:
            for c in group:
                visited.add(c)
            avg_size = sum(chain_sizes[c] for c in group) // len(group)
            similar_groups.append({
                "chains":    group,
                "avg_size":  avg_size,
                "sizes":     {c: chain_sizes[c] for c in group},
            })

    oligomeric_groups = {}
    for g in similar_groups:
        avg = g["avg_size"]
        oligomeric_groups[avg] = g["chains"]

    return {
        "chain_sizes":       chain_sizes,
        "oligomeric_groups": oligomeric_groups,
        "likely_oligomer":   len(similar_groups) > 0,
    }


def _classify_oligomer(n_chains: int) -> str:
    mapping = {
        2: "homodímero",
        3: "homotrímero",
        4: "homotetrámero",
        6: "homohexámero",
        8: "homooctámero",
    }
    return mapping.get(n_chains, f"homo-oligómero ({n_chains} cadenas)")
# ─── Interfaz pública ────────────────────────────────────────────────────────

def validate_protein(path: str | Path) -> ProteinValidationResult:
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"PDB no encontrado: {path}")

    raw = _parse_pdb(path)

    result = ProteinValidationResult(
        source_file    = str(path),
        chains         = raw["chains"],
        total_residues = len(raw["residues"]),
        has_hydrogens  = raw["has_hydrogens"],
        has_hetatm     = raw["hetatm_names"],
        missing_residues = [
            ResidueInfo(
                chain  = m["chain"],
                number = m["number"],
                name   = m["name"],
            )
            for m in raw["missing"]
        ],
        exposed_termini = _detect_exposed_termini(raw["residues"]),
    )
    
    # 2. Detectar oligomerización y enriquecer resultado
    oligo = _detect_oligomeric_state(raw["residues"], raw["chains"])
    result.chain_sizes     = oligo["chain_sizes"]
    result.likely_oligomer = oligo["likely_oligomer"]
    
    if oligo["likely_oligomer"]:
        for avg_size, chains_group in oligo["oligomeric_groups"].items():
            n = len(chains_group)
            classification = _classify_oligomer(n)
            result.oligomeric_state = classification

            # Tamaños reales por cadena para el mensaje
            sizes_str = ", ".join(
                f"{c}:{oligo['chain_sizes'][c]}" for c in chains_group
            )

            result.warnings.append(Warning(
                message  = (
                    f"Posible {classification} detectado: "
                    f"cadenas con tamaños similares [{sizes_str}]"
                ),
                target   = str(path.name),
                severity = Severity.MEDIUM,
            ))
            result.recommendations.append(Recommendation(
                message = (
                    f"Decidir unidad de simulación: "
                    f"monómero funcional (1 cadena) o complejo completo ({n} cadenas)"
                ),
                target  = str(path.name),
                action  = (
                    f"Para monómero: extraer cadena A con 'grep ^ATOM hmg.pdb | awk $5==\"A\"'. "
                    f"Para complejo completo: usar el PDB tal como está. "
                    f"HMG-CoA reductasa es funcionalmente activa como {classification}."
                ),
            ))

    # ─── Inferencia de warnings, risks, recommendations ──────────────────────

    # Residuos faltantes
    if result.missing_residues:
        n = len(result.missing_residues)
        result.warnings.append(Warning(
            message  = f"{n} residuo(s) faltante(s) detectados via REMARK 465",
            target   = str(path.name),
            severity = Severity.HIGH,
        ))
        result.risks.append(Risk(
            message  = "Gaps estructurales pueden generar discontinuidades en la topología",
            target   = str(path.name),
            severity = Severity.HIGH,
        ))
        result.recommendations.append(Recommendation(
            message = "Modelar residuos faltantes con Modeller o SWISS-MODEL antes de continuar",
            target  = str(path.name),
            action  = "Usar loop modeling para gaps internos; ignorar solo si son terminales no funcionales",
        ))

    # HETATM inesperados
    if result.has_hetatm:
        result.warnings.append(Warning(
            message  = f"HETATM no-agua detectados: {result.has_hetatm}",
            target   = str(path.name),
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar si los HETATM son cofactores, ligandos cristalográficos o artefactos",
            target  = str(path.name),
            action  = "Eliminar o parametrizar explícitamente antes de pdb2gmx",
        ))

    # Sin hidrógenos
    if not result.has_hydrogens:
        result.warnings.append(Warning(
            message  = "PDB sin hidrógenos detectados",
            target   = str(path.name),
            severity = Severity.LOW,
        ))
        result.recommendations.append(Recommendation(
            message = "Los hidrógenos serán agregados por pdb2gmx durante la construcción",
            target  = str(path.name),
            action  = "Verificar estado de protonación en residuos His, Asp, Glu antes de pdb2gmx",
        ))

    return result