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
from core.models import Warning, Risk, Recommendation, Severity


# ─── Resultado de validación ─────────────────────────────────────────────────

class ResidueInfo(BaseModel):
    chain:   str
    number:  int
    name:    str


class ProteinValidationResult(BaseModel):
    """
    Contrato de salida del protein_validator.
    Siempre retorna este objeto, independientemente
    de la implementación interna.
    """
    source_file:      str

    # Estructura detectada
    chains:           list[str]         = []
    total_residues:   int               = 0
    has_hydrogens:    bool              = False
    has_hetatm:       list[str]         = []   # nombres de HETATM encontrados

    # Problemas detectados
    missing_residues: list[ResidueInfo] = []
    exposed_termini:  list[str]         = []   # ej: ["A_N", "A_C", "B_N"]

    # Inferencia
    warnings:         list[Warning]         = []
    risks:            list[Risk]            = []
    recommendations:  list[Recommendation] = []


# ─── Parser PDB mínimo (implementación actual) ───────────────────────────────

def _parse_pdb(path: Path) -> dict:
    """
    Lee el PDB línea por línea.
    Retorna datos crudos para que validate_protein() los interprete.

    Reemplazar esta función para integrar MDAnalysis o BioPython
    sin tocar el resto del módulo.
    """
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


# ─── Interfaz pública ────────────────────────────────────────────────────────

def validate_protein(path: str | Path) -> ProteinValidationResult:
    """
    Valida un archivo PDB y retorna ProteinValidationResult.

    Esta es la única función que el resto del sistema debe llamar.
    La implementación interna puede cambiar sin afectar contratos.
    """
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