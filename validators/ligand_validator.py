# validators/ligand_validator.py
"""
Validator estructural y semántico de archivos SDF.

Implementación actual: parser SDF mínimo (línea por línea).
Implementación futura: RDKit (sin cambiar la interfaz pública).

Interfaz pública garantizada:
    validate_ligand(path: str | Path, role: str) -> LigandValidationResult

Agregar RDKit = reemplazar _parse_sdf() únicamente.
El reasoning engine (inferencia semántica) no cambia.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field
from core.models import Warning, Risk, Recommendation, Severity


# ─── Resultado de validación ─────────────────────────────────────────────────

class AtomInfo(BaseModel):
    index:        int
    element:      str
    x:            float
    y:            float
    z:            float
    charge:       int = 0


class BondInfo(BaseModel):
    atom1:      int
    atom2:      int
    bond_type:  int   # 1=single, 2=double, 3=triple, 4=aromatic


class LigandValidationResult(BaseModel):
    """
    Contrato de salida del ligand_validator.
    Siempre retorna este objeto, independientemente
    de la implementación interna.
    """
    source_file:    str
    role:           str

    # Estructura detectada
    n_atoms:        int                 = 0
    n_bonds:        int                 = 0
    atom_elements:  list[str]           = []
    atoms:          list[AtomInfo]      = []
    bonds:          list[BondInfo]      = []
    formal_charges: list[int]           = []
    net_charge:     int                 = 0
    is_complete:    bool                = False

    # Características químicas inferidas
    has_aromatic:       bool            = False
    aromatic_atoms:     int             = 0
    has_sulfur:         bool            = False
    has_phosphorus:     bool            = False
    has_nitrogen:       bool            = False
    has_oxygen:         bool            = False
    has_halogens:       bool            = False
    halogen_types:      list[str]       = []
    n_rotatable_bonds:  int             = 0
    estimated_flexibility: str          = "unknown"  # rigid / moderate / flexible
    estimated_polarity:    str          = "unknown"  # nonpolar / polar / charged
    parametrization_difficulty: str     = "unknown"  # low / medium / high

    # Inferencia
    warnings:         list[Warning]         = []
    risks:            list[Risk]            = []
    recommendations:  list[Recommendation] = []


# ─── Parser SDF mínimo (implementación actual) ───────────────────────────────

def _parse_sdf(path: Path) -> dict:
    

    lines = path.read_text().splitlines()

    if len(lines) < 4:
        return {"error": "Archivo SDF demasiado corto o vacío"}

    # ─── Header ──────────────────────────────────────────────────────────────
    mol_name = lines[0].strip()

    # ─── Counts line (línea 4, índice 3) ─────────────────────────────────────
    counts_line = lines[3]
    try:
        n_atoms = int(counts_line[0:3].strip())
        n_bonds = int(counts_line[3:6].strip())
    except (ValueError, IndexError):
        return {"error": f"Counts line inválida: '{counts_line}'"}

    # ─── Atom block ──────────────────────────────────────────────────────────
    atoms = []
    for i in range(n_atoms):
        line_idx = 4 + i
        if line_idx >= len(lines):
            return {"error": f"Atom block incompleto en línea {line_idx}"}
        line = lines[line_idx]
        try:
            x       = float(line[0:10].strip())
            y       = float(line[10:20].strip())
            z       = float(line[20:30].strip())
            element = line[31:34].strip()
            # charge code: campo opcional en posición 36-39
            charge_code = 0
            try:
                charge_code = int(line[36:39].strip()) if len(line) > 36 else 0
            except ValueError:
                charge_code = 0

            # SDF charge codes → carga formal real
            charge_map = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
            formal_charge = charge_map.get(charge_code, 0)

            atoms.append({
                "index":   i + 1,
                "element": element,
                "x": x, "y": y, "z": z,
                "charge":  formal_charge,
            })
        except (ValueError, IndexError) as e:
            return {"error": f"Error parseando átomo {i+1}: {e}"}

    # ─── Bond block ──────────────────────────────────────────────────────────
    bonds = []
    for i in range(n_bonds):
        line_idx = 4 + n_atoms + i
        if line_idx >= len(lines):
            return {"error": f"Bond block incompleto en línea {line_idx}"}
        line = lines[line_idx]
        try:
            atom1     = int(line[0:3].strip())
            atom2     = int(line[3:6].strip())
            bond_type = int(line[6:9].strip())
            bonds.append({
                "atom1": atom1,
                "atom2": atom2,
                "bond_type": bond_type,
            })
        except (ValueError, IndexError) as e:
            return {"error": f"Error parseando enlace {i+1}: {e}"}

    return {
        "mol_name":  mol_name,
        "n_atoms":   n_atoms,
        "n_bonds":   n_bonds,
        "atoms":     atoms,
        "bonds":     bonds,
        "error":     None,
    }


# ─── Análisis químico básico ──────────────────────────────────────────────────

def _analyze_chemistry(atoms: list[dict], bonds: list[dict]) -> dict:
    
    elements = [a["element"] for a in atoms]
    element_set = set(elements)

    formal_charges = [a["charge"] for a in atoms]
    net_charge = sum(formal_charges)

    # ─── Elementos presentes ─────────────────────────────────────────────────
    has_sulfur     = "S" in element_set
    has_phosphorus = "P" in element_set
    has_nitrogen   = "N" in element_set
    has_oxygen     = "O" in element_set
    halogens       = {"F", "Cl", "Br", "I"}
    halogen_types  = sorted(element_set & halogens)
    has_halogens   = len(halogen_types) > 0

    # ─── Aromaticidad (heurística) ────────────────────────────────────────────
    # Contar átomos en anillos aromáticos via bond_type == 4
    aromatic_bond_atoms = set()
    for b in bonds:
        if b["bond_type"] == 4:
            aromatic_bond_atoms.add(b["atom1"])
            aromatic_bond_atoms.add(b["atom2"])
    has_aromatic  = len(aromatic_bond_atoms) > 0
    aromatic_atoms = len(aromatic_bond_atoms)

    # ─── Rotatable bonds (heurística) ────────────────────────────────────────
    # Single bonds entre átomos no-H que no son parte de anillo aromático
    rotatable = 0
    for b in bonds:
        if b["bond_type"] == 1:
            a1 = next((a for a in atoms if a["index"] == b["atom1"]), None)
            a2 = next((a for a in atoms if a["index"] == b["atom2"]), None)
            if a1 and a2:
                if (a1["element"] != "H" and a2["element"] != "H"
                        and b["atom1"] not in aromatic_bond_atoms
                        and b["atom2"] not in aromatic_bond_atoms):
                    rotatable += 1
    n_rotatable = rotatable

    # ─── Flexibilidad estimada ───────────────────────────────────────────────
    if n_rotatable <= 3:
        flexibility = "rigid"
    elif n_rotatable <= 7:
        flexibility = "moderate"
    else:
        flexibility = "flexible"

    # ─── Polaridad estimada ──────────────────────────────────────────────────
    if net_charge != 0:
        polarity = "charged"
    elif has_nitrogen or has_oxygen or has_sulfur:
        polarity = "polar"
    else:
        polarity = "nonpolar"

    # ─── Dificultad de parametrización ───────────────────────────────────────
    difficulty_score = 0
    if has_sulfur:        difficulty_score += 2
    if has_phosphorus:    difficulty_score += 3
    if has_halogens:      difficulty_score += 1
    if net_charge != 0:   difficulty_score += 2
    if has_aromatic:      difficulty_score += 1
    if len(elements) > 50: difficulty_score += 2

    if difficulty_score <= 2:
        difficulty = "low"
    elif difficulty_score <= 5:
        difficulty = "medium"
    else:
        difficulty = "high"

    return {
        "formal_charges":    formal_charges,
        "net_charge":        net_charge,
        "has_aromatic":      has_aromatic,
        "aromatic_atoms":    aromatic_atoms,
        "has_sulfur":        has_sulfur,
        "has_phosphorus":    has_phosphorus,
        "has_nitrogen":      has_nitrogen,
        "has_oxygen":        has_oxygen,
        "has_halogens":      has_halogens,
        "halogen_types":     halogen_types,
        "n_rotatable_bonds": n_rotatable,
        "flexibility":       flexibility,
        "polarity":          polarity,
        "difficulty":        difficulty,
    }


# ─── Reasoning engine ────────────────────────────────────────────────────────

def _run_reasoning(
    result: LigandValidationResult,
    chem:   dict,
    path:   Path,
) -> LigandValidationResult:
    """
    Inferencia semántica y contextual.
    Este engine NO cambia cuando se reemplace el parser por RDKit.
    Razona sobre características químicas y rol del ligando.
    """

    # ─── Fosfato / CoA (sustrato complejo) ───────────────────────────────────
    if chem["has_phosphorus"] and result.role == "substrate":
        result.warnings.append(Warning(
            message  = "Sustrato con fósforo detectado — posible CoA o nucleótido",
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.risks.append(Risk(
            message  = "Moléculas con grupos fosfato son difíciles de parametrizar con CGenFF/GAFF",
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.recommendations.append(Recommendation(
            message = "Considerar parámetros CHARMM específicos para CoA o usar fragmentación",
            target  = path.name,
            action  = (
                "Verificar disponibilidad de parámetros en CHARMM36 para CoA. "
                "Alternativa: usar ParamChem con penalización y revisión manual."
            ),
        ))

    # ─── Carga neta ───────────────────────────────────────────────────────────
    if chem["net_charge"] != 0:
        result.warnings.append(Warning(
            message  = f"Ligando con carga neta {chem['net_charge']:+d} detectado",
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar estado de protonación a pH fisiológico (7.4)",
            target  = path.name,
            action  = (
                "Usar Epik o propKa para estimar pKa y estado de protonación correcto. "
                "La carga afecta directamente la parametrización y las interacciones electrostáticas."
            ),
        ))

    # ─── Aromaticidad (inhibidor competitivo) ────────────────────────────────
    if chem["has_aromatic"] and result.role == "competitive_ligand":
        result.recommendations.append(Recommendation(
            message = (
                f"Inhibidor aromático ({chem['aromatic_atoms']} átomos aromáticos): "
                f"verificar interacciones π-stacking con residuos del sitio activo"
            ),
            target  = path.name,
            action  = (
                "Inspeccionar sitio activo de HMG-CoA reductasa para Phe, Tyr, Trp, His. "
                "Las xantonas son planares y pueden formar π-stacking con residuos aromáticos."
            ),
        ))

    # ─── Flexibilidad ─────────────────────────────────────────────────────────
    if chem["flexibility"] == "flexible":
        result.warnings.append(Warning(
            message  = (
                f"Ligando flexible: {chem['n_rotatable_bonds']} enlaces rotables detectados"
            ),
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.risks.append(Risk(
            message  = (
                "Alta flexibilidad puede requerir muestreo conformacional extendido "
                "para encontrar la pose de unión correcta"
            ),
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Considerar docking previo para obtener pose inicial razonable",
            target  = path.name,
            action  = (
                "Usar AutoDock Vina o Glide para pose inicial antes de MD. "
                "Un ligando flexible mal posicionado puede no converger en el sitio activo."
            ),
        ))

    # ─── Halogenos ────────────────────────────────────────────────────────────
    if chem["has_halogens"]:
        result.warnings.append(Warning(
            message  = f"Halógenos detectados: {chem['halogen_types']}",
            target   = path.name,
            severity = Severity.LOW,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar parámetros de halógenos en CGenFF — pueden requerir ajuste manual",
            target  = path.name,
            action  = (
                "Los halógenos pesados (Br, I) tienen parámetros limitados en CGenFF. "
                "Revisar penalización en ParamChem y validar geometría post-parametrización."
            ),
        ))

    # ─── Dificultad alta de parametrización ──────────────────────────────────
    if chem["difficulty"] == "high":
        result.risks.append(Risk(
            message  = "Dificultad de parametrización alta — requiere validación manual",
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.recommendations.append(Recommendation(
            message = "Ejecutar ParamChem y revisar penalización antes de continuar",
            target  = path.name,
            action  = (
                "Penalización > 10 en ParamChem indica parámetros no confiables. "
                "Considerar optimización QM (Gaussian/ORCA) para cargas y constantes de fuerza."
            ),
        ))

    return result


# ─── Interfaz pública ────────────────────────────────────────────────────────

def validate_ligand(path: str | Path, role: str) -> LigandValidationResult:
    
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"SDF no encontrado: {path}")

    if path.suffix.lower() not in (".sdf", ".mol"):
        raise ValueError(f"El archivo debe ser SDF o MOL: {path}")

    raw = _parse_sdf(path)

    # Error de parseo
    if raw.get("error"):
        result = LigandValidationResult(
            source_file  = str(path),
            role         = role,
            is_complete  = False,
        )
        result.warnings.append(Warning(
            message  = f"Error al parsear SDF: {raw['error']}",
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.risks.append(Risk(
            message  = "Archivo SDF inválido — no se puede parametrizar",
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar y corregir el archivo SDF antes de continuar",
            target  = path.name,
            action  = (
                "Abrir en Avogadro o RDKit para verificar estructura. "
                "Re-exportar desde ChemDraw, MarvinSketch o PubChem."
            ),
        ))
        return result

    # ─── Análisis químico ─────────────────────────────────────────────────────
    chem = _analyze_chemistry(raw["atoms"], raw["bonds"])

    # ─── Construir resultado ──────────────────────────────────────────────────
    result = LigandValidationResult(
        source_file  = str(path),
        role         = role,
        is_complete  = True,
        n_atoms      = raw["n_atoms"],
        n_bonds      = raw["n_bonds"],
        atom_elements = list({a["element"] for a in raw["atoms"]}),
        atoms        = [AtomInfo(**a) for a in raw["atoms"]],
        bonds        = [BondInfo(**b) for b in raw["bonds"]],
        formal_charges          = chem["formal_charges"],
        net_charge              = chem["net_charge"],
        has_aromatic            = chem["has_aromatic"],
        aromatic_atoms          = chem["aromatic_atoms"],
        has_sulfur              = chem["has_sulfur"],
        has_phosphorus          = chem["has_phosphorus"],
        has_nitrogen            = chem["has_nitrogen"],
        has_oxygen              = chem["has_oxygen"],
        has_halogens            = chem["has_halogens"],
        halogen_types           = chem["halogen_types"],
        n_rotatable_bonds       = chem["n_rotatable_bonds"],
        estimated_flexibility   = chem["flexibility"],
        estimated_polarity      = chem["polarity"],
        parametrization_difficulty = chem["difficulty"],
    )

    # ─── Reasoning engine ────────────────────────────────────────────────────
    result = _run_reasoning(result, chem, path)

    return result
