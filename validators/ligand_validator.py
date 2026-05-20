# validators/ligand_validator.py
"""
Validator estructural y semántico de archivos de ligandos (SDF y PDB).

Arquitectura:
    ligand_validator.py        ← este archivo: interfaz pública + reasoning engine
    ligand_parsers/
        sdf_parser.py          ← _parse_sdf() extraído, reemplazable por RDKit
        pdb_parser.py          ← conectividad por distancia, reemplazable por RDKit

Interfaz pública garantizada:
    validate_ligand(path: str | Path, role: str) -> LigandValidationResult

Para agregar RDKit = reemplazar el parser en ligand_parsers/, no tocar este archivo.
El reasoning engine (inferencia semántica) nunca cambia con el parser.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field
from core.models import Warning, Risk, Recommendation, Severity

from validators.ligand_parsers.sdf_parser import parse_sdf
from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

# Tipos de descriptors — solo para anotaciones en LigandValidationResult
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from descriptors.topology    import TopologyDescriptor
    from descriptors.aromaticity import AromaticityDescriptor
    from descriptors.flexibility import FlexibilityDescriptor
    from descriptors.geometry    import GeometryDescriptor
    from descriptors.polarity    import PolarityDescriptor


# ─── Resultado de validación ─────────────────────────────────────────────────

class AtomInfo(BaseModel):
    index:       int
    element:     str
    x:           float
    y:           float
    z:           float
    charge:      int = 0


class BondInfo(BaseModel):
    atom1:      int
    atom2:      int
    bond_type:  int   # 1=single, 2=double, 3=triple, 4=aromatic


class LigandValidationResult(BaseModel):
    """
    Contrato de salida del ligand_validator.
    Siempre retorna este objeto, independientemente de la implementación interna.
    """
    source_file:    str
    role:           str
    parser_used:    str = "unknown"   # "sdf" | "pdb"

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
    has_aromatic:              bool       = False
    aromatic_atoms:            int        = 0
    has_sulfur:                bool       = False
    has_phosphorus:            bool       = False
    has_nitrogen:              bool       = False
    has_oxygen:                bool       = False
    has_halogens:              bool       = False
    halogen_types:             list[str]  = []
    n_rotatable_bonds:         int        = 0
    estimated_flexibility:     str        = "unknown"  # rigid / moderate / flexible
    estimated_polarity:        str        = "unknown"  # nonpolar / polar / charged
    parametrization_difficulty: str       = "unknown"  # low / medium / high

    # Descriptors completos (disponibles para reasoning avanzado)
    # Se usan Optional[Any] para evitar dependencia circular en imports
    shape_class:           str   = "unknown"   # flat / elongated / globular
    is_globally_planar:    bool  = False
    global_planarity_rms:  float = 0.0
    radius_of_gyration:    float = 0.0
    n_chiral_centers:      int   = 0
    logp_estimate:         float = 0.0
    logp_class:            str   = "unknown"
    hbd_count:             int   = 0
    hba_count:             int   = 0
    lipinski_compliant:    bool  = True
    functional_groups:     list[str] = []      # nombres de grupos detectados
    estimated_solubility:  str   = "unknown"
    n_aromatic_rings:      int   = 0
    n_fused_aromatic:      int   = 0
    scaffold_rigidity:     str   = "unknown"
    sampling_recommendation: str = ""

    # Inferencia
    warnings:        list[Warning]        = []
    risks:           list[Risk]           = []
    recommendations: list[Recommendation] = []


# ─── Análisis químico — descriptor engine ────────────────────────────────────

def _analyze_chemistry(atoms: list[dict], bonds: list[dict]) -> dict:
    from descriptors.topology    import compute_topology
    from descriptors.aromaticity import compute_aromaticity
    from descriptors.flexibility import compute_flexibility
    from descriptors.geometry    import compute_geometry
    from descriptors.polarity    import compute_polarity

    elements    = [a["element"] for a in atoms]
    element_set = set(elements)

    formal_charges = [a["charge"] for a in atoms]
    net_charge     = sum(formal_charges)

    has_sulfur     = "S"  in element_set
    has_phosphorus = "P"  in element_set
    has_nitrogen   = "N"  in element_set
    has_oxygen     = "O"  in element_set
    halogens       = {"F", "Cl", "Br", "I"}
    halogen_types  = sorted(element_set & halogens)
    has_halogens   = len(halogen_types) > 0

    # ── Descriptor engine ─────────────────────────────────────────────────────
    topo = compute_topology(atoms, bonds)
    arom = compute_aromaticity(atoms, bonds, topo)
    flex = compute_flexibility(atoms, bonds, topo, arom)
    geom = compute_geometry(atoms, bonds, topo, arom)
    pol  = compute_polarity(atoms, bonds, topo)

    has_aromatic   = arom.aromatic_atom_count > 0
    aromatic_atoms = arom.aromatic_atom_count
    n_rotatable    = flex.n_rotatable_bonds
    flexibility    = flex.flexibility_class
    polarity       = pol.polarity_class

    # ── Dificultad de parametrización ─────────────────────────────────────────
    difficulty_score = 0
    if has_sulfur:          difficulty_score += 2
    if has_phosphorus:      difficulty_score += 3
    if has_halogens:        difficulty_score += 1
    if net_charge != 0:     difficulty_score += 2
    if has_aromatic:        difficulty_score += 1
    if len(elements) > 50:  difficulty_score += 2

    if difficulty_score <= 2:   difficulty = "low"
    elif difficulty_score <= 5: difficulty = "medium"
    else:                       difficulty = "high"

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
        # Descriptors completos — disponibles para reasoning engine
        "_topo": topo,
        "_arom": arom,
        "_flex": flex,
        "_geom": geom,
        "_pol":  pol,
    }


# ─── Reasoning engine ────────────────────────────────────────────────────────

def _run_reasoning(
    result: LigandValidationResult,
    chem:   dict,
    path:   Path,
) -> LigandValidationResult:
    """
    Inferencia semántica y contextual.
    No cambia cuando se reemplace el parser.
    """

    # ─── Advertencia específica para PDB (conectividad inferida) ─────────────
    if result.parser_used == "pdb":
        result.warnings.append(Warning(
            message  = (
                "Conectividad inferida por distancia (PDB sin CONECT). "
                "Tipos de enlace y aromaticidad son aproximaciones heurísticas."
            ),
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Convertir a SDF con RDKit o OpenBabel para conectividad exacta",
            target  = path.name,
            action  = (
                "obabel ligand.pdb -O ligand.sdf  "
                "o bien: python -c \"from rdkit import Chem; "
                "Chem.MolToMolFile(Chem.MolFromPDBFile('ligand.pdb'), 'ligand.sdf')\""
            ),
        ))

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

    # ─── Halógenos ────────────────────────────────────────────────────────────
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

    # ─── Flexibilidad very_flexible ──────────────────────────────────────────
    if chem["flexibility"] == "very_flexible":
        result.risks.append(Risk(
            message  = (
                f"Ligando muy flexible ({chem['n_rotatable_bonds']} rot. bonds): "
                "MD estándar probablemente insuficiente"
            ),
            target   = path.name,
            severity = Severity.HIGH,
        ))
        result.recommendations.append(Recommendation(
            message = "Considerar REST2 o metadinámica para muestreo conformacional",
            target  = path.name,
            action  = (
                "Con >14 enlaces rotables el espacio conformacional es muy grande. "
                "Opciones: REST2 en GROMACS, metadinámica con PLUMED, "
                "o clustering de poses de docking antes de MD."
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

    # ─── Geometry: planaridad ─────────────────────────────────────────────────
    geom = chem["_geom"]
    if geom.is_globally_planar and result.role == "competitive_ligand":
        result.recommendations.append(Recommendation(
            message = (
                f"Molécula plana (RMS={geom.global_planarity_rms:.3f}Å, "
                f"shape={geom.shape_class}): alta probabilidad de π-stacking e intercalación"
            ),
            target  = path.name,
            action  = (
                "Verificar orientación en sitio activo: las moléculas planas tienden a "
                "intercalarse entre residuos aromáticos. Analizar con distance_analysis "
                "y contact_map en la simulación."
            ),
        ))

    if geom.n_chiral_centers > 0:
        result.warnings.append(Warning(
            message  = f"{geom.n_chiral_centers} centro(s) quiral(es) detectado(s)",
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar estereoquímica antes de parametrizar",
            target  = path.name,
            action  = (
                "Confirmar configuración R/S en cada centro quiral. "
                "CGenFF y GAFF son sensibles a la geometría 3D del input — "
                "un enantiómero incorrecto puede producir parámetros erróneos."
            ),
        ))

    # ─── Polarity: Lipinski y solubilidad ─────────────────────────────────────
    pol = chem["_pol"]
    if not pol.lipinski_compliant:
        result.warnings.append(Warning(
            message  = (
                f"Molécula fuera de regla de Lipinski "
                f"(logP={pol.logp_estimate}, HBD={pol.hbd_count}, HBA={pol.hba_count})"
            ),
            target   = path.name,
            severity = Severity.MEDIUM,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar drug-likeness y biodisponibilidad oral estimada",
            target  = path.name,
            action  = (
                "La violación de Ro5 no invalida la simulación, pero indica que "
                "la molécula puede tener problemas de solubilidad o permeabilidad. "
                "Considerar co-solventes (DMSO) si la solubilidad en agua es baja."
            ),
        ))

    if pol.estimated_solubility == "poor":
        result.warnings.append(Warning(
            message  = f"Solubilidad estimada baja (logP={pol.logp_estimate:.2f})",
            target   = path.name,
            severity = Severity.LOW,
        ))
        result.recommendations.append(Recommendation(
            message = "Verificar concentración y condiciones de solvatación",
            target  = path.name,
            action  = (
                "Moléculas lipofílicas pueden agregar en agua pura. "
                "Considerar añadir co-solvente o usar caja de simulación más grande."
            ),
        ))

    if pol.functional_groups:
        fg_str = ", ".join(
        sorted({fg.name for fg in pol.functional_groups}))
        result.recommendations.append(Recommendation(
            message = f"Grupos funcionales detectados: {fg_str}",
            target  = path.name,
            action  = (
                "Verificar que cada grupo funcional tiene parámetros correctos. "
                "Grupos como fosfato, carboxilo y aminas primarias requieren "
                "atención especial en la asignación de cargas."
            ),
        ))

    return result


# ─── Helpers de error uniforme ────────────────────────────────────────────────

def _error_result(path: Path, role: str, parser: str, message: str) -> LigandValidationResult:
    result = LigandValidationResult(
        source_file = str(path),
        role        = role,
        parser_used = parser,
        is_complete = False,
    )
    result.warnings.append(Warning(
        message  = f"Error al parsear archivo: {message}",
        target   = path.name,
        severity = Severity.HIGH,
    ))
    result.risks.append(Risk(
        message  = "Archivo inválido — no se puede parametrizar",
        target   = path.name,
        severity = Severity.HIGH,
    ))
    result.recommendations.append(Recommendation(
        message = "Verificar y corregir el archivo antes de continuar",
        target  = path.name,
        action  = (
            "Abrir en Avogadro o RDKit para verificar estructura. "
            "Re-exportar desde ChemDraw, MarvinSketch o PubChem."
        ),
    ))
    return result


# ─── Interfaz pública ────────────────────────────────────────────────────────

def validate_ligand(path: str | Path, role: str) -> LigandValidationResult:
    """
    Valida un archivo de ligando (SDF, MOL, o PDB).

    Selección de parser:
        .sdf / .mol → sdf_parser.parse_sdf()
        .pdb        → pdb_parser.parse_pdb_ligand()

    Retorna siempre LigandValidationResult.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Archivo no encontrado: {path}")

    ext = path.suffix.lower()

    # ─── Selección de parser ──────────────────────────────────────────────────
    if ext in (".sdf", ".mol"):
        parser_name = "sdf"
        raw = parse_sdf(path)
    elif ext == ".pdb":
        parser_name = "pdb"
        raw = parse_pdb_ligand(path)
    else:
        raise ValueError(
            f"Formato no soportado: {ext}. "
            f"Soportados: .sdf, .mol, .pdb"
        )

    # ─── Error de parseo ──────────────────────────────────────────────────────
    if raw.get("error"):
        return _error_result(path, role, parser_name, raw["error"])

    # ─── Análisis químico ─────────────────────────────────────────────────────
    chem = _analyze_chemistry(raw["atoms"], raw["bonds"])

    # ─── Construir resultado ──────────────────────────────────────────────────
    geom = chem["_geom"]
    pol  = chem["_pol"]
    arom = chem["_arom"]
    flex = chem["_flex"]

    result = LigandValidationResult(
        source_file  = str(path),
        role         = role,
        parser_used  = parser_name,
        is_complete  = True,
        n_atoms      = raw["n_atoms"],
        n_bonds      = raw["n_bonds"],
        atom_elements = list({a["element"] for a in raw["atoms"]}),
        atoms         = [AtomInfo(**a) for a in raw["atoms"]],
        bonds         = [BondInfo(**b) for b in raw["bonds"]],
        formal_charges            = chem["formal_charges"],
        net_charge                = chem["net_charge"],
        has_aromatic              = chem["has_aromatic"],
        aromatic_atoms            = chem["aromatic_atoms"],
        has_sulfur                = chem["has_sulfur"],
        has_phosphorus            = chem["has_phosphorus"],
        has_nitrogen              = chem["has_nitrogen"],
        has_oxygen                = chem["has_oxygen"],
        has_halogens              = chem["has_halogens"],
        halogen_types             = chem["halogen_types"],
        n_rotatable_bonds         = chem["n_rotatable_bonds"],
        estimated_flexibility     = chem["flexibility"],
        estimated_polarity        = chem["polarity"],
        parametrization_difficulty = chem["difficulty"],
        # Geometry
        shape_class              = geom.shape_class,
        is_globally_planar       = geom.is_globally_planar,
        global_planarity_rms     = geom.global_planarity_rms,
        radius_of_gyration       = geom.radius_of_gyration,
        n_chiral_centers         = geom.n_chiral_centers,
        # Polarity
        logp_estimate            = pol.logp_estimate,
        logp_class               = pol.logp_class,
        hbd_count                = pol.hbd_count,
        hba_count                = pol.hba_count,
        lipinski_compliant       = pol.lipinski_compliant,
        functional_groups        = pol.group_names,
        estimated_solubility     = pol.estimated_solubility,
        # Aromaticity / flexibility summary
        n_aromatic_rings         = arom.n_aromatic_rings,
        n_fused_aromatic         = arom.n_fused_aromatic,
        scaffold_rigidity        = flex.scaffold_rigidity,
        sampling_recommendation  = flex.sampling_recommendation,
    )

    # ─── Reasoning engine ────────────────────────────────────────────────────
    result = _run_reasoning(result, chem, path)

    return result
