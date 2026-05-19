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

from descriptors.topology import (
    compute_topology,
    TopologyDescriptor,
)

from descriptors.aromaticity import (
    compute_aromaticity,
    AromaticityDescriptor,
)

from descriptors.flexibility import (
    compute_flexibility,
    FlexibilityDescriptor,
)


# ─── Resultado de validación ─────────────────────────────────────────────────


class AtomInfo(BaseModel):
    index:   int
    element: str
    x:       float
    y:       float
    z:       float
    charge:  int = 0


class BondInfo(BaseModel):
    atom1:     int
    atom2:     int
    bond_type: int   # 1=single, 2=double, 3=triple, 4=aromatic


class LigandValidationResult(BaseModel):
    """
    Contrato de salida del ligand_validator.
    Siempre retorna este objeto, independientemente de la implementación interna.
    """

    source_file: str
    role:        str
    parser_used: str = "unknown"   # "sdf" | "pdb"

    # Estructura detectada
    n_atoms:        int            = 0
    n_bonds:        int            = 0
    atom_elements:  list[str]      = []
    atoms:          list[AtomInfo] = []
    bonds:          list[BondInfo] = []
    formal_charges: list[int]      = []
    net_charge:     int            = 0
    is_complete:    bool           = False

    # Características químicas inferidas
    has_aromatic:               bool      = False
    aromatic_atoms:             int       = 0
    has_sulfur:                 bool      = False
    has_phosphorus:             bool      = False
    has_nitrogen:               bool      = False
    has_oxygen:                 bool      = False
    has_halogens:               bool      = False
    halogen_types:              list[str] = []
    n_rotatable_bonds:          int       = 0
    estimated_flexibility:      str       = "unknown"
    estimated_polarity:         str       = "unknown"
    parametrization_difficulty: str       = "unknown"

    # Inferencia
    warnings:        list[Warning]        = []
    risks:           list[Risk]           = []
    recommendations: list[Recommendation] = []


# ─── Análisis químico básico ─────────────────────────────────────────────────


def _analyze_chemistry(
    atoms: list[dict],
    bonds: list[dict],
) -> dict:

    elements       = [a["element"] for a in atoms]
    element_set    = set(elements)

    formal_charges = [a["charge"] for a in atoms]
    net_charge     = sum(formal_charges)

    has_sulfur     = "S" in element_set
    has_phosphorus = "P" in element_set
    has_nitrogen   = "N" in element_set
    has_oxygen     = "O" in element_set

    halogens       = {"F", "Cl", "Br", "I"}
    halogen_types  = sorted(element_set & halogens)
    has_halogens   = len(halogen_types) > 0

    # ─── Descriptor engine ──────────────────────────────────────────────────

    topo = compute_topology(atoms, bonds)

    arom = compute_aromaticity(
        atoms,
        bonds,
        topo,
    )

    flex = compute_flexibility(
        atoms,
        bonds,
        topo,
        arom,
    )

    has_aromatic   = arom.aromatic_atom_count > 0
    aromatic_atoms = arom.aromatic_atom_count

    n_rotatable    = flex.n_rotatable_bonds
    flexibility    = flex.flexibility_class

    # ─── Polaridad ──────────────────────────────────────────────────────────

    if net_charge != 0:
        polarity = "charged"

    elif has_nitrogen or has_oxygen or has_sulfur:
        polarity = "polar"

    else:
        polarity = "nonpolar"

    # ─── Dificultad de parametrización ─────────────────────────────────────

    difficulty_score = 0

    if has_sulfur:
        difficulty_score += 2

    if has_phosphorus:
        difficulty_score += 3

    if has_halogens:
        difficulty_score += 1

    if net_charge != 0:
        difficulty_score += 2

    if has_aromatic:
        difficulty_score += 1

    if len(elements) > 50:
        difficulty_score += 2

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

        # Objetos completos para reasoning
        "_topo": topo,
        "_arom": arom,
        "_flex": flex,
    }


# ─── Reasoning engine ───────────────────────────────────────────────────────


def _run_reasoning(
    result: LigandValidationResult,
    chem: dict,
    path: Path,
) -> LigandValidationResult:
    """
    Inferencia semántica y contextual.
    No cambia cuando se reemplace el parser.
    """

    # ─── Advertencia específica para PDB ───────────────────────────────────

    if result.parser_used == "pdb":

        result.warnings.append(Warning(
            message=(
                "Conectividad inferida por distancia "
                "(PDB sin CONECT). "
                "Tipos de enlace y aromaticidad "
                "son aproximaciones heurísticas."
            ),
            target=path.name,
            severity=Severity.MEDIUM,
        ))

        result.recommendations.append(Recommendation(
            message=(
                "Convertir a SDF con RDKit u OpenBabel "
                "para conectividad exacta"
            ),
            target=path.name,
            action=(
                "obabel ligand.pdb -O ligand.sdf "
                "o usar RDKit para reconstrucción de conectividad."
            ),
        ))

    # ─── Fosfato / CoA ─────────────────────────────────────────────────────

    if chem["has_phosphorus"] and result.role == "substrate":

        result.warnings.append(Warning(
            message=(
                "Sustrato con fósforo detectado — "
                "posible CoA o nucleótido"
            ),
            target=path.name,
            severity=Severity.HIGH,
        ))

        result.risks.append(Risk(
            message=(
                "Moléculas con grupos fosfato son "
                "difíciles de parametrizar"
            ),
            target=path.name,
            severity=Severity.HIGH,
        ))

    # ─── Carga neta ────────────────────────────────────────────────────────

    if chem["net_charge"] != 0:

        result.warnings.append(Warning(
            message=(
                f"Ligando con carga neta "
                f"{chem['net_charge']:+d} detectado"
            ),
            target=path.name,
            severity=Severity.MEDIUM,
        ))

    # ─── Aromaticidad ──────────────────────────────────────────────────────

    if chem["has_aromatic"] and result.role == "competitive_ligand":

        result.recommendations.append(Recommendation(
            message=(
                f"Inhibidor aromático "
                f"({chem['aromatic_atoms']} átomos aromáticos)"
            ),
            target=path.name,
            action=(
                "Verificar posibles interacciones "
                "π-stacking en el sitio activo."
            ),
        ))

    # ─── Flexibilidad ──────────────────────────────────────────────────────

    if chem["flexibility"] in ("flexible", "very_flexible"):

        severity = (
            Severity.HIGH
            if chem["flexibility"] == "very_flexible"
            else Severity.MEDIUM
        )

        result.warnings.append(Warning(
            message=(
                f"Ligando {chem['flexibility']}: "
                f"{chem['n_rotatable_bonds']} enlaces rotables"
            ),
            target=path.name,
            severity=severity,
        ))

        if chem["flexibility"] == "very_flexible":

            result.risks.append(Risk(
                message=(
                    "MD estándar insuficiente — "
                    "requiere REST2 o metadinámica"
                ),
                target=path.name,
                severity=Severity.HIGH,
            ))

        result.recommendations.append(Recommendation(
            message=chem["_flex"].sampling_recommendation,
            target=path.name,
            action=(
                "Ver flexibility descriptor "
                "para detalle de enlaces rotables"
            ),
        ))

    # ─── Halógenos ─────────────────────────────────────────────────────────

    if chem["has_halogens"]:

        result.warnings.append(Warning(
            message=(
                f"Halógenos detectados: "
                f"{chem['halogen_types']}"
            ),
            target=path.name,
            severity=Severity.LOW,
        ))

    # ─── Dificultad alta ───────────────────────────────────────────────────

    if chem["difficulty"] == "high":

        result.risks.append(Risk(
            message=(
                "Dificultad de parametrización alta — "
                "requiere validación manual"
            ),
            target=path.name,
            severity=Severity.HIGH,
        ))

    return result


# ─── Helpers de error uniforme ──────────────────────────────────────────────


def _error_result(
    path: Path,
    role: str,
    parser: str,
    message: str,
) -> LigandValidationResult:

    result = LigandValidationResult(
        source_file=str(path),
        role=role,
        parser_used=parser,
        is_complete=False,
    )

    result.warnings.append(Warning(
        message=f"Error al parsear archivo: {message}",
        target=path.name,
        severity=Severity.HIGH,
    ))

    result.risks.append(Risk(
        message="Archivo inválido — no se puede parametrizar",
        target=path.name,
        severity=Severity.HIGH,
    ))

    return result


# ─── Interfaz pública ───────────────────────────────────────────────────────


def validate_ligand(
    path: str | Path,
    role: str,
) -> LigandValidationResult:
    """
    Valida un archivo de ligando (SDF, MOL, o PDB).
    """

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Archivo no encontrado: {path}")

    ext = path.suffix.lower()

    # ─── Selección de parser ───────────────────────────────────────────────

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

    # ─── Error de parseo ───────────────────────────────────────────────────

    if raw.get("error"):

        return _error_result(
            path,
            role,
            parser_name,
            raw["error"],
        )

    # ─── Análisis químico ──────────────────────────────────────────────────

    chem = _analyze_chemistry(
        raw["atoms"],
        raw["bonds"],
    )

    # ─── Construir resultado ───────────────────────────────────────────────

    result = LigandValidationResult(
        source_file=str(path),
        role=role,
        parser_used=parser_name,

        is_complete=True,

        n_atoms=raw["n_atoms"],
        n_bonds=raw["n_bonds"],

        atom_elements=list({
            a["element"]
            for a in raw["atoms"]
        }),

        atoms=[
            AtomInfo(**a)
            for a in raw["atoms"]
        ],

        bonds=[
            BondInfo(**b)
            for b in raw["bonds"]
        ],

        formal_charges=chem["formal_charges"],
        net_charge=chem["net_charge"],

        has_aromatic=chem["has_aromatic"],
        aromatic_atoms=chem["aromatic_atoms"],

        has_sulfur=chem["has_sulfur"],
        has_phosphorus=chem["has_phosphorus"],
        has_nitrogen=chem["has_nitrogen"],
        has_oxygen=chem["has_oxygen"],

        has_halogens=chem["has_halogens"],
        halogen_types=chem["halogen_types"],

        n_rotatable_bonds=chem["n_rotatable_bonds"],

        estimated_flexibility=chem["flexibility"],
        estimated_polarity=chem["polarity"],

        parametrization_difficulty=chem["difficulty"],
    )

    # ─── Reasoning engine ─────────────────────────────────────────────────

    result = _run_reasoning(
        result,
        chem,
        path,
    )

    return result