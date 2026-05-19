# core/models.py

from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, field_validator

from core.ontology import (
    COMPONENT_ROLES,
    FORCEFIELDS,
    WATER_MODELS,
    BIOLOGICAL_CONTEXT,
    RESTRAINT_TYPES,
    ANALYSIS_TYPES,
    OUTPUT_FORMATS,
    SIMULATION_GOALS,
)

from enum import Enum

class Severity(str, Enum):
    LOW    = "low"
    MEDIUM = "medium"
    HIGH   = "high"


class Warning(BaseModel):
    message:  str
    target:   Optional[str] = None
    severity: Severity = Severity.MEDIUM


class Risk(BaseModel):
    message:  str
    target:   Optional[str] = None
    severity: Severity = Severity.HIGH


class Recommendation(BaseModel):
    message:  str
    target:   Optional[str] = None
    action:   Optional[str] = None

class ComponentModel(BaseModel):
    id: str
    role: str
    file: str
    biological_context: list[str] = []

    @field_validator("role")
    @classmethod
    def validate_role(cls, v):
        if v not in COMPONENT_ROLES:
            raise ValueError(f"Role '{v}' no reconocido. Válidos: {COMPONENT_ROLES}")
        return v

    @field_validator("biological_context")
    @classmethod
    def validate_bio_context(cls, v):
        for ctx in v:
            if ctx not in BIOLOGICAL_CONTEXT:
                raise ValueError(f"Contexto biológico '{ctx}' no reconocido. Válidos: {BIOLOGICAL_CONTEXT}")
        return v


class MembraneConfig(BaseModel):
    enabled: bool = False
    type: Optional[str] = None
    thickness: Optional[float] = None


class SolventConfig(BaseModel):
    water_model: str = "tip3p"

    @field_validator("water_model")
    @classmethod
    def validate_water(cls, v):
        if v not in WATER_MODELS:
            raise ValueError(f"Modelo de agua '{v}' no reconocido. Válidos: {WATER_MODELS}")
        return v


class IonsConfig(BaseModel):
    concentration: float = 0.15
    positive: str = "NA"
    negative: str = "CL"


class EnvironmentModel(BaseModel):
    membrane: MembraneConfig = MembraneConfig()
    solvent: SolventConfig = SolventConfig()
    ions: IonsConfig = IonsConfig()


class ForcefieldsModel(BaseModel):
    protein: str
    ligands: Optional[str] = None
    coarse_grained: Optional[str] = None

    @field_validator("protein")
    @classmethod
    def validate_protein_ff(cls, v):
        valid = FORCEFIELDS["atomistic"]
        if v not in valid:
            raise ValueError(f"FF de proteína '{v}' no reconocido. Válidos: {valid}")
        return v

    @field_validator("ligands")
    @classmethod
    def validate_ligand_ff(cls, v):
        if v is None:
            return v
        valid = FORCEFIELDS["ligands"]
        if v not in valid:
            raise ValueError(f"FF de ligando '{v}' no reconocido. Válidos: {valid}")
        return v


class RestraintModel(BaseModel):
    type: str
    target: str
    force_constant: Optional[float] = None

    @field_validator("type")
    @classmethod
    def validate_restraint(cls, v):
        if v not in RESTRAINT_TYPES:
            raise ValueError(f"Tipo de restraint '{v}' no reconocido. Válidos: {RESTRAINT_TYPES}")
        return v


class AnalysisModel(BaseModel):
    type: str
    selection: Optional[dict] = None

    @field_validator("type")
    @classmethod
    def validate_analysis(cls, v):
        if v not in ANALYSIS_TYPES:
            raise ValueError(f"Análisis '{v}' no reconocido. Válidos: {ANALYSIS_TYPES}")
        return v


class FiguresConfig(BaseModel):
    format: str = "svg"
    style: str = "publication"

    @field_validator("format")
    @classmethod
    def validate_format(cls, v):
        if v not in OUTPUT_FORMATS:
            raise ValueError(f"Formato '{v}' no reconocido. Válidos: {OUTPUT_FORMATS}")
        return v


class ReportsConfig(BaseModel):
    format: str = "pdf"


class OutputModel(BaseModel):
    figures: FiguresConfig = FiguresConfig()
    reports: ReportsConfig = ReportsConfig()


class ProjectModel(BaseModel):
    name: str
    description: Optional[str] = None


class SystemState(BaseModel):
    """
    Representación computacional del experimento molecular.
    Contrato entre todos los módulos de la plataforma.
    Solo estructura y validación de campos.
    La inferencia biológica vive en core/inference.py
    """
    project: ProjectModel
    components: list[ComponentModel]
    environment: EnvironmentModel = EnvironmentModel()
    forcefields: ForcefieldsModel
    simulation_objectives: list[str] = []
    restraints: list[RestraintModel] = []
    analysis: list[AnalysisModel] = []
    output: OutputModel = OutputModel()

    # Estado interno — enriquecido por inference.py
    warnings: list[Warning] = []
    risks: list[Risk] = []
    recommendations: list[Recommendation] = []
    inferred_system_type: Optional[str] = None

    @field_validator("simulation_objectives")
    @classmethod
    def validate_objectives(cls, v):
        for obj in v:
            if obj not in SIMULATION_GOALS:
                raise ValueError(f"Objetivo '{obj}' no reconocido. Válidos: {SIMULATION_GOALS}")
        return v

    # ─── Helpers de consulta ──────────────────────────────────────────────────
    def get_components_by_role(self, role: str) -> list[ComponentModel]:
        return [c for c in self.components if c.role == role]

    def has_membrane(self) -> bool:
        return self.environment.membrane.enabled

    def has_biological_context(self, context: str) -> bool:
        return any(context in c.biological_context for c in self.components)

    def component_ids(self) -> list[str]:
        return [c.id for c in self.components]