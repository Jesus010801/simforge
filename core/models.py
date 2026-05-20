# core/models.py
"""
Modelos de datos del SystemState.

Refactor arquitectónico: los resultados de validación y descriptores
ahora viven DENTRO del SystemState, en cada ComponentModel.

Árbol de datos:

    SystemState
    │
    ├── project           ProjectModel
    ├── forcefields       ForcefieldsModel
    ├── environment       EnvironmentModel
    ├── simulation_objectives
    ├── restraints
    ├── analysis
    ├── output
    │
    ├── components        list[ComponentModel]
    │    └── ComponentModel
    │         ├── config       (id, role, file, biological_context)  ← del YAML
    │         ├── validation   ComponentValidation | None            ← validators/
    │         ├── descriptors  ComponentDescriptors | None           ← descriptors/
    │         └── reasoning    ComponentReasoning | None             ← inference/
    │
    └── global_reasoning  GlobalReasoning                            ← inference.py

Principios:
    - Pydantic solo valida estructura, NO infiere
    - Inferencia y validación viven en sus módulos propios
    - Los validators escriben en component.validation
    - Los descriptors escriben en component.descriptors
    - El reasoning engine escribe en component.reasoning y global_reasoning
    - El decision engine lee TODO desde SystemState — sin acceso a archivos
"""

from __future__ import annotations
from typing import Optional, Any
from pydantic import BaseModel, field_validator, Field

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


# ═══════════════════════════════════════════════════════════════════════════════
# Primitivos compartidos
# ═══════════════════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════════════════
# Resultado de validación por componente
# ═══════════════════════════════════════════════════════════════════════════════

class ComponentValidation(BaseModel):
    """
    Envelope genérico que almacena el resultado de cualquier validator.

    Cada validator (protein_validator, ligand_validator) escribe aquí
    su resultado completo. El campo `data` acepta el dict serializado
    del resultado Pydantic del validator (model_dump()).

    El decision engine lee warnings/risks/recommendations directamente.
    Para acceder al resultado tipado completo, usa el campo `data`
    y reconstruye con el modelo específico si lo necesita.
    """
    validator_used:   str              # "protein_validator" | "ligand_validator" | ...
    is_valid:         bool   = False   # el archivo es parseable y completo
    validation_error: Optional[str] = None   # mensaje si el parser falló

    # Issues agregados (subset de los que reportó el validator)
    warnings:         list[Warning]        = []
    risks:            list[Risk]           = []
    recommendations:  list[Recommendation] = []

    # Resultado completo serializado (dict) — preserva todos los campos del validator
    # sin crear dependencia circular en models.py hacia validators/
    data:             dict = Field(default_factory=dict)

    class Config:
        arbitrary_types_allowed = True


# ═══════════════════════════════════════════════════════════════════════════════
# Descriptores moleculares por componente
# ═══════════════════════════════════════════════════════════════════════════════

class ComponentDescriptors(BaseModel):
    """
    Envelope genérico para descriptores moleculares calculados.

    Para ligandos: topology, aromaticity, flexibility, geometry, polarity.
    Para proteínas: reservado (MDAnalysis futuro).

    Igual que ComponentValidation, usa `data` para no crear dependencias
    circulares. Los campos de resumen permiten al decision engine operar
    sin deserializar el dict completo.
    """
    # Resumen ejecutivo — campos planos para el decision engine
    n_heavy_atoms:        int   = 0
    n_rings:              int   = 0
    n_aromatic_rings:     int   = 0
    n_fused_aromatic:     int   = 0
    flexibility_class:    str   = "unknown"
    n_rotatable_bonds:    int   = 0
    effective_rot_bonds:  float = 0.0
    scaffold_rigidity:    str   = "unknown"
    shape_class:          str   = "unknown"
    is_globally_planar:   bool  = False
    global_planarity_rms: float = 0.0
    radius_of_gyration:   float = 0.0
    polarity_class:       str   = "unknown"
    polarity_score:       float = 0.0
    hb_donors:            int   = 0
    hb_acceptors:         int   = 0
    net_charge:           int   = 0
    is_zwitterion:        bool  = False
    amphipathic_class:    str   = "unknown"
    n_functional_groups:  int   = 0
    lipinski_hbd:         int   = 0
    lipinski_hba:         int   = 0
    passes_lipinski:      bool  = True
    sampling_recommendation: str = ""

    # Resultados completos serializados por módulo
    topology:    dict = Field(default_factory=dict)
    aromaticity: dict = Field(default_factory=dict)
    flexibility: dict = Field(default_factory=dict)
    geometry:    dict = Field(default_factory=dict)
    polarity:    dict = Field(default_factory=dict)

    class Config:
        arbitrary_types_allowed = True


# ═══════════════════════════════════════════════════════════════════════════════
# Reasoning por componente
# ═══════════════════════════════════════════════════════════════════════════════

class ComponentReasoning(BaseModel):
    """
    Interpretación contextual de un componente dentro del sistema.

    El reasoning engine escribe aquí las conclusiones sobre el componente
    en el contexto del sistema completo (tipo de sistema, forcefields, rol).
    Esto es diferente de los warnings del validator — esos son estructurales;
    estos son contextuales (ej: "ligando flexible en sistema competitive-inhibition
    requiere docking exhaustivo antes de MD").
    """
    component_id:    str
    component_role:  str

    # Conclusiones contextuales
    warnings:        list[Warning]        = []
    risks:           list[Risk]           = []
    recommendations: list[Recommendation] = []

    # Flags semánticos de alto nivel (leídos por decision engine)
    needs_special_sampling:  bool = False
    needs_parametrization_review: bool = False
    needs_protonation_check: bool = False
    needs_pose_validation:   bool = False
    parametrization_difficulty: str = "unknown"   # low / medium / high

    # Notas libres del reasoning engine
    notes: list[str] = []


# ═══════════════════════════════════════════════════════════════════════════════
# ComponentModel — config + resultados
# ═══════════════════════════════════════════════════════════════════════════════

class ComponentModel(BaseModel):
    """
    Representa un componente del sistema molecular.

    ┌── Config (del YAML) ──────────────────────────────────────────┐
    │  id, role, file, biological_context                           │
    └───────────────────────────────────────────────────────────────┘
    ┌── Validation (escrito por validators/) ───────────────────────┐
    │  ComponentValidation | None                                    │
    └───────────────────────────────────────────────────────────────┘
    ┌── Descriptors (escrito por descriptors/) ─────────────────────┐
    │  ComponentDescriptors | None                                   │
    └───────────────────────────────────────────────────────────────┘
    ┌── Reasoning (escrito por inference.py) ───────────────────────┐
    │  ComponentReasoning | None                                     │
    └───────────────────────────────────────────────────────────────┘
    """

    # ── Config (del YAML, inmutable después de parse) ─────────────────────────
    id:                 str
    role:               str
    file:               str
    biological_context: list[str] = []

    # ── Resultados (escritos por el pipeline, opcionales hasta que se calculen) ─
    validation:  Optional[ComponentValidation]  = None
    descriptors: Optional[ComponentDescriptors] = None
    reasoning:   Optional[ComponentReasoning]   = None

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
                raise ValueError(
                    f"Contexto biológico '{ctx}' no reconocido. Válidos: {BIOLOGICAL_CONTEXT}"
                )
        return v

    # ── Helpers de acceso rápido ──────────────────────────────────────────────

    @property
    def is_validated(self) -> bool:
        return self.validation is not None and self.validation.is_valid

    @property
    def has_descriptors(self) -> bool:
        return self.descriptors is not None and bool(self.descriptors.topology)

    @property
    def all_warnings(self) -> list[Warning]:
        """Todos los warnings del componente (validación + reasoning)."""
        ws = []
        if self.validation:
            ws.extend(self.validation.warnings)
        if self.reasoning:
            ws.extend(self.reasoning.warnings)
        return ws

    @property
    def all_risks(self) -> list[Risk]:
        rs = []
        if self.validation:
            rs.extend(self.validation.risks)
        if self.reasoning:
            rs.extend(self.reasoning.risks)
        return rs

    @property
    def all_recommendations(self) -> list[Recommendation]:
        recs = []
        if self.validation:
            recs.extend(self.validation.recommendations)
        if self.reasoning:
            recs.extend(self.reasoning.recommendations)
        return recs

    class Config:
        arbitrary_types_allowed = True


# ═══════════════════════════════════════════════════════════════════════════════
# Reasoning global del sistema
# ═══════════════════════════════════════════════════════════════════════════════

class GlobalReasoning(BaseModel):
    """
    Razonamiento a nivel de sistema completo.

    Contiene las conclusiones del pipeline de inferencia que emergen
    de la combinación de todos los componentes (ej: compatibilidad de
    forcefields, riesgos de sistema competitive-inhibition, gaps de análisis).

    El decision engine lee este objeto para generar el plan de simulación.
    """
    inferred_system_type: Optional[str] = None

    warnings:        list[Warning]        = []
    risks:           list[Risk]           = []
    recommendations: list[Recommendation] = []

    # Flags de decisión para el decision engine
    system_is_ready:        bool = False   # True cuando todos los componentes validan
    has_blocking_errors:    bool = False   # True cuando hay riesgos HIGH que bloquean
    needs_membrane:         bool = False
    needs_special_sampling: bool = False   # algún ligando muy flexible
    n_components_validated: int  = 0
    n_components_with_errors: int = 0

    notes: list[str] = []


# ═══════════════════════════════════════════════════════════════════════════════
# Modelos de configuración (sin cambios respecto al original)
# ═══════════════════════════════════════════════════════════════════════════════

class MembraneConfig(BaseModel):
    enabled:   bool           = False
    type:      Optional[str]  = None
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
    positive:      str   = "NA"
    negative:      str   = "CL"


class EnvironmentModel(BaseModel):
    membrane: MembraneConfig = MembraneConfig()
    solvent:  SolventConfig  = SolventConfig()
    ions:     IonsConfig     = IonsConfig()


class ForcefieldsModel(BaseModel):
    protein:        str
    ligands:        Optional[str] = None
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
    type:           str
    target:         str
    force_constant: Optional[float] = None

    @field_validator("type")
    @classmethod
    def validate_restraint(cls, v):
        if v not in RESTRAINT_TYPES:
            raise ValueError(f"Tipo de restraint '{v}' no reconocido. Válidos: {RESTRAINT_TYPES}")
        return v


class AnalysisModel(BaseModel):
    type:      str
    selection: Optional[dict] = None

    @field_validator("type")
    @classmethod
    def validate_analysis(cls, v):
        if v not in ANALYSIS_TYPES:
            raise ValueError(f"Análisis '{v}' no reconocido. Válidos: {ANALYSIS_TYPES}")
        return v


class FiguresConfig(BaseModel):
    format: str = "svg"
    style:  str = "publication"

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
    name:        str
    description: Optional[str] = None


# ═══════════════════════════════════════════════════════════════════════════════
# SystemState — contrato global
# ═══════════════════════════════════════════════════════════════════════════════

class SystemState(BaseModel):
    """
    Representación computacional del experimento molecular.

    Contrato entre todos los módulos de SimForge.
    Solo estructura y validación de campos — la inferencia vive en inference.py,
    la validación en validators/, los descriptores en descriptors/.

    El pipeline de parse_yaml() construye este objeto en etapas:
        1. YAML → config (project, components[config], forcefields, ...)
        2. run_inference() → global_reasoning.inferred_system_type
        3. run_validation() → component.validation para cada componente
        4. run_descriptors() → component.descriptors para ligandos
        5. run_component_reasoning() → component.reasoning
        6. run_global_reasoning() → global_reasoning (flags finales)

    Invariante: el decision engine solo lee desde SystemState,
    nunca llama a validators ni a descriptores directamente.
    """

    # ── Config del YAML ───────────────────────────────────────────────────────
    project:               ProjectModel
    components:            list[ComponentModel]
    environment:           EnvironmentModel  = EnvironmentModel()
    forcefields:           ForcefieldsModel
    simulation_objectives: list[str]         = []
    restraints:            list[RestraintModel] = []
    analysis:              list[AnalysisModel]  = []
    output:                OutputModel       = OutputModel()

    # ── Reasoning global (escrito por inference.py) ───────────────────────────
    global_reasoning: GlobalReasoning = Field(default_factory=GlobalReasoning)

    # ── Compatibilidad hacia atrás con código existente ───────────────────────
    # Estos campos se mantienen para no romper inference.py ni test_parser.py.
    # Son aliases que reflejan global_reasoning.* pero se escriben directamente.
    # En el futuro se pueden eliminar cuando inference.py use global_reasoning.
    warnings:             list[Warning]        = []
    risks:                list[Risk]           = []
    recommendations:      list[Recommendation] = []
    inferred_system_type: Optional[str]        = None

    @field_validator("simulation_objectives")
    @classmethod
    def validate_objectives(cls, v):
        for obj in v:
            if obj not in SIMULATION_GOALS:
                raise ValueError(f"Objetivo '{obj}' no reconocido. Válidos: {SIMULATION_GOALS}")
        return v

    # ── Helpers de consulta ───────────────────────────────────────────────────

    def get_components_by_role(self, role: str) -> list[ComponentModel]:
        return [c for c in self.components if c.role == role]

    def get_component(self, component_id: str) -> Optional[ComponentModel]:
        for c in self.components:
            if c.id == component_id:
                return c
        return None

    def has_membrane(self) -> bool:
        return self.environment.membrane.enabled

    def has_biological_context(self, context: str) -> bool:
        return any(context in c.biological_context for c in self.components)

    def component_ids(self) -> list[str]:
        return [c.id for c in self.components]

    # ── Helpers de estado del pipeline ───────────────────────────────────────

    def all_validated(self) -> bool:
        """True si todos los componentes con archivo tienen validación completa."""
        return all(
            c.validation is not None and c.validation.is_valid
            for c in self.components
            if c.file
        )

    def has_high_risks(self) -> bool:
        """True si existe algún riesgo HIGH en cualquier componente o global."""
        for r in self.risks:
            if r.severity == Severity.HIGH:
                return True
        for c in self.components:
            for r in c.all_risks:
                if r.severity == Severity.HIGH:
                    return True
        return False

    def collect_all_warnings(self) -> list[tuple[str, Warning]]:
        """
        Colección plana de todos los warnings del sistema.
        Retorna lista de (source_id, Warning).
        """
        result: list[tuple[str, Warning]] = []
        for w in self.warnings:
            result.append(("system", w))
        for c in self.components:
            for w in c.all_warnings:
                result.append((c.id, w))
        return result

    def collect_all_risks(self) -> list[tuple[str, Risk]]:
        result: list[tuple[str, Risk]] = []
        for r in self.risks:
            result.append(("system", r))
        for c in self.components:
            for r in c.all_risks:
                result.append((c.id, r))
        return result

    def collect_all_recommendations(self) -> list[tuple[str, Recommendation]]:
        result: list[tuple[str, Recommendation]] = []
        for rec in self.recommendations:
            result.append(("system", rec))
        for c in self.components:
            for rec in c.all_recommendations:
                result.append((c.id, rec))
        return result

    class Config:
        arbitrary_types_allowed = True
