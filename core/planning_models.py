# core/planning_models.py
"""
Modelos de datos para Interactive Scientific Planning.

Árbol:
    PlanningQuestion
        kind:             QuestionKind
        priority:         QuestionPriority
        reasoning_trigger: str   ← qué razonamiento científico generó esta pregunta
        context:          str   ← párrafo educativo: situación + consecuencias
        question:         str
        options:          list[PlanningOption]
            key, label, description, is_default, state_patch

    PlanningAnswer          ← respuesta seleccionada para una pregunta
    PlanningSessionRecord   ← provenance científica completa (persistida en JSON)

Principios:
    - Las respuestas producen state_patches, no modifican el plan directamente
    - Los patches se aplican al SystemState ANTES de que el decision_engine corra
    - El decision_engine sigue siendo la única fuente de verdad del plan
    - No hay input libre en v1 — solo opciones numeradas y deterministas
"""

from __future__ import annotations

from enum import Enum
from typing import Optional
from pydantic import BaseModel


# ═══════════════════════════════════════════════════════════════════════════════
# Enumeraciones
# ═══════════════════════════════════════════════════════════════════════════════

class QuestionKind(str, Enum):
    MISSING_PARAMETER   = "missing_parameter"
    POLICY_SELECTION    = "policy_selection"
    RISK_CONFIRMATION   = "risk_confirmation"
    SCIENTIFIC_TRADEOFF = "scientific_tradeoff"
    PROTOCOL_SELECTION  = "protocol_selection"


class QuestionPriority(int, Enum):
    BLOCKING  = 0   # debe responderse — el compile no continúa sin respuesta
    IMPORTANT = 1   # afecta el plan científico significativamente
    ADVISORY  = 2   # útil, pero puede saltarse con --no-plan


# ═══════════════════════════════════════════════════════════════════════════════
# Opción de respuesta
# ═══════════════════════════════════════════════════════════════════════════════

class PlanningOption(BaseModel):
    key:         str    # identificador interno: "quick_validation", "rest2", etc.
    label:       str    # texto breve para la UI: "Quick validation (0.1 ns)"
    description: str    # rationale científico del tradeoff
    is_default:  bool  = False  # recomendación del sistema
    state_patch: dict  = {}     # qué modificar en SystemState si se selecciona


# ═══════════════════════════════════════════════════════════════════════════════
# Pregunta de planning
# ═══════════════════════════════════════════════════════════════════════════════

class PlanningQuestion(BaseModel):
    id:                str
    kind:              QuestionKind
    priority:          QuestionPriority
    reasoning_trigger: str             # cadena descriptiva del razonamiento que generó esto
    context:           str             # párrafo educativo mostrado al usuario
    question:          str             # pregunta concreta
    options:           list[PlanningOption]
    applies_to:        Optional[str] = None   # component_id o "system"


# ═══════════════════════════════════════════════════════════════════════════════
# Respuesta del usuario
# ═══════════════════════════════════════════════════════════════════════════════

class PlanningAnswer(BaseModel):
    question_id:    str
    selected_key:   str
    selected_label: str
    state_patch:    dict


# ═══════════════════════════════════════════════════════════════════════════════
# Registro de sesión — provenance científica
# ═══════════════════════════════════════════════════════════════════════════════

class PlanningSessionRecord(BaseModel):
    """
    Registro completo de una sesión de planning.

    Persistido en metadata/planning_session.json dentro del workspace.
    Permite reproducir exactamente las decisiones tomadas durante la compilación.
    """
    created_at:      str
    config_path:     str
    skipped:         bool              = False   # True si se usó --no-plan
    questions_shown: list[dict]        = []      # preguntas serializadas
    answers:         list[PlanningAnswer] = []
    patches_applied: list[dict]        = []      # union de todos los state_patches
