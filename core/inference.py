# core/inference.py
"""
Pipeline de inferencia biológica del SystemState.

Orden de ejecución garantizado:
  1. infer_system_type      — qué tipo de sistema es
  2. infer_biological_risks — qué riesgos biológicos existen
  3. infer_analysis_gaps    — qué análisis faltan dado el objetivo

Cada función recibe y retorna SystemState.
Agregar nuevas reglas = agregar función + registrarla en run_inference().
"""

from core.models import SystemState


# ─── 1. Inferencia de tipo de sistema ────────────────────────────────────────

def infer_system_type(state: SystemState) -> SystemState:
    roles = {c.role for c in state.components}

    has_protein     = "protein" in roles
    has_substrate   = "substrate" in roles
    has_competitive = "competitive_ligand" in roles
    has_allosteric  = "allosteric_ligand" in roles
    has_membrane    = state.environment.membrane.enabled

    if has_protein and has_substrate and has_competitive:
        state.inferred_system_type = "competitive-inhibition"

    elif has_protein and has_allosteric:
        state.inferred_system_type = "allosteric-modulation"

    elif has_protein and has_membrane and (has_substrate or has_competitive):
        state.inferred_system_type = "protein-membrane-ligand"

    elif has_protein and has_membrane:
        state.inferred_system_type = "protein-membrane"

    elif has_protein and (has_substrate or has_competitive):
        state.inferred_system_type = "protein-ligand"

    elif has_protein:
        state.inferred_system_type = "protein"

    else:
        state.inferred_system_type = "multicomponent-system"

    return state


# ─── 2. Inferencia de riesgos biológicos ─────────────────────────────────────

def infer_biological_risks(state: SystemState) -> SystemState:
    restraint_types = {r.type for r in state.restraints}

    # membrane_associated sin membrana habilitada
    if state.has_biological_context("membrane_associated") and not state.has_membrane():
        state.warnings.append(
            "membrane_associated detectado pero membrane.enabled=False. "
            "La proteína puede mostrar inestabilidad en terminales transmembranales. "
            "Se recomienda aplicar terminal_restraints."
        )

    # membrane_associated sin terminal_restraints
    if (state.has_biological_context("membrane_associated")
            and "terminal_restraints" not in restraint_types):
        state.warnings.append(
            "Proteína membrane_associated sin terminal_restraints definidos. "
            "Los extremos transmembranales pueden desenrollarse durante equilibración."
        )

    # partially_truncated
    if state.has_biological_context("partially_truncated"):
        state.warnings.append(
            "Proteína parcialmente truncada detectada. "
            "Verificar que los terminales artificiales no generen artefactos. "
            "Considerar capping (ACE/NME) en los extremos expuestos."
        )

    return state


# ─── 3. Inferencia de gaps en análisis ───────────────────────────────────────

def infer_analysis_gaps(state: SystemState) -> SystemState:
    analysis_types = {a.type for a in state.analysis}

    # competitive-inhibition requiere distance_analysis
    if (state.inferred_system_type == "competitive-inhibition"
            and "distance_analysis" not in analysis_types):
        state.warnings.append(
            "Sistema competitive-inhibition sin distance_analysis definido. "
            "Se recomienda monitorear distancias substrate↔active_site y ligand↔active_site."
        )

    # distance_analysis sin selection explícita
    for da in state.analysis:
        if da.type == "distance_analysis" and da.selection is None:
            state.warnings.append(
                "distance_analysis definido sin selection explícita. "
                "Especificar group1 y group2 para substrate y competitive_ligand."
            )

    return state


# ─── Pipeline principal ───────────────────────────────────────────────────────

def run_inference(state: SystemState) -> SystemState:
    """
    Ejecuta el pipeline de inferencia en orden garantizado.
    Cada paso depende del anterior.
    """
    state = infer_system_type(state)       # paso 1: qué es
    state = infer_biological_risks(state)  # paso 2: qué riesgos tiene
    state = infer_analysis_gaps(state)     # paso 3: qué falta analizar
    return state
