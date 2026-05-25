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

from core.models import SystemState, Warning, Risk, Recommendation, Severity


# ─── 1. Inferencia de tipo de sistema ────────────────────────────────────────

def infer_system_type(state: SystemState) -> SystemState:
    roles = {c.role for c in state.components}

    has_protein     = "protein" in roles
    has_substrate   = "substrate" in roles
    has_competitive = "competitive_ligand" in roles
    has_allosteric  = "allosteric_ligand" in roles
    has_membrane    = state.environment.membrane.enabled

    # Secondary membrane signals (set by run_semantic_normalization before us):
    #   - workflow_hints.membrane_required: set when membrane objectives detected
    #   - objectives already normalized: membrane_perturbation/membrane_insertion/permeability
    _MEM_OBJ = {"membrane_perturbation", "membrane_insertion", "permeability"}
    has_membrane_hints = getattr(state.workflow_hints, "membrane_required", False)
    has_membrane_obj   = bool(set(state.simulation_objectives) & _MEM_OBJ)

    # implicit_membrane: treat as membrane system when objectives say so but
    # explicit flag is missing — covers old YAMLs and manual edits.
    # Guard: only when no ligand/substrate that would classify differently.
    implicit_membrane = (
        (has_membrane_hints or has_membrane_obj)
        and not has_membrane
        and not (has_substrate or has_competitive)
        and not has_allosteric
    )

    if has_protein and has_substrate and has_competitive:
        state.inferred_system_type = "competitive-inhibition"

    elif has_protein and has_allosteric:
        state.inferred_system_type = "allosteric-modulation"

    elif has_protein and has_membrane and (has_substrate or has_competitive):
        state.inferred_system_type = "protein-membrane-ligand"

    elif has_protein and has_membrane:
        state.inferred_system_type = "protein-membrane"

    elif has_protein and implicit_membrane:
        # Membrane inferred from objectives — enable bilayer flag so pipeline
        # and downstream hints work correctly.
        state.environment.membrane.enabled = True
        state.inferred_system_type = "protein-membrane"
        state.global_reasoning.notes.append(
            "[inference] membrane.enabled was False but membrane objectives "
            f"({sorted(set(state.simulation_objectives) & _MEM_OBJ)}) detected — "
            "routing to protein-membrane workflow. "
            "Add 'environment.membrane.enabled: true' to suppress this fallback."
        )

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
        state.warnings.append(Warning(
            message  = "membrane_associated detectado pero membrane.enabled=False",
            target   = "environment.membrane",
            severity = Severity.MEDIUM,
        ))
        state.risks.append(Risk(
            message  = "Inestabilidad en terminales transmembranales durante equilibración",
            target   = "protein_1",
            severity = Severity.HIGH,
        ))
        state.recommendations.append(Recommendation(
            message = "Habilitar membrana o aplicar restraints en terminales expuestos",
            target  = "protein_1",
            action  = "Agregar terminal_restraints o establecer membrane.enabled=True",
        ))

    # membrane_associated sin terminal_restraints
    if (state.has_biological_context("membrane_associated")
            and "terminal_restraints" not in restraint_types):
        state.warnings.append(Warning(
            message  = "Proteína membrane_associated sin terminal_restraints definidos",
            target   = "restraints",
            severity = Severity.HIGH,
        ))
        state.risks.append(Risk(
            message  = "Extremos transmembranales pueden desenrollarse durante equilibración",
            target   = "protein_1",
            severity = Severity.HIGH,
        ))
        state.recommendations.append(Recommendation(
            message = "Definir terminal_restraints para protein_1",
            target  = "restraints",
            action  = "Agregar type: terminal_restraints con target: protein_1 en el YAML",
        ))

    # partially_truncated
    if state.has_biological_context("partially_truncated"):
        state.warnings.append(Warning(
            message  = "Proteína parcialmente truncada detectada",
            target   = "protein_1",
            severity = Severity.MEDIUM,
        ))
        state.risks.append(Risk(
            message  = "Terminales artificiales pueden generar artefactos estructurales",
            target   = "protein_1",
            severity = Severity.MEDIUM,
        ))
        state.recommendations.append(Recommendation(
            message = "Verificar y aplicar capping en terminales artificiales",
            target  = "protein_1",
            action  = "Aplicar ACE/NME en extremos expuestos durante protein_builder",
        ))

    return state


# ─── 3. Inferencia de gaps en análisis ───────────────────────────────────────

def infer_analysis_gaps(state: SystemState) -> SystemState:
    analysis_types = {a.type for a in state.analysis}

    # competitive-inhibition requiere distance_analysis
    if (state.inferred_system_type == "competitive-inhibition"
            and "distance_analysis" not in analysis_types):
        state.warnings.append(Warning(
            message  = "Sistema competitive-inhibition sin distance_analysis definido",
            target   = "analysis",
            severity = Severity.MEDIUM,
        ))
        state.recommendations.append(Recommendation(
            message = "Agregar distance_analysis para monitorear competencia en sitio activo",
            target  = "analysis",
            action  = "Definir distance_analysis con group1: substrate_1 y group2: ligand_1",
        ))

    # distance_analysis sin selection explícita
    for da in state.analysis:
        if da.type == "distance_analysis" and da.selection is None:
            state.warnings.append(Warning(
                message  = "distance_analysis sin selection explícita",
                target   = "analysis.distance_analysis",
                severity = Severity.LOW,
            ))
            state.recommendations.append(Recommendation(
                message = "Especificar group1 y group2 en distance_analysis",
                target  = "analysis.distance_analysis",
                action  = "Agregar selection con group1: substrate_1 y group2: ligand_1",
            ))

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