# core/scientific_planner.py
"""
Scientific planner de SimForge.

Lee un SystemState completamente enriquecido y genera PlanningQuestions
para decisiones científicas ambiguas, riesgosas o incompletas.

Filosofía:
    - NO preguntar parámetros ya definidos explícitamente por el usuario.
    - SÍ preguntar parámetros importantes ausentes (MISSING_PARAMETER / SCIENTIFIC_TRADEOFF).
    - SÍ confirmar parámetros explícitos que parecen riesgosos o inusuales (RISK_CONFIRMATION).
    - Las preguntas enseñan implícitamente buenas prácticas científicas.
    - Cada pregunta tiene contexto, consecuencias y recomendación explícita.

Flujo:
    SystemState (con reasoning completo)
        ↓
    detect_questions(state) → list[PlanningQuestion]   ordenadas por prioridad
        ↓
    CLI dialogue → PlanningAnswer por cada pregunta
        ↓
    apply_patches(state, answers) → SystemState modificado
        ↓
    decision_engine → SimulationPlan (con política ya resuelta)
"""

from __future__ import annotations

from core.models import SystemState
from core.planning_models import (
    QuestionKind,
    QuestionPriority,
    PlanningOption,
    PlanningQuestion,
    PlanningAnswer,
)


# ═══════════════════════════════════════════════════════════════════════════════
# API pública
# ═══════════════════════════════════════════════════════════════════════════════

def detect_questions(state: SystemState) -> list[PlanningQuestion]:
    """
    Genera preguntas de planning para un SystemState enriquecido.

    Solo genera preguntas cuando:
        a) Hay un parámetro importante faltante sin default razonable, o
        b) Un parámetro explícito es estadísticamente riesgoso o inconsistente.

    Retorna la lista ordenada por QuestionPriority (BLOCKING primero).
    """
    questions: list[PlanningQuestion] = []

    # ── Risk confirmations (parámetros explícitos riesgosos) ──────────────────

    if state.environment.temperature_K > 380:
        questions.append(_q_high_temperature(state))

    if state.environment.temperature_K < 270:
        questions.append(_q_low_temperature(state))

    if (
        state.environment.duration_ns is not None
        and state.environment.duration_ns < 0.005
    ):
        questions.append(_q_very_short_duration(state))

    if state.environment.ions.concentration == 0:
        questions.append(_q_no_ions(state))

    if (
        state.has_biological_context("membrane_associated")
        and not state.has_membrane()
    ):
        questions.append(_q_membrane_approximation(state))

    # ── Scientific tradeoffs / missing parameters ─────────────────────────────

    if state.environment.duration_ns is None:
        questions.append(_q_production_duration(state))

    if state.global_reasoning.needs_special_sampling:
        questions.append(_q_sampling_strategy(state))

    # Ordenar: BLOCKING(0) → IMPORTANT(1) → ADVISORY(2)
    questions.sort(key=lambda q: q.priority.value)

    return questions


def apply_patches(
    state:   SystemState,
    answers: list[PlanningAnswer],
) -> SystemState:
    """
    Aplica los state_patches de las respuestas al SystemState.

    Se llama ANTES de que el decision_engine construya el plan,
    para que la WorkflowPolicy refleje las decisiones del usuario.

    Claves de patch soportadas:
        duration_ns        → state.environment.duration_ns
        temperature_K      → state.environment.temperature_K
        ion_concentration  → state.environment.ions.concentration
        membrane_enabled   → state.environment.membrane.enabled
        enhanced_sampling  → state.global_reasoning.needs_special_sampling
    """
    for answer in answers:
        for key, value in answer.state_patch.items():

            if key == "duration_ns":
                state.environment.duration_ns = value

            elif key == "temperature_K":
                state.environment.temperature_K = value

            elif key == "ion_concentration":
                state.environment.ions.concentration = value

            elif key == "membrane_enabled":
                state.environment.membrane.enabled = value

            elif key == "enhanced_sampling":
                state.global_reasoning.needs_special_sampling = value

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# Builders de preguntas
# ═══════════════════════════════════════════════════════════════════════════════

def _q_high_temperature(state: SystemState) -> PlanningQuestion:
    t = state.environment.temperature_K
    return PlanningQuestion(
        id               = "high_temperature",
        kind             = QuestionKind.RISK_CONFIRMATION,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"temperature_K={t} > 380 K — inusualmente alto para MD biomolecular estándar",
        context          = (
            f"temperature_K: {t} K\n\n"
            f"El rango estándar para MD biomolecular es 270–330 K. "
            f"{t} K produce desnaturalización proteica acelerada y puede "
            f"causar inestabilidad numérica en sistemas con agua explícita. "
            f"Esta temperatura es válida en protocolos específicos como "
            f"muestreo conformacional acelerado o recocido simulado, "
            f"pero es inusual para simulaciones de producción estándar."
        ),
        question         = f"¿Por qué se definió temperature_K = {t} K?",
        applies_to       = "system",
        options          = [
            PlanningOption(
                key         = "intended_sampling",
                label       = "Muestreo conformacional acelerado",
                description = "Alta temperatura para explorar el espacio conformacional más rápido. Válido si es intencional.",
                is_default  = False,
                state_patch = {},
            ),
            PlanningOption(
                key         = "intended_annealing",
                label       = "Recocido simulado (simulated annealing)",
                description = "Protocolo de minimización de energía por enfriamiento gradual. Requiere configuración adicional.",
                is_default  = False,
                state_patch = {},
            ),
            PlanningOption(
                key         = "correct_to_300",
                label       = "Corregir a 300 K (estándar)",
                description = "Temperatura estándar para MD biomolecular en agua. Correcto para la mayoría de estudios.",
                is_default  = True,
                state_patch = {"temperature_K": 300.0},
            ),
            PlanningOption(
                key         = "correct_to_physiological",
                label       = "Corregir a 309.65 K (temperatura fisiológica humana)",
                description = "36.5 °C — apropiado para estudios de proteínas humanas en condiciones fisiológicas.",
                is_default  = False,
                state_patch = {"temperature_K": 309.65},
            ),
        ],
    )


def _q_low_temperature(state: SystemState) -> PlanningQuestion:
    t = state.environment.temperature_K
    return PlanningQuestion(
        id               = "low_temperature",
        kind             = QuestionKind.RISK_CONFIRMATION,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"temperature_K={t} < 270 K — por debajo del punto de congelación del agua",
        context          = (
            f"temperature_K: {t} K\n\n"
            f"El agua explícita (TIP3P, SPC/E, TIP4P) se congela por debajo de ~270 K. "
            f"Simular por debajo de este umbral con agua líquida produce "
            f"comportamientos físicamente incorrectos: la caja puede cristalizarse "
            f"parcialmente o producir artefactos de densidad. "
            f"Esta temperatura solo es apropiada para estudios de cryo-MD "
            f"o simulaciones con cosolventes especiales."
        ),
        question         = f"¿Por qué se definió temperature_K = {t} K?",
        applies_to       = "system",
        options          = [
            PlanningOption(
                key         = "intended_cryo",
                label       = "Cryo-MD / simulación de baja temperatura",
                description = "Intencional. Asegúrate de usar un modelo de agua apropiado para este rango.",
                is_default  = False,
                state_patch = {},
            ),
            PlanningOption(
                key         = "correct_to_277",
                label       = "Corregir a 277 K (agua fría, densidad máxima)",
                description = "4 °C — temperatura del agua a densidad máxima. Límite inferior práctico con agua líquida.",
                is_default  = False,
                state_patch = {"temperature_K": 277.0},
            ),
            PlanningOption(
                key         = "correct_to_300",
                label       = "Corregir a 300 K (estándar)",
                description = "Temperatura estándar para MD biomolecular. Correcto para la mayoría de estudios.",
                is_default  = True,
                state_patch = {"temperature_K": 300.0},
            ),
        ],
    )


def _q_very_short_duration(state: SystemState) -> PlanningQuestion:
    d = state.environment.duration_ns
    return PlanningQuestion(
        id               = "very_short_duration",
        kind             = QuestionKind.RISK_CONFIRMATION,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"duration_ns={d} < 5 ps — demasiado corto para equilibrio real",
        context          = (
            f"duration_ns: {d} ns ({d*1000:.1f} ps)\n\n"
            f"Una simulación de {d*1000:.1f} ps es insuficiente para que el sistema "
            f"alcance equilibrio termodinámico después de la equilibración. "
            f"El mínimo práctico para obtener datos científicamente válidos "
            f"es ~100 ps (0.1 ns) — suficiente para verificar el pipeline. "
            f"Para análisis conformacional real, se necesitan al menos 10 ns. "
            f"Una duración tan corta solo tiene sentido para debugging técnico."
        ),
        question         = f"duration_ns = {d} ns parece inusualmente corto. ¿Es intencional?",
        applies_to       = "system",
        options          = [
            PlanningOption(
                key         = "intended_debug",
                label       = "Sí — debugging técnico del pipeline",
                description = "Intencional. El objetivo es verificar que el workflow corra, no producir datos científicos.",
                is_default  = True,
                state_patch = {},
            ),
            PlanningOption(
                key         = "extend_to_100ps",
                label       = "Extender a 0.1 ns (mínimo para validación)",
                description = "100 ps — suficiente para verificar comportamiento físico básico post-equilibración.",
                is_default  = False,
                state_patch = {"duration_ns": 0.1},
            ),
            PlanningOption(
                key         = "extend_to_10ns",
                label       = "Extender a 10 ns (exploratorio)",
                description = "10 ns — suficiente para dinámica inicial y generación de hipótesis.",
                is_default  = False,
                state_patch = {"duration_ns": 10.0},
            ),
        ],
    )


def _q_no_ions(state: SystemState) -> PlanningQuestion:
    return PlanningQuestion(
        id               = "no_ions",
        kind             = QuestionKind.RISK_CONFIRMATION,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= "ions.concentration = 0 — sistema sin sal iónica",
        context          = (
            "ions.concentration: 0.0 M\n\n"
            "Las proteínas biológicas evolucionaron en ambientes con electrolitos "
            "(~150 mM NaCl fisiológico). Simular sin iones produce:"
            "\n  • Cargas superficiales no neutralizadas → artefactos electrostáticos"
            "\n  • Comportamiento conformacional aberrante en loops cargados"
            "\n  • Resultados de binding erróneos si hay sitios iónicos"
            "\n\n"
            "Omitir la sal es válido solo para sistemas específicos como "
            "simulaciones en vacío, peptidos sin carga, o estudios de agua pura."
        ),
        question         = "¿Deseas agregar sal iónica al sistema?",
        applies_to       = "system",
        options          = [
            PlanningOption(
                key         = "physiological",
                label       = "Agregar 0.15 M NaCl (fisiológico)",
                description = "Concentración fisiológica estándar para sistemas biológicos.",
                is_default  = True,
                state_patch = {"ion_concentration": 0.15},
            ),
            PlanningOption(
                key         = "low_salt",
                label       = "Agregar 0.05 M NaCl (baja concentración)",
                description = "Útil para estudios de proteínas en condiciones de baja sal.",
                is_default  = False,
                state_patch = {"ion_concentration": 0.05},
            ),
            PlanningOption(
                key         = "no_salt_intended",
                label       = "Mantener sin sal — es intencional",
                description = "Continuar sin iones. Asegúrate de que el sistema sea eléctricamente neutro.",
                is_default  = False,
                state_patch = {},
            ),
        ],
    )


def _q_membrane_approximation(state: SystemState) -> PlanningQuestion:
    membrane_comps = [
        c.id for c in state.components
        if "membrane_associated" in c.biological_context
    ]
    comp_str = ", ".join(membrane_comps) or "proteína"
    return PlanningQuestion(
        id               = "membrane_approximation",
        kind             = QuestionKind.RISK_CONFIRMATION,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"{comp_str}: biological_context=membrane_associated pero membrane.enabled=False",
        context          = (
            f"{comp_str} está marcado como asociado a membrana, "
            f"pero no se configuró un entorno de membrana.\n\n"
            f"Simular una proteína de membrana en agua explícita sin membrana "
            f"produce artefactos severos: el dominio transmembrana se expone al "
            f"solvente, se colapsa o forma estructuras no nativas. "
            f"Los resultados conformacionales y de binding en estas condiciones "
            f"son típicamente incorrectos.\n\n"
            f"Alternativa válida: usar la proteína truncada o el dominio soluble "
            f"únicamente, si el objetivo científico lo permite."
        ),
        question         = f"¿Cómo manejar la proteína de membrana {comp_str}?",
        applies_to       = membrane_comps[0] if membrane_comps else None,
        options          = [
            PlanningOption(
                key         = "soluble_approx",
                label       = "Aproximación soluble — continuar sin membrana",
                description = "Válido si solo se estudia el dominio soluble o si la membrana no es relevante para el objetivo.",
                is_default  = False,
                state_patch = {},
            ),
            PlanningOption(
                key         = "abort",
                label       = "Abortar — necesito configurar la membrana primero",
                description = "Edita el YAML para agregar membrane: enabled: true y el tipo de membrana.",
                is_default  = True,
                state_patch = {"abort": True},
            ),
        ],
    )


def _q_production_duration(state: SystemState) -> PlanningQuestion:
    objectives   = set(state.simulation_objectives)
    system_type  = state.inferred_system_type or "unknown"

    # Determinar recomendación basada en objetivos
    if "competitive_binding" in objectives:
        recommended_key = "production_grade"
        obj_context = (
            "Sistemas de inhibición competitiva requieren que el ligando y el "
            "sustrato exploren el sitio activo y alcancen estados metaestables. "
            "Corridas cortas (~100 ps) solo capturan relajación inicial; "
            "la dinámica competitiva real emerge a lo largo de decenas de nanosegundos."
        )
    elif "binding" in objectives:
        recommended_key = "production_grade"
        obj_context = (
            "Estudios de binding requieren suficiente muestreo para que el ligando "
            "explore el sitio de unión. La comunidad utiliza 50 ns como estándar "
            "para resultados de calidad de publicación en sistemas típicos."
        )
    elif "stability" in objectives or "active_site_stability" in objectives:
        recommended_key = "exploratory"
        obj_context = (
            "Estudios de estabilidad proteica pueden producir resultados útiles "
            "a partir de 10 ns para dinámica general. 50 ns es preferible "
            "para datos de publicación."
        )
    elif "conformational_sampling" in objectives:
        recommended_key = "production_grade"
        obj_context = (
            "El muestreo conformacional requiere simulaciones largas para explorar "
            "el espacio conformacional de forma representativa."
        )
    else:
        recommended_key = "exploratory"
        obj_context = (
            "Sin un objetivo de simulación específico, 10 ns es un punto "
            "de partida razonable para dinámica exploratoria inicial."
        )

    options = [
        PlanningOption(
            key         = "quick_validation",
            label       = "Quick validation (0.1 ns)",
            description = "Verificación del pipeline. No produce datos científicamente interpretables.",
            is_default  = recommended_key == "quick_validation",
            state_patch = {"duration_ns": 0.1},
        ),
        PlanningOption(
            key         = "exploratory",
            label       = "Exploratorio (10 ns)",
            description = "Dinámica inicial. Útil para generar hipótesis y caracterizar el sistema.",
            is_default  = recommended_key == "exploratory",
            state_patch = {"duration_ns": 10.0},
        ),
        PlanningOption(
            key         = "production_grade",
            label       = "Producción (50 ns)",
            description = "Estándar comunitario para resultados de calidad de publicación en sistemas típicos.",
            is_default  = recommended_key == "production_grade",
            state_patch = {"duration_ns": 50.0},
        ),
        PlanningOption(
            key         = "extended",
            label       = "Extendido (100 ns)",
            description = "Para sistemas complejos, ligandos muy flexibles o estudios de eventos raros.",
            is_default  = False,
            state_patch = {"duration_ns": 100.0},
        ),
    ]

    return PlanningQuestion(
        id               = "production_duration",
        kind             = QuestionKind.SCIENTIFIC_TRADEOFF,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"duration_ns no definido; sistema={system_type}; objetivos={list(objectives)}",
        context          = obj_context,
        question         = "¿Cuánto tiempo debe durar la simulación de producción?",
        applies_to       = "system",
        options          = options,
    )


def _q_sampling_strategy(state: SystemState) -> PlanningQuestion:
    flexible_comps = [
        c.id for c in state.components
        if c.reasoning and c.reasoning.needs_special_sampling
    ]
    comp_str = ", ".join(flexible_comps) or "ligando"

    return PlanningQuestion(
        id               = "sampling_strategy",
        kind             = QuestionKind.SCIENTIFIC_TRADEOFF,
        priority         = QuestionPriority.IMPORTANT,
        reasoning_trigger= f"{comp_str}: needs_special_sampling=True (alta flexibilidad conformacional)",
        context          = (
            f"{comp_str} tiene alta flexibilidad conformacional. "
            f"La MD estándar puede quedar atrapada en mínimos locales "
            f"y no muestrear el espacio conformacional relevante.\n\n"
            f"REST2 (Replica Exchange with Solute Tempering) aplica "
            f"temperatura efectiva solo al soluto, mejorando el muestreo "
            f"sin el costo de escalar todo el sistema."
        ),
        question         = "¿Qué estrategia de muestreo usar?",
        applies_to       = flexible_comps[0] if flexible_comps else None,
        options          = [
            PlanningOption(
                key         = "standard_md",
                label       = "MD estándar",
                description = "Sin muestreo mejorado. Más rápido pero puede subestimar la flexibilidad real.",
                is_default  = False,
                state_patch = {"enhanced_sampling": False},
            ),
            PlanningOption(
                key         = "rest2",
                label       = "REST2 (recomendado)",
                description = "Enhanced sampling sobre el soluto. Mejor muestreo conformacional sin escalar todo el sistema.",
                is_default  = True,
                state_patch = {"enhanced_sampling": True},
            ),
        ],
    )
