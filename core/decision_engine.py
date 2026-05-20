# core/decision_engine.py
"""
Decision engine de SimForge.

Transforma un SystemState completamente enriquecido en un
SimulationPlan ejecutable a nivel semántico.

Principios arquitectónicos:
    - NO parsea archivos
    - NO calcula descriptores
    - NO hace inferencia química
    - SOLO consume SystemState

Entrada:
    SystemState enriquecido por parser.py

Salida:
    SimulationPlan

Pipeline lógico:
    SystemState
        ↓
    decision_engine
        ↓
    SimulationPlan
"""

from __future__ import annotations

from core.models import (
    SystemState,
    Severity,
)

from core.execution_models import (
    SimulationPlan,
    SimulationStep,
    BlockingIssue,
    CheckItem,
    PlanStatus,
    StepStage,
    StepType,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Roles
# ═══════════════════════════════════════════════════════════════════════════════

_LIGAND_ROLES = {
    "substrate",
    "competitive_ligand",
    "allosteric_ligand",
    "cofactor",
    "essential_oil_component",
}

_PROTEIN_ROLES = {
    "protein",
    "peptide",
}


# ═══════════════════════════════════════════════════════════════════════════════
# API pública
# ═══════════════════════════════════════════════════════════════════════════════

def build_simulation_plan(state: SystemState) -> SimulationPlan:
    """
    Construye un plan de simulación desde un SystemState enriquecido.

    El decision engine:
        - interpreta flags globales
        - traduce reasoning → acciones
        - construye workflow semántico
        - detecta issues bloqueantes
        - genera checklist previo a ejecución
    """

    plan = SimulationPlan(
        status=PlanStatus.READY,
        inferred_system_type=state.inferred_system_type,
    )

    _evaluate_global_status(state, plan)

    _collect_blocking_issues(state, plan)

    _deduplicate_blocking_issues(plan)

    _build_preparation_steps(state, plan)

    _build_parametrization_steps(state, plan)

    _build_validation_steps(state, plan)

    _build_assembly_steps(state, plan)

    _build_md_steps(state, plan)

    _build_analysis_steps(state, plan)

    _build_special_protocols(state, plan)

    _build_checklist(state, plan)

    _finalize_plan(state, plan)

    return plan


# ═══════════════════════════════════════════════════════════════════════════════
# Estado global
# ═══════════════════════════════════════════════════════════════════════════════

def _evaluate_global_status(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    gr = state.global_reasoning

    if gr.has_blocking_errors:
        plan.status = PlanStatus.BLOCKED

    elif (
        gr.needs_special_sampling
        or len(gr.warnings) > 0
    ):
        plan.status = PlanStatus.NEEDS_REVIEW

    else:
        plan.status = PlanStatus.READY


# ═══════════════════════════════════════════════════════════════════════════════
# Issues bloqueantes
# ═══════════════════════════════════════════════════════════════════════════════

def _collect_blocking_issues(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for comp in state.components:

        for risk in comp.all_risks:

            if risk.severity == Severity.HIGH:

                plan.blocking_issues.append(
                    BlockingIssue(
                        source=comp.id,
                        message=risk.message,
                        severity=risk.severity,
                    )
                )

    for risk in state.risks:

        if risk.severity == Severity.HIGH:

            plan.blocking_issues.append(
                BlockingIssue(
                    source="system",
                    message=risk.message,
                    severity=risk.severity,
                )
            )

def _deduplicate_blocking_issues(
    plan: SimulationPlan,
) -> None:

    unique = {}

    for issue in plan.blocking_issues:

        key = (
            issue.source,
            issue.message,
        )

        unique[key] = issue

    plan.blocking_issues = list(
        unique.values()
    )
# ═══════════════════════════════════════════════════════════════════════════════
# Preparation
# ═══════════════════════════════════════════════════════════════════════════════


def _build_preparation_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for comp in state.components:

        if comp.role in _PROTEIN_ROLES:

            plan.steps.append(
                SimulationStep(
                    step_id=f"prepare_{comp.id}",
                    title=f"Preparar proteína {comp.id}",
                    stage=StepStage.PREPARATION,
                    engine="gromacs:pdb2gmx",
                    target_components=[comp.id],
                    notes=[
                        "Agregar hidrógenos",
                        "Asignar protonación",
                        "Construir topología",
                    ],
                )
            )

        elif comp.role in _LIGAND_ROLES:

            plan.steps.append(
                SimulationStep(
                    step_id=f"prepare_{comp.id}",
                    title=f"Preparar ligando {comp.id}",
                    stage=StepStage.PREPARATION,
                    engine="ligand_preparation",
                    target_components=[comp.id],
                )
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Parametrización
# ═══════════════════════════════════════════════════════════════════════════════

def _build_parametrization_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    ligand_ff = state.forcefields.ligands or "unknown"

    for comp in state.components:

        if comp.role not in _LIGAND_ROLES:
            continue

        blocking = False

        if (
            comp.reasoning
            and comp.reasoning.needs_parametrization_review
        ):
            blocking = True

        plan.steps.append(
            SimulationStep(
                step_id=f"parametrize_{comp.id}",
                title=f"Parametrizar {comp.id}",
                stage=StepStage.PARAMETRIZATION,
                engine=ligand_ff,
                target_components=[comp.id],
                blocking=blocking,
                depends_on=[f"prepare_{comp.id}"],
                notes=[
                    f"Forcefield seleccionado: {ligand_ff}",
                ],
            )
        )
# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

def _build_validation_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for comp in state.components:

        reas = comp.reasoning

        if reas is None:
            continue

        # ── revisión manual de parametrización ───────────────────────────────

        if reas.needs_parametrization_review:

            plan.steps.append(
                SimulationStep(
                    step_id=f"review_parametrization_{comp.id}",

                    title=f"Revisión manual de parametrización: {comp.id}",

                    stage=StepStage.VALIDATION,

                    step_type=StepType.MANUAL,

                    engine="manual_review",

                    target_components=[comp.id],

                    blocking=True,

                    depends_on=[
                        f"parametrize_{comp.id}"
                    ],

                    notes=[
                        "Revisar penalizaciones ParamChem",
                        "Verificar cargas",
                        "Validar constantes de fuerza",
                    ],
                )
            )

        # ── validación de pose ───────────────────────────────────────────────

        if reas.needs_pose_validation:

            plan.steps.append(
                SimulationStep(
                    step_id=f"validate_pose_{comp.id}",

                    title=f"Validar pose inicial: {comp.id}",

                    stage=StepStage.VALIDATION,

                    step_type=StepType.VALIDATION,

                    engine="pose_inspection",

                    target_components=[comp.id],

                    depends_on=[
                        f"prepare_{comp.id}"
                    ],

                    notes=[
                        "Inspeccionar orientación inicial",
                        "Verificar clashes",
                        "Validar contactos del sitio activo",
                    ],
                )
            )

# ═══════════════════════════════════════════════════════════════════════════════
# Assembly
# ═══════════════════════════════════════════════════════════════════════════════

def _build_assembly_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    dependencies = [
        step.step_id
        for step in plan.steps
        if step.stage in (
            StepStage.PREPARATION,
            StepStage.PARAMETRIZATION,
            StepStage.VALIDATION,
        )
    ]

    plan.steps.append(
        SimulationStep(
            step_id="assemble_system",
            title="Construcción del sistema completo",
            stage=StepStage.ASSEMBLY,
            engine="system_builder",
            depends_on=dependencies,
            notes=[
                "Integrar proteína, ligandos y solvente",
                "Construir caja de simulación",
            ],
        )
    )

    if state.has_membrane():

        plan.steps.append(
            SimulationStep(
                step_id="build_membrane",
                title="Construcción de membrana",
                stage=StepStage.ASSEMBLY,
                engine="charmm-gui",
                depends_on=["assemble_system"],
                notes=[
                    "Insertar proteína en bicapa lipídica",
                ],
            )
        )

    plan.steps.append(
        SimulationStep(
            step_id="solvate_system",
            title="Solvatación",
            stage=StepStage.ASSEMBLY,
            engine="gromacs:solvate",
            depends_on=["assemble_system"],
        )
    )

    plan.steps.append(
        SimulationStep(
            step_id="add_ions",
            title="Adición de iones",
            stage=StepStage.ASSEMBLY,
            engine="gromacs:genion",
            depends_on=["solvate_system"],
        )
    )


# ═══════════════════════════════════════════════════════════════════════════════
# MD core
# ═══════════════════════════════════════════════════════════════════════════════

def _build_md_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    plan.steps.append(
        SimulationStep(
            step_id="energy_minimization",
            title="Minimización energética",
            stage=StepStage.MINIMIZATION,
            engine="gromacs",
            depends_on=["add_ions"],
        )
    )

    equilibration_notes = []

    if state.global_reasoning.needs_special_sampling:
        equilibration_notes.append(
            "Aplicar equilibración extendida"
        )

    plan.steps.append(
        SimulationStep(
            step_id="equilibration",
            title="Equilibración",
            stage=StepStage.EQUILIBRATION,
            engine="gromacs",
            depends_on=["energy_minimization"],
            notes=equilibration_notes,
        )
    )

    production_notes = []

    if state.inferred_system_type == "competitive-inhibition":

        production_notes.append(
            "Monitorear competencia sustrato/inhibidor"
        )


    # ── dependencias dinámicas ─────────────────────────────

    production_dependencies = ["equilibration"]

    if state.global_reasoning.needs_special_sampling:

        production_dependencies.append(
            "rest2_sampling"
        )


    # ── production MD ─────────────────────────────────────

    plan.steps.append(
        SimulationStep(
            step_id="production_md",
            title="Dinámica molecular de producción",
            stage=StepStage.PRODUCTION,
            engine="gromacs",
            depends_on=production_dependencies,
            notes=production_notes,
        )
    )

# ═══════════════════════════════════════════════════════════════════════════════
# Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def _build_analysis_steps(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for analysis in state.analysis:

        analysis_type = analysis.type

        plan.steps.append(
            SimulationStep(
                step_id=f"analysis_{analysis_type}",
                title=f"Análisis {analysis_type}",
                stage=StepStage.ANALYSIS,
                engine="analysis_pipeline",
                depends_on=["production_md"],
                notes=[
                    f"Tipo de análisis: {analysis_type}",
                ],
            )
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Protocolos especiales
# ═══════════════════════════════════════════════════════════════════════════════

def _build_special_protocols(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for comp in state.components:

        reas = comp.reasoning

        if reas is None:
            continue

        if reas.needs_special_sampling:

            if "enhanced_sampling" not in plan.special_protocols:
                plan.special_protocols.append(
                    "enhanced_sampling"
                )

            if "REST2" not in plan.special_protocols:

                plan.special_protocols.append(
                    "REST2"
                )

                plan.steps.append(
                    SimulationStep(
                        step_id="rest2_sampling",
                        title="REST2 enhanced sampling",
                        stage=StepStage.ENHANCED_SAMPLING,
                        step_type=StepType.EXTERNAL,
                        engine="plumed:gromacs",
                        depends_on=["equilibration"],
                        notes=[
                            "Enhanced conformational sampling",
                            "Replica exchange with solute tempering",
                        ],
                    )
                )

        if reas.needs_pose_validation:

            if "pose_validation" not in plan.special_protocols:
                plan.special_protocols.append(
                    "pose_validation"
                )

        if reas.needs_parametrization_review:

            if "manual_parametrization_review" not in plan.special_protocols:
                plan.special_protocols.append(
                    "manual_parametrization_review"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# Checklist
# ═══════════════════════════════════════════════════════════════════════════════

def _build_checklist(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    for comp in state.components:

        reas = comp.reasoning

        if reas is None:
            continue

        if reas.needs_parametrization_review:

            plan.checklist.append(
                CheckItem(
                    message=(
                        f"Revisar parametrización de {comp.id}"
                    ),
                    required=True,
                    component=comp.id,
                )
            )

        if reas.needs_protonation_check:

            plan.checklist.append(
                CheckItem(
                    message=(
                        f"Verificar protonación de {comp.id}"
                    ),
                    required=True,
                    component=comp.id,
                )
            )

        if reas.needs_pose_validation:

            plan.checklist.append(
                CheckItem(
                    message=(
                        f"Validar pose inicial de {comp.id}"
                    ),
                    required=True,
                    component=comp.id,
                )
            )

    if state.has_membrane():

        plan.checklist.append(
            CheckItem(
                message="Verificar orientación proteína-membrana",
                required=True,
            )
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Finalización
# ═══════════════════════════════════════════════════════════════════════════════

def _finalize_plan(
    state: SystemState,
    plan: SimulationPlan,
) -> None:

    plan.notes.append(
        f"Sistema inferido: {state.inferred_system_type}"
    )

    plan.notes.append(
        f"Componentes: {len(state.components)}"
    )

    plan.notes.append(
        f"Steps generados: {len(plan.steps)}"
    )

    if plan.status == PlanStatus.BLOCKED:

        plan.notes.append(
            "El sistema contiene issues bloqueantes "
            "que deben resolverse antes de producción."
        )

    elif plan.status == PlanStatus.NEEDS_REVIEW:

        plan.notes.append(
            "El sistema requiere revisión manual "
            "antes de ejecutar producción."
        )

    else:

        plan.notes.append(
            "Sistema listo para ejecución."
        )