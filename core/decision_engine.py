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
    WorkflowPolicy,
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
        workflow_policy=_build_workflow_policy(state),
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

    _populate_step_params(plan, state)

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

                    condition="needs_parametrization_review",

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

                    condition="needs_pose_validation",

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
                        condition="needs_special_sampling",
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


# ═══════════════════════════════════════════════════════════════════════════════
# Workflow policy
# ═══════════════════════════════════════════════════════════════════════════════

def _build_workflow_policy(state: SystemState) -> WorkflowPolicy:
    """
    Produce WorkflowPolicy desde SystemState.

    Centraliza todas las decisiones científicas globales del pipeline:
    duración, temperatura, presión, timestep, estrategia de sampling.
    Ningún builder hardcodea estos valores — los leen desde step.params.
    """
    policy = WorkflowPolicy()

    policy.temperature_K = state.environment.temperature_K

    # Estrategia de sampling
    if state.global_reasoning.needs_special_sampling:
        policy.enhanced_sampling     = True
        policy.sampling_method       = "REST2"
        policy.equilibration_time_ns = 0.5   # equilibración extendida antes de REST2

    # Duración de producción según objetivos
    objectives = set(state.simulation_objectives)
    if "competitive_binding" in objectives:
        policy.production_time_ns = 50.0
    elif "active_site_stability" in objectives:
        policy.production_time_ns = 20.0
    # else: default 10ns

    return policy


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers internos para step params
# ═══════════════════════════════════════════════════════════════════════════════

# Mapeo de modelo de agua → archivo de referencia de GROMACS
_WATER_GRO: dict[str, str] = {
    "tip3p": "spc216.gro",
    "spce":  "spc216.gro",
    "tip4p": "tip4p.gro",
    "tip5p": "tip5p.gro",
}


def _tc_grps(state: SystemState, temperature_K: float) -> tuple[str, str, str]:
    """
    Devuelve (tc_grps, tau_t, ref_t) coherentes con la composición del sistema.
    Los tres valores deben tener el mismo número de tokens.
    """
    if state.has_membrane():
        groups = ["Protein", "Membrane", "Non-Protein"]
    else:
        groups = ["Protein", "Non-Protein"]

    n       = len(groups)
    tc_grps = " ".join(groups)
    tau_t   = " ".join(["0.1"] * n)
    ref_t   = " ".join([str(temperature_K)] * n)
    return tc_grps, tau_t, ref_t


def _nsteps_from_ns(time_ns: float, dt_ps: float) -> int:
    """Convierte tiempo en ns + timestep en ps → número de pasos."""
    return int(time_ns * 1_000.0 / dt_ps)


# ═══════════════════════════════════════════════════════════════════════════════
# Builders de params por stage
# ═══════════════════════════════════════════════════════════════════════════════

def _build_minimization_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    # emtol más estricto para sistemas con ligandos flexibles
    has_flexible_ligand = any(
        c.descriptors and c.descriptors.flexibility_class in ("flexible", "very_flexible")
        for c in state.components
        if c.role in _LIGAND_ROLES
    )
    emtol = 100.0 if has_flexible_ligand else 1000.0

    return {
        "integrator": "steep",
        "emtol":      emtol,
        "emstep":     0.01,
        "nsteps":     policy.minimization_steps,
    }


def _build_equilibration_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    tc_grps, tau_t, ref_t = _tc_grps(state, policy.temperature_K)
    nsteps = _nsteps_from_ns(policy.equilibration_time_ns, policy.timestep_ps)

    return {
        "dt":          policy.timestep_ps,
        "nvt_nsteps":  nsteps,
        "npt_nsteps":  nsteps,
        "temperature": policy.temperature_K,
        "pressure":    policy.pressure_bar,
        "tc_grps":     tc_grps,
        "tau_t":       tau_t,
        "ref_t":       ref_t,
        "constraints": "h-bonds",
    }


def _build_production_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    tc_grps, tau_t, ref_t = _tc_grps(state, policy.temperature_K)
    nsteps = _nsteps_from_ns(policy.production_time_ns, policy.timestep_ps)

    return {
        "dt":                 policy.timestep_ps,
        "nsteps":             nsteps,
        "temperature":        policy.temperature_K,
        "pressure":           policy.pressure_bar,
        "tc_grps":            tc_grps,
        "tau_t":              tau_t,
        "ref_t":              ref_t,
        "constraints":        "h-bonds",
        "nstxout_compressed": 5000,
        "nstenergy":          1000,
        "nstlog":             1000,
    }


def _build_analysis_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    # Inferir tipo desde step_id ("analysis_rmsd" → "rmsd")
    analysis_type = step.step_id.removeprefix("analysis_")
    config = next((a for a in state.analysis if a.type == analysis_type), None)

    params: dict = {"analysis_type": analysis_type}
    if config and config.selection:
        params["selection"] = config.selection
    return params


def _build_assembly_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    if step.step_id == "solvate_system":
        wm = state.environment.solvent.water_model
        return {
            "box_type":    "dodecahedron",
            "box_distance": 1.2,
            "water_model": wm,
            "water_gro":   _WATER_GRO.get(wm, "spc216.gro"),
        }

    if step.step_id == "add_ions":
        return {
            "concentration": state.environment.ions.concentration,
            "positive_ion":  state.environment.ions.positive,
            "negative_ion":  state.environment.ions.negative,
        }

    return {}


def _build_enhanced_sampling_params(
    step:   SimulationStep,
    state:  SystemState,
    policy: WorkflowPolicy,
) -> dict:
    return {
        "method":     policy.sampling_method,
        "n_replicas": 8,
        "temp_low":   policy.temperature_K,
        "temp_high":  policy.temperature_K * 2.0,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch table y orquestador
# ═══════════════════════════════════════════════════════════════════════════════

_STEP_PARAMS_BUILDERS = {
    StepStage.MINIMIZATION:      _build_minimization_params,
    StepStage.EQUILIBRATION:     _build_equilibration_params,
    StepStage.PRODUCTION:        _build_production_params,
    StepStage.ANALYSIS:          _build_analysis_params,
    StepStage.ASSEMBLY:          _build_assembly_params,
    StepStage.ENHANCED_SAMPLING: _build_enhanced_sampling_params,
}


def _populate_step_params(
    plan:  SimulationPlan,
    state: SystemState,
) -> None:
    """
    Llena step.params para cada step del plan usando el dispatch modular.
    Stages sin builder registrado quedan con params={} (preparation, parametrization, etc.).
    """
    policy = plan.workflow_policy

    for step in plan.steps:
        builder = _STEP_PARAMS_BUILDERS.get(step.stage)
        if builder is not None:
            step.params = builder(step, state, policy)