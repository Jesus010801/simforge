# core/parser.py
"""
Pipeline de construcción del SystemState.

Etapas en orden garantizado:
    1. parse_yaml()          — YAML → SystemState (solo config)
    2. run_inference()       — infiere tipo de sistema y riesgos globales
    3. run_validation()      — valida archivos de cada componente
    4. run_descriptors()     — calcula descriptores moleculares (ligandos)
    5. run_component_reasoning() — razonamiento contextual por componente
    6. run_global_reasoning()    — consolida flags globales para decision engine

Resultado: SystemState completamente enriquecido donde el decision engine
puede leer todo desde state.components[i].* y state.global_reasoning.*
sin acceder a validators ni a descriptores directamente.

Cada etapa es idempotente y puede ejecutarse independientemente.
"""

from __future__ import annotations
from pathlib import Path
import yaml

from core.models import (
    SystemState,
    ComponentModel,
    ComponentValidation,
    ComponentDescriptors,
    ComponentReasoning,
    GlobalReasoning,
    Warning,
    Risk,
    Recommendation,
    Severity,
)
from core.inference import run_inference


# ═══════════════════════════════════════════════════════════════════════════════
# Etapa 1 — YAML → SystemState (config solamente)
# ═══════════════════════════════════════════════════════════════════════════════

def _load_yaml(path: Path) -> dict:
    """Carga y pre-procesa el YAML: resuelve paths relativos de componentes."""
    with open(path, "r") as f:
        raw = yaml.safe_load(f)

    base_dir = path.parent
    if "components" in raw:
        for component in raw["components"]:
            if "file" in component:
                cp = Path(component["file"])
                if not cp.is_absolute():
                    component["file"] = str(base_dir / cp)
    return raw


# ═══════════════════════════════════════════════════════════════════════════════
# Etapa 3 — Validación de componentes
# ═══════════════════════════════════════════════════════════════════════════════

def _validate_protein_component(comp: ComponentModel) -> ComponentValidation:
    """
    Ejecuta protein_validator y empaqueta el resultado en ComponentValidation.
    Si el archivo no existe o falla, retorna una validación con error.
    """
    from validators.protein_validator import validate_protein
    from validators.file_resolver import resolve_structure_file, StructureFileError

    raw_path = Path(comp.file)

    # Smart directory detection before calling validator
    try:
        path, auto_warning = resolve_structure_file(raw_path, component_id=comp.id)
        # Patch the component's file path so downstream stages use the resolved path
        comp.file = str(path)
    except StructureFileError as e:
        return ComponentValidation(
            validator_used   = "protein_validator",
            is_valid         = False,
            validation_error = str(e),
            warnings=[Warning(
                message  = str(e),
                target   = comp.id,
                severity = Severity.HIGH,
            )],
            risks=[Risk(
                message  = "Archivo de estructura no resuelto — pipeline bloqueado",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )
    except FileNotFoundError:
        return ComponentValidation(
            validator_used = "protein_validator",
            is_valid       = False,
            validation_error = f"Archivo no encontrado: {raw_path}",
            warnings = [Warning(
                message  = f"Archivo de proteína no encontrado: {raw_path.name}",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )

    extra_warnings: list[Warning] = []
    if auto_warning:
        extra_warnings.append(Warning(
            message  = auto_warning,
            target   = comp.id,
            severity = Severity.MEDIUM,
        ))

    try:
        pv = validate_protein(path)
        return ComponentValidation(
            validator_used   = "protein_validator",
            is_valid         = True,
            validation_error = None,
            warnings         = extra_warnings + pv.warnings,
            risks            = pv.risks,
            recommendations  = pv.recommendations,
            data             = pv.model_dump(),
        )
    except Exception as e:
        return ComponentValidation(
            validator_used   = "protein_validator",
            is_valid         = False,
            validation_error = str(e),
            warnings = extra_warnings + [Warning(
                message  = f"Error al validar proteína: {e}",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
            risks = [Risk(
                message  = "Proteína no pudo ser validada — pipeline bloqueado",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )


def _validate_ligand_component(comp: ComponentModel) -> ComponentValidation:
    """
    Ejecuta ligand_validator y empaqueta el resultado en ComponentValidation.
    """
    from validators.ligand_validator import validate_ligand
    from validators.file_resolver import resolve_structure_file, StructureFileError

    raw_path = Path(comp.file)

    try:
        path, auto_warning = resolve_structure_file(raw_path, component_id=comp.id)
        comp.file = str(path)
    except StructureFileError as e:
        return ComponentValidation(
            validator_used   = "ligand_validator",
            is_valid         = False,
            validation_error = str(e),
            warnings=[Warning(
                message  = str(e),
                target   = comp.id,
                severity = Severity.HIGH,
            )],
            risks=[Risk(
                message  = "Archivo de ligando no resuelto — parametrización imposible",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )
    except FileNotFoundError:
        return ComponentValidation(
            validator_used   = "ligand_validator",
            is_valid         = False,
            validation_error = f"Archivo no encontrado: {raw_path}",
            warnings = [Warning(
                message  = f"Archivo de ligando no encontrado: {raw_path.name}",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )

    extra_warnings: list[Warning] = []
    if auto_warning:
        extra_warnings.append(Warning(
            message  = auto_warning,
            target   = comp.id,
            severity = Severity.MEDIUM,
        ))

    try:
        lv = validate_ligand(path, role=comp.role)
        return ComponentValidation(
            validator_used   = "ligand_validator",
            is_valid         = lv.is_complete,
            validation_error = None if lv.is_complete else "Parser reportó archivo incompleto",
            warnings         = extra_warnings + lv.warnings,
            risks            = lv.risks,
            recommendations  = lv.recommendations,
            data             = lv.model_dump(),
        )
    except Exception as e:
        return ComponentValidation(
            validator_used   = "ligand_validator",
            is_valid         = False,
            validation_error = str(e),
            warnings = extra_warnings + [Warning(
                message  = f"Error al validar ligando: {e}",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
            risks = [Risk(
                message  = "Ligando no pudo ser validado — parametrización imposible",
                target   = comp.id,
                severity = Severity.HIGH,
            )],
        )


_LIGAND_ROLES = {"substrate", "competitive_ligand", "allosteric_ligand",
                 "cofactor", "essential_oil_component"}
_PROTEIN_ROLES = {"protein", "peptide"}


def run_validation(state: SystemState) -> SystemState:
    """
    Etapa 3: valida los archivos de todos los componentes.
    Escribe el resultado en component.validation.
    No modifica ningún otro campo del estado.
    """
    for comp in state.components:
        if not comp.file:
            continue

        if comp.role in _PROTEIN_ROLES:
            comp.validation = _validate_protein_component(comp)

        elif comp.role in _LIGAND_ROLES:
            comp.validation = _validate_ligand_component(comp)

        # Roles como "membrane", "solvent", "ion" no tienen validator por ahora
        # → se dejan con validation=None

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# Etapa 4 — Descriptores moleculares (ligandos)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_descriptors_from_ligand_validation(
    comp: ComponentModel,
) -> ComponentDescriptors | None:
    """
    Construye ComponentDescriptors desde los resultados del ligand_validator.

    El ligand_validator ya ejecutó el descriptor engine internamente
    (topology, aromaticity, flexibility, geometry, polarity).
    Aquí extraemos los datos del dict serializado en validation.data
    y los mapeamos al modelo ComponentDescriptors.

    Si la validación falló o no hay datos suficientes, retorna None.
    """
    if comp.validation is None or not comp.validation.is_valid:
        return None

    data = comp.validation.data
    if not data:
        return None

    # ── Extraer campos del resumen desde el resultado del ligand_validator ───
    # El ligand_validator ya calcula y expone estos campos directamente.
    def _get(key, default=None):
        return data.get(key, default)

    try:
        desc = ComponentDescriptors(
            # Conteos básicos
            n_heavy_atoms        = _get("n_atoms", 0),

            # Aromaticidad (del resumen del validator)
            n_aromatic_rings     = _get("n_aromatic_rings", 0),
            n_fused_aromatic     = _get("n_fused_aromatic", 0),

            # Flexibilidad
            flexibility_class    = _get("estimated_flexibility", "unknown"),
            n_rotatable_bonds    = _get("n_rotatable_bonds", 0),
            scaffold_rigidity    = _get("scaffold_rigidity", "unknown"),
            sampling_recommendation = _get("sampling_recommendation", ""),

            # Geometría
            shape_class          = _get("shape_class", "unknown"),
            is_globally_planar   = _get("is_globally_planar", False),
            global_planarity_rms = _get("global_planarity_rms", 0.0),
            radius_of_gyration   = _get("radius_of_gyration", 0.0),

            # Polaridad
            polarity_class       = _get("estimated_polarity", "unknown"),
            hb_donors            = _get("hbd_count", 0),
            hb_acceptors         = _get("hba_count", 0),
            net_charge           = _get("net_charge", 0),
            is_zwitterion        = False,   # el validator actual no expone esto
            amphipathic_class    = "unknown",
            n_functional_groups  = len(_get("functional_groups", []) or []),
            lipinski_hbd         = _get("hbd_count", 0),
            lipinski_hba         = _get("hba_count", 0),
            passes_lipinski      = _get("lipinski_compliant", True),
        )
        return desc

    except Exception:
        # Si algo falla en la extracción, retornar descriptores vacíos pero válidos
        return ComponentDescriptors()


def run_descriptors(state: SystemState) -> SystemState:
    """
    Etapa 4: construye ComponentDescriptors para cada componente ligando
    desde los datos ya calculados en component.validation.

    Para proteínas: reservado (MDAnalysis futuro).
    """
    for comp in state.components:
        if comp.role in _LIGAND_ROLES:
            comp.descriptors = _build_descriptors_from_ligand_validation(comp)
        # Para proteínas, comp.descriptors queda None por ahora

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# Etapa 5 — Reasoning contextual por componente
# ═══════════════════════════════════════════════════════════════════════════════

def run_component_reasoning(state: SystemState) -> SystemState:
    """
    Etapa 5: genera ComponentReasoning para cada componente en el
    contexto del sistema completo.

    Lee:
        - comp.descriptors (descriptores del ligando)
        - comp.validation  (issues del validator)
        - comp.role        (rol en el sistema)
        - state.inferred_system_type
        - state.forcefields

    Escribe:
        - comp.reasoning   (flags + warnings contextuales)

    Principio: los warnings del validator son estructurales ("este archivo
    tiene fósforo"). Los warnings de reasoning son contextuales ("en un sistema
    competitive-inhibition, un sustrato con fosfato requiere parámetros especiales
    porque el forcefield es CGenFF").
    """
    ff_ligands = state.forcefields.ligands or "unknown"
    sys_type   = state.inferred_system_type or "unknown"

    for comp in state.components:
        reasoning = ComponentReasoning(
            component_id   = comp.id,
            component_role = comp.role,
        )

        desc = comp.descriptors
        val  = comp.validation

        if comp.role in _LIGAND_ROLES and desc is not None:

            # ── Flexibilidad very_flexible ────────────────────────────────────
            if desc.flexibility_class == "very_flexible":
                reasoning.needs_special_sampling = True
                reasoning.warnings.append(Warning(
                    message  = (
                        f"{comp.id} es muy flexible ({desc.n_rotatable_bonds} rot. bonds) "
                        f"en sistema {sys_type}: MD estándar probablemente insuficiente"
                    ),
                    target   = comp.id,
                    severity = Severity.HIGH,
                ))
                reasoning.risks.append(Risk(
                    message  = (
                        "Muestreo conformacional incompleto puede producir pose inicial "
                        "subóptima y energías de unión incorrectas"
                    ),
                    target   = comp.id,
                    severity = Severity.HIGH,
                ))
                reasoning.recommendations.append(Recommendation(
                    message = f"Aplicar estrategia de muestreo extendido para {comp.id}",
                    target  = comp.id,
                    action  = desc.sampling_recommendation or
                              "Considerar REST2, metadinámica o docking exhaustivo previo.",
                ))

            elif desc.flexibility_class == "flexible":
                reasoning.needs_pose_validation = True
                reasoning.recommendations.append(Recommendation(
                    message = f"{comp.id} flexible — verificar pose inicial antes de MD",
                    target  = comp.id,
                    action  = desc.sampling_recommendation,
                ))

            # ── Dificultad de parametrización ─────────────────────────────────
            if val and val.data:
                pdiff = val.data.get("parametrization_difficulty", "unknown")
                reasoning.parametrization_difficulty = pdiff

                if pdiff == "high":
                    reasoning.needs_parametrization_review = True
                    reasoning.risks.append(Risk(
                        message  = (
                            f"{comp.id}: dificultad de parametrización alta "
                            f"con {ff_ligands} — requiere validación manual"
                        ),
                        target   = comp.id,
                        severity = Severity.HIGH,
                    ))
                    reasoning.recommendations.append(Recommendation(
                        message = f"Revisar penalización de ParamChem para {comp.id}",
                        target  = comp.id,
                        action  = (
                            "Penalización > 10 indica parámetros no confiables. "
                            "Considerar QM (Gaussian/ORCA) para cargas y constantes de fuerza."
                        ),
                    ))
                elif pdiff == "medium":
                    reasoning.recommendations.append(Recommendation(
                        message = f"{comp.id}: verificar parámetros de {ff_ligands} antes de producción",
                        target  = comp.id,
                        action  = "Revisar RMSD de geometría post-minimización vs. estructura de input.",
                    ))

            # ── Carga neta → verificar protonación ────────────────────────────
            if desc.net_charge != 0:
                reasoning.needs_protonation_check = True
                reasoning.warnings.append(Warning(
                    message  = (
                        f"{comp.id} tiene carga neta {desc.net_charge:+d}: "
                        "verificar estado de protonación a pH 7.4"
                    ),
                    target   = comp.id,
                    severity = Severity.MEDIUM,
                ))

            # ── Planitud + rol competitive_ligand → π-stacking ────────────────
            if desc.is_globally_planar and comp.role == "competitive_ligand":
                reasoning.recommendations.append(Recommendation(
                    message = (
                        f"{comp.id} es plano (RMS={desc.global_planarity_rms:.3f}Å) "
                        "— verificar interacciones π-stacking en sitio activo"
                    ),
                    target  = comp.id,
                    action  = (
                        "Analizar residuos Phe/Tyr/Trp/His en sitio activo. "
                        "Usar contact_map y distance_analysis en la simulación."
                    ),
                ))

            # ── Sistema competitive-inhibition: sustrato flexible es riesgo ──
            if sys_type == "competitive-inhibition" and comp.role == "substrate":
                if desc.flexibility_class in ("flexible", "very_flexible"):
                    reasoning.risks.append(Risk(
                        message  = (
                            f"Sustrato flexible ({desc.n_rotatable_bonds} rot. bonds) "
                            "en competitive-inhibition: la pose en sitio activo puede variar"
                        ),
                        target   = comp.id,
                        severity = Severity.MEDIUM,
                    ))

        # ── Proteínas: reasoning desde validation ──────────────────────────────
        elif comp.role in _PROTEIN_ROLES and val is not None and val.data:
            pdata = val.data

            if pdata.get("missing_residues"):
                n_miss = len(pdata["missing_residues"])
                reasoning.risks.append(Risk(
                    message  = (
                        f"{comp.id}: {n_miss} residuo(s) faltante(s) — "
                        "modelado de loops requerido antes de simulación"
                    ),
                    target   = comp.id,
                    severity = Severity.HIGH,
                ))

            if pdata.get("likely_oligomer"):
                reasoning.warnings.append(Warning(
                    message  = (
                        f"{comp.id}: oligomerización probable — "
                        f"definir unidad de simulación ({pdata.get('oligomeric_state', '?')})"
                    ),
                    target   = comp.id,
                    severity = Severity.MEDIUM,
                ))

        comp.reasoning = reasoning

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# Etapa 6 — Reasoning global (consolidación)
# ═══════════════════════════════════════════════════════════════════════════════

def run_global_reasoning(state: SystemState) -> SystemState:
    """
    Etapa 6: consolida el estado del sistema en global_reasoning.

    Lee todos los componentes validados y razonados, y produce:
        - flags de decisión para el decision engine
        - resumen de issues globales
        - notas de compatibilidad entre componentes
    """
    gr = state.global_reasoning
    gr.inferred_system_type = state.inferred_system_type

    # Propagar warnings/risks de inference.py al global_reasoning
    gr.warnings.extend(state.warnings)
    gr.risks.extend(state.risks)
    gr.recommendations.extend(state.recommendations)

    # Contar componentes validados
    validated    = [c for c in state.components if c.validation and c.validation.is_valid]
    with_errors  = [c for c in state.components if c.validation and not c.validation.is_valid]
    gr.n_components_validated   = len(validated)
    gr.n_components_with_errors = len(with_errors)

    # ── Flag: has_blocking_errors ─────────────────────────────────────────────
    gr.has_blocking_errors = (
        len(with_errors) > 0
        or any(r.severity == Severity.HIGH for r in gr.risks)
        or any(
            r.severity == Severity.HIGH
            for c in state.components
            for r in c.all_risks
        )
    )

    # ── Flag: needs_special_sampling ─────────────────────────────────────────
    gr.needs_special_sampling = any(
        c.reasoning and c.reasoning.needs_special_sampling
        for c in state.components
    )

    # ── Flag: system_is_ready ────────────────────────────────────────────────
    gr.system_is_ready = (
        len(validated) == sum(1 for c in state.components if c.file)
        and not gr.has_blocking_errors
    )

    # ── Nota de estado general ────────────────────────────────────────────────
    if gr.system_is_ready:
        gr.notes.append(
            f"Sistema listo: {len(validated)} componente(s) validado(s), "
            f"tipo={gr.inferred_system_type}, sin errores bloqueantes."
        )
    else:
        issues = []
        if with_errors:
            issues.append(f"{len(with_errors)} componente(s) con errores")
        if gr.has_blocking_errors:
            issues.append("riesgos HIGH presentes")
        gr.notes.append(
            f"Sistema NO listo: {'; '.join(issues)}. "
            "Resolver antes de continuar al decision engine."
        )

    state.global_reasoning = gr
    return state


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline público
# ═══════════════════════════════════════════════════════════════════════════════

def parse_yaml(path: str | Path) -> SystemState:
    """
    Pipeline completo: YAML → SystemState enriquecido.

    Etapas:
        1. Carga y valida el YAML → SystemState (solo config)
        2. run_inference()             → inferred_system_type + risks globales
        3. run_validation()            → component.validation
        4. run_descriptors()           → component.descriptors
        5. run_component_reasoning()   → component.reasoning
        6. run_global_reasoning()      → global_reasoning (flags finales)
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Archivo no encontrado: {path}")

    if path.suffix not in (".yaml", ".yml"):
        raise ValueError(f"El archivo debe ser YAML: {path}")

    # ── Etapa 1: YAML → config ────────────────────────────────────────────────
    raw = _load_yaml(path)
    try:
        state = SystemState(**raw)
    except Exception as e:
        raise ValueError(f"Error validando YAML:\n{e}")

    # ── Etapa 2: Inferencia biológica ─────────────────────────────────────────
    state = run_inference(state)

    # ── Etapa 3: Validación de archivos ───────────────────────────────────────
    state = run_validation(state)

    # ── Etapa 4: Descriptores moleculares ────────────────────────────────────
    state = run_descriptors(state)

    # ── Etapa 5: Reasoning por componente ────────────────────────────────────
    state = run_component_reasoning(state)

    # ── Etapa 6: Reasoning global ─────────────────────────────────────────────
    state = run_global_reasoning(state)

    return state
