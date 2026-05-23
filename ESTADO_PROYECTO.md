# Estado del Proyecto — SimForge

## Fecha
[23/05/2026]

---

# Estado actual del proyecto

SimForge ya funciona como un:

## workflow compiler molecular

Capaz de:

```text
YAML
↓
semantic reasoning
↓
workflow planning
↓
DAG orchestration
↓
workspace generation
↓
simulation artifact materialization
Arquitectura actual
simforge/
├── builders/
│   ├── __init__.py
│   ├── builder_registry.py
│   ├── workspace_builder.py
│   ├── test_workspace_builder.py
│   └── step_builders/
│       ├── __init__.py
│       ├── minimization_builder.py
│       ├── equilibration_builder.py
│       ├── production_builder.py
│       └── analysis_builder.py
│
├── configs/
│   └── hmg_competition.yaml
│
├── core/
│   ├── ontology.py
│   ├── models.py
│   ├── execution_models.py
│   ├── inference.py
│   ├── parser.py
│   ├── decision_engine.py
│   ├── compiler.py
│   ├── compiler_models.py
│   ├── test_parser.py
│   ├── test_decision_engine.py
│   └── test_compiler.py
│
├── descriptors/
│   ├── __init__.py
│   ├── topology.py
│   ├── aromaticity.py
│   ├── flexibility.py
│   ├── geometry.py          ← pendiente/refactor futuro
│   └── polarity.py          ← pendiente/refactor futuro
│
├── pipelines/
│   ├── __init__.py
│   ├── base_pipeline.py
│   ├── md_pipeline.py
│   └── inhibition_pipeline.py
│
├── validators/
│   ├── protein_validator.py
│   ├── ligand_validator.py
│   └── ligand_parsers/
│       ├── __init__.py
│       ├── sdf_parser.py
│       └── pdb_parser.py
│
├── workflows/
│   └── workflow_graph.py
│
└── simforge_runs/
    └── competitive-inhibition/
Completado — Core Compiler Architecture
compiler.py

API pública principal:

compiler = SimulationCompiler()

result = compiler.compile(
    "configs/hmg_competition.yaml"
)

Pipeline completo:

YAML
↓
SystemState
↓
SimulationPlan
↓
WorkflowGraph
↓
CompilationResult
Completado — Pipeline System
pipelines/

Separación formal entre:

infraestructura
vs
estrategia científica

Pipelines actuales:

MDPipeline
InhibitionPipeline

Sistema preparado para:

docking workflows
membrane workflows
free energy workflows
QM/MM workflows
Completado — WorkflowGraph
workflows/workflow_graph.py

Funcionalidades:

DAG formal
validación de dependencias
orden topológico
Mermaid export
user workflow view
execution ordering
Completado — WorkspaceBuilder
builders/workspace_builder.py

Generación automática de:

simforge_runs/
└── competitive-inhibition/
    ├── workflow/
    ├── metadata/
    ├── reports/
    └── steps/

Exporta automáticamente:

workflow.mmd
workflow.txt
summary.json
Completado — Builder Registry
builders/builder_registry.py

Sistema dinámico:

SimulationStep
↓
builder registry
↓
step builder
↓
artifact generation

Dispatch automático por:

step.stage.value
Completado — Step Builders
minimization_builder.py

Genera:

em.mdp
run.sh
metadata.json
equilibration_builder.py

Genera:

nvt.mdp
npt.mdp
run_nvt.sh
run_npt.sh
metadata.json
production_builder.py

Genera:

md.mdp
run_md.sh
metadata.json

Incluye:

PME
Parrinello-Rahman
trajectory compression
checkpoint continuity
analysis_builder.py

Genera:

run_analysis.sh
analysis_config.json
outputs/
plots/
tables/

Análisis actuales:

rmsd
hydrogen_bonds
distance_analysis
Estado actual VALIDADO

Workspace generado correctamente:

prepare
↓
parametrization
↓
assembly
↓
minimization
↓
equilibration
↓
REST2
↓
production
↓
analysis

con artefactos físicos reales.

Hallazgos arquitectónicos importantes
SimForge ya NO es:
parser molecular

Ahora es:

workflow compiler platform
Separación formal lograda
core/
→ infraestructura universal

pipelines/
→ estrategias científicas

builders/
→ materialización física

workflows/
→ DAG y ejecución lógica

descriptors/
→ percepción fisicoquímica

validators/
→ integridad estructural
Próximo gran milestone
executors/

Nueva fase:

workspace
↓
executor
↓
runtime state
↓
logging
↓
failure detection
↓
adaptive reasoning

Arquitectura futura:

executors/
├── base_executor.py
├── shell_executor.py
├── gromacs_executor.py
└── execution_state.py
Objetivo siguiente

Primer execution engine:

executor.run_workspace(...)

Inicialmente:

dry-run
logging
state tracking
subprocess orchestration

SIN ejecutar MD real todavía.

Estado conceptual actual

SimForge ya puede:

✅ interpretar workflows científicos
✅ construir DAGs ejecutables
✅ generar workspaces reproducibles
✅ materializar simulaciones GROMACS
✅ materializar análisis científicos
✅ organizar execution order correctamente

---

Completado — Execution Engine
executors/

Capa de ejecución completada con:

base_executor.py
→ ABC con contrato _run_step()
→ loop topológico sobre steps
→ detección de blocking failures
→ serialización de estado a execution_state.json

shell_executor.py
→ ejecución real via subprocess
→ dry-run mode
→ detección de outputs esperados

gromacs_executor.py
→ extiende ShellExecutor con percepción GROMACS
→ GROMACSLogParser: parsea convergencia, crashes, LINCS warnings
→ GROMACSStepDiagnostic: verdict system (converged/crashed/warning)
→ OutputFileStatus: valida archivos de output por tamaño

execution_state.py
→ WorkspaceExecutionState: estado completo serializable
→ StepExecutionRecord: registro por step con timing, stdout/stderr, outputs
→ remediations: historial de remediación embebido por step
→ gromacs_diagnostic: diagnóstico GROMACS serializado como dict (evita dependencia circular)

adaptive_reasoner.py
→ pattern matching sobre stdout/stderr crudo → DiagnosisResult → RemediationPlan
→ dos capas: detección de señales (text) + generación de plan (por categoría)
→ heurísticas: LINCS, NaN, Fmax no convergido, poor equilibration, timeout, missing output

adaptive_reasoning.py
→ razonamiento sobre WorkspaceExecutionState + GROMACSStepDiagnostic estructurado
→ produce AdaptiveReasoningResult con veredicto global y análisis por step
→ niveles: CONTINUE / REVIEW / REMEDIATE / ABORT

remediation_executor.py
→ loop completo: execute → diagnose → plan → apply → retry
→ aplica acciones: PATCH_MDP, SCALE_TIMESTEP, INJECT_RESTRAINTS, RESET_STEP, WRITE_FILE, etc.
→ desbloqueador de steps dependientes post-remediación

Validación realizada:
dry-run end-to-end sobre hmg_competition.yaml
14 steps completados, 2 skipped (manual/external), 0 failed

---

Fix — Persistencia de remediation history [2026-05-22]

Problema resuelto:
_REMEDIATION_HISTORY era un dict global a nivel de módulo.
Volátil entre sesiones, no persistido, no restartable, invisible al estado serializado.

Causa raíz:
El campo remediations en StepExecutionRecord ya existía en execution_state.py,
pero remediation_executor.py nunca lo usó — siguió usando el dict global como
solución provisional que nunca se eliminó.

Cambios aplicados en remediation_executor.py:
- eliminado _REMEDIATION_HISTORY (dict global)
- eliminados _get_remediations() y _add_remediation() (helpers de módulo)
- record.n_remediations() como fuente de verdad para conteo de intentos
- record.add_remediation(rem_record) para persistir cada ciclo
- _save_state(state) inmediatamente después de apply(), antes de retry()
  → garantiza restartability: si retry() crash, la remediación aplicada ya está en disco
- record.remediations[-1] = rem_record.model_dump() para actualizar el resultado del retry
- remediation_history(record) ahora acepta StepExecutionRecord, reconstruye desde dicts
- full_report(state) ahora acepta WorkspaceExecutionState, itera state.steps

Arquitectura validada:
El estado de remediación es ahora completamente local al record y se persiste
con él en execution_state.json. Sin estado global, sin pérdida en reinicios.

---

Fix — Duplicación de steps en dry-run [2026-05-22]

Problema resuelto:
Steps como prepare_ligand_1, parametrize_substrate_1, analysis_rmsd
aparecían ejecutados dos veces en el dry-run.

Causa raíz (dos factores):

1. WorkspaceBuilder acumulaba directorios stale:
   steps_dir.mkdir(parents=True, exist_ok=True) no limpiaba contenido previo.
   Cada llamada a builder.build() añadía nuevos dirs encima de los anteriores.

2. Detonante: el commit 567f4cd introdujo desempate alfabético en _topological_sort().
   Antes, el orden era FIFO (inserción en plan.steps): prepare_protein_1 primero.
   Con desempate alfabético: prepare_ligand_1 primero.
   Resultado: el mismo step_id aparece en dos dirs con distinto prefijo numérico
   (ej: 01_prepare_protein_1 del run anterior + 02_prepare_protein_1 del run nuevo).

3. Víctima: _initialize_state() en base_executor lee todo lo que encuentra en steps/
   sin filtrar ni validar contra el plan actual.

Cambio aplicado en workspace_builder.py:
    if steps_dir.exists():
        shutil.rmtree(steps_dir)
    steps_dir.mkdir(parents=True)

Arquitectura validada:
DAG y planner están sanos — WorkflowGraph._build_graph() lanza ValueError
en step_ids duplicados, no lo lanzó: el plan siempre fue correcto.
El problema vivía exclusivamente en la capa de materialización (builder → filesystem).

Dry-run post-fix: 16 dirs en steps/, 16 steps únicos, 0 duplicados.

---

Decisiones arquitectónicas registradas

[2026-05-22] Dos capas de reasoning con nombres distintos
adaptive_reasoner.py → trabaja sobre texto crudo (stdout/stderr)
adaptive_reasoning.py → trabaja sobre estado estructurado (WorkspaceExecutionState + GROMACSStepDiagnostic)
Son capas distintas y válidas. El naming actual es confuso (pendiente de rename).

[2026-05-22] Pipeline selection hardcodeada en compiler.py
_select_pipeline() usa if/elif sobre strings.
Identificado como deuda que crece con cada nuevo sistema soportado.
Fix propuesto: dict-based registry. No implementado aún.

[2026-05-22] DAG no conectado al executor
_is_blocked() en base_executor usa lógica conservadora secuencial:
cualquier fallo anterior bloquea todos los steps siguientes.
No lee depends_on desde metadata.json ni desde el WorkflowGraph.
La lógica de parallel waves y critical path del DAG es invisible durante ejecución real.

---

Deuda técnica abierta (identificada, no urgente)

ALTA PRIORIDAD:
→ Conectar DAG al executor: añadir depends_on en metadata.json de cada step builder,
  modificar _is_blocked() para leerlo. Sin esto, parallel waves del DAG no se usan.

MEDIA PRIORIDAD:
→ Rename de capas de reasoning para declarar distinción:
  adaptive_reasoner.py → signal_detector.py (o text_signal_reasoner.py)
  adaptive_models.py → execution_reasoning_models.py
→ Pipeline registry: reemplazar if/elif en compiler._select_pipeline() con dict dispatch.
→ RemediationPlan en adaptive_models.py tiene nombre idéntico al de remediation_models.py.
  Renombrar a ProposedRemediationPlan para evitar ambigüedad en imports.

BAJA PRIORIDAD (esperar segundo engine):
→ Abstracción de engine interface (SimulationEngine ABC)
→ GROMACS coupling en remediation_executor → encapsular cuando haya segundo engine

---

Completado — Manifest-driven execution [2026-05-22]

Milestone: el DAG compilado es ahora la fuente de verdad operativa del executor.

Problema resuelto:
El executor descubría steps haciendo filesystem scan sobre steps/.
El orden dependía del sort de directorios, no del DAG.
depends_on no existía en runtime — _is_blocked() era secuencial conservador.

Cambios aplicados:

executors/execution_state.py
→ añadido depends_on: list[str] = [] a StepExecutionRecord
→ campo aditivo, backward compatible: workspaces sin depends_on deserializan con []

builders/workspace_builder.py
→ genera metadata/execution_manifest.json al final de cada build
→ contiene: compiled_at, system_type, n_steps, lista ordenada de steps
→ cada entry: step_id, dir_name, stage, step_type, blocking, depends_on
→ dir_name es el nombre exacto del directorio creado (ej: "04_parametrize_ligand_1")
→ fuente: result.execution_order (orden topológico del DAG compilado)

executors/base_executor.py — _initialize_state()
→ tres paths en orden de preferencia:
  1. execution_state.json → resume (solo non-dry-run, comportamiento previo)
  2. execution_manifest.json → manifest-driven (nuevo, fuente del DAG)
  3. filesystem scan → fallback backward compat + log [WARN]
→ log en primera línea indica cuál path se usó

executors/base_executor.py — _is_blocked()
→ si record.depends_on está poblado: usa dependencias reales del DAG
  (any dep_status == FAILED → bloqueado)
→ si está vacío: fallback conservador secuencial (workspace antiguo sin manifest)

Validación realizada:
dry-run completo — primera línea del log:
  "Manifest cargado → 16 steps desde metadata/execution_manifest.json"

depends_on correctamente propagado en execution_state.json:
  parametrize_ligand_1   → ['prepare_ligand_1']
  assemble_system        → ['prepare_protein_1', 'prepare_substrate_1', ...]
  analysis_rmsd          → ['production_md']
  analysis_distance_analysis → ['production_md']
  analysis_hydrogen_bonds    → ['production_md']

Los tres steps de análisis tienen el mismo predecesor — semánticamente paralelos.
_is_blocked() ya los trataría correctamente si hubiera un executor paralelo.

Arquitectura validada:
compilation → materialization → execution están ahora formalmente separados:
  SimulationCompiler → produce execution_order (DAG topológico)
  WorkspaceBuilder   → materializa + serializa manifest desde execution_order
  BaseExecutor       → reconstruye runtime state desde manifest, no desde filesystem

---

Completado — Rename capas de reasoning [2026-05-22]

executors/adaptive_reasoner.py → executors/signal_detector.py
→ nombre refleja su responsabilidad real: detección de señales en texto crudo
→ clase AdaptiveReasoner conserva el nombre (API pública no cambia)

executors/adaptive_models.py → executors/execution_reasoning_models.py
→ nombre distingue estos modelos de los modelos de compilación/planning
→ RemediationPlan → ProposedRemediationPlan (resuelve clash con remediation_models.py)

Imports actualizados: adaptive_reasoning.py, remediation_executor.py,
test_adaptive_reasoning.py, remediation_models.py (comentario)
Archivos eliminados: adaptive_reasoner.py, adaptive_models.py
Validado: todos los imports OK, dry-run 14/14

---

Deuda técnica abierta (actualizada [2026-05-22])

ALTA PRIORIDAD:
(ninguna)

MEDIA PRIORIDAD:
(ninguna)

Completado — Resume + manifest enrichment [2026-05-22]

executors/base_executor.py — _initialize_state(), path 1 (resume)
→ después de deserializar execution_state.json, si el manifest existe:
  construye deps_map = {step_id: depends_on} desde manifest
  enriquece cada record con depends_on=[] usando el manifest como fuente
→ log "[RESUME] depends_on enriquecido en N records" cuando aplica
→ sin efecto si records ya tienen depends_on (idempotente)

Caso resuelto: workspace ejecutado con versión anterior del executor (sin manifest-driven PR)
reanudado con versión nueva → _is_blocked() usa dependencias reales del DAG, no fallback secuencial.

Validado: borrado manual de depends_on en 3 records con deps reales → todos recuperados desde manifest.

BAJA PRIORIDAD (esperar segundo engine):
→ Abstracción de engine interface (SimulationEngine ABC)
→ GROMACS coupling en remediation_executor → encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()

---

Completado — IR completo en SimulationPlan [2026-05-22]

Milestone: decision_engine.py es ahora el compilador semántico central.
SystemState + GlobalReasoning → SimulationPlan completamente declarativo.

Cambios aplicados:

core/execution_models.py — Step 1
→ WorkflowPolicy: decisiones científicas globales del pipeline
  (production_time_ns, equilibration_time_ns, minimization_steps, temperature_K,
   pressure_bar, timestep_ps, enhanced_sampling, sampling_method)
→ SimulationStep.params: dict[str, Any] = {} para IR engine-specific
→ SimulationStep.condition: Optional[str] para auditoría semántica
→ SimulationPlan.workflow_policy: WorkflowPolicy = WorkflowPolicy()

core/decision_engine.py — Step 2
→ _build_workflow_policy(state): produce WorkflowPolicy desde SystemState
  - competitive_binding → 50ns production
  - needs_special_sampling → REST2, 0.5ns equilibration, enhanced_sampling=True
→ _populate_step_params(plan, state): dispatch modular por StepStage
→ _STEP_PARAMS_BUILDERS: dict[StepStage, Callable] con builders por etapa:
  - MINIMIZATION: emtol=100.0 si ligando flexible, sino 1000.0
  - EQUILIBRATION: nsteps calculado desde policy, tc_grps dinámicos
  - PRODUCTION: nsteps=25M para 50ns@2fs, PME params completos
  - ANALYSIS: analysis_type + selection desde state.analysis
  - ASSEMBLY: solvate (box_type, water_model, water_gro) + add_ions (concentración)
  - ENHANCED_SAMPLING: method, n_replicas, temp range
→ condition añadida a: review_parametrization, validate_pose, rest2_sampling
→ builders NO modificados — siguen funcionando con fallback interno

Validación realizada:
→ energy_minimization.params: emtol=100.0, nsteps=50000 ✓
→ equilibration.params: 250k steps (0.5ns/0.002ps), tc_grps=Protein Non-Protein ✓
→ production_md.params: 25M steps (50ns/0.002ps) ✓
→ analysis_distance_analysis.params: selection propagada desde YAML ✓
→ workflow_policy: enhanced_sampling=True, sampling_method=REST2, production_time_ns=50.0 ✓
→ dry-run completo: 14 done, 2 skipped, 0 failed ✓

Arquitectura validada:
YAML → SystemState → WorkflowPolicy + step.params → SimulationPlan (IR declarativo)
Builders leen params — eliminar hardcodes es el Paso 3 (pendiente).

---

Deuda técnica abierta (actualizada [2026-05-22])

ALTA PRIORIDAD:
→ Paso 3: builders leen step.params con fallback al valor hardcodeado actual.
  Empezar por minimization_builder.py (más simple) para validar el patrón.
→ Rename de capas de reasoning (naming confuso, afecta legibilidad futura):
  adaptive_reasoner.py → signal_detector.py (trabaja sobre texto crudo)
  adaptive_models.py → execution_reasoning_models.py
→ Pipeline registry: reemplazar if/elif en compiler._select_pipeline() con dict dispatch.

MEDIA PRIORIDAD:
→ Paso 4: step_dir_map en builders para rutas inter-step correctas.
→ RemediationPlan en adaptive_models.py tiene nombre idéntico al de remediation_models.py.
  Renombrar a ProposedRemediationPlan.
→ Resume + manifest: al reanudar desde execution_state.json, enriquecer records con
  depends_on del manifest para correctness en _is_blocked().

BAJA PRIORIDAD (esperar segundo engine):
→ Abstracción de engine interface (SimulationEngine ABC)
→ GROMACS coupling en remediation_executor → encapsular cuando haya segundo engine

---

Completado — Paso 3: builders leen step.params [2026-05-22]

Milestone: los builders leen desde el IR en lugar de hardcodear valores.
La cadena YAML → decision_engine → step.params → artefactos MDP está completa.

Cambios aplicados:

minimization_builder.py
→ lee: integrator, emtol, emstep, nsteps desde step.params
→ fallback: valores anteriores si params vacío (backward compat)
→ metadata.json incluye params efectivos

equilibration_builder.py
→ lee: dt, nvt_nsteps, npt_nsteps, temperature, pressure, tc_grps, tau_t, ref_t, constraints
→ tc_grps/tau_t/ref_t dinámicos: len(tc_grps.split()) determina número de grupos
→ nvt.mdp y npt.mdp generados con valores del IR
→ metadata.json incluye params efectivos

production_builder.py
→ lee: dt, nsteps, temperature, pressure, tc_grps, tau_t, ref_t, constraints,
        nstxout_compressed, nstenergy, nstlog
→ metadata.json incluye params efectivos

Validación realizada:
→ em.mdp: emtol=100.0 (ligando flexible, viene del IR), nsteps=50000 ✓
→ nvt.mdp: nsteps=250000 (0.5ns policy), tc-grps desde IR ✓
→ md.mdp: nsteps=25000000 (50ns competitive_binding), ref_p=1.0 ✓
→ dry-run completo: 14 done, 0 failed ✓

Arquitectura validada:
YAML → SystemState → decision_engine → step.params → builder → MDP/artifacts
Cadena end-to-end completa. Builders son ahora renderizadores de templates, no fuentes de verdad.

---

Completado — Paso 5: validación end-to-end de propagación IR [2026-05-22]

Test arquitectónico: cambio de temperatura en YAML → artefactos físicos.

Cambio habilitante:
→ core/models.py: EnvironmentModel.temperature_K: float = 300.0
→ core/decision_engine.py: _build_workflow_policy lee state.environment.temperature_K
→ configs/hmg_competition.yaml: temperature_K como campo explícito bajo environment

Resultado con temperature_K: 310:
1. workflow_policy.temperature_K = 310.0                      ✓
2. equilibration.params["temperature"] = 310.0                ✓
3. nvt.mdp: ref_t = 310.0 310.0                               ✓
4. npt.mdp: ref_t = 310.0 310.0                               ✓
5. md.mdp:  ref_t = 310.0 310.0                               ✓
6. metadata.json serializa params efectivos por step           ✓
7. dry-run: 14 done, 0 failed                                  ✓

Conclusión arquitectónica:
El planner es ahora la única fuente de verdad científica.
Los builders no contienen parámetros ocultos.
La cadena IR → materialization está cerrada correctamente.

---

Completado — AssemblyBuilder + step_dir_map [2026-05-22]

Milestone: arquitectura declarativa consolidada end-to-end.
Todos los builders son ahora renderizadores de templates puros.

Cambios aplicados:

builders/workspace_builder.py — dos pasadas
→ Pass 1: crear todos los step_dirs + construir step_dir_map: dict[str, Path]
  (step_id → directorio absoluto del step)
→ Pass 2: materializar builders con step_dir_map completo
→ firma de llamada: builder.build(step, step_dir, step_dir_map)

builders/step_builders/assembly_builder.py — refactor completo
→ _build_assemble: itera step.depends_on, filtra via _component_gro(),
  genera variables bash PROTEIN_1, SUBSTRATE_1, LIGAND_1 con rutas relativas dinámicas
→ _build_solvate: lee box_type, box_distance, water_model, water_gro desde step.params;
  resuelve assemble_system dir via step_dir_map
→ _build_ions: lee concentration, positive_ion, negative_ion desde step.params;
  resuelve solvate_system y assemble_system (para topol.top) via step_dir_map
→ metadata.json serializa params efectivos en todos los sub-builders
→ helper _rel(from_dir, to_dir): os.path.relpath para rutas relativas entre steps
→ helper _component_gro(step_id): filtra deps relevantes para assembly
  (prepare_protein_* → *_processed.gro, parametrize_* → *.gro, resto → None)

Todos los demás builders — firma actualizada
→ step_dir_map: dict = {} añadido como kwarg (forward-compatible, ignorado por ahora)
→ builders afectados: minimization, equilibration, production, analysis,
  preparation, parametrization, validation, enhanced_sampling

Validación realizada:
→ assemble_system/run.sh: rutas generadas desde DAG, no hardcodeadas ✓
  PROTEIN_1="../02_prepare_protein_1/protein_1_processed.gro"
  SUBSTRATE_1="../05_parametrize_substrate_1/substrate_1.gro"
  LIGAND_1="../04_parametrize_ligand_1/ligand_1.gro"
→ solvate_system/run.sh: water_model=tip3p, spc216.gro, box=dodecahedron 1.2nm ✓
→ add_ions/run.sh: concentration=0.15M, NA/CL, topol.top desde assemble_system ✓
→ metadata.json serializa params efectivos en todos los assembly steps ✓
→ dry-run: 14 done, 0 failed ✓

Propiedad clave: si el DAG reordena steps, ninguna ruta se rompe.
Los paths en run.sh son calculados en tiempo de build desde el step_dir_map real,
no desde prefijos numéricos hardcodeados.

---

Estado arquitectónico actual

La cadena declarativa está completa end-to-end:

YAML
→ SystemState (parser)
→ WorkflowPolicy + step.params (decision_engine)
→ SimulationPlan IR
→ WorkflowGraph DAG
→ step_dir_map (WorkspaceBuilder pass 1)
→ builders como templates puros (WorkspaceBuilder pass 2)
→ artefactos físicos reproducibles

El planner es la única fuente de verdad científica.
Ningún builder contiene decisiones científicas ocultas.

---

Deuda técnica abierta (actualizada [2026-05-22])

ALTA PRIORIDAD:
→ Rename de capas de reasoning:
  adaptive_reasoner.py → signal_detector.py
  adaptive_models.py → execution_reasoning_models.py
→ Pipeline registry: reemplazar if/elif en compiler._select_pipeline() con dict dispatch.

MEDIA PRIORIDAD:
→ RemediationPlan en adaptive_models.py → renombrar a ProposedRemediationPlan.
→ Resume + manifest: enriquecer records de execution_state.json con depends_on del manifest.

BAJA PRIORIDAD:
→ Abstracción de engine interface (SimulationEngine ABC)
→ GROMACS coupling en remediation_executor → encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()

---

Completado — Pipeline registry [2026-05-22]

core/compiler.py
→ _PIPELINE_REGISTRY: dict[str, type[BasePipeline]] — dispatch por system_type
→ _select_pipeline() reducido a 2 líneas: lookup + fallback a MDPipeline
→ eliminado if/elif hardcodeado
→ comentarios explicitan los system_types futuros pendientes de pipeline propio
→ import limpiado: BasePipeline importado explícitamente (era implícito)
→ indentación rota del método original también corregida

Añadir un nuevo pipeline = 1 línea en el registry. No tocar compiler.py.

Validado:
→ competitive-inhibition → InhibitionPipeline ✓
→ unknown-system → MDPipeline (fallback) ✓
→ dry-run: 14 done, 0 failed ✓

---

Deuda técnica abierta (actualizada [2026-05-22])

ALTA PRIORIDAD:
→ Rename de capas de reasoning:
  adaptive_reasoner.py → signal_detector.py
  adaptive_models.py → execution_reasoning_models.py

MEDIA PRIORIDAD:
→ RemediationPlan en adaptive_models.py → renombrar a ProposedRemediationPlan.
→ Resume + manifest: enriquecer records de execution_state.json con depends_on del manifest.

BAJA PRIORIDAD:
→ Abstracción de engine interface (SimulationEngine ABC)
→ GROMACS coupling en remediation_executor → encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()

---

Completado — Rutas inter-step en run scripts [2026-05-22]

Problema resuelto:
Los scripts de equilibration y production referenciaban em.gro, npt.gro, topol.top
sin path — rotos si el usuario los ejecutaba desde el directorio del step.

Cambios aplicados:

builders/step_builders/_utils.py (nuevo)
→ rel(from_dir, to_dir): helper compartido para rutas relativas entre steps
→ reemplaza las copias locales en assembly_builder.py

builders/step_builders/assembly_builder.py
→ ahora importa _rel desde _utils (eliminada copia local)

builders/step_builders/equilibration_builder.py
→ resuelve EM_DIR desde step.depends_on (energy_minimization)
→ resuelve TOPOL_DIR desde step_dir_map["assemble_system"]
→ run_nvt.sh: -c "$EM_DIR/em.gro" -p "$TOPOL_DIR/topol.top"
→ run_npt.sh: -p "$TOPOL_DIR/topol.top" (nvt.gro es local → OK)

builders/step_builders/production_builder.py
→ resuelve EQ_DIR desde step.depends_on (equilibration)
→ resuelve TOPOL_DIR desde step_dir_map["assemble_system"]
→ run_md.sh: -c "$EQ_DIR/npt.gro" -t "$EQ_DIR/npt.cpt" -p "$TOPOL_DIR/topol.top"

Validado:
→ run_nvt.sh: EM_DIR="../10_energy_minimization", TOPOL_DIR="../07_assemble_system" ✓
→ run_npt.sh: TOPOL_DIR="../07_assemble_system", nvt.gro local ✓
→ run_md.sh: EQ_DIR="../11_equilibration", TOPOL_DIR="../07_assemble_system" ✓
→ dry-run: 14 done, 0 failed ✓

Propiedad: si el DAG reordena steps, los números en los paths se recalculan automáticamente.

---

Completado — CLI [2026-05-22]

cli.py — punto de entrada con Typer + Rich

simforge compile <yaml> [--output-dir <dir>] [--no-build]
→ compila YAML → SimulationPlan + DAG
→ muestra panel con system_type, steps, policy, blocking issues
→ tabla con execution_order completo
→ materializa workspace y muestra path
→ sugiere next command

simforge run <workspace> [--dry-run/--real] [--executor shell|gromacs]
→ dry-run por defecto (--real requiere confirmación interactiva)
→ ejecuta executor completo con output de logs en tiempo real
→ resumen: done/failed/skipped al terminar

simforge status <workspace> [-v]
→ lee execution_state.json
→ enriquece stage desde manifest (no almacenado en state)
→ tabla con iconos (✓ ✗ – ⊘ ▶) por step
→ barra de progreso visual [█████░░░░░] N/M steps (X%)
→ verbose: timing + nota de error por step

Validado:
→ compile: panel + tabla + workspace materializados ✓
→ run --dry-run: 14 done, 2 skipped, 0 failed ✓
→ status: tabla con stage, iconos, barra de progreso ✓
→ status --verbose: timing por step ✓

---

Completado — Analysis builders con comandos GROMACS reales [2026-05-22]

builders/step_builders/analysis_builder.py — reescrito completo

Dispatch por analysis_type desde step.params:
→ rmsd:               gmx rms  — Backbone (4) reference, Protein (1) RMSD
→ rmsf:               gmx rmsf — per-residuo, Backbone + C-alpha
→ hydrogen_bonds:     gmx hbond — Protein-Protein (1 1) + Protein-SOL (1 13)
→ distance:           gmx distance + make_ndx.sh scaffold de guía
→ energy:             gmx energy — Potential, Kinetic, Total, Temp, Pressure
→ radius_of_gyration: gmx gyrate — Protein (1)
→ genérico:           template vacío con paths de referencia

Propiedades arquitectónicas:
→ Grupos: GROMACS built-in exclusivamente (sin make_ndx para análisis estándar)
→ PROD_DIR resuelto desde step_dir_map["production_md"] — relativo y DAG-correcto
→ metadata.json incluye params + gromacs_groups + outputs esperados
→ Fallback a _build_generic para tipos no implementados (sin crash)

builders/step_builders/_utils.py (nuevo)
→ rel(from_dir, to_dir): helper compartido (importado desde todos los builders)

Validado:
→ analysis_rmsd/run_analysis.sh:        echo "4 4" | gmx rms ... PROD_DIR="../13_production_md" ✓
→ analysis_hydrogen_bonds/run_analysis.sh: echo "1 1" | gmx hbond ... ✓
→ analysis_distance_analysis/run_analysis.sh: gmx distance + make_ndx.sh ✓
→ dry-run: 14 done, 2 skipped (manual), 0 failed ✓
→ simforge status: tabla completa con stage, iconos, barra de progreso ✓

---

Estado arquitectónico final [2026-05-22]

La cadena declarativa está completa end-to-end:

YAML
→ SystemState (parser)
→ WorkflowPolicy + step.params (decision_engine — única fuente de verdad científica)
→ SimulationPlan IR
→ WorkflowGraph DAG
→ step_dir_map (WorkspaceBuilder pass 1)
→ builders como templates puros con rutas DAG-dinámicas (WorkspaceBuilder pass 2)
→ execution_manifest.json (fuente de verdad operativa del executor)
→ BaseExecutor con manifest-driven ordering + DAG-aware _is_blocked()
→ artefactos físicos reproducibles: MDP, run scripts, analysis scripts

CLI operativo: simforge compile | run | status

---

Fix — Protocolo de solvatación y hardware-aware mdrun [2026-05-22]

Problema raíz:
Péptido de ~80 aa generaba caja enorme (~500k moléculas de agua).

Causa: editconf sin -princ — la caja se calcula sobre la orientación original del péptido,
no sobre sus ejes principales. Para péptidos elongados esto infla la caja masivamente.

Cambios aplicados:

builders/step_builders/assembly_builder.py
→ editconf incluye -princ: orienta el péptido en sus ejes principales antes de calcular la caja
→ box_type: dodecahedron → triclinic (del protocolo de referencia Dinámica_molecular)

core/decision_engine.py
→ solvate_system: box_type="triclinic", box_distance=1.2nm
→ add_ions: concentration default 0.154M (referencia fisiológica estándar)

---

Completado — Hardware-aware mdrun [2026-05-22]

Problema: scripts generaban mdrun con flags GPU hardcodeados.
Usuarios sin GPU no podían ejecutar los scripts sin edición manual.

Solución: hybrid auto-detection — parámetro + detección en runtime.

builders/step_builders/_utils.py — mdrun_block(deffnm, hardware)
→ hardware="gpu":  flags GPU completos del protocolo de referencia
    -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu -ntomp 10 -nstlist 150
    -pin on -tunepme no -pmefft gpu -dlb auto
→ hardware="cpu":  fallback CPU con nproc dinámico
    -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on
→ hardware="auto": bloque bash con detección nvidia-smi en runtime (default)

core/execution_models.py
→ WorkflowPolicy.hardware: str = "auto"  — el usuario declara "gpu"|"cpu"|"auto" en el YAML

core/decision_engine.py
→ policy.hardware propagado a params de: minimization, equilibration, production

builders: minimization_builder, equilibration_builder, production_builder
→ todos importan mdrun_block desde _utils
→ leen hardware desde step.params con default "auto"
→ grompp incluye -maxwarn 1 en NVT, NPT y producción

equilibration_builder.py — refactorizado (soporte membrana)
→ pcoupltype: isotropic (proteína) | semiisotropic (membrana) — controlado desde params
→ pcoupl_npt, ref_p_xy, ref_p_z, tau_p, rcoulomb, rvdw, disp_corr desde step.params
→ builders no contienen lógica condicional: el decision_engine decide, el builder renderiza

Validado: 9 tests, 0 failed.

---

Completado — simforge init (YAML wizard) [2026-05-22]

cli.py — nueva subcomanda `simforge init`

Wizard interactivo CLI con Typer prompts:
→ Nombre del proyecto
→ Tipo de sistema: protein-in-water | protein-ligand | protein-membrane
→ Archivos PDB (proteína + ligandos opcionales, con rol)
→ Forcefield: opls-aa | charmm36 | amber99sb
→ Water model: spce | tip3p | tip4p
→ Temperatura: 300K | 309.65K | custom
→ Duración: 1 | 10 | 50 | 100 ns | custom
→ Análisis: multi-select (rmsd, rmsf, energy, hbonds, gyration, distance)
→ Hardware: auto | gpu | cpu
→ Genera configs/<nombre>.yaml listo para simforge compile

Helpers internos: _ask(), _ask_multi(), _ask_float(), _ask_str()
Objetivos auto-inferidos por system_type (competitive_binding, stability, etc.)

---

Completado — Packaging moderno PEP517/518 [2026-05-22]

Problema: pyproject.toml usaba setuptools.backends.legacy:build (no existe).
Fix: migrado a setuptools.build_meta.

Cambios:
pyproject.toml — build-backend correcto + include explícito (evita escanear venv/)
simforge/__init__.py — paquete instalable, __version__ = "0.1.0"
simforge/__main__.py — habilita python -m simforge
simforge/cli.py — entry point del paquete (importa cli.py raíz como `app`)
core/__init__.py, builders/__init__.py, validators/__init__.py — faltaban
requirements.txt — dependencias runtime + dev
.gitignore — expandido con Python estándar + artefactos GROMACS

Entry point: simforge = "simforge.cli:app"
Instalación: pip install -e .
Disponible desde cualquier directorio: simforge compile | run | dry-run | status | validate | inspect | init

---

Completado — CLI profesional: compile mejorado + validate + inspect + dry-run [2026-05-22]

cli.py — reescritura completa del comando compile + 3 nuevos comandos

compile — pipeline visual por etapas:
  [1/6] Parsing YAML          → timing + n_components + system_type
  [2/6] Scientific planning   → answers applied o skipped
  [3/6] Building plan         → n_steps planned
  [4/6] Computing DAG         → n_nodes, n_warnings, n_blocking
  [5/6] Materializing         → workspace path
  [6/6] Writing reports       → compile_report.md + execution_manifest.json
  → Panel System Summary: tipo, temperatura, duración, steps, componentes, runtime estimate
  → Panel Scientific Warnings colorizados (⚠ yellow, ✗ red, ◆ cyan)
  → Panel Workflow DAG: Rich Tree con colores por stage, ramas paralelas visibles
  → Panel ✓ Done: paths, tiempo total, next command
  → Errores: Panel rojo con causa + suggested fix (sin tracebacks)

validate — parse + compile sin escribir nada:
  → ✓/✗ por etapa, warnings científicos, blocking issues
  → Exit 1 si hay blocking issues (apto para CI)

inspect — muestra IR completo sin escribir archivos:
  → System Summary, Component Descriptors (heavy atoms, flexibility, charge, H-bond D/A)
  → Workflow DAG tree, Step Parameters IR snapshot (emtol, nsteps, temperature, box_type, hardware)

dry-run — alias de run con dry_run=True forzado
run — refactorizado, lógica compartida en _execute_workspace()

Helpers nuevos: _stage_ok(), _collect_scientific_warnings(), _build_dag_tree(),
_show_system_summary(), _show_scientific_warnings(), _write_compile_report()

_write_compile_report() genera workspace/metadata/compile_report.md:
  → tabla de workflow, parámetros científicos, warnings, metadata de reproducibilidad

---

Próximo milestone — simforge init (YAML wizard)

Estado actual: el usuario escribe el YAML manualmente.
Siguiente paso: `simforge init` — wizard interactivo CLI que genera el YAML.

Fases propuestas:

Fase 1 — wizard CLI (Typer prompts):
→ simforge init → preguntas sobre: tipo de sistema, componentes, temperatura, duración,
  hardware disponible, análisis deseados
→ genera configs/<nombre>.yaml listo para simforge compile

Fase 2 — templates predefinidos:
→ simforge init --template peptide-in-water
→ simforge init --template protein-ligand
→ simforge init --template protein-membrane

Fase 3 (futuro):
→ UI web simple para configuración visual
→ validación estructural integrada (verificar PDB antes de compilar)

Precondición: la arquitectura del YAML y SystemState están estabilizadas.

---

Deuda técnica abierta (actualizada [2026-05-22])

BAJA PRIORIDAD:
→ Abstracción de engine interface (SimulationEngine ABC) — esperar segundo engine
→ GROMACS coupling en remediation_executor — encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()
  (Precondición cumplida: _is_blocked() ya usa depends_on reales.
   Los tres análisis tienen depends_on=['production_md'] — listos para paralelo.)
→ production_builder.py: MDP usa valores hardcodeados pese a leer rcoulomb/rvdw/pcoupltype
  desde params — completar la conexión cuando se implemente soporte membrana en producción

---

Completado — Pipeline membrana end-to-end [2026-05-22]

Milestone: SimForge compila y materializa un workspace de 12 pasos para sistemas
proteína-membrana (DPPC + OPLS-AA). Todos los artefactos físicos son correctos.

## Arquitectura nueva

### Capas de conocimiento y adaptadores

core/membrane_knowledge.py
→ Fuente de verdad para constantes físicas: APL targets, lipid residue names,
  atom names, bilayer geometry, equilibration defaults
→ API pública: apl_target(), lipid_residue_name(), lipid_atom_names(),
  inflation_factor(), bilayer_for_box()
→ MEMBRANE_EQUILIBRATION_DEFAULTS: NVT 100ps, NPT Berendsen 1ns, prod dt=0.001
→ SHRINK_LOOP_EMTOL=1000.0, STRONG_POSRES_FC=(100000,100000,100000)

docs/architecture/membrane_runtime_design.md
→ Documento arquitectónico: 11 secciones, decisiones científicas por step,
  señales de fallo, módulos nuevos, MDP de referencia

adapters/base.py — contrato unificado para herramientas externas
→ AdapterError, ToolNotAvailableError, PreconditionError
→ AvailabilityResult, PreconditionViolation, AdapterResult (Pydantic)
→ ExternalToolAdapter ABC: check_availability, validate_preconditions, run,
  _which, _make_result

adapters/inflategro_adapter.py (implementación completa)
→ Wrappea inflategro-Jorge.pl vía subprocess
→ Metadata garantizada: apl_nm2, apl_ang2, scale_applied, lipid_name, n_lines_gro

adapters/water_deletor_adapter.py (reimplementación Python pura)
→ Reimplementa water_deletor.pl sin dependencia de Perl
→ Algoritmo: calcula fronteras Z de headgroup, elimina SOL con OW en interior
→ Metadata: waters_removed, atoms_in, atoms_out, z_top_nm, z_bot_nm
→ check_availability() siempre True (sin deps externas)

adapters/movememb_adapter.py (reimplementación Python pura)
→ Reimplementa MoveMemb.f sin gfortran
→ compute_z_shift(): alinea midplane de bicapa con Z-centro de proteína
→ run(): aplica shift, combina GROs, renumera átomos, escribe sistema combinado
→ Metadata: z_shift_nm, protein/bilayer Z extents, atoms_total

### Validators

validators/membrane_validators.py — 4 validadores, solo percepción (sin remediación)
→ validate_apl_convergence(area_dat, target, tolerance) → APLConvergenceResult
→ validate_no_overlap(gro, lipid_name, cutoff=0.2nm) → OverlapResult
→ validate_tm_orientation(gro, tm_ranges, threshold=30°) → OrientationResult
   (eje end-to-end CA, no power iteration — robusto con hélices alineadas con X/Y)
→ validate_no_water_in_bilayer(gro, headgroup_atom, tail_atom) → WaterConsistencyResult

### Builders extendidos

builders/step_builders/embedding_builder.py (nuevo)
→ Meta-step shrink loop: inflate(4.0) → minimize → deflate(0.95) loop
→ Genera: minim_shrink.mdp, run.sh auto-contenido, shrink_telemetry.json
→ APL parseado en Python desde area_2.dat (sin Fortran AperR)
→ converged.gro: output canónico; exit 1 si no converge
→ StepStage.MEMBRANE_EMBEDDING añadido a execution_models.py

builders/step_builders/minimization_builder.py — extendido para membrana
→ lee define ("-DPOSRES -DSTRONG_POSRES") desde step.params
→ lee rcoulomb, rvdw, disp_corr desde step.params (membrana: 1.2nm, EnerPres)
→ topol_ref: prefiere generate_topology > assemble_system

builders/step_builders/equilibration_builder.py — extendido para membrana
→ pcoupltype semiisotropic: ref_p = xy  z, compressibility = 4.5e-5  4.5e-5
→ pcoupl_npt, ref_p_xy, ref_p_z, tau_p, rcoulomb, rvdw, disp_corr desde params
→ Bugs corregidos: gen_vel=yes en NVT, continuation=yes en NPT
→ topol_ref: prefiere generate_topology > assemble_system

builders/step_builders/production_builder.py — extendido para membrana
→ tcoupl (Nose-Hoover), pcoupltype semiisotropic, tau_p, ref_p_xy, ref_p_z
→ rcoulomb, rvdw, disp_corr desde params
→ topol_ref: prefiere generate_topology > assemble_system

builders/step_builders/assembly_builder.py — extendido para membrana
→ solvate_membrane: copia topol.top localmente, corre gmx solvate desde converged.gro
→ clean_water: genera run_clean_water.py (Python puro), actualiza conteo SOL en topol.top
→ add_ions: resuelve input_gro desde clean_water > solvate step;
  topol.top desde clean_water > solvate_membrane > assemble_system
→ _FF_GROMACS_NAME: añadida entrada "opls-aa-membrane" → "oplsaa_membrane"

### Pipeline y compiler

pipelines/membrane_pipeline.py — MembraneWorkflowOPLSAA
→ pipeline_type = "protein-membrane"
→ 12 steps: 3 MANUAL + 9 AUTO/EXTERNAL + análisis
→ Steps 1-3: orient_protein, match_box_to_bilayer, embed_in_bilayer (MANUAL)
→ Steps 4-5: generate_topology (pdb2gmx oplsaa_membrane), membrane_embedding (shrink loop)
→ Steps 6-8: solvate_membrane, clean_water (EXTERNAL), add_ions
→ Steps 9-11: energy_minimization, equilibration, production_md
→ Parámetros físicos de membrane_knowledge.py exclusivamente

core/compiler.py — _PIPELINE_REGISTRY actualizado
→ "protein-membrane": MembraneWorkflowOPLSAA

### Benchmark

benchmarks/membrane_dppc_oplsaa/test_pipeline.py — 26 tests
→ Pipeline selection, step count, stages, steps manuales
→ 11 tests de parámetros MDP críticos: gen_vel, continuation, Berendsen,
  semiisotropic, dt=0.001, Nose-Hoover, Parrinello-Rahman, DispCorr
→ Broken case: YAML sin membrana no activa MembraneWorkflowOPLSAA
→ Smoke tests de adapters (always available)
→ Tests funcionales WaterDeletor y MoveMembAdapter con GROs construidos en memoria

## Validación end-to-end

python3 cli.py compile configs/membrane_test.yaml
→ 12 steps, workspace generado en simforge_runs/protein-membrane/
→ NVT: gen_vel=yes, tc-grps=system, rcoulomb=1.2, DispCorr=EnerPres ✓
→ NPT: continuation=yes, Berendsen, semiisotropic, ref_p=0.5 0.5 ✓
→ Producción: dt=0.001, Nose-Hoover, Parrinello-Rahman, semiisotropic, ref_p=1.0 1.0 ✓
→ EM: define=-DPOSRES -DSTRONG_POSRES, rcoulomb=1.2, DispCorr=EnerPres ✓
→ solvate_membrane/run.sh: copia topol.top, gmx solvate desde converged.gro ✓
→ clean_water/run_clean_water.py: Python puro, actualiza SOL en topol.top ✓
→ add_ions/run.sh: INPUT_DIR=../07_clean_water, TOPOL_DIR=../07_clean_water ✓
→ 94 tests, 0 failed ✓

## Cadena de topología (membrana)

generate_topology/topol.top (original)
  → copia en solvate_membrane/topol.top (SOL count añadido por gmx solvate)
  → modifica clean_water/topol.top (SOL count corregido por WaterDeletorAdapter)
  → modifica add_ions/run.sh (ion count añadido por gmx genion)
  → usa energy_minimization (TOPOL_DIR=../08_add_ions)

---

Deuda técnica abierta (actualizada [2026-05-22])

BAJA PRIORIDAD:
→ Abstracción de engine interface (SimulationEngine ABC) — esperar segundo engine
→ GROMACS coupling en remediation_executor — encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()
→ InflateGRO shrink loop: ejercitar contra GROMACS real con DPPC bilayer
→ WaterDeletorAdapter: probar con GRO real de sistema membrana solvado
→ MoveMembAdapter: probar con par proteína/bicapa reales
→ EmbeddingBuilder: path a inflategro-Jorge.pl configurable (actualmente hardcodeado)
→ scientific_planner.py: preguntas específicas para sistemas de membrana
  (selección de lípido, verificación de orientación, protocolo)

---

Completado — Runtime observable y resiliente [2026-05-23]

Milestone: el runtime ya puede ejecutar pipelines GROMACS completos (lysozyme en agua,
hasta production_md) de forma autónoma, con cache inteligente y UX de recovery.

## runtime/ package — 8 módulos nuevos

runtime/events.py
→ ExecutionEvent, EventType (21 tipos), EventSeverity, EventBus, BoundPublisher
→ Pub/sub estructurado para todos los eventos de ejecución

runtime/journal.py
→ JournalWriter: log append-only JSONL en metadata/execution_journal.jsonl
→ Skips HEARTBEAT + STDOUT para no inflar el journal con ruido

runtime/stream.py — AsyncProcessRunner (reescrito completo)
→ asyncio.create_subprocess_exec con limit=10MB (evita LimitOverrunError en logs GROMACS)
→ _iter_lines(): atrapa asyncio.LimitOverrunError, lee bytes y continúa (never crash)
→ Heartbeat como cancellable task (NO en gather — evita deadlock con proc.wait())
→ Execution correctness (returncode) completamente desacoplada de observability (stdout)
→ drain_results = await asyncio.gather(..., return_exceptions=True) → streaming falla = WARNING, no crash

runtime/metrics.py
→ SystemMetricsCollector: psutil + nvidia-smi, background async task
→ SystemSnapshot, GROMACSPerformance (ns/day parseado de stdout)

runtime/artifacts.py
→ ArtifactRef (SHA-256, semantic_role, step_id), ArtifactRegistry → metadata/artifact_registry.json

runtime/cache.py
→ StepCacheManager: fingerprint = SHA-256(params + input checksums)
→ metadata/step_cache.json, is_cached(), record(), invalidate()

runtime/executor.py — RuntimeExecutor(BaseExecutor)
→ ResumePlan: reusable | resume_point | pending | invalidated
→ _pre_run_hook(): _plan_resume() antes del loop → muestra summary de resume plan
→ _should_skip(): verifica _cached_steps + emite CACHE_HIT event
→ Cache integrity: cache hit SOLO si fingerprint matches AND todos los expected_outputs existen físicamente
→ Auto-invalidación de entries stale (fingerprint OK pero artifacts desaparecidos)

## Fixes críticos

### Cache integrity (crítico)
Problema: "⚡ cache hit" aunque add_ions/aaions.gro no existía en disco.
Fix: después de is_cached(), verificar todos los expected_outputs; si falta alguno → invalidate() + rerun.

### Streaming resilience (crítico)
Problema: production_md crash tras ~6 min con "Separator is not found, chunk exceed limit"
Causa: asyncio readline() con buffer 64KB vs logs GROMACS de >64KB por línea.
Fix: limit=10MB + LimitOverrunError caught en _iter_lines() + heartbeat desacoplado.

### DAG blocking (crítico)
Problema: _is_blocked() solo chequeaba == FAILED, no BLOCKED → bloqueo transitivo roto.
Fix: check in {FAILED, BLOCKED} en frozenset.

## Resume / recovery UX

Antes de ejecutar: _pre_run_hook() clasifica todos los steps:
  ✓ reusable   — fingerprint match + artifacts presentes → silently skipped
  ↻ resume_point — primer step que debe correr
  ○ pending    — steps downstream del resume point

CLI muestra panel estructurado en lugar de ⚡ por step.

---

Completado — Persistent project/workspace semantics [2026-05-23]

Milestone: cada compile crea un run timestamped inmutable. Los runs anteriores NUNCA se borran.

## Estructura de directorios

{output_dir}/{project_name}/
  project.json                      ← registro de todos los runs
  runs/
    2026-05-23_14-30-00/            ← un directorio por compile
      metadata/
        execution_manifest.json
        run_info.json               ← provenance del run
        compile_report.md
        planning_session.json
      steps/
      inputs/

## Módulos nuevos/modificados

core/project_manager.py (nuevo)
→ ProjectManager.project_dir(), create_run_dir(), update_project_registry(), get_run_history()
→ run_dir = {output_dir}/{project_name}/runs/{YYYY-MM-DD_HH-MM-SS}/

builders/workspace_builder.py
→ añadido workspace_path parameter (override de output_dir+workspace_name)
→ ELIMINADO shutil.rmtree(steps_dir) → reemplazado con exist_ok=True
→ backward compat: output_dir+workspace_name siguen funcionando

cli.py compile
→ crea project_dir y run_dir antes de build
→ pasa workspace_path=run_dir a WorkspaceBuilder
→ escribe metadata/run_info.json + actualiza project.json
→ panel final muestra "Project →" y "Run →" separados con run count

---

Completado — Semantic objective system [2026-05-23]

Milestone: SimForge ya NO falla con "Objetivo no reconocido". En cambio normaliza,
infiere, sugiere y aplica presets científicos.

## Módulos nuevos

core/semantic_objectives.py
→ CANONICAL_OBJECTIVES: 11 objetivos autorizados con descripciones
→ _ALIASES: 45+ strings informales → lista de canonicals
  (espacios, guiones, acrónimos, sinónimos, multi-canonical todos funcionan)
→ SIMULATION_PRESETS: 8 perfiles named con objectives + workflow hints:
  membrane_protein, soluble_protein, protein_ligand_binding, idr_sampling,
  competitive_inhibition, membrane_lipid_dynamics, allosteric_modulation, ion_channel
→ normalize_objective(text) → (canonicals, note) — alias → fuzzy → [] (sin excepción)
→ suggest_objectives(text) → top 3 canonical para mensajes de error/autocomplete

core/semantic_inference.py
→ detect_membrane_protein_signals(state) → lista de señales detectadas
  (biological_context=transmembrane/membrane_associated, membrane.enabled, lipid roles, mem objectives)
→ run_semantic_normalization(state) → stage 1.5 del pipeline:
  1. Aplica simulation_profile preset (si especificado)
  2. Normaliza simulation_objectives (aliases → canonical)
  3. Auto-inyecta membrane_perturbation + stability si señales de membrana detectadas
  4. Registra notas en global_reasoning + emite Warning con suggestions para unknowns

## Archivos modificados

core/models.py
→ SystemState.simulation_profile: Optional[str] = None (nuevo campo YAML)
→ validate_objectives: normaliza inline con aliases (sin ValueError) — unknowns pasan a semantic stage

core/ontology.py
→ SIMULATION_GOALS derivado de CANONICAL_OBJECTIVES (single source of truth, backward compat)

core/parser.py
→ stage 1.5: run_semantic_normalization() insertado entre YAML load e inferencia biológica

cli.py
→ _show_parse_error(): detecta errores de objetivo y añade suggestions + preset hint
→ _show_semantic_normalization_notes(): panel "Semantic Normalization" después de parse

## Ejemplos de normalización

'membrane_protein_dynamics'  → ['membrane_perturbation', 'stability']
'protein stability'          → ['stability']
'binding affinity'           → ['binding']
'idr'                        → ['conformational_sampling']
'pore_formation'             → ['membrane_perturbation', 'permeability']
'completely_unknown'         → []  suggestions: ['competitive_binding', 'conformational_sampling', 'stability']

## Tests

core/test_semantic_inference.py — 38 tests
→ TestNormalizeObjective, TestSuggestObjectives, TestPresets
→ TestMembraneSignalDetection, TestSemanticNormalizationStage
→ TestSystemStateObjectiveValidator, TestSimulationGoalsCompat

---

Estado de tests [2026-05-23]

300 tests totales, 0 fallos:
→ core/test_semantic_inference.py    — 38 tests (semantic objectives + inference)
→ core/test_geometry_advisor.py      — 21 tests (GeometryAdvisor: aspect ratio, advisories)
→ runtime/test_runtime.py            — 29 tests (cache integrity, stream robustness, events)
→ builders/step_builders/test_mdrun_stage_compat.py — 29 tests (GPU flags por stage)
→ benchmarks/membrane_dppc_oplsaa/test_pipeline.py — 26 tests (pipeline membrana)
→ executors/test_dag_blocking.py     — 10 tests (bloqueo transitivo)
→ builders/, core/, executors/ otros — 147 tests restantes

---

Deuda técnica abierta [2026-05-23]

CHECKPOINT-AWARE RECOVERY (pendiente, alta prioridad científica):
→ Nuevos estados: INTERRUPTED, RECOVERABLE, PARTIAL en StepStatus
→ Detección de .cpt en _plan_resume() → classify como RECOVERABLE
→ run_md_resume.sh (sin grompp, usa -cpi md.cpt -append)
→ run_nvt_resume.sh, run_npt_resume.sh equivalentes
→ _print_resume_plan() muestra categoría RECOVERABLE con "↻ checkpoint detected"
→ _find_script() prioriza scripts de resume cuando step en _resumable_steps

BAJA PRIORIDAD:
→ Abstracción engine interface (SimulationEngine ABC) — esperar segundo engine
→ GROMACS coupling en remediation_executor — encapsular cuando haya segundo engine
→ Parallel wave execution: ThreadPoolExecutor sobre to_parallel_waves()
→ InflateGRO shrink loop: ejercitar contra GROMACS real con DPPC bilayer
→ EmbeddingBuilder: path a inflategro-Jorge.pl configurable (actualmente hardcodeado)
→ scientific_planner.py: preguntas específicas para sistemas de membrana
→ Decision engine: leer workflow hints de simulation_profile presets
  (semiisotropic_coupling, conservative_timestep, etc. actualmente guardados en global_reasoning.notes
   pero el decision_engine aún no los lee)