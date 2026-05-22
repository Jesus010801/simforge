# Estado del Proyecto — SimForge

## Fecha
[20/05/2026]

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
→ Manifest-driven execution: _initialize_state() debería leer orden de steps desde
  un manifest.json generado por WorkspaceBuilder, no desde scan de filesystem.
  Previene recurrencia de bugs como el de stale dirs.

---

Próximo milestone lógico

Conectar DAG → executor vía depends_on en metadata.json

Pasos concretos:
1. Añadir depends_on: list[str] a metadata.json en cada step builder
   (el campo blocking ya existe, depends_on es el que falta)
2. Modificar _is_blocked() en base_executor para leer depends_on desde metadata
3. Test: verificar que un step bloqueado solo por su predecesor directo
   no bloquea steps de ramas paralelas

Beneficio: el DAG pasa de ser una estructura de planificación a una estructura
que dirige la ejecución real. Desbloquea parallel waves en el futuro.