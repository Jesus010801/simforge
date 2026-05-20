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