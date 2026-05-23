# Estado del Proyecto — SimForge

## Fecha
[23/05/2026]

---

# Estado actual del proyecto

SimForge funciona como un:

## comparative MD study analyzer + workflow compiler

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
↓
multi-system comparative analysis    ← NUEVO
```

---

# Evolución arquitectónica completada hoy

## Study Layer — análisis comparativo multi-sistema

SimForge ya no asume `1 directorio = 1 simulación`.
Ahora soporta `1 directorio = estudio completo` con múltiples sistemas, réplicas y observables.

### Archivos nuevos

```
runtime/
├── study_models.py       ← modelos de datos: Study, SystemGroup, Replica,
│                             ObservableSeries, AggregateMetrics, ComparativeSummary
├── observable_resolver.py ← clasificación robusta de XVG por nombre, alias,
│                             regex y matching semántico
└── study_analyzer.py     ← motor principal: auto-discovery, grouping,
                              aggregate stats, outlier detection, findings
```

### Convención de nombres soportada

```
AA-A1rmsd_protein.xvg     → system=AA, replica=A1, observable=protein_rmsd
LP-A4_rmsd-ligand.xvg     → system=LP, replica=A4, observable=ligand_rmsd
HMG-A5mindist_lig_active.xvg → system=HMG, replica=A5, observable=mindist
```

El `ObservableResolver` tolera variantes heterogéneas:
`rmsd_protein`, `protein_rmsd`, `rmsd-protein`, `rmsdProtein`, etc.

### Nuevo comando CLI

```bash
simforge study [PATH]    # PATH opcional, default="."
```

Salida:
- Sistemas detectados (N réplicas, observables)
- Tabla de estadísticas agregadas (mean ± inter-replica std)
- Outliers detectados (umbral de Grubbs adaptado a n pequeño)
- Análisis comparativo textual
- Ranking de estabilidad por sistema

### CLI ergonomics mejorado

- `simforge analyze`  → PATH="." por defecto
- `simforge status`   → PATH="." por defecto
- `simforge summary`  → PATH="." por defecto
- `simforge study`    → PATH="." por defecto
- Shell completion habilitado (`--install-completion`)

### Scientific UX

- Eliminado lenguaje agresivo: "system may have exploded" → 
  "possible large conformational change, aggregation, or incomplete equilibration"
- Evidencia probabilística y contextual

---

# Arquitectura completa

```
simforge/
├── cli.py                        ← CLI con study + defaults mejorados
│
├── core/
│   ├── compiler.py               ← compilador IR → plan
│   ├── parser.py                 ← parse YAML → SystemState
│   ├── decision_engine.py        ← razonamiento semántico
│   ├── ontology.py               ← tipos biológicos
│   ├── models.py                 ← modelos Pydantic
│   ├── execution_models.py       ← IR pasos
│   ├── variant_compiler.py       ← compilación multi-variante
│   ├── semantic_inference.py     ← inferencia automática
│   ├── scientific_planner.py     ← diálogo de planning
│   ├── geometry_advisor.py       ← advisory geométrico
│   ├── project_manager.py        ← gestión de runs
│   └── md_knowledge/             ← base de conocimiento MD
│
├── builders/
│   ├── workspace_builder.py
│   └── step_builders/
│       ├── preparation_builder.py
│       ├── minimization_builder.py
│       ├── equilibration_builder.py
│       ├── production_builder.py
│       └── analysis_builder.py
│
├── runtime/
│   ├── xvg_parser.py             ← parse GROMACS .xvg
│   ├── convergence_analyzer.py   ← análisis convergencia/estabilidad
│   ├── quality_classifier.py     ← clasificación calidad simulación
│   ├── trajectory_ingestor.py    ← discovery automático de trayectorias
│   ├── scientific_summary.py     ← resumen científico por workspace
│   ├── study_models.py           ← modelos Study Layer [NUEVO]
│   ├── observable_resolver.py    ← resolver heterogeneous naming [NUEVO]
│   ├── study_analyzer.py         ← motor análisis comparativo [NUEVO]
│   └── executor.py               ← ejecución de workspaces
│
├── executors/
│   ├── shell_executor.py
│   ├── remediation_executor.py
│   └── signal_detector.py
│
└── adapters/
    ├── inflategro_adapter.py
    ├── movememb_adapter.py
    └── water_deletor_adapter.py
```

---

# Comandos disponibles

```bash
simforge compile   <yaml>          # compilar config → workspace
simforge run       [workspace]     # ejecutar workspace compilado
simforge dry-run   [workspace]     # simular sin ejecutar
simforge validate  <yaml>          # validar config
simforge inspect   <yaml>          # inspeccionar IR y DAG
simforge status    [workspace]     # estado de ejecución (default: .)
simforge summary   [workspace]     # resumen científico (default: .)
simforge analyze   [path]          # análisis calidad simulación (default: .)
simforge study     [path]          # análisis estudio multi-sistema (default: .)  ← NUEVO
simforge init                      # wizard interactivo de configuración
simforge clean     <workspace>     # limpiar workspace
simforge recompile <yaml>          # recompilar
```

---

# Capacidades actuales

## Workflow compiler
- Parse YAML → SystemState → SimulationPlan → ExecutionDAG
- Semantic inference: tipo de sistema, forcefield, objectives
- Scientific planning dialogue
- Geometry advisory
- Multi-variante (variants YAML)
- Workspace materialization con scripts GROMACS

## Analysis layer (single simulation)
- Parse XVG: RMSD, RMSF, energy, pressure, temperature
- Convergence analysis (last-20% plateau detection)
- Quality classification (5 tiers: converged → insufficient_data)
- Context-aware classification (globular_protein, membrane_protein, etc.)
- Scientific summary generation

## Study layer (multi-system comparative)  ← NUEVO
- Auto-discovery de XVG con naming heterogéneo
- Inferencia sistema/réplica/observable desde nombre de archivo
- ObservableResolver: aliases, regex, semántica
- Estadísticas agregadas por sistema × observable
- Detección de outliers (umbral de Grubbs n-dependiente)
- Análisis comparativo generado automáticamente
- Ranking de estabilidad por sistema

---

# Tests

- 101 tests pasando (pytest)
- Cobertura: core compiler, parser, decision engine, quality classifier,
  XVG parser, convergence analyzer, scientific summary, trajectory ingestor

---

# Pendiente / Próximos pasos

## Prioridad alta
- Tests para study_models, observable_resolver, study_analyzer
- Refinar findings comparativos (más contexto biológico)
- Plateu stability score más preciso para series largas

## Prioridad media
- Exportación de resumen de estudio a Markdown/PDF
- Comparación temporal entre runs del mismo sistema
- Integración con simforge compile → simforge study pipeline

## No priorizado todavía
- Plotting / dashboards
- Web UI
- ML models
- Análisis FEP
