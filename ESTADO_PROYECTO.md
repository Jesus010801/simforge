# Estado del Proyecto — SimForge
**Actualizado:** 2026-05-25 (sesión 5)
**Líneas de código:** ~29,000 Python en 143 archivos  
**Tests:** 708 pasando, 0 fallos

---

## ¿Qué es SimForge hoy?

Dos productos en un repositorio:

```
┌─────────────────────────────────────────────────────────────────┐
│  PRODUCTO A — Workflow Compiler                                  │
│  YAML → IR → DAG → scripts GROMACS → ejecución                  │
│  Estado: funcional en papel, nunca ejecutado contra GROMACS real │
├─────────────────────────────────────────────────────────────────┤
│  PRODUCTO B — Comparative MD Study Analyzer      ← VALOR REAL   │
│  directorio con XVGs → clasificación científica + ranking        │
│  Estado: funcional y útil HOY con datos reales                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Arquitectura completa (estado real)

```
simforge/
│
├── cli.py                         2127 líneas  ← MONOLITO (problema)
│
├── core/                          pipeline de compilación
│   ├── compiler.py
│   ├── parser.py                  640 líneas
│   ├── decision_engine.py         970 líneas
│   ├── models.py                  600 líneas (Pydantic + MembraneOrientation)
│   ├── execution_models.py        AutomationLevel + SimulationStep.automation_level
│   ├── membrane_geometry.py       rotación TM desde EC/IC residuos
│   ├── structural_annotation.py   NEW — StructuralAnnotation, MembraneTopologyAnnotation,
│   │                                    OrientationAnnotation; overlap_warnings()
│   ├── semantic_inference.py      normalización objetivos
│   ├── semantic_objectives.py     aliases y presets
│   ├── scientific_planner.py      482 líneas
│   ├── variant_compiler.py        360 líneas
│   ├── geometry_advisor.py        advisory geométrico
│   ├── workflow_hints.py          bridge semántica → policy
│   ├── workspace_fingerprint.py   SHA-256 invalidación
│   ├── project_manager.py         gestión timestamped runs
│   └── md_knowledge/              base de conocimiento MD
│       ├── states.py              9 estados simulación
│       ├── patterns.py            6 patrones temporales
│       ├── contexts.py            9 contextos de sistema
│       ├── heuristics.py          10 observables × 9 contextos
│       ├── evidence.py            acumulación de evidencia
│       └── interpreter.py        interpretación single-sim
│
├── builders/                      materializadores de workspace
│   ├── workspace_builder.py
│   └── step_builders/
│       ├── preparation_builder.py  dispatch → membrane_orient_builder
│       ├── membrane_orient_builder.py  NEW — orient_protein AUTOMATED/GUIDED
│       ├── assembly_builder.py
│       ├── minimization_builder.py
│       ├── equilibration_builder.py
│       ├── production_builder.py
│       ├── analysis_builder.py
│       ├── embedding_builder.py   shrink loop membrana
│       ├── parametrization_builder.py
│       ├── validation_builder.py
│       ├── enhanced_sampling_builder.py
│       └── _utils.py              mdrun_block / mdrun_resume_block
│
├── runtime/                       capa de análisis (el valor real)
│   ├── xvg_parser.py              parse GROMACS .xvg
│   ├── convergence_analyzer.py    plateau detection, drift
│   ├── trajectory_ingestor.py     auto-discovery de trayectorias
│   ├── quality_classifier.py      533 líneas — clasifica calidad
│   ├── scientific_summary.py      resumen por workspace
│   │
│   ├── study_models.py            Study, SystemGroup, Replica, AggregateMetrics
│   ├── observable_resolver.py     aliases heterogéneos → canonical
│   ├── study_analyzer.py          parse_study() — motor principal
│   │
│   ├── synthesis_models.py        modelos Scientific Synthesis
│   ├── interaction_interpreter.py 446 líneas — 9 reglas ponderadas
│   ├── consensus_engine.py        replica consensus via median
│   ├── event_detector.py          eventos temporales (3 detectores)
│   ├── scientific_synthesis.py    orchestrator + narrative gen
│   │
│   ├── executor.py                642 líneas — RuntimeExecutor (async)
│   ├── execution_backend.py       ABC ExecutionBackend
│   ├── stream.py                  AsyncProcessRunner
│   ├── events.py                  EventBus (21 tipos)
│   ├── journal.py                 append-only JSONL
│   ├── metrics.py                 psutil + nvidia-smi
│   ├── artifacts.py               SHA-256 registry
│   └── cache.py                   fingerprint cache
│
├── executors/                     ← FRAGMENTACIÓN (problema)
│   ├── shell_executor.py          legacy — ya no es el default
│   ├── gromacs_executor.py        747 líneas — percepción GROMACS
│   ├── base_executor.py           457 líneas
│   ├── adaptive_reasoning.py      697 líneas — parcialmente cableado
│   ├── remediation_executor.py    783 líneas — loop de remediación
│   ├── signal_detector.py         876 líneas — el más grande
│   └── execution_state.py         estado por step
│
├── adapters/                      herramientas externas
│   ├── inflategro_adapter.py      334 líneas
│   ├── movememb_adapter.py
│   └── water_deletor_adapter.py
│
├── pipelines/                     definiciones de pipeline
│   └── membrane_pipeline.py       385 líneas — 12 steps DPPC/OPLS-AA
│
├── validators/
│   ├── protein_validator.py
│   ├── ligand_validator.py        553 líneas
│   └── membrane_validators.py     568 líneas
│
└── descriptors/                   percepción fisicoquímica
    ├── flexibility.py
    ├── polarity.py                 421 líneas
    ├── geometry.py
    ├── aromaticity.py
    └── topology.py
```

---

## Comandos disponibles

```bash
# Workflow compiler
simforge init                      # wizard interactivo → YAML
simforge compile   <yaml>          # YAML → workspace con scripts GROMACS
simforge validate  <yaml>          # validar sin generar artefactos
simforge inspect   <yaml>          # IR + DAG + descriptores
simforge run       <workspace>     # ejecutar con RuntimeExecutor
simforge dry-run   <workspace>     # simular sin ejecutar
simforge clean     <workspace>     # limpiar steps/ + estado
simforge recompile <yaml>          # clean + recompile

# Analysis layer (valor real hoy)
simforge analyze   [path]                     # clasifica calidad de una simulación
simforge study     [path]                     # análisis comparativo multi-sistema ← POTENTE
simforge study     [path] --report out.md     # ídem + exporta reporte Markdown
simforge study     [path] --output out.json   # ídem + exporta JSON
simforge summary   [workspace]                # resumen científico de workspace
simforge status    [workspace]                # estado de ejecución

# Structural annotation (protein_membrane)
simforge annotate-structure <config.yaml>    # wizard EC/IC/TM → escribe structural_annotation al YAML
```

---

## Capacidades actuales en detalle

### Workflow Compiler
- YAML → SystemState → SimulationPlan → DAG → manifest
- Inferencia semántica: tipo sistema, forcefield, objectives
- Scientific planning dialogue (preguntas interactivas)
- Geometry advisory (box size, PBC warnings)
- Multi-variante (variants YAML)
- Workspace timestamped con SHA-256 fingerprinting
- Checkpoint recovery (resume desde .cpt)
- Membrane pipeline: 12 steps DPPC + OPLS-AA

### Analysis Layer (single simulation — `simforge analyze`)
- Parse XVG: RMSD, RMSF, energía, presión, temperatura
- Convergence: last-20% plateau detection, drift linear fit
- Quality classification: 5 tiers (converged → insufficient_data)
- Context-aware: 9 contextos de sistema, 10 observables, soft ranges
- QualityReport: quality, confidence [0-1], evidence[], warnings[], recommendations[]

### Study Layer (multi-sistema — `simforge study`)
- Auto-discovery XVG con naming heterogéneo
- Inferencia sistema/réplica/observable desde nombre de archivo
- ObservableResolver: aliases, regex, semántica (tolera 20+ variantes)
- Estadísticas agregadas (mean ± inter-replica std) por sistema × observable
- Grubbs outlier detection (n-dependiente: n=3→1.15σ, n=4→1.48σ…)
- Análisis comparativo textual generado automáticamente

### Scientific Synthesis Layer (`simforge study` — síntesis profunda)
- Normalización sigmoid centrada en media del estudio
- 9 reglas ponderadas con `min_data_fraction` guard
- Binding states: stable_binding, weak_binding, transient_binding,
  ligand_destabilization, possible_dissociation, interaction_persistent
- Structural states: structurally_stable, flexible_but_stable, conformational_rearrangement
- Composite score: binding×0.5 + stability×0.3 + convergence×0.2
- Replica consensus via mediana del estudio
- Eventos temporales: late_destabilization, abrupt_transition, contact_loss, ligand_drift
- Narrativa científica auto-generada (5 párrafos)

---

## Tests

| Módulo | Tests | Estado |
|--------|-------|--------|
| core/compiler | ✓ | cubierto |
| core/parser | ✓ | cubierto |
| core/decision_engine | ✓ | cubierto |
| core/semantic_inference | 38 | cubierto |
| core/md_knowledge | 59 | cubierto |
| core/geometry_advisor | 21 | cubierto |
| runtime/quality_classifier | 57 | cubierto |
| runtime/trajectory_ingestor | ✓ | cubierto |
| runtime/executor | 29 | cubierto |
| runtime/checkpoint_recovery | 25 | cubierto |
| **runtime/study_analyzer** | **0** | **SIN TESTS** |
| **runtime/observable_resolver** | **0** | **SIN TESTS** |
| **runtime/interaction_interpreter** | **0** | **SIN TESTS** |
| **runtime/consensus_engine** | **0** | **SIN TESTS** |
| **runtime/event_detector** | **0** | **SIN TESTS** |
| **runtime/scientific_synthesis** | **0** | **SIN TESTS** |
| core/variant_compiler | 0 | sin tests |
| core/workspace_fingerprint | 0 | sin tests |
| core/scientific_planner | 0 | sin tests |

**Total: 708 tests, 0 fallos** (sesión 4: +1 test → test_embed_in_bilayer_is_automated)

---

# Análisis de virtudes y problemas

## Virtudes reales

**1. Arquitectura declarativa limpia (core)**
El pipeline YAML → SystemState → SimulationPlan → DAG está bien diseñado.
Cada capa tiene una responsabilidad clara. El `decision_engine` es
sofisticado sin ser opaco.

**2. Study Layer + Synthesis — potencia real**
`simforge study` hace algo que no existe en herramientas standard:
clasifica automáticamente sistemas desde XVGs crudos con razonamiento
explícito. Un investigador con 24 XVGs de 3 sistemas × 4 réplicas × 2
observables obtiene en segundos: estado de interacción, ranking,
narrativa científica. Eso tiene valor real hoy.

**3. Sin dependencias científicas pesadas**
Pure Python. Sin numpy, sin MDAnalysis, sin scipy. Corre en cualquier
máquina con Python 3.10+. Eso es una ventaja de distribución enorme.

**4. Test suite sólida para el compilador**
548 tests para los módulos del compilador es impressive. La arquitectura
está protegida de regresiones en esa capa.

**5. UX rico con Rich**
Los paneles con colores, barras de progreso, tablas comparativas y
narrativa hacen que la salida sea legible para investigadores no-técnicos.
Es fácil de usar desde el terminal.

---

## Problemas reales (sin filtro)

### P1. DEUDA TÉCNICA CRÍTICA: 8 módulos sin tests
Todo el Study Layer y Scientific Synthesis Layer (los módulos más nuevos
y más complejos) tienen cobertura CERO. Eso incluye la lógica de
normalización sigmoid, las 9 reglas de interacción, el detector de
eventos, y el motor de consenso. Si cualquier cosa se rompe aquí, no hay
red de seguridad.

**Impacto:** Alto — cualquier refactor futuro puede romper silenciosamente
el output científico.

### P2. DOS CAMINOS DE INTERPRETACIÓN PARALELOS Y DESCONECTADOS
```
core/md_knowledge/interpreter.py   ← para simulación individual
runtime/interaction_interpreter.py ← para estudio multi-sistema
```
Ambos hacen: normalizar → clasificar con evidencia → generar estado.
El primero usa soft ranges contextuales. El segundo usa sigmoid + reglas.
No comparten código ni modelos. Tienen estados distintos
(`SimulationState` enum vs strings literales).

**Impacto:** Mantenimiento duplicado. Cuando se mejora uno, el otro no
se beneficia. El usuario obtiene respuestas distintas de `analyze` y
`study` para el mismo sistema.

### P3. FRAGMENTACIÓN DE EXECUTORS
Hay tres caminos de ejecución:
```
executors/shell_executor.py     ← legacy, ya no es default
executors/gromacs_executor.py   ← percepción GROMACS, 747 líneas
runtime/executor.py             ← RuntimeExecutor, 642 líneas (el actual default)
```
El CLI selecciona entre `gromacs` y `shell/runtime` con un flag string.
`adaptive_reasoning.py` (697 líneas) existe pero la conexión real entre
"step falla → llama adaptive reasoning → decide remediación" no está
completamente cerrada en un loop de producción.

**Impacto:** Confusión sobre qué executor usar. Código de remediación
sofisticado que puede nunca ejecutarse en la práctica.

### P4. cli.py CON 2127 LÍNEAS ES UN MONOLITO
El archivo CLI contiene: helpers de display, lógica de planning, 
lógica de geometría, compilación, variants, ejecución, análisis, study.
Mezcla UI con orquestación con formateo.

**Impacto:** Difícil de leer, navegar y mantener. Imposible de testear
unitariamente.

### P5. EL COMPILER NUNCA HA VISTO GROMACS REAL
Todo el workflow compiler genera scripts `.sh` que invocan comandos
GROMACS, pero jamás se ha ejecutado una simulación completa end-to-end
con GROMACS real. Los tests son unitarios (sin subprocess real).

**Impacto:** No sabemos si los scripts generados funcionan. La membrana
pipeline de 12 pasos es arquitecturalmente correcta pero no validada.

### P6. ESTADO_PROYECTO.md SIEMPRE DESACTUALIZADO
El archivo de estado no refleja la Scientific Synthesis Layer que se acaba
de construir (5 módulos nuevos). Cada sesión empieza con contexto stale.

---

# Evaluación de utilidad como usuario

## ¿Qué funciona y da valor HOY?

```
simforge study /path/to/xvgs/
```

Esto es el producto real. Un investigador con datos de MD comparativa
(múltiples sistemas, réplicas, observables) puede:
1. Soltar sus XVGs en un directorio
2. Ejecutar `simforge study .`
3. Obtener en segundos:
   - Qué sistema tiene mejor binding
   - Estado de interacción (stable/weak/transient) con evidencia
   - Réplicas outlier con razón
   - Eventos temporales detectados
   - Narrativa científica lista para copiar a un paper

Eso no existe en GROMACS tools. Es útil.

```
simforge analyze /path/to/simulation/
```

También útil: clasifica calidad de cualquier simulación existente.

## ¿Qué NO está listo para uso real?

```
simforge compile → simforge run
```

El loop completo de compilación → ejecución no está validado con GROMACS
real. No se puede poner en manos de un investigador y decir "compila esto
y corre tu simulación". Puede generar scripts incorrectos silenciosamente.

---

# Trabajo inteligente: cómo mejorar la eficiencia

## Problema actual de las sesiones

Cada conversación:
1. Reconstruye contexto desde cero (costoso en tokens)
2. Mezcla dominos (compiler + análisis + synthesis en la misma sesión)
3. El ESTADO_PROYECTO.md es la única fuente de verdad pero siempre está
   desactualizado

## Solución propuesta: sesiones por dominio + ESTADO como contrato

**Regla 1: ESTADO_PROYECTO.md como primer mensaje implícito**
Al inicio de cada sesión, leer este archivo primero. Evita preguntas
de "¿dónde estábamos?"

**Regla 2: Una sesión = un dominio**
- Sesión A: solo Study Layer / Synthesis (runtime/)
- Sesión B: solo Compiler (core/ + builders/)
- Sesión C: solo CLI / UX

**Regla 3: Tests antes de features**
Antes de agregar features a Study/Synthesis, escribir los tests que
faltan. Esto permite sesiones más cortas y seguras.

**Regla 4: Scope declarations al inicio**
"Esta sesión: solo escribir test_study.py. No tocar más nada."
Eso evita la tentación de arreglar cosas adyacentes y gastar contexto.

---

# El salto: qué implementar para dar el próximo nivel

## Nivel 1 — Consolidación (1-2 sesiones, alta prioridad)

### 1a. test_study.py — cerrar la deuda técnica
Cubrir los 8 módulos sin tests:
```python
# test_study.py
- test_observable_resolver: aliases, prioridad ligand>protein, regex
- test_study_analyzer: filename parsing, aggregate stats, outlier detection
- test_interaction_interpreter: sigmoid norm, reglas, tiebreaker
- test_consensus_engine: median split, labels
- test_event_detector: late_destab, abrupt, ligand_drift
- test_scientific_synthesis: pipeline completo con datos sintéticos
```
Esto toma ~1 sesión y bloquea regresiones futuras.

### 1b. Unificar los dos intérpretes
`core/md_knowledge/interpreter.py` y `runtime/interaction_interpreter.py`
resuelven el mismo problema. La solución: hacer que `simforge analyze`
use el mismo backend que `simforge study` para sistemas single-replica.
Un único `InteractionInterpreter` que funciona con 1 o N réplicas.

---

## Nivel 2 — El salto real: report generado (1 sesión)

### `simforge study --report report.md`

Exportar el resultado de `study` a un Markdown/PDF listo para usar en
un paper o tesis. Incluye:
- Tabla de estadísticas agregadas
- Clasificación de estados con evidencia
- Ranking con barras
- Narrativa científica completa
- Sección de outliers y eventos

**Por qué es el salto:** Convierte SimForge de "tool para el terminal" a
"herramienta que genera un artefacto entregable". Un investigador puede
ejecutar `simforge study . --report resultados.md` y tener el draft de
la sección de Resultados de su paper en segundos.

---

## Nivel 3 — Validación del compiler (requiere GROMACS)

### End-to-end test con GROMACS real
El compiler necesita al menos un test de integración real:
```bash
simforge compile configs/lysozyme_simple.yaml
simforge run simforge_runs/lysozyme/
# Verificar: minimización converge, equilibración completa
```

Sin esto, el compiler sigue siendo software de demostración.

---

## Nivel 4 — Multi-observable study enriquecido

### Catalytic distance + H-bonds como observables first-class
Los estudios de inhibición usan `catalytic_distance` y `hydrogen_bonds`
como observables clave. Actualmente están en las reglas pero raramente
aparecen en los XVGs porque el ObservableResolver los detecta poco.

### Temporal profiles por réplica
En vez de solo mostrar mean±std, mostrar si hay divergencia temporal:
"réplica A3 converge, réplica A4 deriva". Esto requiere extender
`event_detector` para reportar a nivel réplica en la tabla principal.

---

## Cambios recientes (2026-05-25) — sesión 5

### Paridad `init` ↔ `annotate-structure` en validación de anotación EC/IC/TM

**Contexto:** `simforge init` integra el subwizard de anotación estructural
(EC/IC/TM) cuando el sistema es `protein_membrane`. El comando independiente
`simforge annotate-structure` hace lo mismo sobre un YAML ya existente.
Los dos caminos mostraban diferente información al usuario tras el wizard.

**Cambios en `cli.py`:**

- **Default del prompt EC/IC/TM → "No"** (`default=2`): antes el usuario que
  pulsaba Enter entraba al subwizard. Ahora Enter salta la anotación, consistente
  con la filosofía "No path adds zero overhead" y compatibilidad de flujo.

- **`overlap_warnings()` en `init`**: tras mostrar `validation_warnings()`, el
  wizard de init ahora también llama `_topology.overlap_warnings()` y muestra
  los solapamientos en rojo (`[red]✗[/red]`), igual que hace `annotate-structure`
  en las líneas 2865–2869. Ambos caminos son ahora metodológicamente equivalentes.

**No hay cambios de tests:** la lógica nueva está en el wizard interactivo (rama
`if want_ann == "yes"`) que no está cubierta por los tests automatizados actuales.

---

## Cambios recientes (2026-05-24) — sesión 4

### `embed_in_bilayer` promovido de MANUAL → AUTOMATED

**Problema:** `embed_in_bilayer` caía en `_build_generic()` (README vacío). Era el único
step de la ruta crítica sin builder real. `system.gro` no se generaba automáticamente,
bloqueando `generate_topology` → `membrane_embedding` → todo lo demás.

**Cambios:**

- **`builders/step_builders/assembly_builder.py`**:
  - Nuevo dispatch `elif sid == "embed_in_bilayer":` en `build()`
  - `_build_embed_in_bilayer()`: genera `run_embed.py` (self-contained) + `run.sh` wrapper.
    El script usa `MoveMembAdapter` (Python puro, sin gfortran) para centrar el
    midplano del bilayer en el Z-centro de la proteína → `system.gro`.
    Luego `gmx genrestr -fc 100000 100000 100000` → `strong_posre.itp`.
    Luego `gmx editconf -resnr 1` para renumerar residuos desde 1.
  - Resolución del bilayer GRO: primero en CWD, luego en `Prot-Memb_FILES/` (fallback).
  - `protein_boxed.gro` leído desde `match_box_to_bilayer/` via `_rel()`.

- **`pipelines/membrane_pipeline.py`**:
  - `step_type`: `MANUAL` → `AUTOMATIC`
  - `automation_level`: `AUTOMATED` añadido
  - `params`: eliminado el bloque `note` con instrucciones manuales (ya no necesario)
  - Docstring DAG actualizado

- **`benchmarks/membrane_dppc_oplsaa/test_pipeline.py`**:
  - `test_manual_steps`: eliminada aserción `"embed_in_bilayer" in step_ids`; añadida
    aserción `"embed_in_bilayer" not in step_ids` con mensaje claro
  - `test_embed_in_bilayer_is_automated`: nuevo test — verifica `step_type == "automatic"`
    y `automation_level == AUTOMATED`
  - Importado `AutomationLevel` desde `core.execution_models`

**Resultado:** 27 → 28 tests en test_pipeline.py | 707 → 708 tests totales | 0 fallos

---

## Bugfixes recientes (2026-05-23)

### Sesión 1 — membrane runtime blocking
- **`PreparationBuilder`**: steps MANUAL con engine != pdb2gmx generan solo README.md
- **`RuntimeExecutor` + `ShellExecutor`**: step_type verificado antes que `_find_script()` — MANUAL siempre salta
- **`_BLOCKING_STATUSES`**: `SKIPPED` añadido — pasos manuales pendientes bloquean downstream

### Sesión 2 — semantic inference de membrana
- **`simforge init`**: serializa `environment.membrane.enabled: true` + pregunta tipo lípido
- **`core/inference.py`**: fallback por objetivos — `membrane_perturbation` infiere "protein-membrane" con warning
- **Trazabilidad**: panel "Semantic Normalization" muestra ✓ pipeline + lípido + hints

### Sesión 3 — orient_protein automation + AutomationLevel

**Diagnóstico del bug (`orient_protein` siempre SKIPPED):**
El runtime leía `step_type` de `metadata.json` (string "manual"/"automatic"). El modelo binario
no podía expresar "automático cuando orientation está definido, guiado sin ella". Cualquier
workspace compilado sin orientación (o compilado con código viejo) quedaba SKIPPED aunque
tuviera `run.sh`.

**Cambios implementados:**

- **`core/execution_models.py`**: `AutomationLevel` enum formal (MANUAL/GUIDED/SEMI_AUTOMATED/AUTOMATED)
  con `from_step_type()`, `needs_user` property, y `effective_automation_level()` en `SimulationStep`
- **`core/models.py`**: `MembraneOrientation` model + campo `orientation` en `MembraneConfig`
  — YAML acepta `extracellular_residues`, `intracellular_residues`, `tm_segments`
- **`core/membrane_geometry.py`** (NEW): biblioteca pura — `parse_residue_range()`,
  `compute_orient_rotation(gro, ec_resids, ic_resids) → OrientRotation(rx,ry,rz)`
  Tabla de rotación: TM±X→270/90°Y, TM±Y→∓90°X, TM+Z→0°, TM-Z→180°X
- **`builders/step_builders/membrane_orient_builder.py`** (NEW): genera `run.sh` + `orient_helper.py`
  (self-contained, sin imports SimForge) cuando orientation está presente; README.md sin script si no
- **`pipelines/membrane_pipeline.py`**: `orient_protein` recibe `automation_level=AUTOMATED` o
  `GUIDED` según `membrane.orientation` — ya no binario
- **`builders/workspace_builder.py`**: `execution_manifest.json` incluye campo `automation_level`
- **`runtime/executor.py` + `executors/shell_executor.py`**: skip logic lee `automation_level`
  primero, fallback a `step_type` (backward compat con workspaces viejos)
- **`builders/step_builders/preparation_builder.py`**: `_build_manual_readme` escribe `automation_level: guided`
- **`configs/membrane_orient_test.yaml`** (NEW): config de test con `membrane.orientation`
- **`builders/test_orient_protein_flow.py`** (NEW): 15 tests integración (pipeline→workspace→metadata→dry-run)
- **`core/test_membrane_geometry.py`** (NEW): 17 tests unitarios de geometría

**Invariante clave**: workspaces compilados sin `automation_level` en metadata.json siguen
funcionando — el runtime cae al campo `step_type` legado.

---

# Próximos pasos concretos

## Prioridad 1 (bloquea todo lo demás)
- [ ] `runtime/test_study.py` — cubrir 8 módulos sin tests
- [ ] Actualizar ESTADO_PROYECTO.md al final de cada sesión (regla)

## Prioridad 2 — membrana (secuencial, uno por uno)
Estado: `orient_protein` ✓ | `match_box_to_bilayer` ✓ | `embed_in_bilayer` ✓

- [x] `orient_protein` — AUTOMATED con `membrane.orientation` en YAML
- [x] `match_box_to_bilayer` — MatchBoxBuilder: protein_oriented.gro → box_match_report.json + protein_boxed.gro
- [x] `embed_in_bilayer` — AssemblyBuilder: MoveMembAdapter + gmx genrestr → system.gro + strong_posre.itp
- [ ] `clean_water` — script Python → `WaterDeletorAdapter` + actualizar SOL en topol.top

Cada step: impl + tests antes de continuar al siguiente.

**Próximo desbloqueado:** `clean_water` es el único step de assembly que sigue como EXTERNAL.
Completarlo cierra el loop assembly → ions → minimización, dejando el pipeline compilable
y ejecutable end-to-end.

## Prioridad 3 (alto valor, baja complejidad)
- [ ] `simforge study --report <file>` — exportar Markdown
- [ ] Unificar `analyze` y `study` en un solo backend de interpretación

## Prioridad 4 (validación real)
- [ ] Test end-to-end con GROMACS real (lysozyme, 1ns)
- [ ] Calibrar thresholds de Grubbs con datos reales

## No priorizado
- Plotting / dashboards
- Web UI
- ML models
- FEP
- SLURM / cloud backend
