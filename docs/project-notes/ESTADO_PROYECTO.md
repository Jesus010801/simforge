# Estado del Proyecto — SimForge
**Actualizado:** 2026-06-11 (sesión 8)
**Líneas de código:** ~33,500 Python en ~162 archivos
**Tests:** 1036 pasando, 0 fallos, 6 skipped (benchmarks excluidos; 3 skipped requieren rdkit_env)

---

## ¿Qué es SimForge hoy?

Tres productos en un repositorio:

```
┌─────────────────────────────────────────────────────────────────┐
│  PRODUCTO A — Workflow Compiler                                  │
│  YAML → IR → DAG → scripts GROMACS → ejecución                  │
│  Estado: DAG correcto, pipeline membrana completa, nunca         │
│          ejecutado contra GROMACS real end-to-end                │
├─────────────────────────────────────────────────────────────────┤
│  PRODUCTO B — Comparative MD Study Analyzer      ← VALOR REAL   │
│  directorio con XVGs → clasificación científica + ranking        │
│  Estado: funcional y útil HOY con datos reales                   │
├─────────────────────────────────────────────────────────────────┤
│  PRODUCTO C — Ligand Parameterization Toolkit                    │
│  SDF/PDB → LigParGen-ready export → normalización → assembly     │
│  Estado: export validado experimentalmente con LigParGen real    │
│          Charge reporting completo. CLI maduro.                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Arquitectura completa (estado real)

```
simforge/
│
├── cli.py                         ~3100 líneas  ← MONOLITO (deuda conocida)
│   ├── simforge compile/run/validate/inspect/...
│   ├── simforge study/analyze/summary/status
│   ├── simforge annotate-structure
│   └── simforge ligand export-ligpargen  (sesión 7-8: --legacy --smiles --charge)
│
├── core/
│   ├── compiler.py
│   ├── parser.py
│   ├── decision_engine.py
│   ├── models.py
│   ├── execution_models.py        AutomationLevel + SimulationStep.automation_level
│   ├── membrane_geometry.py
│   ├── structural_annotation.py   StructuralAnnotation, MembraneTopologyAnnotation
│   ├── semantic_inference.py      normalización objetivos
│   ├── semantic_objectives.py     aliases y presets
│   ├── scientific_planner.py
│   ├── variant_compiler.py
│   ├── geometry_advisor.py
│   ├── workflow_hints.py          bridge semántica → policy
│   ├── workspace_fingerprint.py   SHA-256 invalidación
│   ├── project_manager.py         gestión timestamped runs
│   ├── ligand_workflow_models.py  Pydantic models fases 1-4 + formal_charge ← sesión 8
│   └── md_knowledge/
│       ├── states.py / patterns.py / contexts.py
│       ├── heuristics.py / evidence.py / interpreter.py
│
├── ligand/                        toolkit parameterización ligandos
│   ├── rdkit_reader.py            carga mol (lazy RDKit)
│   ├── export.py                  legacy/smiles + charge helpers  ← sesión 8
│   ├── preparation.py             validate_ligand_for_parameterization
│   ├── legacy_writer.py           LigParGenLegacyWriter (ATOM-record PDB)
│   ├── normalization.py           normalize_ligpargen_outputs
│   ├── ligpargen_import_validator.py
│   ├── pose_rewriter.py           PoseRewriter (Kabsch + coord transfer)
│   └── test_*.py                  53 tests (normalization, validator, pose)
│
├── builders/
│   ├── workspace_builder.py
│   └── step_builders/
│       ├── preparation_builder.py
│       ├── membrane_orient_builder.py  orient_protein AUTOMATED/GUIDED
│       ├── assembly_builder.py         clean_water AUTOMATED
│       ├── match_box_builder.py        match_box_to_bilayer AUTOMATED
│       ├── minimization/equilibration/production/analysis builders
│       ├── embedding_builder.py   shrink loop membrana
│       ├── parametrization_builder.py
│       └── _utils.py              mdrun_block / mdrun_resume_block
│
├── runtime/                       capa de análisis (el valor real)
│   ├── xvg_parser.py / convergence_analyzer.py / trajectory_ingestor.py
│   ├── quality_classifier.py / scientific_summary.py
│   ├── study_models.py / observable_resolver.py / study_analyzer.py
│   ├── synthesis_models.py / interaction_interpreter.py
│   ├── consensus_engine.py / event_detector.py / scientific_synthesis.py
│   ├── executor.py                RuntimeExecutor (async)
│   ├── execution_backend.py
│   └── gates/
│       ├── gate_runner.py / water_gate.py / apl_gate.py
│       ├── overlap_gate.py / topology_gate.py
│       ├── box_match_gate.py / orientation_gate.py
│
├── executors/
│   ├── shell_executor.py / gromacs_executor.py / base_executor.py
│   ├── adaptive_reasoning.py / remediation_executor.py
│   ├── signal_detector.py / execution_state.py
│
├── adapters/
│   ├── inflategro_adapter.py / movememb_adapter.py
│   └── water_deletor_adapter.py   WaterDeletorAdapter (Python puro, sin Perl)
│
├── pipelines/
│   └── membrane_pipeline.py       12 steps DPPC/OPLS-AA
│
└── validators/
    ├── protein_validator.py / ligand_validator.py / membrane_validators.py
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
simforge study     [path]                     # análisis comparativo multi-sistema
simforge study     [path] --report out.md     # ídem + exporta reporte Markdown
simforge study     [path] --output out.json   # ídem + exporta JSON
simforge summary   [workspace]                # resumen científico de workspace
simforge status    [workspace]                # estado de ejecución

# Structural annotation (protein_membrane)
simforge annotate-structure <config.yaml>     # wizard EC/IC/TM → escribe al YAML

# Ligand toolkit (sesiones 7-8 — maduro)
simforge ligand export-ligpargen <lig.sdf> --legacy
    # PDB validado experimentalmente con servidor LigParGen real
    # Escribe 4 archivos: LIG_ligpargen_legacy.pdb + LIG_ligpargen.smi
    #                     + LIG_meta.json + LIG_charge.txt

simforge ligand export-ligpargen <lig.sdf> --smiles
    # SMILES canónico + metadata JSON + charge advisory
    # Escribe 3 archivos: LIG.smi + LIG_meta.json + LIG_charge.txt

simforge ligand export-ligpargen <lig.sdf> --legacy --mol-name A001 --output-dir ./out/
    # Opciones completas
```

---

## Estado del pipeline membrana (DAG verificado)

| Step | automation_level | dry-run | downstream input |
|---|---|---|---|
| orient_protein | automated¹ | DONE | → match_box_to_bilayer |
| match_box_to_bilayer | automated | DONE | → embed_in_bilayer/protein_boxed.gro |
| embed_in_bilayer | automated | DONE | → generate_topology/system.gro |
| generate_topology | automatic² | DONE | → solvate_membrane/topol.top |
| membrane_embedding | automatic² | DONE | → solvate_membrane/converged.gro |
| solvate_membrane | automatic² | DONE | → clean_water/solvated.gro |
| clean_water | automated | DONE | → add_ions/system_clean.gro |
| add_ions | automatic² | DONE | → energy_minimization/aaions.gro |
| energy_minimization | automatic² | DONE | → equilibration |
| equilibration | automatic² | DONE | → production_md |
| production_md | automatic² | DONE | → analysis |

¹ Requiere `structural_annotation` completa en YAML; sin ella → GUIDED (SKIPPED).
² `automatic` = legacy `step_type`; executor lo trata idéntico a `automated`.

**Invariantes topología:**
- `clean_water` escribe `topol.top` (SOL decrementado) y `system_clean.gro`
- `add_ions` lee ambos de `../07_clean_water/`
- `energy_minimization` y posteriores leen `topol.top` de `../08_add_ions/`

---

## Toolkit de ligandos — arquitectura completa

```
Fase 1 — Validación estructural
  ligand_validator.py:validate_ligand()           ← hookeado en compiler

Fase 2 — Exportación para LigParGen  ← CLI maduro (sesiones 7-8)
  preparation.py:validate_ligand_for_parameterization()   (pre-flight check)
  export.py:export_for_ligpargen_legacy()          (PDB ATOM, validado con servidor real)
  export.py:export_for_ligpargen_smiles()          (SMILES canónico)
  export.py:export_for_ligpargen()                 (PDB moderno HETATM, referencia)
  legacy_writer.py:LigParGenLegacyWriter
  export._compute_formal_charge() / _charge_label() / _write_charge_txt()

  ► CLI: simforge ligand export-ligpargen <sdf> --legacy|--smiles [opciones]

Fase 3 — Importación de outputs LigParGen
  normalization.py:normalize_ligpargen_outputs()
  ligpargen_import_validator.py:LigParGenImportValidator

Fase 4 — Reescritura de pose
  pose_rewriter.py:PoseRewriter.rewrite()          (transfiere coords docked → GRO)

Pendiente:
  - Fases 3-4 sin exposición CLI
  - Assembly integration (ligand + protein → system.gro) no implementado
  - parametrize_ligand step en compiler todavía GUIDED (manual)
```

### Outputs generados por `--legacy` (validados con LigParGen real)

```
ligpargen_export/
├── LIG_ligpargen_legacy.pdb   ← subir al servidor web LigParGen
├── LIG_ligpargen.smi          ← SMILES canónico companion
├── LIG_meta.json              ← nombre, SMILES, heavy_atom_count, formal_charge
└── LIG_charge.txt             ← advisory de carga (SIEMPRE revisar antes de submit)
```

`LIG_charge.txt` contiene:
```
Molecule: LIG
Formal charge: +1
Recommended LigParGen charge selection: +1
```

**Lección aprendida (sesión 8):** el formato PDB legacy generado por SimForge SÍ es
aceptado por el servidor LigParGen. El fallo anterior fue un **charge mismatch**:
la molécula tenía carga formal +1 y se sometió con charge=0. SimForge ahora reporta
la carga formal explícitamente y emite un WARNING amarillo cuando es ≠ 0.

**Dependencia RDKit:**
- Fases 2-4 requieren RDKit — activar `rdkit_env` antes de usar `simforge ligand export-ligpargen`
- El CLI importa RDKit de forma lazy (solo al ejecutar el comando, no en `--help`)
- El resto de SimForge funciona sin RDKit (pure Python)

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
- Membrane pipeline: 12 steps DPPC + OPLS-AA (DAG verificado en dry-run)
- Gates físicos: APL, overlap, topology consistency, water-in-bilayer, box match, orientation

### Membrane Pipeline — clean_water (completado sesión 6)
- `run_clean_water.py` usa `WaterDeletorAdapter` (Python puro, sin Perl)
- Cuenta SOL antes/después; escribe `clean_water_report.json` con:
  `input_water_count`, `removed_water_count`, `final_water_count`,
  `cutoff_used` (z_bot/z_top nm), `output_gro_path`, `topology_updated`
- `water_gate.py` lee `clean_water_report.json` primero, fallback a `water_report.json`

### Analysis Layer (single simulation — `simforge analyze`)
- Parse XVG: RMSD, RMSF, energía, presión, temperatura
- Convergence: last-20% plateau detection, drift linear fit
- Quality classification: 5 tiers (converged → insufficient_data)
- Context-aware: 9 contextos, 10 observables, soft ranges

### Study Layer (multi-sistema — `simforge study`)
- Auto-discovery XVG con naming heterogéneo
- ObservableResolver: aliases, regex, semántica (20+ variantes)
- Estadísticas agregadas (mean ± inter-replica std) por sistema × observable
- Grubbs outlier detection (n-dependiente)
- Scientific Synthesis: 9 reglas ponderadas, composite score, narrativa 5 párrafos
- Invariante verificado (sesión 7): parse_study → synthesize_study siempre encadenados

### Ligand Export (sesiones 7-8 — CLI maduro)
- **`--legacy`**: PDB ATOM-record + CONECT, validado con servidor LigParGen real
- **`--smiles`**: SMILES canónico para submit directo al campo de texto de LigParGen
- Ambos modos computan `formal_charge` con `Chem.GetFormalCharge(mol)`
- Ambos escriben `_meta.json` (con `formal_charge`) y `_charge.txt` (advisory)
- CLI muestra `Formal charge: +1` + `LigParGen charge: +1`
- WARNING panel amarillo cuando `formal_charge ≠ 0`: "In LigParGen select charge +1 instead of 0"
- Valida RMSD heavy-atom pre/post (Kabsch) — warn si > 0.05 Å
- Mensaje claro si RDKit ausente: "Activate the rdkit_env environment"
- Instalable desde cualquier directorio (`ligand*` y `utils*` en pyproject.toml)

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
| core/structural_annotation | 52 | cubierto |
| runtime/quality_classifier | 57 | cubierto |
| runtime/trajectory_ingestor | ✓ | cubierto |
| runtime/executor | 29 | cubierto |
| runtime/checkpoint_recovery | 25 | cubierto |
| runtime/study_analyzer + observable_resolver + study_layer | 132 | cubierto (sesión 7) |
| runtime/interaction_interpreter | ✓ | cubierto |
| runtime/consensus_engine | ✓ | cubierto |
| runtime/event_detector | ✓ | cubierto |
| runtime/scientific_synthesis | ✓ | cubierto |
| runtime/membrane_gates | 36 | cubierto |
| builders/clean_water | 31 | cubierto (sesión 6) |
| builders/orient_protein | 15 | cubierto |
| builders/topology_chain | ✓ | cubierto |
| ligand/normalization | 13 | cubierto |
| ligand/ligpargen_import_validator | 19 | cubierto |
| ligand/pose_rewriter | 21 | cubierto |
| ligand/export + preparation + legacy | skipped¹ | requiere rdkit_env |
| CLI/ligand export-ligpargen | 58 | cubierto (sesión 8, mocked) ← era 22 |
| benchmarks/membrane_dppc_oplsaa | 26 | pre-existing failures |

¹ Los tests de export/preparation/legacy_writer usan `pytest.importorskip("rdkit.Chem")`.
  Correr con: `conda run -n rdkit_env python -m pytest ligand/ -v`

**Total: 1036 tests, 0 fallos** (6 skipped por RDKit; benchmarks excluidos)

### Nuevos tests sesión 8 (test_cli_ligand.py: 22 → 58)
- `TestCLIExternalDirectory` — regresión subprocess desde `/tmp` (import fix)
- `TestExportSmiles` — modo `--smiles`: dispatch, SMILES en output, error cases
- `TestHelpText` — help documenta `--legacy` como validado, menciona carga
- `TestPDBAdvisory` — legacy lista companion files (smi, meta, charge.txt)
- `TestLegacyBackwardCompat` — regresión backward compat
- `TestChargeReporting` (18 tests) — neutro/positivo/negativo en legacy y smiles,
  `_charge_label`, `_write_charge_txt`, WARNING panel

---

# Análisis de virtudes y problemas

## Virtudes reales

**1. Arquitectura declarativa limpia (core)**
Pipeline YAML → SystemState → SimulationPlan → DAG bien diseñado.

**2. Study Layer + Synthesis — potencia real**
`simforge study` clasifica sistemas desde XVGs crudos con razonamiento
explícito. Ranking, estados de interacción, narrativa científica en segundos.

**3. Sin dependencias científicas pesadas (base)**
Pure Python en el compilador y el estudio. RDKit es opt-in (solo ligandos).

**4. Test suite sólida**
1036 tests con 0 fallos. CLI ligand maduró de 22 a 58 tests en esta sesión.

**5. Pipeline membrana DAG-completo**
12 steps verificados en dry-run con automation_level correcto.

**6. Ligand export validado con servidor real**
`--legacy` PDB verificado experimentalmente con LigParGen web.
Charge reporting explícito elimina el error más común (charge mismatch).

---

## Problemas actuales

### P1. cli.py con ~3100 líneas es un monolito
Crecer con más subcomandos ligand sin refactorizar aumenta la deuda.

**Split recomendado:**
- `cli_compile.py` → compile, validate, inspect, init, recompile, clean
- `cli_runtime.py` → run, dry-run, status
- `cli_study.py` → study, analyze, summary
- `cli_annotation.py` → annotate-structure
- `cli_ligand.py` → ligand subapp (ya coherente, 150 líneas)
- `cli.py` → router fino + helpers compartidos

### P2. DOS CAMINOS DE INTERPRETACIÓN — CERRADO
`parse_study` → `synthesize_study` encadenados en producción.
Test de regresión en `runtime/test_study.py::TestStudyCLIPipeline`.

### P3. FRAGMENTACIÓN DE EXECUTORS
Tres caminos de ejecución:
- `executors/shell_executor.py` — legacy
- `executors/gromacs_executor.py` — percepción GROMACS
- `runtime/executor.py` — RuntimeExecutor (default actual)

### P4. EL COMPILER NUNCA HA VISTO GROMACS REAL
Todo el workflow genera scripts que nunca se han ejecutado end-to-end.

### P5. LIGAND CLI PARCIALMENTE EXPUESTO
Solo `export-ligpargen` tiene CLI. Las fases 3-4 (importación, pose
rewriting, assembly) son solo biblioteca.

### P6. parametrize_ligand TODAVÍA ES GUIDED
El step de parametrización en el compiler sigue siendo GUIDED (manual).
La automatización requiere integrar fases 3-4 con el builder.

---

## Historial de sesiones

| Sesión | Foco principal |
|--------|---------------|
| 1-5 | Core compiler, study layer, membrane pipeline |
| 6 | clean_water AUTOMATED + WaterDeletorAdapter |
| 7 | Ligand toolkit CLI (export-ligpargen), study layer fixes |
| 8 | Import fix (pyproject.toml), --smiles mode, charge reporting, 36 tests nuevos |
