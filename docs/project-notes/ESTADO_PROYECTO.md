# Estado del Proyecto — SimForge
**Actualizado:** 2026-06-20 (sesión 10)
**Líneas de código:** ~38,000 Python en ~166 archivos
**Tests:** 1317 pasando, 0 fallos, 6 skipped (benchmarks excluidos; skipped requieren rdkit_env)

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
│  PRODUCTO C — Ligand Parameterization Toolkit    ← COMPLETO      │
│  PDB complejo → prepare → LigParGen → integrate → complex.gro   │
│  Estado: pipeline completo y funcionando. Assembly validado.     │
│          prepare + integrate con CLI maduro.                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Arquitectura completa (estado real)

```
simforge/
│
├── cli.py                         ~3500 líneas  ← MONOLITO (deuda conocida)
│   ├── simforge compile/run/validate/inspect/...
│   ├── simforge study/analyze/summary/status
│   ├── simforge annotate-structure
│   ├── simforge ligand export-ligpargen  (sesión 7-8)
│   ├── simforge ligand prepare           (sesión 9-10) ← NUEVO
│   └── simforge ligand integrate         (sesión 9-10) ← NUEVO
│
├── core/
│   ├── compiler.py / parser.py / decision_engine.py / models.py
│   ├── execution_models.py / membrane_geometry.py
│   ├── structural_annotation.py / semantic_inference.py / semantic_objectives.py
│   ├── scientific_planner.py / variant_compiler.py / geometry_advisor.py
│   ├── workflow_hints.py / workspace_fingerprint.py / project_manager.py
│   ├── ligand_workflow_models.py
│   └── md_knowledge/
│
├── ligand/                        toolkit parameterización ligandos
│   ├── rdkit_reader.py            carga mol (lazy RDKit)
│   ├── export.py                  legacy/smiles + charge helpers
│   ├── preparation.py             validate_ligand_for_parameterization
│   ├── legacy_writer.py           LigParGenLegacyWriter
│   ├── normalization.py           normalize_ligpargen_outputs
│   ├── ligpargen_import_validator.py
│   ├── pose_rewriter.py           PoseRewriter (Kabsch + coord transfer)
│   ├── hydrogenation.py           backends H opcionales  ← NUEVO (sesión 9)
│   ├── prepare.py                 extracción + evaluación H desde PDB  ← NUEVO (sesión 9-10)
│   ├── integrate.py               assembly system GROMACS + pdb2gmx  ← NUEVO (sesión 10)
│   └── test_integrate.py          195 tests (prepare + integrate)  ← NUEVO
│
├── builders/ / runtime/ / executors/ / adapters/ / pipelines/ / validators/
│   (sin cambios respecto sesión 8)
```

---

## Comandos disponibles

```bash
# Workflow compiler
simforge init / compile / validate / inspect / run / dry-run / clean / recompile / status

# Analysis layer
simforge analyze / study / summary

# Structural annotation
simforge annotate-structure <config.yaml>

# Ligand toolkit — pipeline completo
simforge ligand export-ligpargen <lig.sdf> --legacy|--smiles
    # PDB/SMILES para upload a LigParGen (sesiones 7-8)

simforge ligand prepare <complex.pdb>
    # Extrae ligando de complejo docked → evalúa protonación → hydrogenación automática
    # Outputs: ligand_for_ligpargen.pdb, protein_only.pdb, ligand_report.yaml
    # --ligand-resname E20 (o auto-detect)
    # --hydrogenation auto|rdkit|openbabel|none

simforge ligand integrate \
    --protein protein_only.pdb \
    --ligand-itp LIG.itp \
    --ligand-gro LIG.gro \
    --out system/
    # Ensambla sistema GROMACS completo desde outputs de LigParGen
    # Outputs: protein.gro, protein.itp, ligand_atomtypes.itp,
    #          LIG.itp, topol.top, complex.gro, assembly_report.yaml
    # --ff oplsaa --water spce --no-grompp

# Flujo completo proteína-ligando:
#   1. simforge ligand prepare docked.pdb --ligand-resname E20
#   2. Subir ligand_for_ligpargen_H.pdb a LigParGen (charge=N)
#   3. Descargar LIG.itp + LIG.gro
#   4. simforge ligand integrate --protein protein_only.pdb \
#        --ligand-itp LIG.itp --ligand-gro LIG.gro
```

---

## Ligand Toolkit — Arquitectura completa (sesiones 7-10)

```
Fase 1 — Validación estructural
  ligand_validator.py:validate_ligand()

Fase 2a — Export para LigParGen (SDF como input)         ← sesiones 7-8
  export.py: export_for_ligpargen_legacy() / export_for_ligpargen_smiles()
  CLI: simforge ligand export-ligpargen <sdf> --legacy|--smiles

Fase 2b — Prepare desde complejo docked (PDB como input) ← sesión 9-10
  prepare.py: extract_ligand_from_complex()
    - detecta residuos HETATM automáticamente
    - evalúa hydrogenation_status: "complete"|"incomplete"|"missing"|"unknown"
    - hydrogenación automática: RDKit → Open Babel → ManualRequired
    - estimación de carga formal (RDKit, lazy)
    - escribe ligand_report.yaml con todos los campos
  hydrogenation.py: backends opcionales (ninguno hard-required)
    - RDKitHydrogenationBackend    (confidence=1.0)
    - OpenBabelHydrogenationBackend (confidence=0.9, subprocess obabel -h)
    - ManualRequiredBackend         (siempre disponible, reporta error claro)
    - select_backend(mode) → None|Backend
  CLI: simforge ligand prepare

Fase 3 — Importación outputs LigParGen
  normalization.py + ligpargen_import_validator.py
  (sin exposición CLI aún)

Fase 4 — Reescritura de pose
  pose_rewriter.py:PoseRewriter (Kabsch + coord transfer)
  (sin exposición CLI aún)

Fase 5 — Assembly sistema GROMACS                         ← sesión 10
  integrate.py: assemble_system()
    - gmx pdb2gmx con fallback adaptativo -ignh
    - split topol.top → protein.itp
    - extrae [ atomtypes ] del .itp del ligando
    - genera master topol.top (include order correcto OPLS-AA)
    - merge protein.gro + ligand.gro → complex.gro
    - validación: nombre mol, atom count, include order, grompp dry-run
    - escribe assembly_report.yaml (con metadata completo de pdb2gmx)
  CLI: simforge ligand integrate
```

---

## Novedades sesión 9-10

### `simforge ligand prepare` — extracción + evaluación H

**Problema resuelto:** el toolkit anterior solo sabía exportar desde SDF (requería RDKit).
Ahora se puede trabajar directamente con el PDB del complejo docked.

**Evaluación de protonación:**
- `has_hydrogens = n_hydrogen_atoms > 0` (preservado)
- `hydrogenation_status`: semántica en 4 niveles
  - `"complete"` — RDKit verifica que el conteo es correcto
  - `"incomplete"` — ratio H/heavy < 0.15 (ej. 28 heavy + 1 H → ratio 0.036)
  - `"missing"` — sin H explícitos
  - `"unknown"` — H presentes, ratio plausible, pero sin bond orders (PDB puro)
- Disparo automático de hydrogenación si status es `"missing"` o `"incomplete"`
- `--hydrogenation none` bloquea si status no es `"complete"` (protege contra submit incorrecto a LigParGen)

**Regression test clave:** 28 heavy atoms + 1 H → `"incomplete"` (no `"present"`).

**Umbral heurístico:** `_INCOMPLETE_RATIO_THRESHOLD = 0.15`
- 1/28 = 0.036 → `"incomplete"` ← E20
- 1/2 = 0.50  → `"unknown"` (sin RDKit no se puede probar)

**Backends de hydrogenación** (`ligand/hydrogenation.py`):
- Sin dependencia hard en Open Babel ni RDKit (todo lazy/optional)
- `select_backend("auto")` → RDKit → obabel → ManualRequired
- Mock paths limpios para tests: `ligand.hydrogenation._rdkit_available` / `_obabel_available`

**Report YAML campos nuevos:**
```yaml
has_hydrogens: true
n_hydrogen_atoms: 1
n_heavy_atoms: 28
hydrogenation_status: incomplete
hydrogenation_complete: false
hydrogenation_required: true
hydrogenation_backend: openbabel
hydrogenation_performed: true
hydrogenation_confidence: 0.9
recommended_ligpargen_input: /path/ligand_for_ligpargen_H.pdb
```

---

### `simforge ligand integrate` — assembly sistema GROMACS

**Funcionando en producción.** Probado con complejo real.

**Pipeline interno:**
1. `gmx pdb2gmx` con fallback adaptativo `-ignh`
2. Split `topol.top` → `protein.itp`
3. Extrae `[ atomtypes ]` del `.itp` del ligando → `ligand_atomtypes.itp`
4. Limpia `.itp` del ligando (sin `[ atomtypes ]`)
5. Genera `topol.top` con include order correcto (atomtypes antes que ligand.itp)
6. Merge `protein.gro + ligand.gro` → `complex.gro` (renumeración secuencial)
7. Validación (mol name, atom count, include order, grompp dry-run opcional)
8. `assembly_report.yaml`

**pdb2gmx adaptive -ignh fallback:**
- Intento 1: sin `-ignh`
- Si falla y GROMACS sugiere explícitamente `-ignh` en la salida:
  - `"Option -ignh will ignore all hydrogens in the input"`
  - `"For a hydrogen, this can be a different protonation state"`
- → Intento 2 con `-ignh` (GROMACS re-añade H desde la plantilla del forcefield)
- Errores no relacionados con H (residuo desconocido, ff ausente, PDB mal formado) fallan inmediatamente, sin retry

**CLI print:**
```
⚠  pdb2gmx failed due to protein hydrogen mismatch; retried with -ignh
   (input hydrogens ignored, GROMACS will re-add them from the force-field template)
```

**`assembly_report.yaml` incluye metadata de ambos intentos:**
```yaml
pdb2gmx_initial_attempt_exit_code: 1
pdb2gmx_initial_command: "gmx pdb2gmx -f protein.pdb ..."
pdb2gmx_initial_error_summary: "Atom HD1 in residue HIS 212 ..."
pdb2gmx_fallback_used: true
pdb2gmx_fallback_reason: "Option -ignh will ignore all hydrogens ..."
pdb2gmx_fallback_command: "gmx pdb2gmx -f protein.pdb ... -ignh"
pdb2gmx_final_exit_code: 0
```

---

## Estado del pipeline membrana (sin cambios respecto sesión 8)

| Step | automation_level | Estado |
|---|---|---|
| orient_protein | automated¹ | DONE |
| match_box_to_bilayer | automated | DONE |
| embed_in_bilayer | automated | DONE |
| generate_topology | automatic² | DONE |
| solvate_membrane | automatic² | DONE |
| clean_water | automated | DONE |
| add_ions | automatic² | DONE |
| energy_minimization → production | automatic² | DONE |

¹ Requiere `structural_annotation` completa; sin ella → GUIDED.
² `automatic` = legacy; executor lo trata idéntico a `automated`.

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
| runtime/trajectory_ingestor | ✓ | cubierto (XVG time units sesión 9) |
| runtime/executor | 29 | cubierto |
| runtime/checkpoint_recovery | 25 | cubierto |
| runtime/study_analyzer + observable_resolver | 132 | cubierto |
| runtime/interaction_interpreter + consensus + event + synthesis | ✓ | cubierto |
| runtime/membrane_gates | 36 | cubierto |
| builders/clean_water | 31 | cubierto |
| builders/orient_protein | 15 | cubierto |
| builders/topology_chain | ✓ | cubierto |
| ligand/normalization | 13 | cubierto |
| ligand/ligpargen_import_validator | 19 | cubierto |
| ligand/pose_rewriter | 21 | cubierto |
| ligand/export + preparation + legacy | skipped¹ | requiere rdkit_env |
| CLI/ligand export-ligpargen | 58 | cubierto (mocked) |
| ligand/test_integrate.py | **195** | NUEVO (sesión 9-10) |
| benchmarks/membrane_dppc_oplsaa | 26 | pre-existing failures |

¹ Correr con: `conda run -n rdkit_env python -m pytest ligand/ -v`

**Total (excl. benchmarks): 1317 tests, 0 fallos** (6 skipped por RDKit)
**Anterior sesión 8: 1036 tests → +281 nuevos**

### Nuevos tests sesión 9-10 (ligand/test_integrate.py — 195 tests)

**TestHasHydrogens** (7) — `_has_hydrogens` por columna elemento y nombre átomo

**TestHydrogenHandling** (39) — `extract_ligand_from_complex`:
- ligando con H (happy path): 8 tests
- no H + obabel éxito: 8 tests
- no H + obabel falla: 3 tests
- ambos backends no disponibles: 7 tests
- `hydrogenation_mode="none"`: 3 tests
- campos YAML report: 10 tests

**TestHydrationStatus** (25) — evaluación semántica:
- `_count_heavy_atoms`, `_determine_hydrogenation_status` unit tests
- **Regression test E20:** 28 heavy + 1 H → `"incomplete"`
- Trigger automático de hydrogenación en "incomplete"
- CLI muestra "incomplete (1 H; 28 heavy)" no "present"
- Campos YAML: `n_heavy_atoms`, `hydrogenation_status`, `hydrogenation_complete`

**TestHydrogenationCLIOption** (5) — `--hydrogenation none/openbabel/rdkit`

**TestPrepareCLI** (8) — CLI `simforge ligand prepare`

**TestHydrogenationBackendModule** (17) — `select_backend`, backends individuales

**TestPdb2gmxFallback** (25) — `simforge ligand integrate`:
- `_should_retry_with_ignh`: 7 unit tests
- `run_pdb2gmx` con mocks: 9 tests (éxito, mismatch H, error no-H, doble falla)
- `assemble_system` integración: 9 tests (fallback, reporte con 7 campos)

---

## Problemas actuales

### P1. cli.py con ~3500 líneas es un monolito (sigue creciendo)
Ver split recomendado en sesión 8. Ahora suma `ligand prepare` e `integrate`.

### P3. FRAGMENTACIÓN DE EXECUTORS
Tres caminos de ejecución (shell, gromacs, runtime). Sin cambios.

### P4. EL COMPILER NUNCA HA VISTO GROMACS REAL
El compiler genera scripts que nunca se han ejecutado end-to-end.
El toolkit `prepare`/`integrate` SÍ habla con GROMACS real.

### P6. parametrize_ligand EN EL COMPILER SIGUE SIENDO GUIDED
El step en el DAG del compiler es GUIDED. El workflow real (prepare → LigParGen →
integrate) funciona como CLI directo pero no está conectado al compiler DAG.

---

## Virtudes reales

**1. Pipeline proteína-ligando completo y funcionando**
`prepare` → LigParGen (manual) → `integrate` produce `complex.gro` + `topol.top`
listos para `gmx solvate`. Validado contra GROMACS real.

**2. Detección de protonación insuficiente**
`hydrogenation_status = "incomplete"` detecta casos como 28 heavy + 1 H (ratio 0.036).
Ya no se puede pasar un ligando mal protonado a LigParGen sin advertencia bloqueante.

**3. pdb2gmx resiliente**
El fallback `-ignh` evita que histidinas con protonación no-estándar rompan el assembly.
El reporte YAML registra ambos intentos para debugging.

**4. Backends de hydrogenación sin dependencias hard**
RDKit y Open Babel son opcionales. El sistema falla con un mensaje claro si ninguno
está disponible, en vez de silenciosamente.

**5. Test suite sólida — 1317 tests**
+281 tests desde sesión 8, cubriendo todos los nuevos flujos con mocks limpios.

---

## Historial de sesiones

| Sesión | Foco principal |
|--------|---------------|
| 1-5 | Core compiler, study layer, membrane pipeline |
| 6 | clean_water AUTOMATED + WaterDeletorAdapter |
| 7 | Ligand toolkit CLI (export-ligpargen), study layer fixes |
| 8 | Import fix (pyproject.toml), --smiles mode, charge reporting, 36 tests |
| 9 | XVG time units fix; hydrogenation backend system (ligand/hydrogenation.py + prepare.py) |
| 10 | H completeness semantics (incomplete/unknown/missing); pdb2gmx -ignh fallback; integrate CLI; assembly validado con GROMACS real |
