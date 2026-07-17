# Estado del Proyecto — SimForge
**Actualizado:** 2026-07-17 — cierre de validación geométrica Protein–Membrane Builder v0.1-beta
**Líneas de código:** ~52,000 Python en ~220+ archivos
**Tests al cierre:** 36 pasando y 1 skipped en regresiones focalizadas; suite completa: 1955 pasando, 9 fallos legacy en `runtime/test_membrane_gates.py`, 9 skipped

---

## Protein–Membrane Builder v0.1-beta — fase geométrica cerrada

**Estado:** validación beta cualitativa completada. La generación automática produce ahora sistemas iniciales proteína–membrana visualmente razonables tanto para proteínas tipo receptor como para proteínas tipo canal. El pipeline todavía no es completamente desatendido ni está listo para producción larga, pero sí es adecuado para una validación controlada mediante minimización y equilibración corta. El control de calidad visual y por métricas sigue siendo obligatorio antes de iniciar MD de producción prolongada.

La validación cualitativa cubrió dos arquitecturas estructuralmente distintas:

- **Receptor tipo GLP-1R:** se eliminó el gran vacío anular artificial observado en versiones anteriores. El empaquetamiento actual es aceptable para validar equilibración corta, aunque continúa algo menos compacto que la referencia curada manualmente en algunas regiones del receptor.
- **Canal grande:** se conserva una envolvente lipídica anular externa razonable y el poro/lumen central permanece abierto, como corresponde a un canal. No se observa un cráter anular catastrófico y el sistema es apto para validación corta de minimización/equilibración.

No se continuará ajustando geometría estática en esta fase salvo que aparezca un defecto bloqueante.

### Correcciones completadas

1. **Paridad del cutoff de InflateGRO.** El cutoff inicial de exclusión lipídica se corrigió a **1.4 nm**, en paridad con el método original de protmemfiles. El cutoff del shrink loop se corrigió a **0.0 nm**, evitando la eliminación repetida de lípidos durante la compactación.
2. **Sincronización de topología.** Se añadió sincronización de topología después de operaciones sobre coordenadas que cambian el número de lípidos, evitando incompatibilidades de conteo de átomos entre coordenadas y topología.
3. **Refinamiento de exclusión pre-shrink.** Se sustituyó la eliminación agresiva por solapamiento contra toda la proteína por una exclusión consciente de la región TM/slab, reduciendo la eliminación errónea de lípidos anulares.
4. **Refinamiento del detector de lípidos atrapados.** Se separaron los lípidos realmente atrapados en cavidades de los contactos anulares normales con la superficie proteica.
5. **Infraestructura conservadora de reparación anular.** La reparación permanece desactivada por defecto y la inserción también permanece desactivada. Se corrigió su integración de compilación y ejecución: la ruta desactivada escribe un reporte sin modificar el sistema ni fallar; la ruta activada utiliza residuos TM expandidos durante el build y no variables `TM_RESIDUES` sin resolver.
6. **Interpretación consciente de canales.** Los poros y lúmenes centrales legítimos no deben penalizarse como vacíos de membrana. El QC futuro debe distinguir un vacío anular patológico de un poro funcional válido.

### Limitaciones actuales

- La geometría estática ha mejorado sustancialmente, pero no es perfecta.
- El empaquetamiento generado para GLP-1R sigue siendo ligeramente menos compacto que la referencia manual.
- Receptores y canales requieren lógica de QC diferente.
- La MD de producción larga todavía no está validada.
- La siguiente fase debe medir estabilidad dinámica, no continuar acumulando parches visuales sobre geometría estática.

### Siguiente fase — validación dinámica corta

Ejecutar minimización, NVT restringido, NPT restringido y una equilibración corta sin restricciones o con restricciones débiles. Monitorizar:

- convergencia de energía;
- estabilidad de presión;
- área por lípido (APL);
- espesor de membrana;
- deriva e inclinación de la proteína;
- penetración de lípidos en cavidades de receptores o canales;
- cobertura de lípidos anulares;
- retención o eliminación de agua dentro de poros legítimos;
- consistencia entre topología y coordenadas.

---

## ¿Qué es SimForge hoy?

Tres productos en un repositorio:

```
┌─────────────────────────────────────────────────────────────────┐
│  PRODUCTO A — Workflow Compiler                                  │
│  YAML → IR → DAG → scripts GROMACS → ejecución                  │
│  Estado: Protein–Membrane Builder v0.1-beta validado             │
│          cualitativamente; listo para validación dinámica corta  │
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
├── cli.py                         ~3800 líneas  ← MONOLITO (deuda conocida)
│   ├── simforge compile/run/validate/inspect/...
│   ├── simforge study/analyze/summary/status
│   ├── simforge annotate-structure
│   ├── simforge ligand export-ligpargen  (sesión 7-8)
│   ├── simforge ligand prepare           (sesión 9-10)
│   ├── simforge ligand integrate         (sesión 9-10)
│   ├── simforge analyze-md               (sesión 13)  ← NUEVO
│   ├── simforge fel run                  (sesión 13)  ← NUEVO
│   └── simforge fel extract-minima       (sesión 14)  ← NUEVO
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
│   (membrane pipeline + topology refactor + spatial classifier)
│
├── analysis/                      paquete análisis post-producción
│   ├── md/                        MD run discovery + observable planning (sesión 13)
│   │   ├── models.py / discovery.py / observables.py / provenance.py / report.py / cli.py
│   └── fel/                       Free Energy Landscape pipeline (sesiones 13-14)
│       ├── models.py              FELWarning, XVGSeries, FELConfig, FELSurface, FELMinimum, ExtractConfig…
│       ├── xvg.py                 load_xvg_series, default_unit_for_label
│       ├── features.py            build_feature_table (time alignment, ps↔ns auto-correct)
│       ├── surface.py             compute_fel_surface (numpy + fallback puro Python)
│       ├── workflow.py            orchestrador FEL run
│       ├── minima.py              detect_local_minima, deduplicación greedy
│       ├── frame_map.py           map_frames_to_minima (nearest-frame determinístico)
│       ├── extract.py             run_extract_minima, filtro ΔG, enriquecimiento counts
│       ├── extract_report.py      report.md con tabla minima + rejected
│       ├── extract_cli.py         fel_extract_minima_fn (Typer)
│       ├── report.py / provenance.py / cli.py
│       └── tests/analysis/fel/test_fel.py   92 tests
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

# MD run discovery
simforge analyze-md <run_dir> [--out DIR] [--dry-run]
    # Descubre archivos en un directorio MD existente
    # Outputs: analysis/metadata/ + report.md con plan de observables

# Free Energy Landscape
simforge fel run --xvg-x rmsd.xvg --xvg-y rg.xvg \
    --temperature 300 --bins 50,50 --out fel_analysis/
    # Genera 2D FEL por inversión Boltzmann
    # Outputs: data/free_energy_surface.csv + histogram_counts.csv + features.csv + plots/

simforge fel run <run_dir> --features rmsd,rg
    # Modo discovery: busca XVGs en run_dir

simforge fel extract-minima FEL_DIR \
    --trajectory md.xtc --topology md.tpr \
    --dry-run                        # planifica sin llamar gmx
    # ó
    --execute                        # extrae estructuras con gmx trjconv
    --n-minima 5                     # max minima (default 5)
    --max-delta-g-kj 2.5             # filtro ΔG (none = sin filtro)
    --min-count 1                    # excluir bins con pocos frames
    --system-name HMGCoA_R           # nombre del sistema para filename
    --compat-receptor-name           # crea alias receptor.pdb
    --force                          # sobrescribe directorio existente
    # Outputs: representative_states/state_001/{name}__FEL-M01__t44p7ns.pdb
    #          + metadata.json + frame_mapping.csv + report.md
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

## Novedades sesión 11 — Topology Refactor (2026-06-26)

**Problema resuelto:** `generate_topology` llamaba pdb2gmx sobre `embed_in_bilayer/system.gro`
(GRO mixto proteína+lípidos), perdiendo la semántica de cadena/TER/terminus → error
`Atom OXT in residue LYS 439 was not found in rtp entry LYSH`.

**Nueva arquitectura:**
```
inputs/protein_1.pdb
  → generate_protein_topology (pdb2gmx solo sobre PDB original)
  → orient / match_box / embed_in_bilayer
  → generate_topology (topology:assemble_embed — copia system.gro, ensambla topología)
  → membrane_embedding (shrink loop)
  → assemble_system_topology (lee converged.gro, counts finales)
  → solvate_membrane
```

**Nuevos archivos:** `core/topology_models.py`, `validators/topology_guardrail.py`,
`tests/test_topology_assembly.py` (76 tests), `tests/test_gromacs_error_classification.py`.

**Fix crítico:** `preparation_builder.py` — `source_file` de `protein.file` es ruta absoluta
(el parser la resuelve). Ahora usa `os.path.relpath(source_file, step_dir)`.

---

## Novedades sesión 12 — DPP ITP + Atom Count + GLP-1R (2026-07-07)

### Fix 1: DPP ITP desde template RTP (sin pdb2gmx en lípidos)

**Regla:** nunca llamar pdb2gmx sobre GRO mixtos o solo-lípidos.

- `core/lipid_itp_builder.py` (NUEVO) — parsea `aminoacids.rtp`, genera ITP completo:
  50 átomos, 49 enlaces, 57 ángulos, 53 dihedros (Ryckaert-Bellemans func=3), 53 pares, 2 impropios.
- `docs/Prot-Memb_FILES/oplsaa_membrane.ff/dpp.itp` — asset estático pre-construido.
- `assembly_builder.py` — usa `_LIPID_ITP_MAP = {'DPP': 'dpp.itp', 'DPPC': 'dpp.itp'}` en lugar de pdb2gmx.
- `tests/test_membrane_embedding_topology_bootstrap.py` (NUEVO) — 24 tests.

### Fix 2: Atom count mismatch proteína (540 heavy ≠ 1079 con-H)

**Causa:** `membrane_orient_builder.py` usaba el PDB original (solo heavy atoms) como input de editconf.
`protein_processed.gro` de pdb2gmx tiene hidrógenos añadidos → 1079 átomos. El sistema embebido
tenía solo 540 átomos de proteína, causando error fatal de grompp (`work.gro 26140 ≠ topology 26679`).

**Fix:** `membrane_orient_builder.py` ahora usa `protein_processed.gro` del step
`generate_protein_topology` como input de editconf. Fallback al PDB original solo si ese step no está en el DAG.

### Fix 3: Topology lipid count sync post-InflateGRO

InflateGRO elimina lípidos que solapan con la proteína en la inflación inicial (cutoff=1.4Å).
La topología tenía 512 DPP pero el GRO inflado tenía menos → mismatch grompp.

`embedding_builder.py` — lee `ORIG_LIPID_COUNT` ANTES de InflateGRO y llama
`update_topology_lipid_count` inmediatamente después de la inflación inicial si los counts difieren.

### Fix 4: Shrink loop MDP usa reaction-field

`coulombtype = PME` → `reaction-field` + `epsilon-rf = 0`.
**Por qué:** GROMACS 2025.2 trata PME + carga neta como FATAL ERROR. Reaction-field es suficiente
para minimización steepest descent (eliminación de clashes).

### GLP-1R run — clean_water gate false positives (2026-06-24)

**Fix 1:** `validate_no_water_in_bilayer` usaba `min/max(Z headgroups)` como boundaries del core
→ abarca superficies externas de ambas leaflets, causando falsos positivos masivos (95759+ aguas
flaggeadas). Ahora usa midplane por leaflet:
- `z_midplane` = mean Z de átomos de cola
- `z_core_bot` = mean Z de headgroups leaflet inferior (Z ≤ midplane)
- `z_core_top` = mean Z de headgroups leaflet superior (Z > midplane)

**Fix 2:** `water_gate` leía `final_water_count` (total aguas en sistema) como "aguas en core bilayer".
Nuevo campo: `n_water_oxygens_remaining_in_core`. Reports legacy sin este campo: warning, nunca bloquea.

**Nuevos campos en water_report.json:** `bilayer_midplane_z`, `core_z_min`, `core_z_max`,
`n_water_molecules_removed`, `n_water_oxygens_remaining_in_core`, `cleanup_passed`.

### Tests (sesión 12)
- `tests/test_topology_assembly.py` + `test_membrane_embedding_topology_bootstrap.py`: **142 passed**
- `test_membrane_integration.py`: **66 passed**
- builders/adapters/pipelines/runtime/core/tests: **1215 passed, 3 skipped**

---

## Novedades sesión 12B — Phase 10B Spatial Classifier (2026-07-08)

Clasificador 3D de ocupación espacial membrana-proteína para `clean_water` pore-aware.

**Algoritmo:** Grid numpy 3D (spacing 0.15 nm) + flood-fill por dilatación iterativa desde
semillas bulk superior/inferior. Agua dentro del footprint XY del bundle TM que está conectada
al bulk superior o inferior → `pore_water` (preservar). Agua en core membrana fuera del footprint
TM → `membrane_core_water` (eliminar). Lípido COM dentro del footprint TM → `forbidden_pore_lipid`.

**No usa scipy** — dilatación iterativa pura Python+numpy.

**Archivos nuevos:**
- `validators/membrane_protein_spatial_classifier.py`
- `tests/test_membrane_protein_spatial_classifier.py` — 12 tests
- `tests/test_clean_water.py` — tests pore-aware

**Integración:** `validators/lipid_refill.py` — acepta `spatial_classifier_report` para
saltar gap_clusters en la región de poro.

**Tests:** 1236 (previos) + 21 nuevos = **1257 total, 0 fallos**.

---

## Novedades sesión 13 — Analysis/FEL Foundation (2026-07-08)

### MD Run Discovery (`analysis/md/`)

Nuevo paquete para descubrimiento de runs MD existentes y planificación de observables.

```
analysis/md/
  models.py      — DiscoveredFile, AnalysisWarning, DiscoveredMDRun
  discovery.py   — discover_md_run(run_dir): scan recursivo .tpr/.xtc/.trr/.gro/.xvg/.edr/…
  observables.py — plan_observables(run): 7 observables (rmsd/rmsf/rg/sasa/hbonds/energy_qc/existing_xvg)
  provenance.py  — Provenance, AnalysisArtifact, build_provenance()
  report.py      — write_report(): Markdown con archivos detectados + comandos gmx propuestos
  cli.py         — analyze_md_fn (Typer plain function)
```

CLI: `simforge analyze-md <run_dir> [--out DIR] [--dry-run]`

**Nota:** `simforge analyze` (calidad XVG) ya existe como `@cli.command`. Typer 0.17.4 no permite
el mismo nombre que una sub-app → `analyze-md` evita rotura de backward compat.

### FEL Pipeline (`analysis/fel/`)

Free Energy Landscape nativo Python desde XVGs de GROMACS.

```
Flujo: XVG parsing → alineación temporal → histograma 2D → inversión Boltzmann
```

**Módulos:**
- `xvg.py` — `load_xvg_series(path, column)`: selección de columna 1-based, detección de unidades, errores estrictos. `default_unit_for_label(label)`: word-boundary-safe (evita que "Energy" devuelva "nm").
- `features.py` — `build_feature_table(sx, sy)`: alineación temporal, auto-corrección ps↔ns.
- `surface.py` — `compute_fel_surface()`: histograma numpy 2D + fallback puro Python; k_B = 0.00831446261815324 kJ/mol/K.
- `workflow.py` — orquesta el run completo, escribe todos los archivos de salida.

**Estructura de salida:**
```
fel_analysis/
  config.json, report.md
  metadata/  provenance.json, warnings.json
  data/  features.csv, histogram_counts.csv, probability.csv, free_energy_surface.csv
  minima/  README.md
  plots/  fel_contour.png  (si matplotlib disponible)
```

**Diseño:**
- Corrección automática de unidades cuando header dice "ps" pero valores son de rango ns (o viceversa)
- `time_tolerance=1e-3 ns` para tolerancia fp en datos GROMACS reales
- `tests/__init__.py` removido de `tests/analysis/` — evita colisión con el paquete `analysis`

**Tests:** 28 FEL tests. **Total: 1871 passed, 9 skipped**.

---

## Novedades sesión 14 — FEL Phase 3/3B/3C: Extracción de Minima (2026-07-08)

Pipeline completo de extracción de estructuras representativas desde minima FEL.

### Phase 3 — Extract-Minima Core

**Nuevos módulos:**

`analysis/fel/minima.py` — Detección de minima locales:
- `detect_local_minima(fel_dir, n_minima, min_distance_bins)`
- 8-vecinos: bin es mínimo local si TODOS los vecinos finitos tienen G estrictamente mayor (`g2 <= g` descalifica → superficies planas sin mínimos)
- Deduplicación greedy por distancia euclidiana en espacio de bins
- Fallback `_global_minimum` cuando no hay mínimos locales → warning `no_local_minima`

`analysis/fel/frame_map.py` — Mapping frames→minima:
- `map_frames_to_minima(minima, features_csv, x_label, y_label, x_unit, y_unit)`
- Frame más cercano en espacio (x,y); tie-breaking determinístico: dist → time_ns → row_index

`analysis/fel/extract.py` — Orquestador:
- `run_extract_minima(config)` → lee config.json del FEL dir, detecta mínimos, mapea frames, extrae con `gmx trjconv -dump {time_ps:.3f}`
- `_process_state(minimum, mapping, config, state_dir, x_label, y_label, system_name)`

`analysis/fel/extract_report.py` — Reporte:
- Tabla: `Display Name | X | Y | ΔG (kJ/mol) | Time (ns) | Output File | Status`
- Sub-tabla "Rejected by ΔG filter" cuando hay minima rechazados

`analysis/fel/extract_cli.py` — CLI `simforge fel extract-minima`

### Phase 3B — Fixes de Producción (HMG-CoA-R real run)

**A. Unidades nulas cuando el ylabel XVG no tiene anotación `(nm)`:**
`default_unit_for_label(label)` en `xvg.py` — fallback por label semántico:
- Exact set: `{"rmsd", "rg", "rmsf"}` → "nm"
- Starts-with: `("rmsd", "rmsf", "radius")` → "nm"
- Contains: `("gyration",)` → "nm"
- Regex `\brg\b` para "rg" standalone (evita matchear "ene**rg**y")

**B. `--force` no limpiaba dirs stale (state_004, state_005):**
`shutil.rmtree(out_dir)` ejecutado ANTES de cualquier `mkdir` cuando `config.force=True`.

**C. Filtro energético — rechazar mínimos matemáticos con alto ΔG:**
`_apply_energy_filter(minima, max_delta_g_kj)` → `(kept, rejected)`.
- Default: `max_delta_g_kj = 2.5` kJ/mol
- Preserva siempre el mínimo global (minima[0]) aunque todos excedan el umbral
- `rejected_minima` escrito en `metadata/rejected_minima.json`

**D. Count/probability desde `histogram_counts.csv` + filtro `--min-count`:**
`_enrich_with_counts(minima, counts_csv, min_count, warnings)` — lee counts CSV, adjunta
`count` y `probability` a cada `FELMinimum`, filtra bins con pocos frames.

**E. Mensaje de error claro cuando output dir existe:**
`"Output directory already exists. Use --force to overwrite or choose a new --out directory."`

### Phase 3C — Filenames Descriptivos

Patrón: `{system_name}__FEL-M{id:02d}__t{time_safe}ns.{ext}`

Ejemplos reales (HMG-CoA-R):
```
HMGCoA_R__FEL-M01__t44p7ns.pdb
HMGCoA_R__FEL-M02__t51p6ns.pdb
HMGCoA_R__FEL-M03__t63p1ns.pdb
```

**Funciones:**
- `sanitize_system_name(name)` — quita chars fuera de `[A-Za-z0-9_-]`, colapsa `__`
- `format_time_ns(t, decimals=1)` — `f"{t:.1f}".replace(".", "p")` (sin puntos en el stem)
- `build_output_filename(system_name, minimum_id, time_ns, fmt)` — construye el nombre completo
- `resolve_system_name(config)` — cadena de inferencia: `config.system_name` → `fel_dir.name` → `topology.stem` → `trajectory.stem` → `"system"`

**Metadata.json incluye:** `display_name`, `system_name`, `filename_pattern`, `output_filenames`,
`output_paths`, `count`, `probability` (cuando disponibles).

**`--compat-receptor-name`:** crea symlink `receptor.{ext}` → nombre descriptivo (no por defecto).

**No hay ΔG en el filename** — toda la información científica queda en metadata.json y report.md.

### Tests finales
- 92 tests FEL en `tests/analysis/fel/test_fel.py`
- Clases: `TestMinimaDetection`, `TestFrameMapping`, `TestExtractMinima`, `TestUnitMetadata`,
  `TestPhase3BUnitDefaults`, `TestPhase3BForceCleanup`, `TestPhase3BEnergyFilter`,
  `TestPhase3BCountFilter`, `TestPhase3CFilenames`
- **Total: 1919 passed, 9 skipped, 0 fallos**

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
| ligand/test_integrate.py | **195** | cubierto (sesión 9-10) |
| tests/test_topology_assembly.py | 76 | cubierto (sesión 11) |
| tests/test_membrane_embedding_topology_bootstrap.py | 24 | cubierto (sesión 12) |
| tests/test_membrane_protein_spatial_classifier.py | 12 | cubierto (sesión 12B) |
| tests/test_clean_water.py | ~9 | cubierto (sesión 12B) |
| tests/analysis/test_md_discovery.py | 9 | cubierto (sesión 13) |
| tests/analysis/fel/test_fel.py | **92** | cubierto (sesiones 13-14) |
| benchmarks/membrane_dppc_oplsaa | 26 | pre-existing failures |

¹ Correr con: `conda run -n rdkit_env python -m pytest ligand/ -v`

**Total (excl. benchmarks): 1919 tests, 0 fallos** (9 skipped por RDKit)
**Anterior sesión 10: 1317 tests → +602 nuevos (sesiones 11-14)**

### Nuevos tests sesiones 11-14

| Sesión | Módulo | Tests |
|--------|--------|-------|
| 11 | `tests/test_topology_assembly.py` + gromacs error classification | 76 |
| 12 | `tests/test_membrane_embedding_topology_bootstrap.py` | 24 |
| 12 | `runtime/test_membrane_gates.py` — clean_water gate regression | +4 |
| 12B | `tests/test_membrane_protein_spatial_classifier.py` | 12 |
| 12B | `tests/test_clean_water.py` (pore-aware) | 9 |
| 13 | `tests/analysis/test_md_discovery.py` | 9 |
| 13 | `tests/analysis/fel/test_fel.py` (Phase 2 FEL) | 28 |
| 14 | `tests/analysis/fel/test_fel.py` (Phase 3/3B/3C) | +64 |

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

**5. Free Energy Landscape pipeline completo**
`simforge fel run` + `simforge fel extract-minima` producen estructuras representativas de
mínimos energéticos con nombres descriptivos, filtro ΔG, enriquecimiento de counts, y report.md.
Verificado con datos reales de HMG-CoA-R.

**6. Test suite sólida — 1919 tests**
+602 tests desde sesión 10, cubriendo topology refactor, DPP ITP, spatial classifier, FEL pipeline completo.

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
| 11 | Topology refactor — pdb2gmx solo sobre PDB original; generate_protein_topology + assemble_system_topology |
| 12 | DPP ITP desde RTP (sin pdb2gmx en lípidos); atom count fix (protein_processed.gro); topology lipid sync post-InflateGRO; reaction-field MDP; GLP-1R run + clean_water gate fixes |
| 12B | Phase 10B spatial classifier 3D flood-fill pore/lipid; pore-aware clean_water; lipid_refill integration |
| 13 | analysis/md/ discovery; analysis/fel/ Phase 2 — XVG parsing + 2D histogram + Boltzmann + FEL surface; simforge fel run CLI; 1871 tests |
| 14 | FEL Phase 3 (minima detection + frame mapping + gmx trjconv extraction); Phase 3B (unit defaults, --force rmtree, ΔG filter, count enrichment); Phase 3C (filenames descriptivos HMGCoA_R__FEL-M01__t44p7ns.pdb); 1919 tests |
