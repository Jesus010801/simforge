# Estado del Proyecto — SimForge

## Fecha
[18/05/2026]

## Lo que funciona
- core/ontology.py     — categorías formales completas
- core/models.py       — SystemState con Pydantic, solo estructura
- core/inference.py    — pipeline explícito: infer_system_type → infer_biological_risks → infer_analysis_gaps
- core/parser.py       — YAML → SystemState → run_inference()
- configs/hmg_competition.yaml — caso de referencia funcionando

## Decisiones de diseño tomadas
- Pydantic solo valida estructura, NO inferencia
- Inferencia vive en inference.py con orden explícito garantizado
- inferred_system_type = "competitive-inhibition" cuando hay protein + substrate + competitive_ligand
- warnings/risks/recommendations se van a separar en tres listas distintas (PENDIENTE)

## Próximo paso inmediato
Refactorizar inference.py para separar warnings, risks y recommendations
como tres listas en SystemState con un modelo InferenceResult propio.

## Caso de referencia
hmg_competition: proteína membrane_associated + partially_truncated,
sustrato HMG-CoA, inhibidor xantona, sin membrana, charmm36 + cgenff

## Para retomar con Claude
Pega este archivo al inicio de la conversación y di:
"Continuamos SimForge. Lee el estado y seguimos desde el próximo paso."
## Completado
- core/ontology.py
- core/models.py      — SystemState, Warning, Risk, Recommendation, Severity
- core/inference.py   — pipeline: infer_system_type → infer_biological_risks → infer_analysis_gaps
- core/parser.py      — YAML → SystemState → run_inference()

## Decisiones tomadas
- Pydantic solo estructura, inferencia separada en inference.py
- Tres categorías distintas: warnings / risks / recommendations
- Cada item tiene: message, target, severity (warnings/risks) o action (recommendations)
- competitive-inhibition inferido por roles: protein + substrate + competitive_ligand

## Próximo paso
Semana 3-4: protein_validator.py
- Detectar residuos faltantes
- Verificar terminales
- Detectar clashes básicos
- Integrar output al SystemState como ValidationResult

## Completado
- core/ontology.py
- core/models.py          — SystemState, Warning, Risk, Recommendation, Severity
- core/inference.py       — pipeline: infer_system_type → infer_biological_risks → infer_analysis_gaps
- core/parser.py          — YAML → SystemState → run_inference(), resolución de paths
- validators/protein_validator.py
    — parser PDB mínimo (reemplazable por MDAnalysis sin romper interfaz)
    — detección de residuos faltantes via REMARK 465
    — detección de HETATM inesperados
    — detección de terminales expuestos
    — detección de oligomerización con tolerancia 5%
    — ProteinValidationResult como contrato de salida garantizado

## Decisiones de diseño tomadas
- Pydantic solo estructura, inferencia separada en inference.py
- warnings / risks / recommendations como tres categorías distintas
- protein_validator: interfaz pública fija, implementación interna reemplazable
- Oligomerización por similitud de tamaño (tolerancia 5%), no igualdad exacta

## Próximo paso
validators/ligand_validator.py
- Leer SDF de sustrato y ligando competitivo
- Detectar: número de átomos, carga formal, presencia de anillos aromáticos
- Verificar que el archivo SDF es parseable
- Generar LigandValidationResult con mismo patrón de contrato
## Completado
- core/ontology.py
- core/models.py
- core/inference.py
- core/parser.py
- validators/protein_validator.py  — completo y funcionando
- validators/ligand_validator.py   — reasoning engine completo, parsers pendientes

## Próximo paso INMEDIATO — primera tarea de la siguiente sesión
Refactorizar ligand_validator para arquitectura de parsers modulares:

validators/
├── ligand_validator.py        ← interfaz pública + reasoning engine
└── ligand_parsers/
    ├── __init__.py
    ├── sdf_parser.py          ← extraer _parse_sdf() del validator actual
    └── pdb_parser.py          ← construir nuevo, sin CONECT

## Decisión de diseño pendiente
Los PDB de ligandos NO tienen sección CONECT.
Sin CONECT: no hay conectividad explícita.
El pdb_parser inferirá conectividad por distancia entre átomos
(si dos átomos están a <1.9Å se asume enlace).
Aromaticidad se infiere por geometría plana + elementos típicos (C,N,O en anillo).
Esto es heurística, no química exacta — RDKit lo resolverá correctamente después.

Estado del Proyecto — SimForge
Fecha
[19/05/2026]
Arquitectura actual
simforge/
├── core/
│   ├── ontology.py
│   ├── models.py
│   ├── inference.py
│   └── parser.py
├── validators/
│   ├── protein_validator.py
│   ├── ligand_validator.py        ← interfaz pública + reasoning engine
│   └── ligand_parsers/
│       ├── __init__.py
│       ├── sdf_parser.py
│       └── pdb_parser.py          ← parser PDB con inferencia de conectividad
├── descriptors/                   ← NUEVO — Molecular Descriptor Engine
│   ├── __init__.py
│   ├── topology.py                ← COMPLETO
│   ├── aromaticity.py             ← COMPLETO
│   ├── flexibility.py             ← COMPLETO
│   ├── geometry.py                ← PENDIENTE
│   └── polarity.py                ← PENDIENTE
└── configs/
    └── hmg_competition.yaml
Completado
core/

ontology.py          — categorías formales completas
models.py            — SystemState, Warning, Risk, Recommendation, Severity
inference.py         — pipeline: infer_system_type → infer_biological_risks → infer_analysis_gaps
parser.py            — YAML → SystemState → run_inference(), resolución de paths

validators/

protein_validator.py
— parser PDB mínimo (reemplazable por MDAnalysis sin romper interfaz)
— detección de residuos faltantes via REMARK 465
— detección de HETATM inesperados
— detección de terminales expuestos
— detección de oligomerización con tolerancia 5%
— ProteinValidationResult como contrato de salida garantizado
ligand_validator.py
— interfaz pública: validate_ligand(path, role) → LigandValidationResult
— reasoning engine semántico (no cambia con el parser)
— selección automática de parser por extensión (.sdf/.mol → sdf_parser, .pdb → pdb_parser)
— _analyze_chemistry() ahora delega en el descriptor engine
— campo parser_used en resultado
ligand_parsers/sdf_parser.py
— parser SDF mínimo línea por línea
— contrato de salida: {mol_name, n_atoms, n_bonds, atoms, bonds, error}
ligand_parsers/pdb_parser.py
— inferencia de conectividad por distancia con umbrales por par de elementos
— límite de valencia por elemento (evita bonds inflados)
— ordenamiento de candidatos por distancia antes de aplicar valencia
— detección de residuo objetivo (primer HETATM no-agua)
— filtro por residuo: solo átomos del mismo resname/chain/resseq
— extracción de elemento por columna 12 del formato PDB:
· line[12] == ' ' → elemento de 1 letra garantizado (C, N, O, H, S, P)
· line[12] != ' ' → intentar metal bio-relevante (Fe, Zn, Mg...)
· resuelve CA→C, CB→C, CD→C, SC→S, PF→P, OC→O, FE→Fe correctamente

descriptors/

topology.py
— compute_topology(atoms, bonds) → TopologyDescriptor
— grafo de adyacencia solo átomos pesados
— ring detection por DFS (ciclos 3-8 miembros)
— detección de anillos fusionados (≥2 átomos compartidos)
— clasificación: terminal_atoms, branching_atoms
— conteo de elementos, element_counts
aromaticity.py
— compute_aromaticity(atoms, bonds, topology) → AromaticityDescriptor
— tres niveles de evidencia: bond_type_4 > planar_geometry > composition_heuristic
— cálculo de planitud RMS por anillo
— umbrales de planitud por tamaño (5: 0.10Å, 6: 0.08Å, 7: 0.12Å)
— detección de sistemas fusionados aromáticos (Union-Find)
— enriquece RingInfo.is_aromatic in-place
flexibility.py
— compute_flexibility(atoms, bonds, topology, aromaticity) → FlexibilityDescriptor
— exclusión correcta: H-bonds, aromáticos, dobles/triples, mismo anillo, terminales
— flexibilidad efectiva: rigid/moderate/flexible/very_flexible
— flexibility_score continuo (rot_bonds / heavy_atoms)
— scaffold_rigidity: acyclic/rigid_core/flexible_core/mixed
— sampling_recommendation por clase de flexibilidad

Decisiones de diseño tomadas

Pydantic solo estructura, inferencia separada
Tres categorías distintas: warnings / risks / recommendations
Interfaz pública fija, implementación interna reemplazable (parsers → RDKit futuro)
Descriptor engine desacoplado de validators: descriptors/ es independiente
_analyze_chemistry() delega en descriptor engine y expone _topo/_arom/_flex en el dict
PDB parser: conectividad por distancia + límite de valencia + filtro por residuo
Elemento desde nombre de átomo: regla de columna 12 del formato PDB (no heurística de letras)
very_flexible añadido como cuarta clase (≥15 rot. bonds)

Casos de referencia validados

hmg_competition.yaml — sistema competitive-inhibition funcionando end-to-end
hmgcoa.pdb (GROMACS) → 80 átomos, elementos [C,H,N,O,P,S], flexible (29 rot. bonds),
aromático (9 átomos — anillo adenina del CoA), param. dif. high ✓
A1.pdb (xantona, GROMACS) → 35 átomos, elementos [C,H,O], rigid (2 rot. bonds),
aromático (14 átomos — sistema tricíclico fusionado), param. dif. low ✓

Hallazgos arquitectónicos importantes (sesión 19/05)

PDB de GROMACS contienen geometría pero poca semántica química explícita
Flujo correcto: raw atom names → chemical normalization → validated chemistry → reasoning engine
RDKit debe entrar como engine de percepción química, no solo como parser de archivos
El parser mínimo manual tiene valor: transparencia, reasoning explícito, sin dependencias pesadas
aromaticidad ≠ solo ciclos; flexibilidad ≠ solo conteo de enlaces
conectividad no es suficiente para describir comportamiento conformacional

Próximo paso inmediato
Completar el descriptor engine:
geometry.py — PENDIENTE

compute_geometry(atoms, bonds, topology, aromaticity) → GeometryDescriptor
planaridad global de la molécula (no solo por anillo)
detección de centros quirales (átomos C con 4 sustituyentes distintos)
ángulos de torsión de los enlaces rotables
distancia máxima átomo a átomo (extensión molecular)

polarity.py — PENDIENTE

compute_polarity(atoms, bonds, topology) → PolarityDescriptor
momento dipolar estimado (suma vectorial de electronegatividades × posición)
grupos funcionales detectados: OH, NH, C=O, COOH, SO3, PO4...
donors/acceptors de H-bond (reglas Lipinski)
logP estimado (método fragmentario simple, sin RDKit)

Para retomar con Claude
Pega este archivo al inicio de la conversación y di:
"Continuamos SimForge. Lee el estado y seguimos desde el próximo paso."

Estado del Proyecto — SimForge
Fecha
[19/05/2026]
Comando de ejecución
bashcd /home/jesusxd/Escritorio/simforge
python -m core.test_parser
Arquitectura actual
simforge/
├── core/
│   ├── ontology.py
│   ├── models.py
│   ├── inference.py
│   ├── parser.py
│   └── test_parser.py
├── validators/
│   ├── protein_validator.py
│   ├── ligand_validator.py
│   └── ligand_parsers/
│       ├── __init__.py
│       ├── sdf_parser.py
│       └── pdb_parser.py
├── descriptors/
│   ├── __init__.py            ← solo expone topology/aromaticity/flexibility
│   ├── topology.py            ← COMPLETO Y FUNCIONANDO
│   ├── aromaticity.py         ← COMPLETO Y FUNCIONANDO
│   ├── flexibility.py         ← COMPLETO Y FUNCIONANDO
│   ├── geometry.py            ← PRESENTE PERO ROTO — reescribir próxima sesión
│   └── polarity.py            ← PRESENTE PERO ROTO — reescribir próxima sesión
└── configs/
    └── hmg_competition.yaml
Estado de cada módulo
COMPLETADO Y FUNCIONANDO
core/ontology.py — categorías formales completas
core/models.py — SystemState, Warning, Risk, Recommendation, Severity
core/inference.py — pipeline: infer_system_type → infer_biological_risks → infer_analysis_gaps
core/parser.py — YAML → SystemState → run_inference(), resolución de paths
validators/protein_validator.py

Parser PDB mínimo, detección residuos faltantes (REMARK 465), HETATM,
terminales expuestos, oligomerización con tolerancia 5%

validators/ligand_validator.py

validate_ligand(path, role) → LigandValidationResult
Selección automática de parser por extensión
_analyze_chemistry() delega en descriptor engine
Reasoning engine semántico contextual

validators/ligand_parsers/pdb_parser.py

Conectividad por distancia + límite de valencia + filtro por residuo
Elemento desde columna 12 del formato PDB:
line[12]==' ' → 1 letra (C,N,O,H,S,P) | line[12]!=' ' → metal bio-relevante
Resuelve CA→C, CB→C, SC→S, PF→P, FE→Fe correctamente

descriptors/topology.py

compute_topology(atoms, bonds) → TopologyDescriptor
Grafo pesado, ring detection DFS, anillos fusionados, terminal/branching atoms

descriptors/aromaticity.py

compute_aromaticity(atoms, bonds, topology) → AromaticityDescriptor
Evidencia: bond_type_4 > planar_geometry > composition_heuristic
Planitud RMS por anillo, sistemas fusionados (Union-Find)

descriptors/flexibility.py

compute_flexibility(atoms, bonds, topology, aromaticity) → FlexibilityDescriptor
Clases: rigid / moderate / flexible / very_flexible
scaffold_rigidity, flexibility_score continuo, sampling_recommendation

PENDIENTE — próxima sesión
descriptors/geometry.py — reescribir completo

Función pública: compute_geometry(atoms, bonds, topology, aromaticity) → GeometryDescriptor
Planaridad global, centros quirales, extensión molecular,
shape descriptors (elongated/globular/flat), anisotropía, simetría

descriptors/polarity.py — reescribir completo

Función pública: compute_polarity(atoms, bonds, topology) → PolarityDescriptor
Dipolo estimado, grupos funcionales, H-bond donors/acceptors, logP fragmentario

Decisiones de diseño

Pydantic solo estructura, inferencia separada
Interfaz pública fija, implementación interna reemplazable
Descriptor engine desacoplado de validators
very_flexible = cuarta clase (≥15 rot. bonds)
RDKit entra después como engine de percepción, no como parser simple

Casos de referencia validados

hmgcoa.pdb → 80 átomos, [C,H,N,O,P,S], flexible (29 rot.), aromático (9), dif. high ✓
A1.pdb (xantona) → 35 átomos, [C,H,O], rigid (2 rot.), aromático (14), dif. low ✓

Git — commit sugerido
bashgit add descriptors/topology.py descriptors/aromaticity.py descriptors/flexibility.py
git add descriptors/__init__.py
git add validators/ligand_validator.py validators/ligand_parsers/pdb_parser.py
git add ESTADO_PROYECTO.md
git commit -m "feat: molecular descriptor engine (topology, aromaticity, flexibility)

- topology: grafo pesado, ring detection DFS, fusión de anillos
- aromaticity: 3 niveles de evidencia, planitud RMS, sistemas fusionados
- flexibility: exclusión semántica de bonds, 4 clases, scaffold rigidity
- pdb_parser: elemento por col12 PDB, filtro por residuo, límite de valencia
- ligand_validator: _analyze_chemistry delega en descriptor engine"
Para retomar con Claude
Pega este archivo al inicio de la conversación y di:
"Continuamos SimForge. Lee el estado y seguimos desde el próximo paso."
