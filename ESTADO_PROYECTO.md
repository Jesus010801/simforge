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
