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
