#!/bin/bash
# ─── Parametrización CGenFF: substrate_1 ────────────────────────────────────────
# ⚠  BLOCKING: revisión manual requerida antes de continuar

# Opción A: ParamChem online (recomendado)
#   1. Ir a https://cgenff.umaryland.edu
#   2. Subir substrate_1.mol2 o substrate_1.sdf
#   3. Descargar substrate_1.str
#   4. Revisar penalizaciones (penalty score)
#      - Score < 10  → parámetros confiables
#      - Score 10-50 → revisar manualmente
#      - Score > 50  → requiere QM

# Opción B: CHARMM local (si disponible)
#   cgenff substrate_1.mol2 > substrate_1.str

# Post-parametrización: convertir a formato GROMACS
python cgenff_charmm2gmx.py substrate_1 substrate_1.mol2 substrate_1.str charmm36.ff
