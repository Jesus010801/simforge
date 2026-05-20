#!/bin/bash
# ─── Parametrización CGenFF: ligand_1 ────────────────────────────────────────


# Opción A: ParamChem online (recomendado)
#   1. Ir a https://cgenff.umaryland.edu
#   2. Subir ligand_1.mol2 o ligand_1.sdf
#   3. Descargar ligand_1.str
#   4. Revisar penalizaciones (penalty score)
#      - Score < 10  → parámetros confiables
#      - Score 10-50 → revisar manualmente
#      - Score > 50  → requiere QM

# Opción B: CHARMM local (si disponible)
#   cgenff ligand_1.mol2 > ligand_1.str

# Post-parametrización: convertir a formato GROMACS
python cgenff_charmm2gmx.py ligand_1 ligand_1.mol2 ligand_1.str charmm36.ff
