#!/bin/bash
# ─── Preparación de ligando: ligand_1 ────────────────────────────────────────
# Engine: ligand_preparation
# Ejecutar manualmente

# Opción A: si tienes SDF limpio
#   obabel ligand_1.pdb -O ligand_1.sdf --gen3d

# Opción B: si ya tienes SDF
#   cp ligand_1.sdf .

# Verificar estructura en Avogadro o PyMOL antes de parametrizar
echo "Verificar ligand_1.sdf antes de continuar a parametrización"
