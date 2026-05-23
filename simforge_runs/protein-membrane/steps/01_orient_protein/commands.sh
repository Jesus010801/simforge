#!/bin/bash
# ─── Preparación de ligando: ligand ────────────────────────────────────────
# Engine: ligand_preparation
# Ejecutar manualmente

# Opción A: si tienes SDF limpio
#   obabel ligand.pdb -O ligand.sdf --gen3d

# Opción B: si ya tienes SDF
#   cp ligand.sdf .

# Verificar estructura en Avogadro o PyMOL antes de parametrizar
echo "Verificar ligand.sdf antes de continuar a parametrización"
