#!/bin/bash
# ─── Preparación de ligando: substrate_1 ────────────────────────────────────────
# Engine: ligand_preparation
# Ejecutar manualmente

# Opción A: si tienes SDF limpio
#   obabel substrate_1.pdb -O substrate_1.sdf --gen3d

# Opción B: si ya tienes SDF
#   cp substrate_1.sdf .

# Verificar estructura en Avogadro o PyMOL antes de parametrizar
echo "Verificar substrate_1.sdf antes de continuar a parametrización"
