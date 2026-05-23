#!/bin/bash
# ─── Preparación de ligando: substrate_1 ────────────────────────────────────────
# El PDB de entrada vive en workspace/inputs/ — workspace auto-contenido.
# Ejecutar manualmente.

INPUTS_DIR="../../inputs"

# Opción A: convertir desde inputs (recomendado)
#   obabel "$INPUTS_DIR/substrate_1.pdb" -O substrate_1.sdf --gen3d

# Opción B: si ya tienes SDF listo
#   cp /path/to/substrate_1.sdf .

# Verificar estructura en Avogadro o PyMOL antes de parametrizar
echo "Verificar substrate_1.sdf antes de continuar a parametrización"
