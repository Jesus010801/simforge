#!/bin/bash
# ─── Preparación de proteína: peptide_1 ───────────────────────────────────────
# El PDB de entrada vive en workspace/inputs/ — workspace auto-contenido.

INPUTS_DIR="../../inputs"

gmx pdb2gmx \
    -f "$INPUTS_DIR/peptide_1.gro" \
    -o peptide_1_processed.gro \
    -p topol.top \
    -ff oplsaa_membrane \
    -water spce \
    -ignh
