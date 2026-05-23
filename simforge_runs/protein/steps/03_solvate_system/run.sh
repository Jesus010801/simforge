#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────
# water_model=spce  box_type=triclinic  d=1.2nm
# Paths resueltos desde DAG

ASSEMBLE_DIR="../02_assemble_system"

# Definir caja de simulación
gmx editconf \
    -f "$ASSEMBLE_DIR/complex.gro" \
    -o box.gro \
    -c \
    -d 1.2 \
    -bt triclinic

# Agregar agua
gmx solvate \
    -cp box.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p "$ASSEMBLE_DIR/topol.top"
