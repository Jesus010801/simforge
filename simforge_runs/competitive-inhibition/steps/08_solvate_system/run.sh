#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────
# water_model=tip3p  box_type=dodecahedron  d=1.2nm
# Paths resueltos desde DAG

ASSEMBLE_DIR="../07_assemble_system"

# Definir caja de simulación
gmx editconf \
    -f "$ASSEMBLE_DIR/complex.gro" \
    -o box.gro \
    -c \
    -d 1.2 \
    -bt dodecahedron

# Agregar agua
gmx solvate \
    -cp box.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p "$ASSEMBLE_DIR/topol.top"
