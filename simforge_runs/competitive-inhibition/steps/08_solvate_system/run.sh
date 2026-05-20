#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────

# Definir caja de simulación
gmx editconf \
    -f ../assemble_system/complex.gro \
    -o box.gro \
    -c \
    -d 1.2 \
    -bt dodecahedron

# Agregar agua TIP3P
gmx solvate \
    -cp box.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p ../assemble_system/topol.top
