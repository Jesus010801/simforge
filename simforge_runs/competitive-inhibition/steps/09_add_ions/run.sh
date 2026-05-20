#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────

# Preparar archivo de entrada para genion
gmx grompp \
    -f ions.mdp \
    -c ../solvate_system/solvated.gro \
    -p ../assemble_system/topol.top \
    -o ions.tpr \
    -maxwarn 2

# Agregar iones (neutralizar + 0.15M NaCl)
echo "SOL" | gmx genion \
    -s ions.tpr \
    -o aaions.gro \
    -p ../assemble_system/topol.top \
    -pname NA \
    -nname CL \
    -neutral \
    -conc 0.15
