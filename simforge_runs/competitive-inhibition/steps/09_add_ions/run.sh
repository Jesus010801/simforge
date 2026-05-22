#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────
# concentration=0.15M  +=NA  -=CL
# Paths resueltos desde DAG

SOLVATE_DIR="../08_solvate_system"
ASSEMBLE_DIR="../07_assemble_system"

gmx grompp \
    -f ions.mdp \
    -c "$SOLVATE_DIR/solvated.gro" \
    -p "$ASSEMBLE_DIR/topol.top" \
    -o ions.tpr \
    -maxwarn 2

echo "SOL" | gmx genion \
    -s ions.tpr \
    -o aaions.gro \
    -p "$ASSEMBLE_DIR/topol.top" \
    -pname NA \
    -nname CL \
    -neutral \
    -conc 0.15
