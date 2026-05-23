#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────
# concentration=0.154M  +=NA  -=CL
INPUT_DIR="../03_solvate_system"
TOPOL_DIR="../03_solvate_system"

gmx grompp \
    -f ions.mdp \
    -c "$INPUT_DIR/solvated.gro" \
    -p "$TOPOL_DIR/topol.top" \
    -o ions.tpr \
    -maxwarn 2

echo "SOL" | gmx genion \
    -s ions.tpr \
    -o aaions.gro \
    -p "$TOPOL_DIR/topol.top" \
    -pname NA \
    -nname CL \
    -neutral \
    -conc 0.154
