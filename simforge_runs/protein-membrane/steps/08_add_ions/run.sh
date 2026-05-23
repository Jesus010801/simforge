#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────
# concentration=0.15M  +=NA  -=CL
INPUT_DIR="../07_clean_water"
TOPOL_DIR="../07_clean_water"

gmx grompp \
    -f ions.mdp \
    -c "$INPUT_DIR/system_clean.gro" \
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
    -conc 0.15
