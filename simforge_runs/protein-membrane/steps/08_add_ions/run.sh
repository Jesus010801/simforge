#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────
# concentration=0.15M  +=NA  -=CL
INPUT_DIR="../07_clean_water"
TOPOL_SRC="../07_clean_water"

# Copiar topología — gmx genion la modifica in-place (reemplaza SOL por iones).
# La copia local garantiza que el paso anterior quede sin modificar.
cp "$TOPOL_SRC/topol.top" topol.top

gmx grompp \
    -f ions.mdp \
    -c "$INPUT_DIR/system_clean.gro" \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 2

echo "SOL" | gmx genion \
    -s ions.tpr \
    -o aaions.gro \
    -p topol.top \
    -pname NA \
    -nname CL \
    -neutral \
    -conc 0.15
