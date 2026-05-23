#!/bin/bash
# ─── Production MD ───────────────────────────────────────────────────────────
# Paths resueltos desde DAG

EQ_DIR="../11_equilibration"
TOPOL_DIR="../07_assemble_system"

gmx grompp \
    -f md.mdp \
    -c "$EQ_DIR/npt.gro" \
    -t "$EQ_DIR/npt.cpt" \
    -p "$TOPOL_DIR/topol.top" \
    -o md.tpr

gmx mdrun \
    -v \
    -deffnm md \
    -nb gpu