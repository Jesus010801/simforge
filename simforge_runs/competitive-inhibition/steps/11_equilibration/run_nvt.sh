#!/bin/bash
# ─── NVT equilibration ───────────────────────────────────────────────────────
# Paths resueltos desde DAG

EM_DIR="../10_energy_minimization"
TOPOL_DIR="../07_assemble_system"

gmx grompp \
    -f nvt.mdp \
    -c "$EM_DIR/em.gro" \
    -r "$EM_DIR/em.gro" \
    -p "$TOPOL_DIR/topol.top" \
    -o nvt.tpr

gmx mdrun \
    -v \
    -deffnm nvt \
    -nb gpu