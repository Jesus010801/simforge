#!/bin/bash
# ─── NPT equilibration ───────────────────────────────────────────────────────
# nvt.gro y nvt.cpt son outputs locales del step NVT anterior

TOPOL_DIR="../07_assemble_system"

gmx grompp \
    -f npt.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p "$TOPOL_DIR/topol.top" \
    -o npt.tpr

gmx mdrun \
    -v \
    -deffnm npt \
    -nb gpu