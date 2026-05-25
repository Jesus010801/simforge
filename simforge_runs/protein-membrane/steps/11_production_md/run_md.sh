#!/bin/bash
# ─── Production MD ───────────────────────────────────────────────────────────
# Paths resueltos desde DAG

EQ_DIR="../10_equilibration"
TOPOL_DIR="../08_add_ions"

gmx grompp \
    -f md.mdp \
    -c "$EQ_DIR/npt.gro" \
    -t "$EQ_DIR/npt.cpt" \
    -p "$TOPOL_DIR/topol.top" \
    -o md.tpr \
    -maxwarn 1

if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null 2>&1; then
    gmx mdrun -v -deffnm md -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu -ntmpi 1 -ntomp $(nproc) -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto
else
    gmx mdrun -v -deffnm md -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on
fi