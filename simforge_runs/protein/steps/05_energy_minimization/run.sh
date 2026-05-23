#!/bin/bash
# ─── Energy minimization ─────────────────────────────────────────────────────
IONS_DIR="../04_add_ions"
TOPOL_DIR="../02_assemble_system"

gmx grompp \
    -f em.mdp \
    -c "$IONS_DIR/aaions.gro" \
    -p "$TOPOL_DIR/topol.top" \
    -o em.tpr \
    -maxwarn 1

if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null 2>&1; then
    gmx mdrun -v -deffnm em -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu -ntmpi 1 -ntomp $(nproc) -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto
else
    gmx mdrun -v -deffnm em -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on
fi