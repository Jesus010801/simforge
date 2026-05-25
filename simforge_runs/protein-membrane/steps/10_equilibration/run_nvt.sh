#!/bin/bash
# ─── NVT equilibration ───────────────────────────────────────────────────────
EM_DIR="../09_energy_minimization"
TOPOL_DIR="../08_add_ions"

gmx grompp \
    -f nvt.mdp \
    -c "$EM_DIR/em.gro" \
    -r "$EM_DIR/em.gro" \
    -p "$TOPOL_DIR/topol.top" \
    -o nvt.tpr \
    -maxwarn 1

if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null 2>&1; then
    gmx mdrun -v -deffnm nvt -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu -ntmpi 1 -ntomp $(nproc) -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto
else
    gmx mdrun -v -deffnm nvt -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on
fi
