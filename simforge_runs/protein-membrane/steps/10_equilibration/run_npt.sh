#!/bin/bash
# ─── NPT equilibration ───────────────────────────────────────────────────────
TOPOL_DIR="../04_generate_topology"

gmx grompp \
    -f npt.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p "$TOPOL_DIR/topol.top" \
    -o npt.tpr \
    -maxwarn 1

if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null 2>&1; then
    gmx mdrun -v -deffnm npt -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu -ntomp 10 -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto
else
    gmx mdrun -v -deffnm npt -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on
fi
