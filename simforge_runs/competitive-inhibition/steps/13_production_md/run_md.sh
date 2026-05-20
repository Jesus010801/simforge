gmx grompp     -f md.mdp     -c npt.gro     -t npt.cpt     -p topol.top     -o md.tpr

gmx mdrun     -v     -deffnm md     -nb gpu