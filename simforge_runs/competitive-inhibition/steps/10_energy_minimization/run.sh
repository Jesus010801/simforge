gmx grompp     -f em.mdp     -c aaions.gro     -p topol.top     -o em.tpr     -maxwarn 1

gmx mdrun     -v     -deffnm em     -nb gpu