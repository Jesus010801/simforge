#! /bin/bash

echo "Este script minimiza y usa otro script para comprimir lípidos."

gfortran -o AperR AperR.f
for file in $(ls system_shrink_1.gro)
do new=$(echo $file | sed 's/.gro//g')

#####################  1  #######################
echo "Minimizando $file" #'system_shrink_1.gro':

gmx grompp -f minim.mdp -c $file -r $file -p topol.top -o $new.tpr -v -maxwarn 2
gmx mdrun -s $new.tpr -deffnm $new -nice 0 -v -nb gpu
rm \#*

echo "Comprimiendo $new" #'system_shrink_.gro':

perl inflategro-Jorge.pl $new.gro 0.95 DPP 0 $new.gro 5 area_2.dat

./AperR

AperR=$(cat areaAng2.dat)

echo El área por residuo es "$AperR" A^2

################################################
Aexp=62
while [ $AperR -ge $Aexp ]
do
  echo "Minimizando $new" #'system_shrink.gro':

  gmx grompp -f minim.mdp -c $new.gro -r $new.gro -p topol.top -o $new.tpr -v -maxwarn 2
  gmx mdrun -s $new.tpr -deffnm $new -nice 0 -v -nb gpu
  rm \#*

  echo "Comprimiendo $new" #'system_shrink_.gro':

  perl inflategro.pl $new.gro 0.95 DPP 0 $new.gro 5 area_2.dat

  ./AperR

  AperR=$(cat areaAng2.dat)

  echo El área por residuo es "$AperR" A^2

done

################################################
echo "Minimizando $new" #'system_shrink.gro':

gmx grompp -f minim.mdp -c $new.gro -r $new.gro -p topol.top -o $new.tpr -v  -maxwarn 2
gmx mdrun -s $new.tpr -deffnm final_systm_shrnk -nice 0 -v -nb gpu
rm \#*


echo "
La optimización ha terminado.

El área por lípido aproximadamente es de $AperR A^2.

Autor de este script: Camilo - JAAP

"

done
