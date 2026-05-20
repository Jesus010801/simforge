# REST2 Enhanced Sampling

## Engine
`plumed:gromacs` — Replica Exchange with Solute Tempering (REST2)

## Notas
- Enhanced conformational sampling
- Replica exchange with solute tempering

## Requisitos
- GROMACS compilado con soporte MPI
- PLUMED instalado y patcheado en GROMACS

## Configuración básica REST2
```bash
# 1. Definir réplicas (temperatura efectiva del soluto)
#    Típico: 4-8 réplicas entre 300K y 450K

# 2. Preparar MDP para cada réplica
#    Ver rest2_template.mdp

# 3. Configurar PLUMED
#    Ver plumed.dat

# 4. Lanzar
mpirun -np 4 gmx_mpi mdrun \
    -v \
    -deffnm rest2 \
    -multidir rep0 rep1 rep2 rep3 \
    -replex 1000 \
    -hrex \
    -plumed plumed.dat
```

## Referencias
- REST2: Terakawa et al. JCTC 2011
- PLUMED: https://www.plumed.org
