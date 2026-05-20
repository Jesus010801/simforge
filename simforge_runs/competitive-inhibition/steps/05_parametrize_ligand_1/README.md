# Parametrización CGenFF: ligand_1

## Engine
`cgenff` — CHARMM General Force Field

## Estado
Automático con revisión recomendada

## Outputs esperados
- `ligand_1.str`      — parámetros CGenFF
- `ligand_1.itp`      — topología GROMACS
- `ligand_1.prm`      — parámetros en formato CHARMM

## Criterio de aceptación
Penalty score < 10 en ParamChem para todos los átomos.
