# Parametrización CGenFF: substrate_1

## Engine
`cgenff` — CHARMM General Force Field

## Estado
⚠ **BLOCKING** — revisión manual requerida antes de producción

## Outputs esperados
- `substrate_1.str`      — parámetros CGenFF
- `substrate_1.itp`      — topología GROMACS
- `substrate_1.prm`      — parámetros en formato CHARMM

## Criterio de aceptación
Penalty score < 10 en ParamChem para todos los átomos.
