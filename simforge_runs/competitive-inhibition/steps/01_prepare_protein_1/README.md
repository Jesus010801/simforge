# Preparación: protein_1

## Qué hace este step
Convierte el PDB de la proteína a formato GROMACS con topología completa.

## Engine
`gromacs:pdb2gmx`

## Notas
- Agregar hidrógenos
- Asignar protonación
- Construir topología

## Outputs esperados
- `protein_1_processed.gro`
- `topol.top`
- `posre.itp`

## Cómo ejecutar
```bash
bash commands.sh
```
