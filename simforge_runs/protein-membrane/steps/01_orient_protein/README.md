# Preparación: ligand

## Qué hace este step
Prepara el ligando para parametrización.
Convierte de PDB a SDF con conectividad explícita.

## Recomendación
Usar OpenBabel o RDKit para conversión limpia.
Verificar visualmente en Avogadro o PyMOL.

## Notas
- Rotation angle depends on structure — cannot be inferred automatically in v1

## Outputs esperados
- `ligand.sdf`
