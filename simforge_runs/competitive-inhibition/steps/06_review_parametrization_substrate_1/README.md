# Revisión Manual de Parametrización: substrate_1

## ⚠ Este step requiere intervención manual

## Qué revisar
- Revisar penalizaciones ParamChem
- Verificar cargas
- Validar constantes de fuerza

## Criterio de aceptación
- Penalty score < 10 en ParamChem para todos los átomos
- Geometría post-minimización similar al input (RMSD < 0.5Å)
- Cargas parciales razonables (sin valores > ±1.5e para átomos orgánicos típicos)

## Cuando esté listo
Marcar como completado y continuar con el siguiente step.

## Herramientas sugeridas
- ParamChem: https://cgenff.umaryland.edu
- VMD para inspección visual de cargas
- CHARMM o GROMACS para minimización de prueba
