#!/bin/bash
# ─── Preparación de proteína: protein_1 ───────────────────────────────────────
# Engine: gromacs:pdb2gmx
# Ejecutar manualmente — requiere selección interactiva de forcefield

# 1. Verificar archivo de entrada
#    El archivo PDB debe estar limpio (sin HETATM inesperados, sin cadenas rotas)

# 2. Generar topología con pdb2gmx
gmx pdb2gmx \
    -f protein_1.pdb \
    -o protein_1_processed.gro \
    -p topol.top \
    -ignh \
    -ter

# Flags importantes:
#   -ignh     → ignorar hidrógenos existentes, agregar nuevos
#   -ter      → modo interactivo para terminales (N y C)
#
# Seleccionar en modo interactivo:
#   Forcefield → según configs/hmg_competition.yaml (charmm36)
#   Water      → tip3p

# 3. Verificar outputs
#   protein_1_processed.gro  → estructura procesada
#   topol.top               → topología
#   posre.itp               → restraints de posición (generado automáticamente)
