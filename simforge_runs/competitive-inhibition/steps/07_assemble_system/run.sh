#!/bin/bash
# ─── Assembly: combinar proteína y ligandos ───────────────────────────────────

# Combinar GRO de proteína + ligandos
# Ajustar según número de ligandos del sistema

# Proteína procesada (del step de preparation)
PROTEIN="../../01_prepare_protein_1/protein_1_processed.gro"

# Ligandos parametrizados
LIGAND_1="../../03_parametrize_substrate_1/substrate_1.gro"
LIGAND_2="../../04_parametrize_ligand_1/ligand_1.gro"

# Combinar estructuras
cat $PROTEIN $LIGAND_1 $LIGAND_2 > complex_raw.gro

# Actualizar número de átomos en la primera línea (suma total)
# Editar manualmente o usar script de python:
python3 -c "
import sys
lines = open('complex_raw.gro').readlines()
n_atoms = sum(1 for l in lines[2:-1] if l.strip())
lines[1] = f'{n_atoms}\n'
open('complex.gro', 'w').writelines(lines)
print(f'Complex: {n_atoms} atoms')
"

# Combinar topologías (editar topol.top manualmente para incluir ligandos)
echo "Editar topol.top para incluir ligand ITP files"
echo "Agregar al final de topol.top:"
echo '; Ligandos'
echo '#include "substrate_1.itp"'
echo '#include "ligand_1.itp"'
