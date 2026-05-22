#!/bin/bash
# ─── Assembly: combinar proteína y ligandos ───────────────────────────────────
# Paths resueltos desde DAG — no editar manualmente

PROTEIN_1="../02_prepare_protein_1/protein_1_processed.gro"
SUBSTRATE_1="../05_parametrize_substrate_1/substrate_1.gro"
LIGAND_1="../04_parametrize_ligand_1/ligand_1.gro"

# Combinar estructuras
cat $PROTEIN_1 $SUBSTRATE_1 $LIGAND_1 > complex_raw.gro

# Actualizar número de átomos
python3 -c "
lines = open('complex_raw.gro').readlines()
n_atoms = sum(1 for l in lines[2:-1] if l.strip())
lines[1] = f'{n_atoms}\n'
open('complex.gro', 'w').writelines(lines)
print(f'Complex: {n_atoms} atoms')
"

# Combinar topologías
echo "Editar topol.top para incluir ligand ITP files:"
echo '#include "substrate_1.itp"'
echo '#include "ligand_1.itp"'
