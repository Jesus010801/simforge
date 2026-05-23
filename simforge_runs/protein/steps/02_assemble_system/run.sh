#!/bin/bash
# ─── Assembly: combinar proteína y ligandos ───────────────────────────────────
# Paths resueltos desde DAG — no editar manualmente

PROTEIN_1="../01_prepare_protein_1/protein_1_processed.gro"

# Handoff topology from protein prep
PROTEIN_PREP_DIR="../01_prepare_protein_1"
cp "$PROTEIN_PREP_DIR/topol.top" topol.top
cp "$PROTEIN_PREP_DIR/posre.itp" posre.itp

# Combinar estructuras
cat $PROTEIN_1 > complex_raw.gro

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

