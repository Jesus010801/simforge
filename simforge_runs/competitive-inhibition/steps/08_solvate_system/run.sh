#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────
# water_model=tip3p  box_type=triclinic  d=1.2nm
# Paths resueltos desde DAG

ASSEMBLE_DIR="../07_assemble_system"

# Copiar topología — gmx solvate la modifica in-place (añade SOL).
# La copia local asegura que assemble_system/topol.top quede intacta
# y que los pasos siguientes lean la topología actualizada desde aquí.
cp "$ASSEMBLE_DIR/topol.top" topol.top

# Definir caja de simulación
gmx editconf \
    -f "$ASSEMBLE_DIR/complex.gro" \
    -o box.gro \
    -c \
    -d 1.2 \
    -bt triclinic

# Agregar agua
gmx solvate \
    -cp box.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p topol.top
