#!/bin/bash
# ─── Preparación de proteína: protein_1 ───────────────────────────────────────

gmx pdb2gmx \
    -f protein_1.pdb \
    -o protein_1_processed.gro \
    -p topol.top \
    -ff oplsaa \
    -water spce \
    -ignh
