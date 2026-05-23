#!/bin/bash
# ─── Preparación de proteína: peptide_1 ───────────────────────────────────────

gmx pdb2gmx \
    -f peptide_1.pdb \
    -o peptide_1_processed.gro \
    -p topol.top \
    -ff oplsaa_membrane \
    -water spce \
    -ignh
