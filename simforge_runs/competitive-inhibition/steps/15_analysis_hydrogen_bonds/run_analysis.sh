#!/bin/bash
# ─── Hydrogen bond analysis ───────────────────────────────────────────────────
# Protein intra-molecular H-bonds (group 1 → group 1)
# Input groups use GROMACS built-in defaults — no make_ndx required.

PROD_DIR="../13_production_md"

mkdir -p tables plots

# Intra-protein H-bonds
echo "1 1" | gmx hbond \
    -s "$PROD_DIR/md.tpr" \
    -f "$PROD_DIR/md.xtc" \
    -num tables/hbnum_protein.xvg \
    -dist tables/hbdist_protein.xvg

# Protein–solvent H-bonds
echo "1 13" | gmx hbond \
    -s "$PROD_DIR/md.tpr" \
    -f "$PROD_DIR/md.xtc" \
    -num tables/hbnum_protein_sol.xvg

echo "H-bond analysis complete → tables/hbnum_protein.xvg  tables/hbnum_protein_sol.xvg"