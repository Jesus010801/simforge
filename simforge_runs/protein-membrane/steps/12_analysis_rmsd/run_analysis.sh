#!/bin/bash
# ─── RMSD analysis ───────────────────────────────────────────────────────────
# Reference: Backbone (group 4)
# Selection: Backbone
# Input groups use GROMACS built-in defaults — no make_ndx required.

PROD_DIR="../11_production_md"

mkdir -p tables plots

# RMSD vs. initial structure (backbone reference + selection)
# Prompt: group 4 (Backbone) for least-squares fit, then selection for RMSD
echo "4 4" | gmx rms \
    -s "$PROD_DIR/md.tpr" \
    -f "$PROD_DIR/md.xtc" \
    -o tables/rmsd_backbone.xvg \
    -tu ns

# RMSD of selection vs. backbone reference
echo "4 1" | gmx rms \
    -s "$PROD_DIR/md.tpr" \
    -f "$PROD_DIR/md.xtc" \
    -o tables/rmsd_protein.xvg \
    -tu ns

echo "RMSD analysis complete → tables/rmsd_backbone.xvg  tables/rmsd_protein.xvg"