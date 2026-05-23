#!/bin/bash
# ─── Distance analysis ────────────────────────────────────────────────────────
# Requires index.ndx — run make_ndx.sh first.
# Selection from config: {'group1': 'substrate_1', 'group2': 'ligand_1'}

PROD_DIR="../13_production_md"

mkdir -p tables plots

# Replace GROUP_A and GROUP_B with the group names from index.ndx
echo "GROUP_A GROUP_B" | gmx distance \
    -s "$PROD_DIR/md.tpr" \
    -f "$PROD_DIR/md.xtc" \
    -n index.ndx \
    -oav tables/distance_avg.xvg \
    -tu ns

echo "Distance analysis complete → tables/distance_avg.xvg"