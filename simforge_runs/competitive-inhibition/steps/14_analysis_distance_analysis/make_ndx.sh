#!/bin/bash
# ─── Custom index groups for distance analysis ────────────────────────────────
# Edit the 'name' and atom-selection expressions below to match your system.
# Run this BEFORE run_analysis.sh.
# GROMACS default groups (pdb2gmx standard output):
#   0  System          — all atoms
#   1  Protein         — protein atoms
#   2  Protein-H       — protein heavy atoms
#   3  C-alpha         — Cα atoms only
#   4  Backbone        — N, Cα, C
#   5  MainChain       — backbone + Cβ
#   6  MainChain+Cb    — MainChain + Cβ
#   7  MainChain+H     — MainChain + backbone H
#   8  SideChain       — side chain atoms
#   9  SideChain-H     — side chain heavy atoms
#  10  Prot-Masses     — protein with masses
#  11  non-Protein     — water + ions + ligands
#  12  Other           — non-protein non-water
#  13  SOL             — water molecules
#  14  non-Water       — everything except water
#  15+ NA, CL ...      — individual ion species (if present)

# Selection from config: {'group1': 'substrate_1', 'group2': 'ligand_1'}

PROD_DIR="../13_production_md"

# Interactive make_ndx — type selections at the prompt, then 'q' to save.
# Example selections:
#   "r 100 & a CA"    — Cα of residue 100
#   "r 200 & a CA"    — Cα of residue 200
# Then name them:
#   "name 16 ResA_CA"
#   "name 17 ResB_CA"
#   "q"
gmx make_ndx \
    -f "$PROD_DIR/md.tpr" \
    -o index.ndx

echo "Index file created → index.ndx"
echo "Edit run_analysis.sh with the correct group numbers before running."