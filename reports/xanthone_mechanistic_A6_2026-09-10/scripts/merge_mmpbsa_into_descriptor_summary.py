#!/usr/bin/env python3
"""Add 200 ns MM-PBSA columns to mechanistic_descriptor_summary.csv, sourced
from mmpbsa_200ns_summary.json (already reconciled, not rerun here).
Backs up the pre-merge CSV once.
"""
import json, shutil
from pathlib import Path
import pandas as pd

REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")
OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")
CSV = REP / "mechanistic_descriptor_summary.csv"
BAK = REP / "mechanistic_descriptor_summary.pre_mmpbsa200ns.csv.bak"

j = json.loads((OUT / "mmpbsa_200ns" / "mmpbsa_200ns_summary.json").read_text())

# system (as in the CSV's 'system' column) -> (channel-label-in-json, column-prefix)
MAP = {
    "A1":             [("A1", "mmpbsa_200ns")],
    "A3-HMG-R":       [("A3", "mmpbsa_200ns")],
    "HMG-R-200ns-A6": [("A6", "mmpbsa_200ns")],
    "COA":            [("substrate-only (E+S)", "mmpbsa_200ns")],
    "COMP":           [("A1 in ternary (E+S+I, I channel)", "lig_mmpbsa_200ns"),
                       ("CoA in ternary (E+S+I, S channel)", "coa_mmpbsa_200ns"),
                       ("A1<->CoA pairwise (I-S, E+S+I)", "lig_coa_mmpbsa_200ns")],
}
SOURCE_DIR = {"A1": "HMG_CoA_R-A1-R2", "A3": "A3-HMG-R (own mdfit.xtc)",
             "A6": "HMG-R-200ns-A6 (md_final.xtc)",
             "substrate-only (E+S)": "HMG-CoA-R_sustrato",
             "A1 in ternary (E+S+I, I channel)": "Competitive-system-2",
             "CoA in ternary (E+S+I, S channel)": "Competitive-system-2",
             "A1<->CoA pairwise (I-S, E+S+I)": "Competitive-system-2"}

if not BAK.is_file():
    shutil.copy2(CSV, BAK)
    print("backed up ->", BAK)

df = pd.read_csv(CSV)
new_cols = set()
for _, chans in MAP.items():
    for _, prefix in chans:
        new_cols.add(f"{prefix}_kjmol_mean"); new_cols.add(f"{prefix}_kjmol_sd")
        new_cols.add(f"{prefix}_source")
for c in sorted(new_cols):
    if c not in df.columns:
        df[c] = pd.Series([pd.NA] * len(df), dtype="object")

for sys_name, chans in MAP.items():
    mask = df["system"] == sys_name
    if not mask.any():
        print("WARNING: system not found in CSV:", sys_name); continue
    for label, prefix in chans:
        t = j[label]["energy_summary_kjmol"]["Total"]
        df.loc[mask, f"{prefix}_kjmol_mean"] = round(t["mean"], 3)
        df.loc[mask, f"{prefix}_kjmol_sd"] = round(t["sd"], 3)
        df.loc[mask, f"{prefix}_source"] = SOURCE_DIR[label]

df.to_csv(CSV, index=False)
print("updated", CSV)
print(df[["system", "label"] + sorted(new_cols)].to_string(index=False))
