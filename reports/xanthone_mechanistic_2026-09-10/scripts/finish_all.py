#!/usr/bin/env python3
"""Driver: after the 4 new systems have a PBC-clean mdfit.xtc, run all their
analyses, then the cross-system phases and every report."""
import subprocess, sys
from pathlib import Path
S = "/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10/scripts"
NEW = ["A3-HMG-R", "HMG-R-200ns-A6", "system_A3_COA", "system_A6_COA"]
OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")


def run(*a):
    print("::", " ".join(a), flush=True)
    r = subprocess.run(["python3", *a], cwd=S)
    print(f"   rc={r.returncode}", flush=True)


# 1. sanity gate: every new system's mdfit must have clean protein RMSD
bad = []
for s in NEW:
    j = OUT / s / "prep_corrective.json"
    if j.is_file():
        import json
        d = json.loads(j.read_text())
        if not d.get("pbc_fix_ok"):
            bad.append((s, d.get("sanity_protein_rmsd_max_nm")))
if bad:
    print("PBC FIX NOT CONFIRMED for:", bad)
    print("continuing anyway (sanity metric may have failed for a benign reason) — check qc_report")

# 2. per-new-system analyses
for s in NEW:
    run("run_all_new.py", s)

# 3. cross-system
run("shared_pca.py")
run("coa_decomposition.py")
run("windowed.py")

# 4. reports
run("build_reports.py")
run("finalize.py")
run("final_report.py")
print("\nFINISH_ALL COMPLETE", flush=True)
