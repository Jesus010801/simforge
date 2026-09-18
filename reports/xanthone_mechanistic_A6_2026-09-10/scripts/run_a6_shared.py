#!/usr/bin/env python3
"""Phase 2 of the A6 correction: once run_a6_core.py has finished (core
descriptors + PCA/FEL/DCCM + clustering + macro-states for HMG-R-200ns-A6 on
md_final.xtc), rebuild every cross-system product that includes A6:

  1. windowed.py        (unmodified — generic over SYS, now picks up A6's xvgs)
  2. shared_pca_a6.py    (adds A6 to the shared A1/A3 PCA+FEL basis)
  3. build_reports.py   (unmodified — generic; regenerates descriptor_summary,
                          clustering_state_populations, shared_pca_summary.md,
                          fel_state_summary.md, dccm_summary.md [now includes
                          A6 vs A1 / A6 vs A3 diffs], coa_interaction_summary.md)
  4. finalize_a6.py      (QC + execution_summary.md, A6-aware)
"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).parent
PIPE = "/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10/scripts"


def run(script, cwd):
    print(f"\n===== {script} (cwd={cwd}) =====", flush=True)
    r = subprocess.run([sys.executable, script], cwd=cwd)
    print(f"   rc={r.returncode}", flush=True)
    if r.returncode != 0:
        raise SystemExit(f"{script} failed rc={r.returncode}")


run("windowed.py", PIPE)
run(str(HERE / "shared_pca_a6.py"), str(HERE))
run("build_reports.py", PIPE)
run(str(HERE / "finalize_a6.py"), str(HERE))
print("\nRUN_A6_SHARED COMPLETE", flush=True)
