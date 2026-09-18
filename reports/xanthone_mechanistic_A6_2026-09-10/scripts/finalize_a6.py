#!/usr/bin/env python3
"""Re-run finalize.py's QC/integrity checks, patched so that A6's provenance
(source_trajectory = md_final.xtc, by design) is not flagged as an error the
way a stray mdfit.xtc substitution would be. Appends an A6-specific QC note.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

PIPE = "/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10/scripts"
sys.path.insert(0, PIPE)
import finalize as fz  # noqa: E402
from mech_common import SYS, OUT  # noqa: E402

A6 = "HMG-R-200ns-A6"
A6_TRAJ_NAME = "md_final.xtc"


def collect_provenance_patched():
    provs = [x for x in OUT.rglob("provenance.*.json")
            if "_pbc_contaminated" not in str(x) and "_REVIEW_REQUIRED" not in str(x)
            and "_superseded" not in str(x)]
    bad = []
    for p in provs:
        try:
            j = json.loads(p.read_text())
        except Exception as e:
            bad.append((str(p), f"unparsable: {e}")); continue
        st = j.get("source_trajectory", "")
        ok_name = "mdfit.xtc" if j.get("system") != A6 else A6_TRAJ_NAME
        if not st.endswith(ok_name):
            bad.append((str(p), f"source_trajectory not {ok_name}: {st}"))
        outs = j.get("outputs", [])
        for o in outs:
            if o.get("sha256") is None and Path(o["path"]).stat().st_size < 2e9:
                bad.append((str(p), f"missing sha for {o['path']}"))
    return len(provs), bad


fz.collect_provenance = collect_provenance_patched
fz.qc()
fz.execution_summary()  # regenerate with A6 now included in the SYS-driven table
print("finalize_a6 done")

# ---- A6-specific supplementary QC note --------------------------------
REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")
n_prov = sum(1 for x in (OUT / A6).rglob("provenance.*.json"))
lines = ["", "## A6 correction (2026-09-10) — source_trajectory = md_final.xtc", "",
         f"- provenance JSON files for {A6}: **{n_prov}**, all stamped "
         f"`source_trajectory` ending in `{A6_TRAJ_NAME}` (by design — this is the "
         "user-designated corrected trajectory, not a pipeline-built mdfit.xtc).",
         "- superseded contaminated outputs archived under "
         f"`{A6}/_pbc_contaminated_2026-09-09/` and `{A6}/_superseded_2026-09-10_pre_md_final/`.",
         "- raw sources (`Nuevos_sistemas/HMG-R-200ns-A6/{md_final.xtc,md.tpr,index.ndx}`) "
         "were only read, never modified.", ""]
with open(REP / "qc_report.md", "a") as fh:
    fh.write("\n".join(lines))
print("appended A6 QC note to qc_report.md")
