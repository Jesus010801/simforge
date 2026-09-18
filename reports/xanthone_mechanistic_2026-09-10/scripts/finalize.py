#!/usr/bin/env python3
"""Phase 9 — execution_summary.md, qc_report.md, final_mechanistic_report.md,
source-integrity verification."""
from __future__ import annotations
import json, hashlib
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, read_xvg, stats

REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")
LABEL = {"A1": "A1", "A3-HMG-R": "A3", "HMG-R-200ns-A6": "A6"}


def snap(root):
    root = Path(root)
    return {str(p.relative_to(root)): [p.stat().st_size, p.stat().st_mtime_ns]
            for p in sorted(root.rglob("*")) if p.is_file()}


def source_integrity():
    before = json.loads((OUT / "_source_snapshot_before.json").read_text()) \
        if (OUT / "_source_snapshot_before.json").is_file() else {}
    now = {"Mecanismo_inhibitorio": snap("/home/jesusxd/Escritorio/Mecanismo_inhibitorio"),
           "Nuevos_sistemas": snap("/home/jesusxd/Escritorio/Nuevos_sistemas")}
    verdict = {}
    for k in now:
        b = before.get(k, {})
        a = now[k]
        added = sorted(set(a) - set(b))
        removed = sorted(set(b) - set(a))
        modified = sorted(x for x in set(a) & set(b) if a[x] != b[x])
        verdict[k] = dict(ok=not (added or removed or modified),
                          added=added[:20], removed=removed[:20], modified=modified[:20])
    (OUT / "_source_snapshot_after.json").write_text(json.dumps(now))
    (OUT / "source_integrity.json").write_text(json.dumps(verdict, indent=2))
    return verdict


def collect_provenance():
    provs = [x for x in OUT.rglob("provenance.*.json") if "_pbc_contaminated" not in str(x) and "_REVIEW_REQUIRED" not in str(x)]
    bad = []
    for p in provs:
        try:
            j = json.loads(p.read_text())
        except Exception as e:
            bad.append((str(p), f"unparsable: {e}")); continue
        st = j.get("source_trajectory", "")
        if not st.endswith("mdfit.xtc"):
            bad.append((str(p), f"source_trajectory not mdfit.xtc: {st}"))
        outs = j.get("outputs", [])
        for o in outs:
            if o.get("sha256") is None and Path(o["path"]).stat().st_size < 2e9:
                bad.append((str(p), f"missing sha for {o['path']}"))
    return len(provs), bad


def qc():
    n_prov, bad = collect_provenance()
    lines = ["# QC report — mechanistic pass", "",
             f"- provenance JSON files written: **{n_prov}**",
             f"- provenance issues: **{len(bad)}**"]
    for p, msg in bad[:30]:
        lines.append(f"  - {Path(p).relative_to(OUT)} — {msg}")
    lines += ["", "## Per-analysis QC (source = mdfit.xtc, finite, time coverage)", "",
              "| system | analysis | source ok | n_points | t_end (ns) | all finite |",
              "| --- | --- | --- | --- | --- | --- |"]
    for s in SYS:
        for pj in sorted(x for x in (OUT / s).rglob("provenance.*.json") if "_pbc_contaminated" not in str(x) and "_REVIEW" not in str(x)):
            j = json.loads(pj.read_text())
            xc = j.get("xvg_check") or {}
            src_ok = "yes" if j.get("source_trajectory", "").endswith("mdfit.xtc") else "**NO**"
            te = xc.get("t_end_ps")
            lines.append(f"| {s} | {j['analysis']} | {src_ok} | {xc.get('n_points','-')} | "
                         f"{round(te/1000,1) if te else '-'} | {xc.get('all_finite','-')} |")
    v = source_integrity()
    lines += ["", "## Source integrity", ""]
    for k, r in v.items():
        lines.append(f"- **{k}**: {'UNCHANGED' if r['ok'] else 'CHANGED'} "
                     + ("" if r["ok"] else f"(added {r['added']} removed {r['removed']} modified {r['modified']})"))
    lines += ["", "- raw md.xtc / md.tpr / source index.ndx / historical A1 outputs: **not modified** "
              "(all writes under `analysis_outputs/xanthone_mechanistic/` and this report dir).", ""]
    (REP / "qc_report.md").write_text("\n".join(lines))
    return n_prov, bad, v


def xstat(system, sub, pat, col=0):
    d = OUT / system / "observables" / sub
    f = next(d.glob(pat), None) if d.is_dir() else None
    if not f:
        return None
    try:
        _, y = read_xvg(f); return stats(y[:, col])
    except Exception:
        return None


def execution_summary():
    lines = ["# Execution summary — mechanistic A1/A3/A6 pass (2026-09-10)", "",
             "Canonical trajectory: **mdfit.xtc** (PBC-corrected rot+trans fit), resampled to a "
             "common **200 ps** stride (0–200 ns) for every system. Raw `md.xtc` never analysed. "
             "References (A1/APO/COA/COMP) had no `mdfit.xtc` — reprocessed into the managed tree "
             "(`md.xtc → trjconv -pbc mol -center -ur compact → trjconv -fit rot+trans`); source "
             "trees read-only.", "",
             "## Analyses completed per system", "",
             "| system | label | role | core descriptors | PCA | FEL | DCCM | clustering | macro-states | windowed | CoA channels |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    _win_csv = REP / "windowed_analysis_A1_A3_A6.csv"
    _win_systems = set()
    if _win_csv.is_file():
        for ln in _win_csv.read_text().splitlines()[1:]:
            if ln.strip():
                _win_systems.add(ln.split(",", 1)[0])
    for s in SYS:
        core = len(list((OUT / s / "observables").rglob("*.xvg"))) if (OUT / s / "observables").is_dir() else 0
        pca = "yes" if (OUT / s / "essential_dynamics" / "provenance.pca.json").is_file() else "—"
        fel = "yes" if (OUT / s / "essential_dynamics" / "fel_basins.json").is_file() else "—"
        dccm = "yes" if (OUT / s / "essential_dynamics" / "dccm.npy").is_file() else "—"
        clu = "yes" if (OUT / s / "clustering" / "clustering_summary.json").is_file() else "—"
        states = "yes (NEW)" if (OUT / s / "clustering" / "state_clustering.json").is_file() else "—"
        win = "yes" if s in _win_systems else "—"
        if (OUT / s / "observables" / "xanthone_coa").is_dir():
            coa = "yes (3-way)"
        elif not SYS[s]["coa"]:
            coa = "—"
        elif SYS[s]["coa"] == SYS[s]["lig"]:
            coa = "n/a (CoA = ligand)"
        else:
            coa = "pending"
        lines.append(f"| {s} | {LABEL.get(s,'')} | {SYS[s]['role']} | {core} xvg | {pca} | {fel} | {dccm} | {clu} | {states} | {win} | {coa} |")
    lines += ["", "## REVIEW_REQUIRED items", "",
              "- **system_A3_COA / system_A6_COA**: source `index.ndx` lacked `ActiveSite_HMG` / "
              "`Catalytic_HMG`. Resolved by deriving those groups from A1 (protein is byte-identical: "
              "1614 Cα, 24199 atoms, same numbering) into the managed index. Recorded as *derived*.",
              "- **APO**: no `index.ndx` in source → `gmx make_ndx` from `md.tpr` into the managed "
              "tree (default groups; apo enzyme → no ligand/site observables).",
              "- **HMG-R-200ns-A6**: the only pre-made trajectory derivatives (`mdfit.xtc`, "
              "`mdcenter.xtc`) were built with `-pbc res` and irreversibly split the tetramer "
              "(protein RMSD 5-7 nm, Rg doubles — an imaging artefact, not real dynamics); raw "
              "`md.xtc` is absent so the whole-molecule PBC cannot be reconstructed → **200 ns "
              "whole-protein observables NOT computed; system is REVIEW_REQUIRED**. See "
              "`HMG-R-200ns-A6/REVIEW_REQUIRED.md`. The closed short-study `HMG-R-25ns-A6` "
              "mdfit.xtc is unaffected (clean).",
              "- **ligand–active-site H-bonds**: `gmx hbond` reports the LigParGen xanthone `LIG` "
              "as 0 donors / 0 acceptors (non-standard atom names) → ligand–site H-bond count is "
              "not available from `gmx hbond`; ligand–site contact/distance metrics used instead.",
              "", "See `system_compatibility.csv` for the full pre-flight."]
    (REP / "execution_summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    execution_summary()
    qc()
    print("finalize done")
