#!/usr/bin/env python3
"""Phase 8 — assemble every report from the managed outputs.

Reads analysis_outputs/xanthone_mechanistic/**  (xvg + provenance + summary json)
and writes reports/xanthone_mechanistic_2026-09-10/*.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, read_xvg, stats

REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")
SHARED = OUT / "shared_pca"
LABEL = {"A1": "A1", "A3-HMG-R": "A3", "HMG-R-200ns-A6": "A6"}


def xvg_stats(system, subdir, pat, col=0):
    d = OUT / system / "observables" / subdir
    if not d.is_dir():
        return None
    f = next(d.glob(pat), None)
    if not f:
        return None
    try:
        _, y = read_xvg(f)
        return stats(y[:, col])
    except Exception:
        return None


# ── mechanistic_descriptor_summary.csv ──────────────────────────────────────
def descriptor_summary():
    OBS = [
     ("protein_rmsd_nm",        "protein_rmsd",        "rmsd_protein.xvg"),
     ("ligand_rmsd_nm",         "ligand_rmsd",         "rmsd_*.xvg"),
     ("coa_rmsd_nm",            "coa_rmsd",            "rmsd_coa.xvg"),
     ("rmsf_calpha_mean_nm",    "rmsf",                "rmsf_calpha.xvg"),
     ("rg_nm",                  "radius_of_gyration",  "rg_protein.xvg"),
     ("sasa_nm2",               "sasa",                "sasa_protein.xvg"),
     ("hbonds_protein_count",   "hbonds_protein",      "hbnum_protein.xvg"),
     ("hbonds_lig_site_count",  "hbonds_ligand_site",  "hbnum_*_site.xvg"),
     ("active_site_mindist_nm", "active_site",         "mindist_*active.xvg"),
     ("active_site_contacts",   "active_site",         "contacts_*active.xvg"),
     ("catalytic_com_dist_nm",  "catalytic_com_distance", "dist_*catalytic.xvg"),
     ("coa_active_mindist_nm",  "coa_active_site",     "mindist_coa_active.xvg"),
     ("coa_catalytic_com_nm",   "coa_catalytic_com_distance", "dist_coa_catalytic.xvg"),
    ]
    rows = []
    for s in SYS:
        row = {"system": s, "label": LABEL.get(s, ""), "role": SYS[s]["role"]}
        for name, sub, pat in OBS:
            col = 0  # gmx distance -oall writes "time  value" -> read_xvg y[:,0] is the value
            st = xvg_stats(s, sub, pat, col)
            row[f"{name}_mean"] = round(st["mean"], 4) if st and st.get("n") else ""
            row[f"{name}_sd"] = round(st["sd"], 4) if st and st.get("n") else ""
        rows.append(row)
    if not rows:
        return
    cols = ["system", "label", "role"] + [k for k in rows[0] if k not in ("system", "label", "role")]
    with open(REP / "mechanistic_descriptor_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


# ── clustering_state_populations.csv ────────────────────────────────────────
def clustering_csv():
    rows = []
    for s in SYS:
        f = OUT / s / "clustering" / "clustering_summary.json"
        if not f.is_file():
            continue
        j = json.loads(f.read_text())
        for cid, pop in (j.get("populations_overall") or {}).items():
            row = dict(system=s, label=LABEL.get(s, ""), cluster=cid,
                       frames=pop["frames"], fraction_overall=pop["fraction"],
                       middle_frame=(j.get("cluster_middle_frames") or {}).get(str(cid), ""))
            for w, wp in (j.get("populations_by_window") or {}).items():
                row[w] = wp.get(str(cid), wp.get(int(cid), 0.0))
            rows.append(row)
    if not rows:
        return
    _worder = ["0-50ns", "50-100ns", "100-150ns", "150-200ns"]
    wcols = [w for w in _worder if any(w in r for r in rows)]
    cols = ["system", "label", "cluster", "frames", "fraction_overall", "middle_frame"] + wcols
    with open(REP / "clustering_state_populations.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


# ── shared_pca_summary.md ───────────────────────────────────────────────────
def shared_pca_md():
    f = SHARED / "shared_pca_summary.json"
    if not f.is_file():
        (REP / "shared_pca_summary.md").write_text("# Shared PCA — NOT YET RUN\n")
        return
    j = json.loads(f.read_text())
    v = j["shared_variance_explained"]
    lines = ["# Shared PCA basis — A1 / A3 / A6", "",
             j["method"], "",
             f"- Cα per system: **{j['n_calpha']}** (identical topology / residue numbering)",
             f"- frames: {j['frames']}",
             f"- kT = {j['kT_kjmol']} kJ/mol (309.65 K)", "",
             "## Variance explained (shared eigenbasis)", "",
             "| PC | fraction | cumulative |", "| --- | --- | --- |"]
    cum = 0.0
    for i, x in enumerate(v[:6], 1):
        cum += x
        lines.append(f"| PC{i} | {x:.3f} | {cum:.3f} |")
    lines += ["", f"PC1+PC2 = **{j['cum_variance_pc1_2']:.1%}**, PC1–PC5 = {j['cum_variance_pc1_5']:.1%}", "",
              "## PC1/PC2 projections (nm)", "",
              "| system | PC1 centroid | PC2 centroid | PC1 range | PC2 range |",
              "| --- | --- | --- | --- | --- |"]
    for k in j["pc_centroids_nm"]:
        c = j["pc_centroids_nm"][k]
        r = j["pc_ranges_nm"][k]
        lines.append(f"| {k} | {c[0]} | {c[1]} | {r['pc1']} | {r['pc2']} |")
    lines += ["", "## Overlap / separation (shared plane)", "",
              "| pair | Bhattacharyya overlap | centroid separation (nm) |",
              "| --- | --- | --- |"]
    for pair, o in j["pairwise_overlap"].items():
        lines.append(f"| {pair} | {o['bhattacharyya_overlap']} | {o['centroid_separation_nm']} |")
    lines += ["", "## Dominant conformational regions (shared-basis FEL basins, ΔG < 2.5 kJ/mol)", ""]
    for k in j["pc_centroids_nm"]:
        bs = j["shared_fel_basins"].get(k, [])
        lines.append(f"- **{k}**: {len(bs)} basin(s) — "
                     + "; ".join(f"(PC1 {b['pc1']}, PC2 {b['pc2']}) ΔG {b['G']}" for b in bs[:4]))
    lines += ["", "Figures: `analysis_outputs/xanthone_mechanistic/shared_pca/shared_pc12.png`, "
              "`shared_fel.png`.", ""]
    (REP / "shared_pca_summary.md").write_text("\n".join(lines))


# ── fel_state_summary.md ────────────────────────────────────────────────────
def fel_state_md():
    lines = ["# FEL state summary (per-system PCA + shared-basis)", ""]
    for s in SYS:
        f = OUT / s / "essential_dynamics" / "fel_basins.json"
        pca = OUT / s / "essential_dynamics" / "provenance.pca.json"
        if not f.is_file():
            continue
        j = json.loads(f.read_text())
        var = json.loads(pca.read_text()).get("variance_explained_pc1_5", []) if pca.is_file() else []
        lines.append(f"## {s}  ({LABEL.get(s,'reference')})")
        lines.append(f"- per-system PCA variance PC1–PC5: {var}")
        lines.append(f"- FEL basins (ΔG < 3 kJ/mol): **{j['n_basins']}**")
        for b in j["basins"][:5]:
            lines.append(f"  - PC1 {b['pc1']:.2f}, PC2 {b['pc2']:.2f} nm — ΔG {b['G']:.2f} kJ/mol")
        lines.append("")
    sp = SHARED / "shared_pca_summary.json"
    if sp.is_file():
        j = json.loads(sp.read_text())
        lines += ["## Shared-basis FEL (A1/A3/A6, same PC1/PC2 axes)", ""]
        for k in j["pc_centroids_nm"]:
            bs = j["shared_fel_basins"].get(k, [])
            lines.append(f"- **{k}**: {len(bs)} basin(s); global min region "
                         + (f"PC1 {bs[0]['pc1']}, PC2 {bs[0]['pc2']}" if bs else "n/a"))
    (REP / "fel_state_summary.md").write_text("\n".join(lines))


# ── dccm_summary.md ─────────────────────────────────────────────────────────
def dccm_md():
    lines = ["# DCCM summary (Cα–Cα dynamic cross-correlation)", "",
             "Normalised covariance on aligned Cα; matrix `.npy` + `.png` per system in "
             "`analysis_outputs/xanthone_mechanistic/<sys>/essential_dynamics/`.", "",
             "| system | mean |C| | frac strongly anti-corr (< -0.4) | frac strongly corr (> 0.6, off-diag) |",
             "| --- | --- | --- | --- |"]
    mats = {}
    for s in SYS:
        f = OUT / s / "essential_dynamics" / "dccm.npy"
        if not f.is_file():
            continue
        C = np.load(f)
        mats[s] = C
        off = ~np.eye(C.shape[0], dtype=bool)
        lines.append(f"| {s} | {np.abs(C[off]).mean():.3f} | "
                     f"{(C[off] < -0.4).mean():.3f} | {(C[off] > 0.6).mean():.3f} |")
    # A1/A3/A6 difference maps (A6 only if it has a clean DCCM)
    for a, b in [("A3-HMG-R", "A1"), ("HMG-R-200ns-A6", "A1"), ("HMG-R-200ns-A6", "A3-HMG-R")]:
        if a in mats and b in mats:
            D = mats[a] - mats[b]
            np.save(OUT / "shared_pca" / f"dccm_diff_{LABEL[a]}_minus_{LABEL[b]}.npy", D)
            lines.append("")
            lines.append(f"- **ΔDCCM {LABEL[a]} − {LABEL[b]}**: mean|Δ| {np.abs(D).mean():.3f}, "
                         f"max +{D.max():.2f} / {D.min():.2f}; "
                         f"frac |Δ|>0.3 = {(np.abs(D) > 0.3).mean():.3f}")
    (REP / "dccm_summary.md").write_text("\n".join(lines))


# ── coa_interaction_summary.md ──────────────────────────────────────────────
def coa_md():
    lines = ["# CoA-system interaction decomposition", "",
             "Three channels kept separate for COMP / system_A3_COA / system_A6_COA.", ""]
    for s in ["COMP", "system_A3_COA", "system_A6_COA"]:
        f = OUT / s / "observables" / "xanthone_coa" / "coa_channels_summary.json"
        if not f.is_file():
            lines.append(f"## {s} — not yet run\n"); continue
        j = json.loads(f.read_text())
        ch = j["channels"]
        lines.append(f"## {s}")
        def fmt(d, key):
            v = d.get(key) if isinstance(d, dict) else None
            if isinstance(v, dict) and v.get("n"):
                return f"{v['mean']:.3f} ± {v['sd']:.3f}"
            return "n/a"
        lines.append(f"- protein–xanthone: active-site mindist {fmt(ch['protein_xanthone'],'active_site_mindist_nm')} nm, "
                     f"contacts {fmt(ch['protein_xanthone'],'active_site_contacts')}, "
                     f"catalytic COM {fmt(ch['protein_xanthone'],'catalytic_com_distance_nm')} nm")
        pc = ch.get("protein_coa")
        if isinstance(pc, dict):
            lines.append(f"- protein–CoA: active-site mindist {fmt(pc,'active_site_mindist_nm')} nm, "
                         f"contacts {fmt(pc,'active_site_contacts')}, "
                         f"catalytic COM {fmt(pc,'catalytic_com_distance_nm')} nm")
        xc = ch["xanthone_coa"]
        lines.append(f"- xanthone–CoA: mindist {fmt(xc,'mindist_nm')} nm, "
                     f"contacts(<0.6 nm) {fmt(xc,'contacts_lt_0p6nm')}, "
                     f"COM distance {fmt(xc,'com_distance_nm')} nm")
        lines.append("")
    (REP / "coa_interaction_summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    descriptor_summary()
    clustering_csv()
    shared_pca_md()
    fel_state_md()
    dccm_md()
    coa_md()
    print("reports assembled")
