#!/usr/bin/env python3
"""Phase 6 — windowed temporal analysis (0-50 / 50-100 / 100-150 / 150-200 ns).

For every scalar observable already computed, report per-window mean/sd plus a
stationarity assessment (OLS slope over the whole trace + its significance, and
window-to-window drift).  Convergence is NOT claimed from length alone.
Writes windowed_analysis_A1_A3_A6.csv (all 200 ns systems included).
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, read_xvg

WINDOWS = [(0, 50), (50, 100), (100, 150), (150, 200)]

# observable key -> (subdir, filename, unit, column)
OBS = [
 ("protein_rmsd",          "protein_rmsd",        "rmsd_protein.xvg",        "nm",   0),
 ("ligand_rmsd",           "ligand_rmsd",         None,                     "nm",   0),
 ("coa_rmsd",              "coa_rmsd",            "rmsd_coa.xvg",           "nm",   0),
 ("radius_of_gyration",    "radius_of_gyration",  "rg_protein.xvg",         "nm",   0),
 ("sasa",                  "sasa",                "sasa_protein.xvg",       "nm^2", 0),
 ("hbonds_protein",        "hbonds_protein",      "hbnum_protein.xvg",      "count",0),
 ("active_site_mindist",   "active_site",         None,                     "nm",   0),
 ("active_site_contacts",  "active_site",         None,                     "count",0),
 ("catalytic_com_distance","catalytic_com_distance", None,                  "nm",   0),
]


def find_xvg(system, subdir, fname):
    d = OUT / system / "observables" / subdir
    if not d.is_dir():
        return None
    if fname and (d / fname).is_file():
        return d / fname
    cands = sorted(d.glob("*.xvg"))
    if not fname:
        # pick by keyword
        key = {"active_site": "mindist", "ligand_rmsd": "rmsd", "coa_rmsd": "rmsd",
               "catalytic_com_distance": "dist"}.get(subdir)
        for c in cands:
            if subdir == "active_site":
                return None  # handled explicitly below
            if key and key in c.name:
                return c
    return cands[0] if cands else None


def windows_for(system, subdir, fname, col):
    # active_site has two files
    d = OUT / system / "observables" / subdir
    files = []
    if subdir == "active_site":
        files = [("active_site_mindist", next(d.glob("mindist_*active.xvg"), None)),
                 ("active_site_contacts", next(d.glob("contacts_*active.xvg"), None))]
    else:
        f = find_xvg(system, subdir, fname)
        files = [(subdir, f)]
    rows = []
    for label, f in files:
        if not f or not Path(f).is_file():
            continue
        t, y = read_xvg(f)
        v = y[:, col]
        # time to ns
        tn = t / 1000.0 if t.max() > 500 else t
        finite = np.isfinite(v)
        t2, v2 = tn[finite], v[finite]
        if v2.size < 4:
            continue
        # OLS slope over whole trace
        A = np.vstack([t2, np.ones_like(t2)]).T
        slope, intercept = np.linalg.lstsq(A, v2, rcond=None)[0]
        resid = v2 - (A @ [slope, intercept])
        se = np.sqrt((resid**2).sum() / (v2.size - 2)) / np.sqrt(((t2 - t2.mean())**2).sum())
        tstat = slope / se if se else 0.0
        row = dict(system=system, observable=label,
                   overall_mean=round(float(v2.mean()), 4), overall_sd=round(float(v2.std(ddof=1)), 4),
                   ols_slope_per_ns=round(float(slope), 6), slope_tstat=round(float(tstat), 2),
                   drift_over_200ns=round(float(slope * 200), 4),
                   stationary=("yes" if abs(tstat) < 3 else "DRIFT"))
        for w0, w1 in WINDOWS:
            m = (t2 >= w0) & (t2 < w1)
            seg = v2[m]
            row[f"mean_{w0}_{w1}"] = round(float(seg.mean()), 4) if seg.size else ""
            row[f"sd_{w0}_{w1}"] = round(float(seg.std(ddof=1)), 4) if seg.size > 1 else ""
        # window-to-window max shift
        wm = [row[f"mean_{a}_{b}"] for a, b in WINDOWS if row[f"mean_{a}_{b}"] != ""]
        row["max_window_shift"] = round(float(max(wm) - min(wm)), 4) if len(wm) >= 2 else ""
        rows.append(row)
    return rows


def main():
    systems = [s for s in SYS if s not in ()]  # all; 25 ns none here, all 200 ns
    all_rows = []
    for s in systems:
        seen_subdir = set()
        for key, subdir, fname, unit, col in OBS:
            if subdir in seen_subdir:      # active_site yields both metrics in one call
                continue
            seen_subdir.add(subdir)
            all_rows += windows_for(s, subdir, fname, col)
    if not all_rows:
        print("no observables found yet"); return
    cols = ["system", "observable", "overall_mean", "overall_sd",
            "mean_0_50", "sd_0_50", "mean_50_100", "sd_50_100",
            "mean_100_150", "sd_100_150", "mean_150_200", "sd_150_200",
            "max_window_shift", "ols_slope_per_ns", "slope_tstat",
            "drift_over_200ns", "stationary"]
    p = OUT / ".." / ".." / "reports" / "xanthone_mechanistic_2026-09-10" / "windowed_analysis_A1_A3_A6.csv"
    p = p.resolve()
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    (OUT / "_windowed_all.json").write_text(json.dumps(all_rows, indent=2))
    print(f"wrote {p}  ({len(all_rows)} rows)")
    for r in all_rows:
        if r["stationary"] == "DRIFT":
            print(f"  DRIFT: {r['system']}/{r['observable']}  slope {r['ols_slope_per_ns']}/ns "
                  f"(t={r['slope_tstat']}), Δ200ns={r['drift_over_200ns']}")


if __name__ == "__main__":
    main()
