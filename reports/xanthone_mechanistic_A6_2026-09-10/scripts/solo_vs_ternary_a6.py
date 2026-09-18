#!/usr/bin/env python3
"""Single-ligand A6 (HMG-R-200ns-A6, md_final.xtc) vs ternary A6+CoA
(system_A6_COA) comparison, now that A6-solo has real 200 ns data instead of
the 25 ns proxy used in the 2026-09-09 report.

Also assembles the equivalent A1-solo vs COMP (A1+CoA) and A3-solo vs
system_A3_COA rows for a consistent 3-system panel.
"""
from __future__ import annotations
import json
from pathlib import Path

OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")
REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")

PAIRS = [("A1", "A1", "COMP"), ("A3", "A3-HMG-R", "system_A3_COA"),
        ("A6", "HMG-R-200ns-A6", "system_A6_COA")]


def load_json(p):
    return json.loads(Path(p).read_text()) if Path(p).is_file() else None


def clust(sys_name):
    j = load_json(OUT / sys_name / "clustering" / "clustering_summary.json")
    if not j:
        return None, None
    pops = sorted(j["populations_overall"].items(), key=lambda kv: -kv[1]["fraction"])
    top2 = [round(v["fraction"] * 100, 1) for _, v in pops[:2]]
    return j["n_clusters"], top2


def descr(sys_name, key, sub, pat, col=0):
    import sys as _sys
    _sys.path.insert(0, "/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10/scripts")
    from mech_common import read_xvg, stats
    d = OUT / sys_name / "observables" / sub
    f = next(d.glob(pat), None) if d.is_dir() else None
    if not f:
        return None
    _, y = read_xvg(f)
    return stats(y[:, col])


def main():
    lines = ["# Single-ligand vs ternary (E+I+CoA) — A1 / A3 / A6 (2026-09-10 update)", "",
             "A6-solo now has full 200 ns data (`HMG-R-200ns-A6`, `md_final.xtc`), replacing "
             "the 25 ns proxy used previously. Same clustering method (gmx cluster gromos, "
             "Cα/Backbone RMSD, 0.15 nm cutoff) for solo and ternary.", "",
             "| label | solo: catalytic COM (nm) | ternary: xanthone catalytic COM (nm) | "
             "Δ (ternary−solo) | solo: active-site contacts | ternary: active-site contacts | "
             "solo: clusters (top state %) | ternary: clusters (top state %) |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    rows = []
    for label, solo, tern in PAIRS:
        cd_solo = descr(solo, "catalytic", "catalytic_com_distance", "dist_*catalytic.xvg")
        cd_tern = descr(tern, "catalytic", "catalytic_com_distance", "dist_*catalytic.xvg")
        ac_solo = descr(solo, "as", "active_site", "contacts_*active.xvg")
        ac_tern = descr(tern, "as", "active_site", "contacts_*active.xvg")
        nc_solo, top_solo = clust(solo)
        nc_tern, top_tern = clust(tern)
        d = (cd_tern["mean"] - cd_solo["mean"]) if (cd_tern and cd_solo) else None
        rows.append(dict(label=label, solo=solo, ternary=tern,
                         cd_solo=cd_solo, cd_tern=cd_tern, delta=d,
                         ac_solo=ac_solo, ac_tern=ac_tern,
                         nc_solo=nc_solo, top_solo=top_solo, nc_tern=nc_tern, top_tern=top_tern))
        def f(x, d=3):
            return f"{x['mean']:.{d}f}±{x['sd']:.{d}f}" if x and x.get("n") else "n/a"
        lines.append(f"| {label} | {f(cd_solo)} | {f(cd_tern)} | "
                     f"{f'{d:+.3f}' if d is not None else 'n/a'} | "
                     f"{f(ac_solo,0)} | {f(ac_tern,0)} | "
                     f"{nc_solo} ({'/'.join(map(str,top_solo)) if top_solo else '-'}%) | "
                     f"{nc_tern} ({'/'.join(map(str,top_tern)) if top_tern else '-'}%) |")
    (REP / "solo_vs_ternary_A1_A3_A6.md").write_text("\n".join(lines))
    (OUT / "solo_vs_ternary_A1_A3_A6.json").write_text(json.dumps(rows, indent=2, default=str))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
