#!/usr/bin/env python3
"""Phase 7 — CoA-system interaction decomposition (COMP, system_A3_COA, system_A6_COA).

Keep three interaction channels DISTINCT (never a single generic metric):
  * protein-xanthone : LIG-ActiveSite mindist/contacts, LIG-Catalytic COM dist
  * protein-CoA      : COA-ActiveSite mindist/contacts, COA-Catalytic COM dist
  * xanthone-CoA     : LIG-COA mindist, LIG-COA contacts (<0.6 nm), LIG-COA COM dist

The first two channels are already produced by core_descriptors for these
systems; this script adds the xanthone-CoA channel and assembles the summary.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, gmx, mdfit, tpr, ndx, read_xvg, stats, provenance

COA_SYS = ["COMP", "system_A3_COA", "system_A6_COA"]
CUT = "0.6"


def d_for(system):
    d = OUT / system / "observables" / "xanthone_coa"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_system(system):
    if system not in COA_SYS:
        return
    if not Path(mdfit(system)).is_file():
        print(f"  {system}: mdfit not ready"); return
    T, N, fit = tpr(system), ndx(system), mdfit(system)
    d = d_for(system)
    md = d / "mindist_lig_coa.xvg"; cn = d / "contacts_lig_coa.xvg"; cm = d / "dist_lig_coa_com.xvg"
    if not md.is_file():
        gmx(["mindist", "-s", T, "-f", fit, "-n", N, "-d", CUT, "-od", md, "-on", cn],
            stdin="LIG\nCOA\n", log=d / "mindist.log")
    if not cm.is_file():
        gmx(["distance", "-s", T, "-f", fit, "-n", N,
             "-select", 'com of group "LIG" plus com of group "COA"', "-oall", cm],
            stdin="", log=d / "distance.log")
    for a, argv, stdin, outs in [
        ("xanthone_coa_mindist", ["gmx", "mindist", "LIG", "COA"], "LIG\nCOA\n", [md]),
        ("xanthone_coa_contacts", ["gmx", "mindist", "LIG", "COA", "-on"], "LIG\nCOA\n", [cn]),
        ("xanthone_coa_com_distance", ["gmx", "distance", "com LIG + com COA"], "", [cm]),
    ]:
        prov = provenance(system, a, argv=argv, stdin=stdin, outputs=outs)
        (d / f"provenance.{a}.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")

    def s(sub, pat, col=0):
        dd = OUT / system / "observables" / sub
        f = next(dd.glob(pat), None) if dd.is_dir() else None
        if not f:
            return None
        _, y = read_xvg(f)
        return stats(y[:, col])

    summary = dict(
        system=system,
        channels=dict(
            protein_xanthone=dict(
                active_site_mindist_nm=s("active_site", "mindist_*active.xvg"),
                active_site_contacts=s("active_site", "contacts_*active.xvg"),
                catalytic_com_distance_nm=s("catalytic_com_distance", "dist_*catalytic.xvg"),
            ),
            protein_coa=dict(
                active_site_mindist_nm=s("coa_active_site", "mindist_coa_active.xvg"),
                active_site_contacts=s("coa_active_site", "contacts_coa_active.xvg"),
                catalytic_com_distance_nm=s("coa_catalytic_com_distance", "dist_coa_catalytic.xvg"),
            ) if system != "COMP" or (OUT / system / "observables" / "coa_active_site").is_dir() else "core did not run CoA channel — see note",
            xanthone_coa=dict(
                mindist_nm=stats(read_xvg(md)[1][:, 0]) if md.is_file() else None,
                contacts_lt_0p6nm=stats(read_xvg(cn)[1][:, 0]) if cn.is_file() else None,
                com_distance_nm=stats(read_xvg(cm)[1][:, 0]) if cm.is_file() else None,
            ),
        ),
        interpretation_note="three channels kept separate; a low xanthone-CoA "
                            "mindist with sustained contacts indicates direct "
                            "inhibitor-substrate contact in the pocket.",
    )
    (d / "coa_channels_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    xc = summary["channels"]["xanthone_coa"]["mindist_nm"]
    print(f"  {system}: xanthone-CoA mindist mean="
          f"{xc['mean']:.3f} nm" if xc and xc.get('n') else f"  {system}: xanthone-CoA n/a")
    return summary


if __name__ == "__main__":
    for s in (sys.argv[1:] or COA_SYS):
        run_system(s)
