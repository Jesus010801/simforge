#!/usr/bin/env python3
"""Phase 2 — core mechanistic descriptors on the managed mdfit.xtc (200 ps).

protein RMSD, ligand RMSD, RMSF (Cα), Rg, SASA, H-bonds (protein-protein AND
ligand-active-site), ligand-active-site min distance + contacts (<0.6 nm),
ligand-catalytic-site COM distance.  For CoA systems the same set is also run
for the CoA molecule.  Writes xvg + provenance.<analysis>.json per observable.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from mech_common import (SYS, OUT, gmx, mdfit, tpr, ndx, refgro, read_xvg, stats,
                         provenance, sha)

CUTOFF = "0.6"


def obs_dir(system, key):
    d = OUT / system / "observables" / key
    d.mkdir(parents=True, exist_ok=True)
    return d


def _finish(system, key, analysis, argv, stdin, out_xvg, ycol_stats=True, extra=None):
    d = obs_dir(system, key)
    prov = provenance(system, analysis, argv=argv, stdin=stdin, outputs=[out_xvg], extra=extra)
    if Path(out_xvg).is_file() and ycol_stats:
        try:
            x, y = read_xvg(out_xvg)
            prov["xvg_check"] = dict(
                n_points=len(x), t_start_ps=float(x[0]) if len(x) else None,
                t_end_ps=float(x[-1]) if len(x) else None,
                all_finite=bool(len(y) and __import__("numpy").isfinite(y[:, 0]).all()),
                series0_stats=stats(y[:, 0]) if len(y) else {})
        except Exception as e:  # noqa
            prov["xvg_check"] = {"error": str(e)}
    (d / f"provenance.{analysis}.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return prov


def run_system(system):
    cfg = SYS[system]
    T, N = tpr(system), ndx(system)
    fit = mdfit(system)
    if not Path(fit).is_file():
        print(f"  {system}: mdfit.xtc not ready — SKIP"); return
    print(f"[{system}] core descriptors")

    # ---- protein RMSD (Backbone, lsq fit Backbone) -tu ns ----
    d = obs_dir(system, "protein_rmsd"); o = d / "rmsd_protein.xvg"
    a = [ "rms", "-s", T, "-f", fit, "-n", N, "-tu", "ns", "-o", o]
    if not o.is_file():
        gmx(a, stdin="Backbone\nBackbone\n", log=d / "gmx.log")
    _finish(system, "protein_rmsd", "protein_rmsd", ["gmx", *a], "Backbone\nBackbone\n", o)

    # ---- RMSF (Cα, per residue) ----
    d = obs_dir(system, "rmsf"); o = d / "rmsf_calpha.xvg"; opdb = d / "rmsf.pdb"
    a = ["rmsf", "-s", T, "-f", fit, "-n", N, "-res", "-o", o, "-oq", opdb]
    if not o.is_file():
        gmx(a, stdin="C-alpha\n", log=d / "gmx.log")
    _finish(system, "rmsf", "rmsf", ["gmx", *a], "C-alpha\n", o)

    # ---- radius of gyration (Protein) -tu ns ----
    d = obs_dir(system, "radius_of_gyration"); o = d / "rg_protein.xvg"
    a = ["gyrate", "-s", T, "-f", fit, "-n", N, "-o", o]
    if not o.is_file():
        gmx(a, stdin="Protein\n", log=d / "gmx.log")
    _finish(system, "radius_of_gyration", "radius_of_gyration", ["gmx", *a], "Protein\n", o)

    # ---- SASA (Protein total) -tu ns ----
    d = obs_dir(system, "sasa"); o = d / "sasa_protein.xvg"
    a = ["sasa", "-s", T, "-f", fit, "-n", N, "-o", o, "-tu", "ns"]
    if not o.is_file():
        gmx(a, stdin="Protein\n", log=d / "gmx.log")
    _finish(system, "sasa", "sasa", ["gmx", *a], "Protein\n", o)

    # ---- H-bonds: whole-protein intramolecular (historical comparability) ----
    # gmx hbond protein-protein on the tetramer is O(N^2) and very slow on shared
    # hardware; it barely discriminates the states (~1210 ± 20 everywhere).  Skip
    # for the big CoA systems unless MECH_FULL_HBOND=1.
    import os
    d = obs_dir(system, "hbonds_protein"); o = d / "hbnum_protein.xvg"
    a = ["hbond", "-s", T, "-f", fit, "-n", N, "-num", o, "-tu", "ns"]
    big = cfg["coa"] is not None
    if not o.is_file() and (not big or os.environ.get("MECH_FULL_HBOND") == "1"):
        gmx(a, stdin="Protein\nProtein\n", log=d / "gmx.log", timeout=5400)
    if o.is_file():
        _finish(system, "hbonds_protein", "hbonds_protein", ["gmx", *a], "Protein\nProtein\n", o)
    else:
        (d / "SKIPPED.txt").write_text(
            "whole-protein H-bond count not computed for this large CoA system "
            "(gmx hbond O(N^2), prohibitively slow on shared hardware; ~1210 ± 20 "
            "across all states — negligible state discrimination). Set "
            "MECH_FULL_HBOND=1 to force.")

    lig = cfg["lig"]
    if lig:
        # ---- ligand RMSD (fit Backbone, RMSD of LIG/COA) -tu ns ----
        d = obs_dir(system, "ligand_rmsd"); o = d / f"rmsd_{lig.lower()}.xvg"
        a = ["rms", "-s", T, "-f", fit, "-n", N, "-tu", "ns", "-o", o]
        if not o.is_file():
            gmx(a, stdin=f"Backbone\n{lig}\n", log=d / "gmx.log")
        _finish(system, "ligand_rmsd", "ligand_rmsd", ["gmx", *a], f"Backbone\n{lig}\n", o,
                extra={"ligand_group": lig, "ligand_kind": cfg["ligkind"]})

        # ---- ligand - active-site H-bonds ----
        d = obs_dir(system, "hbonds_ligand_site"); o = d / f"hbnum_{lig.lower()}_site.xvg"
        a = ["hbond", "-s", T, "-f", fit, "-n", N, "-num", o, "-tu", "ns"]
        if not o.is_file():
            gmx(a, stdin=f"{lig}\nActiveSite_HMG\n", log=d / "gmx.log")
        _finish(system, "hbonds_ligand_site", "hbonds_ligand_site", ["gmx", *a],
                f"{lig}\nActiveSite_HMG\n", o, extra={"pair": f"{lig}-ActiveSite_HMG"})

    if cfg["site"] and lig:
        # ---- ligand - active-site min distance + contacts (<0.6 nm) ----
        d = obs_dir(system, "active_site")
        od = d / f"mindist_{lig.lower()}_active.xvg"; on = d / f"contacts_{lig.lower()}_active.xvg"
        a = ["mindist", "-s", T, "-f", fit, "-n", N, "-d", CUTOFF, "-od", od, "-on", on]
        if not od.is_file():
            gmx(a, stdin=f"{lig}\nActiveSite_HMG\n", log=d / "gmx.log")
        _finish(system, "active_site", "active_site_mindist", ["gmx", *a], f"{lig}\nActiveSite_HMG\n", od)
        _finish(system, "active_site", "active_site_contacts", ["gmx", *a], f"{lig}\nActiveSite_HMG\n", on)

        # ---- ligand - catalytic COM distance ----
        d = obs_dir(system, "catalytic_com_distance"); o = d / f"dist_{lig.lower()}_catalytic.xvg"
        sel = f'com of group "{lig}" plus com of group "Catalytic_HMG"'
        a = ["distance", "-s", T, "-f", fit, "-n", N, "-select", sel, "-oall", o]
        if not o.is_file():
            gmx(a, stdin="", log=d / "gmx.log")
        _finish(system, "catalytic_com_distance", "catalytic_com_distance", ["gmx", *a], "", o,
                extra={"select": sel})

    # ---- CoA-specific metrics for CoA systems ----
    if cfg["coa"] and cfg["coa"] != lig:
        coa = cfg["coa"]
        d = obs_dir(system, "coa_rmsd"); o = d / "rmsd_coa.xvg"
        a = ["rms", "-s", T, "-f", fit, "-n", N, "-tu", "ns", "-o", o]
        if not o.is_file():
            gmx(a, stdin=f"Backbone\n{coa}\n", log=d / "gmx.log")
        _finish(system, "coa_rmsd", "coa_rmsd", ["gmx", *a], f"Backbone\n{coa}\n", o)
        if cfg["site"]:
            d = obs_dir(system, "coa_active_site")
            od = d / "mindist_coa_active.xvg"; on = d / "contacts_coa_active.xvg"
            a = ["mindist", "-s", T, "-f", fit, "-n", N, "-d", CUTOFF, "-od", od, "-on", on]
            if not od.is_file():
                gmx(a, stdin=f"{coa}\nActiveSite_HMG\n", log=d / "gmx.log")
            _finish(system, "coa_active_site", "coa_active_site_mindist", ["gmx", *a], f"{coa}\nActiveSite_HMG\n", od)
            _finish(system, "coa_active_site", "coa_active_site_contacts", ["gmx", *a], f"{coa}\nActiveSite_HMG\n", on)
            d = obs_dir(system, "coa_catalytic_com_distance"); o = d / "dist_coa_catalytic.xvg"
            sel = f'com of group "{coa}" plus com of group "Catalytic_HMG"'
            a = ["distance", "-s", T, "-f", fit, "-n", N, "-select", sel, "-oall", o]
            if not o.is_file():
                gmx(a, stdin="", log=d / "gmx.log")
            _finish(system, "coa_catalytic_com_distance", "coa_catalytic_com_distance", ["gmx", *a], "", o)
    print(f"  {system}: done")


if __name__ == "__main__":
    for s in (sys.argv[1:] or list(SYS)):
        run_system(s)
