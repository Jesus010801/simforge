#!/usr/bin/env python3
"""Incorporate (NOT rerun) the completed 200 ns MM-PBSA for A1 / A3 / A6, plus
the substrate-only and competitive-ternary channels recovered from the A1
source-reconciliation pass (2026-09-10).

Canonical sources (per user instruction — historical total discrepancy for
A1-alone NOT reopened; HMG_CoA_R-A1-R2 treated as authoritative):

  A1 (alone)      /home/jesusxd/Escritorio/HMG_CoA_R-A1-R2/
                  {energy_summary.csv, residues_energy_summary.csv} — g_mmpbsa
                  on md_final.xtc, -b 150000 -e 200000 -dt 500 (101 configs).
                  ΔG = -236.516 +/- 33.740 kJ/mol.
  A3 (alone)      /home/jesusxd/Escritorio/Nuevos_sistemas/A3-HMG-R/ — own
                  mdfit.xtc, -b 150000 -e 200000 -dt 200 (251 configs).
  A6 (alone)      /home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6/ —
                  md_final.xtc, -b 150000 -e 200000 -dt 500 (101 configs).
  substrate-only  /home/jesusxd/Escritorio/HMG-CoA-R_sustrato/ — md_final.xtc,
                  -b 150000 -e 200000 -dt 500 (101 configs).
                  ΔG = -196.585 +/- 34.952 kJ/mol.
  A1 in ternary   /home/jesusxd/Escritorio/Competitive-system-2/
                  {energy_summary-A1.csv, "A1energy_MM.xvg "} — md_final.xtc,
                  -dt 100 (501 configs). ΔG = -211.691 +/- 26.273 kJ/mol.
  CoA in ternary  /home/jesusxd/Escritorio/Competitive-system-2/
                  {COA-energy_summary.csv, COA-energy_MM.xvg} — -dt 500
                  (101 configs). ΔG = -133.832 +/- 58.683 kJ/mol.
  A1<->CoA pair   /home/jesusxd/Escritorio/Competitive-system-2/
                  {A1-COAenergy_summary.csv, A1-COAenergy_MM.xvg} — -dt 500
                  (101 configs). ΔG = -27.410 +/- 22.365 kJ/mol.

HMG_CoA_R-A1 (no "-R2" suffix, 100 ns, smaller box, no MM-PBSA) is superseded
and NOT used. The A1-alone historical total discrepancy (vs the existing
presentation slide's E+I row) is intentionally NOT re-litigated here — see
a1_source_reconciliation_2026-09-10/a1_long_mmpbsa_recovery.md for the record;
HMG_CoA_R-A1-R2 is used as directed.

Writes analysis_outputs/xanthone_mechanistic/mmpbsa_200ns/ (summary + top
residues per system/channel) and reports/.../mmpbsa_200ns_A1_A3_A6.md.
"""
from __future__ import annotations
import csv, json, re
from pathlib import Path

OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")
REP = Path("/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10")
DST = OUT / "mmpbsa_200ns"
DST.mkdir(exist_ok=True)

# label -> (dir, energy_summary filename, energy_MM filename, residues filename or None,
#           trajectory description)
SRC = {
 "A1":  ("/home/jesusxd/Escritorio/HMG_CoA_R-A1-R2",
         "energy_summary.csv", "energy_MM.xvg", "residues_energy_summary.csv",
         "md_final.xtc (HMG_CoA_R-A1-R2 — canonical A1-alone source per "
         "2026-09-10 reconciliation; supersedes the earlier 'no 200 ns MM-PBSA "
         "on disk' statement, which was true only for Mecanismo_inhibitorio/A1)"),
 "A3":  ("/home/jesusxd/Escritorio/Nuevos_sistemas/A3-HMG-R",
         "energy_summary.csv", "energy_MM.xvg", "residues_energy_summary.csv",
         "mdfit.xtc (A3's own, PBC-clean)"),
 "A6":  ("/home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6",
         "energy_summary.csv", "energy_MM.xvg", "residues_energy_summary.csv",
         "md_final.xtc"),
}
# extra channels, same recovery pass, reported alongside but outside the A1/A3/A6 solo comparison
EXTRA = {
 "substrate-only (E+S)": ("/home/jesusxd/Escritorio/HMG-CoA-R_sustrato",
         "energy_summary.csv", "energy_MM.xvg", "residues_energy_summary.csv",
         "md_final.xtc"),
 "A1 in ternary (E+S+I, I channel)": ("/home/jesusxd/Escritorio/Competitive-system-2",
         "energy_summary-A1.csv", "A1energy_MM.xvg ", None,
         "md_final.xtc (unit1 Protein / unit2 LIG)"),
 "CoA in ternary (E+S+I, S channel)": ("/home/jesusxd/Escritorio/Competitive-system-2",
         "COA-energy_summary.csv", "COA-energy_MM.xvg", "COA-residues_energy_summary.csv",
         "md_final.xtc (unit1 Protein / unit2 COA)"),
 "A1<->CoA pairwise (I-S, E+S+I)": ("/home/jesusxd/Escritorio/Competitive-system-2",
         "A1-COAenergy_summary.csv", "A1-COAenergy_MM.xvg", "A1-COAresidues_energy_summary.csv",
         "md_final.xtc (unit1 COA / unit2 LIG)"),
}


def read_energy_summary(p):
    rows = {}
    with open(p) as fh:
        r = csv.reader(fh)
        next(r)
        for row in r:
            if not row or not row[0].strip():
                continue
            rows[row[0].strip()] = dict(mean=float(row[1]), sd=float(row[2]))
    return rows


def cmdline(energy_mm_xvg):
    for ln in open(energy_mm_xvg, errors="replace"):
        if ln.startswith("#   'g_mmpbsa"):
            return ln.strip("# \n")
    return None


AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
       "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX"}


def read_residues(p, top=15):
    rows, lig_row = [], None
    with open(p) as fh:
        r = csv.DictReader(fh)
        for row in r:
            res = (row.get("Residue") or row.get("Residue ") or "").strip()
            if not res:
                continue
            try:
                total = float(row["total"])
            except (KeyError, ValueError):
                continue
            rec = dict(residue=res, vdw=float(row["vDW"]),
                      elec=float(row["Electrostatic"]), polar=float(row["polar"]),
                      apolar=float(row["apolar"]), total=total,
                      total_sd=float(row.get("total-stddev", 0) or 0))
            if res.split("-")[0].upper() not in AA3:
                lig_row = rec
                continue
            rows.append(rec)
    rows.sort(key=lambda x: x["total"])
    return rows[:top], rows, lig_row


def window_stride_n(cmd):
    if not cmd:
        return "?", "?", "?"
    b = re.search(r"-b (\d+)", cmd); e = re.search(r"-e (\d+)", cmd); dt = re.search(r"-dt (\d+)", cmd)
    win = f"{int(b.group(1))/1000:.0f}-{int(e.group(1))/1000:.0f} ns" if b and e else "?"
    stride = f"{dt.group(1)} ps" if dt else "?"
    ncfg = (int(e.group(1)) - int(b.group(1))) // int(dt.group(1)) + 1 if b and e and dt else "?"
    return win, stride, ncfg


def load_channel(d, es_name, mm_name, res_name):
    d = Path(d)
    es = read_energy_summary(d / es_name)
    cmd = cmdline(d / mm_name)
    top, allres, lig_row = ([], [], None)
    if res_name and (d / res_name).is_file():
        top, allres, lig_row = read_residues(d / res_name)
    return dict(source_dir=str(d), energy_summary_kjmol=es, g_mmpbsa_command=cmd,
               n_residues=len(allres), ligand_self_decomposition=lig_row,
               top15_stabilising_residues=top)


def main():
    summary = {}
    for label, (d, es_name, mm_name, res_name, traj) in {**SRC, **EXTRA}.items():
        ch = load_channel(d, es_name, mm_name, res_name)
        ch["source_trajectory"] = traj
        ch["note"] = "incorporated as-is; NOT rerun by this pass"
        summary[label] = ch
        if ch["top15_stabilising_residues"]:
            safe = label.split()[0].replace("/", "-").replace("(", "").replace(")", "")
            (DST / f"top_residues_{safe}.csv").write_text(
                "residue,vdw,elec,polar,apolar,total,total_sd\n" +
                "\n".join(f"{r['residue']},{r['vdw']},{r['elec']},{r['polar']},{r['apolar']},"
                         f"{r['total']},{r['total_sd']}" for r in ch["top15_stabilising_residues"]))
    (DST / "mmpbsa_200ns_summary.json").write_text(json.dumps(summary, indent=2))

    lines = ["# MM-PBSA, 200 ns window — A1 / A3 / A6 and HMG-CoA-reductase mechanistic controls",
             "",
             "Incorporated from already-completed `g_mmpbsa` runs; **not rerun** by this pass.",
             "**A1-alone now has a canonical 200 ns value** (source:"
             " `HMG_CoA_R-A1-R2`, recovered 2026-09-10 — see"
             " `a1_source_reconciliation_2026-09-10/`). The earlier statement"
             " \"A1 has no 200 ns MM-PBSA run\" is **retracted**: it was true only"
             " for `Mecanismo_inhibitorio/A1`, which has no MM-PBSA output on disk;"
             " the run exists in `HMG_CoA_R-A1-R2`, the production directory"
             " underlying that same A1 system (byte-identical `md.tpr`/`md.cpt`/`md.edr`).",
             "",
             "## A1 / A3 / A6 — single-ligand, 200 ns",
             "",
             "| system | ΔG total (kJ/mol) | vDW | Elec | Polar-solv | Non-polar-solv | window | stride | n cfg | source |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for label in ("A1", "A3", "A6"):
        s = summary[label]
        es = s["energy_summary_kjmol"]
        t = es["Total"]
        win, stride, ncfg = window_stride_n(s["g_mmpbsa_command"])
        src_dir = SRC[label][0]
        lines.append(f"| {label} | **{t['mean']:.1f} ± {t['sd']:.1f}** | "
                     f"{es['vDW']['mean']:.1f} | {es['Electrostatic']['mean']:.1f} | "
                     f"{es['Polar-solvation']['mean']:.1f} | {es['Non-polar-solvation']['mean']:.1f} | "
                     f"{win} | {stride} | {ncfg} | `{src_dir}` ({s['source_trajectory']}) |")

    lines += ["", "## Convergence", ""]
    for label, note in [("A1", ""), ("A3", ""), ("A6", "")]:
        t = summary[label]["energy_summary_kjmol"]["Total"]
        ratio = abs(t["sd"] / t["mean"])
        verdict = "well converged" if ratio < 0.3 else "poorly converged (SD approaching |mean|)"
        lines.append(f"- **{label}**: SD/|mean| = {ratio:.2f} — {verdict}.")
    lines += ["",
             "A1 (0.14) and A6 (0.12) are both well converged and now directly comparable at 200 ns; "
             "A3 (0.87) remains poorly converged — its 200 ns value is reported for completeness but "
             "should not be over-interpreted on its own.", ""]

    lines += ["## HMG-CoA-reductase mechanistic controls — substrate-only and competitive ternary "
             "(200 ns; recovered 2026-09-10, see `a1_source_reconciliation_2026-09-10/`)", "",
             "| channel | ΔG total (kJ/mol) | vDW | Elec | Polar-solv | Non-polar-solv | window | stride | n cfg | source |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for label in EXTRA:
        s = summary[label]
        es = s["energy_summary_kjmol"]; t = es["Total"]
        win, stride, ncfg = window_stride_n(s["g_mmpbsa_command"])
        lines.append(f"| {label} | **{t['mean']:.1f} ± {t['sd']:.1f}** | "
                     f"{es['vDW']['mean']:.1f} | {es['Electrostatic']['mean']:.1f} | "
                     f"{es['Polar-solvation']['mean']:.1f} | {es['Non-polar-solvation']['mean']:.1f} | "
                     f"{win} | {stride} | {ncfg} | `{EXTRA[label][0]}` |")
    lines += ["",
             "These reproduce (4 of 4 checked rows, to 2-3 decimals) the existing presentation's MM-PBSA "
             "decomposition table for E+S, E+S+I(I), E+S+I(S) and I-S — see "
             "`a1_source_reconciliation_2026-09-10/a1_long_mmpbsa_recovery.md` §4 for the row-by-row check.",
             ""]

    for label in ("A1", "A3", "A6", *EXTRA):
        s = summary[label]
        top = s.get("top15_stabilising_residues")
        if not top:
            continue
        lig = s.get("ligand_self_decomposition")
        lines.append(f"## {label} — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)")
        lines.append("")
        if lig:
            lines.append(f"*(ligand's own decomposition row, {lig['residue']}: total "
                         f"{lig['total']:.1f} ± {lig['total_sd']:.1f} — excluded from the ranking below, "
                         f"which is protein residues only)*")
            lines.append("")
        lines.append("| residue | total | vdW | elec | polar | apolar |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for r in top[:10]:
            lines.append(f"| {r['residue']} | {r['total']:.2f} ± {r['total_sd']:.2f} | "
                         f"{r['vdw']:.2f} | {r['elec']:.2f} | {r['polar']:.2f} | {r['apolar']:.2f} |")
        lines.append("")

    lines += ["## A1<->CoA pairwise per-ligand share", "",
             "`A1-COAresidues_energy_summary.csv` decomposes the pairwise ΔG additively between the "
             "two ligand \"residues\" themselves:", ""]
    pw = summary["A1<->CoA pairwise (I-S, E+S+I)"]
    if pw.get("top15_stabilising_residues") is not None:
        # this channel's residues ARE the two ligands - re-read without the AA3 filter for display
        d = Path(EXTRA["A1<->CoA pairwise (I-S, E+S+I)"][0])
        rows = list(csv.DictReader(open(d / "A1-COAresidues_energy_summary.csv")))
        lines.append("| entry | total ΔG (kJ/mol) |")
        lines.append("| --- | --- |")
        for row in rows:
            res = (row.get("Residue") or "").strip()
            if res:
                lines.append(f"| {res} | {float(row['total']):.2f} ± {float(row.get('total-stddev',0) or 0):.2f} |")
        lines.append("")

    (REP / "mmpbsa_200ns_A1_A3_A6.md").write_text("\n".join(lines))
    print("wrote", DST / "mmpbsa_200ns_summary.json", "and", REP / "mmpbsa_200ns_A1_A3_A6.md")
    for label in ("A1", "A3", "A6", *EXTRA):
        t = summary[label]["energy_summary_kjmol"]["Total"]
        print(f"  {label}: {t['mean']:.1f} +/- {t['sd']:.1f} kJ/mol")


if __name__ == "__main__":
    main()
