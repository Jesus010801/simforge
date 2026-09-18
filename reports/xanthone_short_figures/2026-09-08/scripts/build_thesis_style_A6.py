#!/usr/bin/env python3
"""Complete the thesis '25ns-DM' general short-study figures with A6, in the historical Grace style.

Style is reverse-engineered from the historical Grace 5.1.25 SVG exports in
`Resultados-Tesis-maestría/25ns-DM/` (see scratchpad/extract_style.py output and style_audit.md):

  canvas         11 x 8.5 in landscape (792 x 612 pt)
  font           serif (Nimbus Roman / Times), black, normal weight
  data lines     width 1.377 pt, solid
  frame          all four spines black 0.918 pt, ticks inward on all four sides,
                 major + one minor, no grid, white background
  x axis         "Tiempo (ns)", 0-25, majors every 5
  y axis         "RMSD (nm)" / "Distancia (nm)" (data-driven range, not forced to 0 for distances)
  legend         inside the axes, no box, short colour segment + label, order A1..A6
  colours        A1 #FF00FF  A2 #FFFF00  A3 #FFA500  A4 #00FFFF  A5 #7221BC   (recovered, identical
                 across all 12 historical general figures)
                 A6 #00A000  (new, green - the one saturated hue absent from the historical five)

Reads only. Writes only under completadas_A6/. Never touches an original file.
"""
from __future__ import annotations
import json, math, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

HIST = Path("/home/jesusxd/Escritorio/Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM")
A6ROOT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_short")
OUT = HIST / "completadas_A6"
SVGDIR = OUT / "svg"

TARGETS = ["AA", "AG", "HMG", "LP"]
COMPOUNDS = ["A1", "A2", "A3", "A4", "A5", "A6"]

# --- recovered historical colour map (verbatim hex from the Grace SVGs) + new A6
COLOR = {
    "A1": "#FF00FF", "A2": "#FFFF00", "A3": "#FFA500",
    "A4": "#00FFFF", "A5": "#7221BC", "A6": "#00A000",
}
LINE_W = 1.377
FRAME_W = 0.918

# descriptor family -> (out-stem, y-label, historical style template, has_historical_general)
FAMILIES = {
    "protein_rmsd": ("rmsd",     "RMSD (nm)",       {"AA": "AA-rmsd.svg", "AG": "AG-rmsd.svg",
                                                     "HMG": "svg/HMG-rmsd.svg", "LP": "LP-rmsd.svg"}, True),
    "ligand_rmsd":  ("LIG-rmsd", "RMSD (nm)",       {"AA": "AA-LIG-rmsd.svg", "AG": "AG-rmsd-ligand.svg",
                                                     "HMG": "svg/HMG-rmsd-ligand.svg", "LP": "LP-rmsd-ligand.svg"}, True),
    "catalytic_dist": ("dist",   "Distancia (nm)",  {"AA": "AA-dist.svg", "AG": "AG-dist.svg",
                                                     "HMG": "HMG-dist.svg", "LP": "LP-dist.svg"}, True),
    "active_mindist": ("mindist", "Distancia mínima (nm)",
                       {"AA": "AA-dist.svg", "AG": "AG-dist.svg", "HMG": "HMG-dist.svg", "LP": "LP-dist.svg"}, False),
    "active_contacts": ("contacts", "Contactos (< 0.6 nm)",
                        {"AA": "AA-dist.svg", "AG": "AG-dist.svg", "HMG": "HMG-dist.svg", "LP": "LP-dist.svg"}, False),
}

A6_DIR = {"AA": "AA-A6", "AG": "AG-A6", "HMG": "HMG-R-25ns-A6", "LP": "LP-A6"}
A6_REL = {
    "protein_rmsd":    "observables/protein_rmsd/rmsd_protein.xvg",
    "ligand_rmsd":     "observables/ligand_rmsd/rmsd_ligand.xvg",
    "active_mindist":  "observables/active_site/mindist_lig_active.xvg",
    "active_contacts": "observables/active_site/contacts_lig_active.xvg",
    "catalytic_dist":  "observables/catalytic_com_distance/dist_lig_catalytic.xvg",
}


def hist_candidates(target, comp, desc):
    t = target
    if desc == "protein_rmsd":
        return [f"{t}-{comp}rmsd_protein.xvg"]
    if desc == "ligand_rmsd":
        if t == "HMG":
            return [f"HMG_CoA_R-{comp}_rmsd_ligand.xvg"]
        return [f"{t}-{comp}_rmsd_ligand.xvg", f"{t}-{comp}_rmsd-ligand.xvg",
                f"{t}-{comp}rmsd_ligand.xvg", f"{t}-{comp}rmsd-ligand.xvg"]
    if desc == "active_mindist":
        return [f"{t}-{comp}mindist_lig_active.xvg"]
    if desc == "active_contacts":
        return [f"{t}-{comp}contacts_lig_active.xvg"]
    if desc == "catalytic_dist":
        return [f"{t}-{comp}dist_lig_catalytic.xvg"]
    return []


def resolve(target, comp, desc):
    if comp == "A6":
        p = A6ROOT / A6_DIR[target] / A6_REL[desc]
        return (p if p.exists() else None, "A6-simforge")
    for name in hist_candidates(target, comp, desc):
        p = HIST / name
        if p.exists():
            return (p, "historical")
    return (None, "historical")


def parse_xvg(path):
    xs, ys = [], []
    sub = ""
    for line in open(path, "r", errors="replace"):
        s = line.strip()
        if not s or s[0] == "#":
            continue
        if s[0] == "@":
            m = re.match(r'@\s+subtitle\s+"(.*)"', s)
            if m:
                sub = m.group(1)
            continue
        p = s.split()
        try:
            xs.append(float(p[0])); ys.append(float(p[1]))
        except (ValueError, IndexError):
            continue
    t = np.asarray(xs, float); v = np.asarray(ys, float)
    tmax = float(t[-1]) if t.size else 0.0
    t_ns = t / 1000.0 if tmax > 500 else t
    return t_ns, v, dict(subtitle=sub, n=int(v.size),
                         tmin=float(t_ns[0]) if t_ns.size else None,
                         tmax=float(t_ns[-1]) if t_ns.size else None,
                         vmax=float(np.nanmax(v)) if v.size else None)


def apply_grace_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Nimbus Roman", "Liberation Serif", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "figure.figsize": (11.0, 8.5),
        "figure.dpi": 100,
        "savefig.dpi": 200,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.linewidth": FRAME_W,
        "axes.edgecolor": "black",
        "axes.grid": False,
        "axes.labelcolor": "black",
        "axes.titlesize": 20,
        "axes.labelsize": 20,
        "xtick.labelsize": 18, "ytick.labelsize": 18,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 7, "ytick.major.size": 7,
        "xtick.minor.size": 4, "ytick.minor.size": 4,
        "xtick.major.width": FRAME_W, "ytick.major.width": FRAME_W,
        "xtick.minor.width": FRAME_W, "ytick.minor.width": FRAME_W,
        "xtick.color": "black", "ytick.color": "black",
        "legend.fontsize": 17, "legend.frameon": False,
        "legend.handlelength": 1.6, "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.35,
        "svg.fonttype": "none",
        "path.simplify": True, "path.simplify_threshold": 0.0,
    })


def main():
    OUT.mkdir(exist_ok=True)
    SVGDIR.mkdir(exist_ok=True)
    apply_grace_style()

    inventory = []
    qc = []
    audit_rows = []

    for target in TARGETS:
        for desc, (stem, ylabel, template_map, has_hist) in FAMILIES.items():
            traces, comps_present = {}, []
            src_records = []
            for comp in COMPOUNDS:
                path, stype = resolve(target, comp, desc)
                if path is None:
                    src_records.append(dict(compound=comp, found=False, source_type=stype))
                    continue
                t_ns, v, meta = parse_xvg(path)
                traces[comp] = (t_ns, v, meta)
                comps_present.append(comp)
                src_records.append(dict(compound=comp, found=True, source_type=stype,
                                        source_path=str(path), n=meta["n"],
                                        tmin_ns=round(meta["tmin"], 3), tmax_ns=round(meta["tmax"], 3),
                                        vmax=round(meta["vmax"], 4), subtitle=meta["subtitle"]))

            # --- physical-plausibility guard (LP-A6 ligand RMSD PBC artefact)
            artefacts = []
            for comp in comps_present:
                _, v, meta = traces[comp]
                if desc in ("protein_rmsd", "ligand_rmsd") and meta["vmax"] and meta["vmax"] > 3.0:
                    artefacts.append(dict(compound=comp, vmax=round(meta["vmax"], 3)))

            # --- figure
            fig, ax = plt.subplots()
            yhi_clip = None
            if artefacts:
                others = [np.percentile(traces[c][1][np.isfinite(traces[c][1])], 99.5)
                          for c in comps_present
                          if c not in [a["compound"] for a in artefacts]]
                if others:
                    yhi_clip = max(others) * 1.30

            for comp in COMPOUNDS:                      # fixed A1..A6 order = legend order
                if comp not in traces:
                    continue
                t_ns, v, meta = traces[comp]
                vp = v.astype(float).copy()
                lbl = comp
                if yhi_clip is not None and comp in [a["compound"] for a in artefacts]:
                    vp = np.where(v > yhi_clip, np.nan, v)
                    lbl = f"{comp} *"
                ax.plot(t_ns, vp, color=COLOR[comp], lw=LINE_W, solid_capstyle="round",
                        label=lbl, zorder=3 if comp != "A6" else 4)

            ax.set_xlim(0, 25)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            if yhi_clip is not None:
                ax.set_ylim(bottom=min(0.0, ax.get_ylim()[0]), top=yhi_clip)
            ax.set_xlabel("Tiempo (ns)")
            ax.set_ylabel(ylabel)
            for sp in ax.spines.values():
                sp.set_linewidth(FRAME_W)
                sp.set_color("black")
            leg = ax.legend(loc="best")
            for txt in leg.get_texts():
                txt.set_color("black")
            if artefacts:
                note = "; ".join(f"{a['compound']}: excursión a ~{a['vmax']:.1f} nm fuera de escala "
                                 f"(artefacto PBC sospechado; ver qc_report.json)" for a in artefacts)
                fig.text(0.13, 0.02, "* " + note, fontsize=12, color="#7A0000", ha="left", va="bottom")

            fig.tight_layout(rect=(0, 0.045 if artefacts else 0, 1, 1))

            out_stem = f"{target}-{stem}-A1-A6"
            paths_written = []
            for ext in ("svg", "png", "pdf"):
                p = OUT / f"{out_stem}.{ext}"
                fig.savefig(p)
                paths_written.append(str(p))
            fig.savefig(SVGDIR / f"{out_stem}.svg")
            plt.close(fig)

            template = template_map[target]
            template_path = HIST / template
            inventory.append(dict(
                figure=f"{out_stem}.svg", target=target, descriptor=desc,
                compounds=comps_present, complete=(comps_present == COMPOUNDS),
                outputs=[f"{out_stem}.{e}" for e in ("svg", "png", "pdf")],
                style_template=str(template_path) if template_path.exists() else None,
                style_template_is_same_family=has_hist,
                sources=src_records,
            ))
            audit_rows.append(dict(
                figure=f"{out_stem}.svg", target=target, descriptor=desc,
                style_template=(str(template_path) if has_hist else
                                f"{template} (nearest-family template; no historical general "
                                f"figure exists for this descriptor)"),
                a1_a5_colours={c: COLOR[c] for c in ["A1", "A2", "A3", "A4", "A5"]},
                a6_colour=COLOR["A6"],
                x_axis="Tiempo (ns), 0-25  (historical: same; AA-rmsd/LP-rmsd had 'Tiempo'/'TIempo' "
                       "typos - standardised)",
                y_axis=(f"'{ylabel}' - exact match to historical" if has_hist else
                        f"'{ylabel}' - new wording, Spanish, consistent with historical 'Distancia (nm)'"),
                differences=("Rendered with matplotlib rather than Grace: identical palette, line "
                             "width (1.377 pt), serif font (Nimbus Roman), inward ticks on 4 sides, "
                             "no grid, borderless in-axes legend. Legend order forced to A1..A6 "
                             "(historical draw order varied per file). Y-range auto (matches "
                             "historical data-driven ranging)."
                             + (" LP-A6 ligand RMSD off-scale portion masked; raw data unchanged."
                                if artefacts else "")),
            ))
            cov_ok = all(24.0 <= traces[c][2]["tmax"] <= 26.0 and traces[c][2]["tmin"] <= 0.5
                         for c in comps_present)
            qc.append(dict(
                figure=f"{out_stem}.svg", target=target, descriptor=desc,
                compounds_included=comps_present, all_six_present=(comps_present == COMPOUNDS),
                data_sources=[r["source_path"] for r in src_records if r.get("found")],
                time_coverage_ns={c: [round(traces[c][2]["tmin"], 3), round(traces[c][2]["tmax"], 3)]
                                  for c in comps_present},
                n_points={c: traces[c][2]["n"] for c in comps_present},
                coverage_0_25ns_ok=bool(cov_ok),
                suspected_artifacts=artefacts,
                hmg_window="25ns short profile (HMG-R-25ns-A6)" if target == "HMG" else None,
                validation_status=("WARN_ARTIFACT" if artefacts else
                                   "PASS" if (comps_present == COMPOUNDS and cov_ok) else "WARN"),
            ))

    (OUT / "figure_inventory_completed.json").write_text(json.dumps({
        "generated": "2026-09-08",
        "output_dir": str(OUT),
        "historical_style_source": str(HIST),
        "a6_data_source": str(A6ROOT),
        "recovered_colour_map": COLOR,
        "figures": inventory,
    }, indent=2, ensure_ascii=False))

    (OUT / "qc_report.json").write_text(json.dumps({
        "generated": "2026-09-08",
        "sd_convention": "n/a (time-series figures)",
        "lp_a6_ligand_rmsd_caveat": ("LP-A6 rmsd_ligand.xvg reaches ~7.7 nm while the ligand keeps "
                                     "active-site contact (mindist ~0.2 nm, contacts ~1000): suspected "
                                     "periodic-image artefact in that single xvg. Plotted honestly with "
                                     "the off-scale part masked and a caption; raw data NOT modified. "
                                     "Fix later with `gmx trjconv -pbc mol -center` / `-pbc nojump`."),
        "entries": qc,
    }, indent=2, ensure_ascii=False))

    npass = sum(1 for e in qc if e["validation_status"] == "PASS")
    print(f"figures written : {len(inventory)}  (svg+png+pdf each, plus svg/ copies)")
    print(f"QC PASS         : {npass}/{len(qc)}")
    for e in qc:
        if e["validation_status"] != "PASS":
            print("   ", e["validation_status"], e["figure"], e["suspected_artifacts"])
    return audit_rows


if __name__ == "__main__":
    rows = main()
    # dump audit rows for the markdown writer
    Path("/tmp/claude-1000/-home-jesusxd-Escritorio-simforge/"
         "e69d8c09-869c-4cff-acd6-f8b2e374ecae/scratchpad/audit_rows.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False))
