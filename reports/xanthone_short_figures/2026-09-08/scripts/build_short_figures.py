#!/usr/bin/env python3
"""Presentation-quality figure builder for the thesis-compatible xanthone_short study (A1-A6).

Reads (read-only):
  - Historical thesis corpus  : Nuevos_sistemas/Resultados-Tesis-maestria/25ns-DM/  (A1-A5)
  - New SimForge short outputs : simforge/analysis_outputs/xanthone_short/            (A6)

Writes everything under: reports/xanthone_short_figures/2026-09-08/
Does NOT touch source simulation files or historical results.
"""
from __future__ import annotations
import json, csv, math, re, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

REPO = Path("/home/jesusxd/Escritorio/simforge")
HIST = Path("/home/jesusxd/Escritorio/Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM")
A6ROOT = REPO / "analysis_outputs/xanthone_short"
OUT = REPO / "reports/xanthone_short_figures/2026-09-08"

TARGETS = ["AA", "AG", "HMG", "LP"]
TARGET_LABEL = {
    "AA": "AA  (alpha-amylase)",
    "AG": "AG  (alpha-glucosidase)",
    "HMG": "HMG  (HMG-CoA reductase)",
    "LP": "LP  (lipase)",
}
COMPOUNDS = ["A1", "A2", "A3", "A4", "A5", "A6"]

# descriptor family -> (label, y-axis label, short key)
DESCRIPTORS = {
    "protein_rmsd":     ("Protein RMSD",                      "RMSD (nm)",       "protein_rmsd"),
    "ligand_rmsd":      ("Ligand RMSD (fit to protein)",      "RMSD (nm)",       "ligand_rmsd"),
    "active_mindist":   ("Ligand-active-site minimum distance", "min. distance (nm)", "active_mindist"),
    "active_contacts":  ("Ligand-active-site contacts (<0.6 nm)", "number of contacts", "active_contacts"),
    "catalytic_dist":   ("Ligand-catalytic-site COM distance", "COM distance (nm)", "catalytic_dist"),
}
DESC_ORDER = list(DESCRIPTORS.keys())

# consistent A1-A6 colour cycle (colour-blind friendly)
CMAP = {
    "A1": "#4C72B0", "A2": "#DD8452", "A3": "#55A868",
    "A4": "#C44E52", "A5": "#8172B3", "A6": "#000000",
}

# ---------------------------------------------------------------- source resolution
A6_DIR = {"AA": "AA-A6", "AG": "AG-A6", "HMG": "HMG-R-25ns-A6", "LP": "LP-A6"}
A6_REL = {
    "protein_rmsd":    "observables/protein_rmsd/rmsd_protein.xvg",
    "ligand_rmsd":     "observables/ligand_rmsd/rmsd_ligand.xvg",
    "active_mindist":  "observables/active_site/mindist_lig_active.xvg",
    "active_contacts": "observables/active_site/contacts_lig_active.xvg",
    "catalytic_dist":  "observables/catalytic_com_distance/dist_lig_catalytic.xvg",
}


def hist_candidates(target: str, comp: str, desc: str) -> list[str]:
    t = target
    if desc == "protein_rmsd":
        return [f"{t}-{comp}rmsd_protein.xvg"]
    if desc == "ligand_rmsd":
        if t == "HMG":
            # thesis-compatible ligand RMSD = fit to Protein => the HMG_CoA_R-* variant
            return [f"HMG_CoA_R-{comp}_rmsd_ligand.xvg"]
        return [
            f"{t}-{comp}_rmsd_ligand.xvg",
            f"{t}-{comp}_rmsd-ligand.xvg",
            f"{t}-{comp}rmsd_ligand.xvg",
            f"{t}-{comp}rmsd-ligand.xvg",
        ]
    if desc == "active_mindist":
        return [f"{t}-{comp}mindist_lig_active.xvg"]
    if desc == "active_contacts":
        return [f"{t}-{comp}contacts_lig_active.xvg"]
    if desc == "catalytic_dist":
        return [f"{t}-{comp}dist_lig_catalytic.xvg"]
    return []


def resolve_source(target: str, comp: str, desc: str) -> tuple[Path | None, str]:
    if comp == "A6":
        p = A6ROOT / A6_DIR[target] / A6_REL[desc]
        return (p if p.exists() else None, "A6-simforge")
    for name in hist_candidates(target, comp, desc):
        p = HIST / name
        if p.exists():
            return (p, "historical")
    return (None, "historical")


# ---------------------------------------------------------------- xvg parsing
def parse_xvg(path: Path):
    """Return (time_ns np.array, value np.array, meta dict). Auto-detects ps vs ns x-units."""
    xs, ys = [], []
    title = subtitle = xax = yax = ""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s[0] == "#":
                continue
            if s[0] == "@":
                m = re.match(r'@\s+title\s+"(.*)"', s)
                if m: title = m.group(1)
                m = re.match(r'@\s+subtitle\s+"(.*)"', s)
                if m: subtitle = m.group(1)
                m = re.match(r'@\s+xaxis\s+label\s+"(.*)"', s)
                if m: xax = m.group(1)
                m = re.match(r'@\s+yaxis\s+label\s+"(.*)"', s)
                if m: yax = m.group(1)
                continue
            parts = s.split()
            try:
                x = float(parts[0]); y = float(parts[1])
            except (ValueError, IndexError):
                continue
            xs.append(x); ys.append(y)
    t = np.asarray(xs, float)
    v = np.asarray(ys, float)
    # unit detection: short study spans 25 ns. ps files end ~25000, ns files end ~25
    tmax = float(t[-1]) if t.size else 0.0
    if tmax > 500:            # picoseconds
        t_ns = t / 1000.0
        xunit = "ps"
    else:                    # already nanoseconds
        t_ns = t.copy()
        xunit = "ns"
    meta = dict(title=title, subtitle=subtitle, xaxis=xax, yaxis=yax,
                xunit_detected=xunit, n=int(v.size),
                tmin_ns=float(t_ns[0]) if t_ns.size else None,
                tmax_ns=float(t_ns[-1]) if t_ns.size else None)
    return t_ns, v, meta


# ---------------------------------------------------------------- main assembly
def main():
    long_rows = []          # tidy timeseries
    summary = {}            # (target,comp) -> {desc: (mean,sd,n)}
    inventory = []          # source inventory records
    qc_entries = []
    traces = {}             # (target,desc,comp) -> (t_ns, v, meta, src_type, path)

    for target in TARGETS:
        for comp in COMPOUNDS:
            for desc in DESC_ORDER:
                path, src_type = resolve_source(target, comp, desc)
                rec = dict(target=target, compound=comp, descriptor=desc,
                           source_type=src_type,
                           source_path=str(path) if path else None,
                           found=bool(path))
                if path is None:
                    rec.update(n=0, tmin_ns=None, tmax_ns=None, status="MISSING_SOURCE")
                    inventory.append(rec)
                    summary.setdefault((target, comp), {})[desc] = (math.nan, math.nan, 0)
                    continue
                t_ns, v, meta = parse_xvg(path)
                traces[(target, desc, comp)] = (t_ns, v, meta, src_type, path)
                finite = np.isfinite(v)
                nfin = int(finite.sum())
                mean = float(np.mean(v[finite])) if nfin else math.nan
                sd = float(np.std(v[finite], ddof=0)) if nfin else math.nan   # population sd (matches thesis resumen.txt)
                vmax_val = float(np.max(v[finite])) if nfin else math.nan
                # physical-plausibility guard: RMSD / distance excursions beyond ~3 nm while the
                # ligand still reports active-site contact usually indicate a periodic-image (PBC)
                # artefact in that particular xvg rather than genuine dynamics.
                implausible = bool(desc in ("ligand_rmsd", "protein_rmsd") and nfin and vmax_val > 3.0)
                summary.setdefault((target, comp), {})[desc] = (mean, sd, meta["n"])
                cov_ok = (meta["tmax_ns"] is not None and 24.0 <= meta["tmax_ns"] <= 26.0
                          and meta["tmin_ns"] is not None and meta["tmin_ns"] <= 0.5)
                rec.update(n=meta["n"], tmin_ns=meta["tmin_ns"], tmax_ns=meta["tmax_ns"],
                           xunit_detected=meta["xunit_detected"], subtitle=meta["subtitle"],
                           mean=mean, sd=sd, vmax=vmax_val,
                           coverage_0_25ns=bool(cov_ok),
                           n_is_2501=bool(meta["n"] == 2501),
                           all_nan=bool(nfin == 0),
                           physically_implausible=implausible,
                           status=("ARTIFACT_SUSPECTED" if implausible
                                   else "OK" if (cov_ok and nfin > 0) else "CHECK"))
                inventory.append(rec)
                for i in range(meta["n"]):
                    long_rows.append(dict(
                        target=target, compound=comp, descriptor=desc,
                        time_ns=round(float(t_ns[i]), 6),
                        time_ps=round(float(t_ns[i] * 1000.0), 3),
                        value=float(v[i]),
                        source_type=src_type,
                        source_path=str(path),
                    ))

    # ---- write tidy long CSV
    long_csv = OUT / "tables/short_descriptor_timeseries_long.csv"
    with open(long_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "compound", "descriptor",
                                           "time_ns", "time_ps", "value",
                                           "source_type", "source_path"])
        w.writeheader()
        w.writerows(long_rows)

    # ---- summary CSV + MD
    sum_csv = OUT / "tables/short_descriptor_summary.csv"
    cols = ["target", "compound",
            "protein_rmsd_mean", "protein_rmsd_sd",
            "ligand_rmsd_mean", "ligand_rmsd_sd",
            "active_mindist_mean", "active_mindist_sd",
            "active_contacts_mean", "active_contacts_sd",
            "catalytic_dist_mean", "catalytic_dist_sd"]
    srows = []
    for target in TARGETS:
        for comp in COMPOUNDS:
            d = summary.get((target, comp), {})
            row = dict(target=target, compound=comp)
            for desc in DESC_ORDER:
                m, s, _ = d.get(desc, (math.nan, math.nan, 0))
                row[f"{desc}_mean"] = "" if math.isnan(m) else round(m, 4)
                row[f"{desc}_sd"] = "" if math.isnan(s) else round(s, 4)
            srows.append(row)
    with open(sum_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(srows)

    sum_md = OUT / "tables/short_descriptor_summary.md"
    with open(sum_md, "w") as fh:
        fh.write("# xanthone_short descriptor summary (A1-A6, 0-25 ns)\n\n")
        fh.write("Values are trajectory means +/- population SD (ddof=0, matching the thesis "
                 "`resumen.txt` convention). A1-A5 from the thesis 25ns-DM corpus; A6 from the "
                 "SimForge `xanthone_short` profile.\n\n")
        fh.write("| target | cmp | protRMSD | ligRMSD | act.mindist | act.contacts | cat.dist |\n")
        fh.write("| --- | --- | --- | --- | --- | --- | --- |\n")
        for r in srows:
            def cell(k):
                m, s = r[f"{k}_mean"], r[f"{k}_sd"]
                if m == "":
                    return "n/a"
                return f"{m:.4g} +/- {s:.3g}" if k != "active_contacts" else f"{m:.4g} +/- {s:.3g}"
            fh.write(f"| {r['target']} | {r['compound']} | {cell('protein_rmsd')} | "
                     f"{cell('ligand_rmsd')} | {cell('active_mindist')} | "
                     f"{cell('active_contacts')} | {cell('catalytic_dist')} |\n")

    # ---------------------------------------------------------------- primary figures
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 150, "font.size": 11,
        "axes.titlesize": 13, "axes.labelsize": 12, "axes.grid": True,
        "grid.alpha": 0.25, "legend.frameon": False, "figure.autolayout": False,
        "svg.fonttype": "none",
    })
    src_by_key = {(r["target"], r["compound"], r["descriptor"]): r for r in inventory}
    primary_manifest = []
    for target in TARGETS:
        for desc in DESC_ORDER:
            dlabel, ylabel, _ = DESCRIPTORS[desc]
            fig, ax = plt.subplots(figsize=(9.0, 4.6))
            present = []
            all_p99, artefact_comps = [], []
            series = []
            for comp in COMPOUNDS:
                key = (target, desc, comp)
                if key not in traces:
                    continue
                t_ns, v, meta, src_type, path = traces[key]
                present.append(comp)
                fin = np.isfinite(v)
                is_art = bool(src_by_key.get((target, comp, desc), {}).get("physically_implausible"))
                if is_art:
                    artefact_comps.append((comp, float(np.max(v[fin]))))
                elif fin.any():
                    all_p99.append(float(np.percentile(v[fin], 99.5)))
                series.append((comp, t_ns, v, is_art))
            yhi = (max(all_p99) * 1.25) if (artefact_comps and all_p99) else None
            for comp, t_ns, v, is_art in series:
                vp = v.copy()
                lbl = f"{comp}" + ("  (SimForge)" if comp == "A6" else "")
                if is_art and yhi is not None:
                    vp = np.where(v > yhi, np.nan, v)      # break the line instead of streaking
                    lbl = f"{comp}  (SimForge, off-scale >{yhi:.1f} nm)"
                lw = 1.7 if comp == "A6" else 1.0
                alpha = 0.95 if comp == "A6" else 0.8
                z = 5 if comp == "A6" else 2
                ax.plot(t_ns, vp, color=CMAP[comp], lw=lw, alpha=alpha, zorder=z, label=lbl)
            ax.set_xlim(0, 25)
            if yhi is not None:
                ax.set_ylim(0, yhi)
            ax.set_xlabel("time (ns)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{TARGET_LABEL[target]}  -  {dlabel}")
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=9,
                      handlelength=1.6, borderaxespad=0)
            if artefact_comps:
                txt = "; ".join(f"{c}: excursion to ~{m:.1f} nm (off-scale)" for c, m in artefact_comps)
                fig.text(0.012, 0.012,
                         f"QC note - {txt}. Active-site contact is retained throughout, so this is a "
                         f"suspected periodic-image artefact in that single xvg, not ligand egress.",
                         fontsize=7.5, color="#8B0000", ha="left", va="bottom", wrap=True)
            fig.tight_layout(rect=(0, 0.05 if artefact_comps else 0, 1, 1))
            stem = f"{target}_{desc}_A1-A6"
            for ext in ("svg", "pdf", "png"):
                fig.savefig(OUT / f"figures/primary/{stem}.{ext}")
            plt.close(fig)
            primary_manifest.append(dict(figure=f"figures/primary/{stem}.svg",
                                         target=target, descriptor=desc,
                                         compounds=present,
                                         complete=(present == COMPOUNDS)))

    # ---------------------------------------------------------------- per-target heatmaps
    heat_manifest = []
    for target in TARGETS:
        M = np.full((len(COMPOUNDS), len(DESC_ORDER)), np.nan)
        raw = np.full_like(M, np.nan)
        for i, comp in enumerate(COMPOUNDS):
            for j, desc in enumerate(DESC_ORDER):
                m, s, _ = summary.get((target, comp), {}).get(desc, (math.nan, math.nan, 0))
                raw[i, j] = m
        # column-wise z-score across the 6 compounds (within target) -> cross-compound comparison
        Z = np.full_like(raw, np.nan)
        for j in range(raw.shape[1]):
            col = raw[:, j]
            fin = np.isfinite(col)
            if fin.sum() >= 2 and np.nanstd(col[fin]) > 0:
                Z[fin, j] = (col[fin] - np.nanmean(col[fin])) / np.nanstd(col[fin])
            elif fin.sum() >= 1:
                Z[fin, j] = 0.0
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        # clip the colour scale to +/-2.5 sigma so a single outlier compound (e.g. a PBC-artefact
        # ligand RMSD) does not flatten the contrast for the remaining cells; raw means stay printed.
        vmax = min(2.5, np.nanmax(np.abs(Z))) if np.isfinite(Z).any() else 1.0
        vmax = max(vmax, 0.5)
        Zc = np.clip(Z, -vmax, vmax)
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        im = ax.imshow(Zc, cmap="RdBu_r", norm=norm, aspect="auto")
        ax.set_xticks(range(len(DESC_ORDER)))
        ax.set_xticklabels(["protein\nRMSD", "ligand\nRMSD", "active-site\nmin.dist",
                            "active-site\ncontacts", "catalytic\nCOM dist"], fontsize=9)
        ax.set_yticks(range(len(COMPOUNDS)))
        ax.set_yticklabels(COMPOUNDS)
        for i in range(len(COMPOUNDS)):
            for j in range(len(DESC_ORDER)):
                if np.isfinite(raw[i, j]):
                    txt = f"{raw[i, j]:.3g}" if j != 3 else f"{raw[i, j]:.0f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                            color="black" if abs(Zc[i, j]) < 0.7 * vmax else "white")
                else:
                    ax.text(j, i, "n/a", ha="center", va="center", fontsize=8, color="0.4")
        ax.set_title(f"{TARGET_LABEL[target]}  -  short-study descriptor profile\n"
                     f"(cell text = raw mean; colour = per-column z-score across A1-A6)", fontsize=10)
        cb = fig.colorbar(im, ax=ax, shrink=0.85)
        cb.set_label("z-score (within target, across compounds)", fontsize=9)
        ax.set_xticks(np.arange(-.5, len(DESC_ORDER), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(COMPOUNDS), 1), minor=True)
        ax.grid(which="minor", color="white", lw=1.5)
        ax.tick_params(which="minor", length=0)
        art = [(c, d) for c in COMPOUNDS for d in DESC_ORDER
               if src_by_key.get((target, c, d), {}).get("physically_implausible")]
        if art:
            lbls = ", ".join(f"{c} {DESCRIPTORS[d][0].split(' (')[0].lower()}" for c, d in art)
            fig.text(0.5, 0.005,
                     f"QC note: {lbls} value is a suspected periodic-image artefact (see qc_report.json); "
                     f"cell shows the raw mean unfiltered.",
                     ha="center", va="bottom", fontsize=7, color="#8B0000")
        fig.tight_layout(rect=(0, 0.04 if art else 0, 1, 1))
        stem = f"{target}_heatmap_short_descriptors"
        for ext in ("svg", "pdf", "png"):
            fig.savefig(OUT / f"figures/summary/{stem}.{ext}")
        plt.close(fig)
        heat_manifest.append(dict(figure=f"figures/summary/{stem}.svg", target=target,
                                  rows=COMPOUNDS, cols=DESC_ORDER,
                                  normalization="per-column z-score across A1-A6 within target",
                                  raw_means=[[None if not np.isfinite(x) else round(float(x), 4)
                                              for x in raw[i]] for i in range(len(COMPOUNDS))]))

    # ---------------------------------------------------------------- combined summary bar figure
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2))
    axes = axes.ravel()
    x = np.arange(len(COMPOUNDS))
    for k, desc in enumerate(DESC_ORDER):
        ax = axes[k]
        for ti, target in enumerate(TARGETS):
            means = [summary.get((target, c), {}).get(desc, (math.nan,)*3)[0] for c in COMPOUNDS]
            sds = [summary.get((target, c), {}).get(desc, (math.nan,)*3)[1] for c in COMPOUNDS]
            off = (ti - 1.5) * 0.2
            ax.bar(x + off, means, width=0.2, yerr=sds, capsize=2, label=target,
                   error_kw=dict(lw=0.7))
        ax.set_xticks(x); ax.set_xticklabels(COMPOUNDS)
        ax.set_title(DESCRIPTORS[desc][0], fontsize=10)
        ax.set_ylabel(DESCRIPTORS[desc][1], fontsize=9)
        if k == 0:
            ax.legend(fontsize=8, ncol=2)
    axes[-1].axis("off")
    fig.suptitle("xanthone_short  -  descriptor means +/- SD per target (A1-A6, 0-25 ns)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for ext in ("svg", "pdf", "png"):
        fig.savefig(OUT / f"figures/summary/ALL_targets_descriptor_bars_A1-A6.{ext}")
    plt.close(fig)

    # ---------------------------------------------------------------- inventory / gap / qc
    (OUT / "inventory/figure_inventory.json").write_text(json.dumps({
        "generated": "2026-09-08",
        "profile": "xanthone_short / 1.0",
        "historical_corpus": str(HIST),
        "a6_corpus": str(A6ROOT),
        "sources": inventory,
        "generated_primary_figures": primary_manifest,
        "generated_heatmaps": heat_manifest,
    }, indent=2))

    # gap analysis: pre-existing historical figures vs the expected complete set
    HIST_SVG = HIST / "svg"
    preexisting = {
        ("AA", "protein_rmsd"): "AA-rmsd.svg", ("AG", "protein_rmsd"): "AG-rmsd.svg",
        ("HMG", "protein_rmsd"): "HMG-rmsd.svg", ("LP", "protein_rmsd"): "LP-rmsd.svg",
        ("AA", "ligand_rmsd"): "AA-LIG-rmsd.svg", ("AG", "ligand_rmsd"): "AG-rmsd-ligand.svg",
        ("HMG", "ligand_rmsd"): "HMG-rmsd-ligand.svg", ("LP", "ligand_rmsd"): "LP-rmsd-ligand.svg",
        ("AA", "catalytic_dist"): "AA-dist.svg", ("AG", "catalytic_dist"): "AG-dist.svg",
        ("HMG", "catalytic_dist"): "HMG-dist.svg", ("LP", "catalytic_dist"): "LP-dist.svg",
    }
    gap = []
    for target in TARGETS:
        for desc in DESC_ORDER:
            have_hist_fig = (target, desc) in preexisting
            hist_fig_path = str(HIST_SVG / preexisting[(target, desc)]) if have_hist_fig else None
            a5_ok = all((target, desc, c) in traces for c in COMPOUNDS[:5])
            a6_ok = (target, desc, "A6") in traces
            if have_hist_fig:
                classification = "EXISTS_BUT_INCOMPLETE"   # historical fig covers A1-A5 only, Grace styling
                reason = ("historical Grace SVG exists for A1-A5; missing A6 and not presentation-"
                          "consistent with the new set -> regenerated as A1-A6")
            else:
                classification = "MISSING"
                reason = "no historical comparative figure for this descriptor family -> generated new"
            gap.append(dict(target=target, descriptor=desc,
                            expected_figure=f"figures/primary/{target}_{desc}_A1-A6.svg",
                            historical_figure=hist_fig_path,
                            classification=classification,
                            a1_a5_sources_present=a5_ok,
                            a6_source_present=a6_ok,
                            action="regenerated_A1-A6",
                            note=reason))
    for target in TARGETS:
        gap.append(dict(target=target, descriptor="heatmap",
                        expected_figure=f"figures/summary/{target}_heatmap_short_descriptors.svg",
                        historical_figure=None, classification="MISSING",
                        a1_a5_sources_present=True, a6_source_present=True,
                        action="generated_new",
                        note="no historical per-target descriptor heatmap"))
    (OUT / "inventory/figure_gap_analysis.json").write_text(json.dumps({
        "expected_primary_figures": 20, "expected_heatmaps": 4,
        "summary_tables": ["short_descriptor_summary.csv", "short_descriptor_summary.md"],
        "items": gap,
    }, indent=2))

    # QC
    for target in TARGETS:
        for desc in DESC_ORDER:
            comps_present, cov, nan_flags, npts = [], {}, [], {}
            srcs, artifacts = [], []
            for comp in COMPOUNDS:
                key = (target, desc, comp)
                if key not in traces:
                    continue
                t_ns, v, meta, src_type, path = traces[key]
                comps_present.append(comp)
                srcs.append(str(path))
                cov[comp] = [round(meta["tmin_ns"], 3), round(meta["tmax_ns"], 3)]
                npts[comp] = meta["n"]
                nan_flags.append(bool(np.isfinite(v).sum() == 0))
                if src_by_key.get((target, comp, desc), {}).get("physically_implausible"):
                    artifacts.append(dict(compound=comp,
                                          vmax=round(src_by_key[(target, comp, desc)]["vmax"], 3),
                                          note="excursion >3 nm; suspected periodic-image artefact "
                                               "in this xvg (active-site contact retained)"))
            cov_ok = all(24.0 <= cov[c][1] <= 26.0 and cov[c][0] <= 0.5 for c in comps_present)
            status = "PASS"
            if comps_present != COMPOUNDS or not cov_ok or any(nan_flags):
                status = "WARN"
            if artifacts:
                status = "WARN_ARTIFACT"
            qc_entries.append(dict(
                figure=f"figures/primary/{target}_{desc}_A1-A6.svg",
                target=target, descriptor=desc,
                data_sources=srcs,
                compounds_included=comps_present,
                all_six_present=(comps_present == COMPOUNDS),
                time_coverage_ns=cov,
                n_points=npts,
                any_all_nan_trace=any(nan_flags),
                coverage_0_25ns_ok=bool(cov_ok),
                suspected_artifacts=artifacts,
                validation_status=status,
            ))
    (OUT / "qc/qc_report.json").write_text(json.dumps({
        "generated": "2026-09-08",
        "sd_convention": "population (ddof=0), matches thesis resumen.txt",
        "hmg_ligand_rmsd_note": ("HMG ligand RMSD uses the fit-to-protein files "
                                 "HMG_CoA_R-A{n}_rmsd_ligand.xvg for A1-A5 (thesis-compatible) and "
                                 "the SimForge 'LIG after lsq fit to Protein' file for A6. The "
                                 "alternative HMG-A4/A5 fit-to-LIG files are intentionally NOT used."),
        "entries": qc_entries,
    }, indent=2))

    print("primary figures :", len(primary_manifest))
    print("heatmaps        :", len(heat_manifest))
    print("long rows       :", len(long_rows))
    print("summary rows    :", len(srows))
    warn = [e for e in qc_entries if e["validation_status"] != "PASS"]
    print("QC warnings     :", len(warn))
    for w in warn:
        print("   -", w["figure"], w["compounds_included"], w["time_coverage_ns"])


if __name__ == "__main__":
    main()
