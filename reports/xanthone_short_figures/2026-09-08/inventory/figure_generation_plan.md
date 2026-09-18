# xanthone_short (A1–A6) — figure generation plan

Date: 2026-09-08
Profile: `xanthone_short` / `1.0` (thesis-compatible)
Scope: presentation figure set for the **short** xanthone series only. Mechanistic-profile
figures are explicitly **out of scope** for this phase.

## Data sources (read-only)

| role | path | compounds |
| --- | --- | --- |
| historical thesis corpus | `Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM/` | A1–A5 |
| new SimForge short outputs | `simforge/analysis_outputs/xanthone_short/{AA-A6,AG-A6,HMG-R-25ns-A6,LP-A6}/observables/` | A6 |

No scientific MD analysis was rerun. Figures and statistics are derived only by reading the
existing `.xvg` files. No source or historical file is modified or overwritten.

## Targets × descriptors

Targets: **AA** (α-amylase), **AG** (α-glucosidase), **HMG** (HMG-CoA reductase), **LP** (lipase).
Compounds: **A1–A6**.
Descriptor families (5):

| key | descriptor | y-unit | historical file stem | A6 file |
| --- | --- | --- | --- | --- |
| `protein_rmsd` | protein RMSD (fit to protein) | nm | `{T}-{C}rmsd_protein.xvg` | `protein_rmsd/rmsd_protein.xvg` |
| `ligand_rmsd` | ligand RMSD, **fit to protein** | nm | `{T}-{C}_rmsd_ligand.xvg` / `HMG_CoA_R-{C}_rmsd_ligand.xvg` | `ligand_rmsd/rmsd_ligand.xvg` |
| `active_mindist` | ligand–active-site minimum distance | nm | `{T}-{C}mindist_lig_active.xvg` | `active_site/mindist_lig_active.xvg` |
| `active_contacts` | ligand–active-site contacts (<0.6 nm) | count | `{T}-{C}contacts_lig_active.xvg` | `active_site/contacts_lig_active.xvg` |
| `catalytic_dist` | ligand–catalytic-site COM distance | nm | `{T}-{C}dist_lig_catalytic.xvg` | `catalytic_com_distance/dist_lig_catalytic.xvg` |

Filename normalisation handled by the resolver (`build_short_figures.py :: hist_candidates`):
hyphen/underscore variants (`LP-A2_rmsd-ligand.xvg`), the `AG-A4rmsd_ligand.xvg` form, and the
HMG-specific `HMG_CoA_R-*` ligand files. XVG identity is validated from the header
`subtitle` / axis labels, and x-units (ps vs ns) are auto-detected per file.

## Expected complete set

* 20 primary comparative figures = 4 targets × 5 descriptors, each an A1–A6 time series (x = time ns).
* 4 per-target summary heatmaps (rows A1–A6, columns = 5 descriptors, colour = per-column
  z-score across the six compounds within that target, cell text = raw mean).
* 1 combined summary table (`short_descriptor_summary.csv` + `.md`).
* 1 optional combined bar figure (means ± SD per target).

## Classification (vs. gap analysis)

| descriptor family | historical figure | classification | action |
| --- | --- | --- | --- |
| protein RMSD (×4) | `svg/{AA,AG,HMG,LP}-rmsd.svg` (Grace, A1–A5) | EXISTS_BUT_INCOMPLETE | regenerate A1–A6, unified style |
| ligand RMSD (×4) | `svg/{AA-LIG-rmsd, AG-rmsd-ligand, HMG-rmsd-ligand, LP-rmsd-ligand}.svg` (A1–A5) | EXISTS_BUT_INCOMPLETE | regenerate A1–A6 |
| catalytic COM distance (×4) | `svg/{AA,AG,HMG,LP}-dist.svg` (A1–A5) | EXISTS_BUT_INCOMPLETE | regenerate A1–A6 |
| active-site min distance (×4) | none | MISSING | generate new |
| active-site contacts (×4) | none | MISSING | generate new |
| per-target heatmaps (×4) | none | MISSING | generate new |

No historical figure is presentation-consistent with the new set (Grace exports, A1–A5 only,
mixed axis units), so all 12 incomplete figures are regenerated rather than edited in place.
The historical `.svg` files are left untouched.

## Output layout

```
reports/xanthone_short_figures/2026-09-08/
├── figures/primary/    20 × {svg,pdf,png}   {TARGET}_{descriptor}_A1-A6
├── figures/summary/    4 heatmaps + 1 bar chart × {svg,pdf,png}
├── tables/             short_descriptor_timeseries_long.csv, short_descriptor_summary.{csv,md}
├── inventory/          figure_inventory.json, figure_gap_analysis.json, figure_generation_plan.md
├── qc/                 qc_report.json
├── scripts/            build_short_figures.py
└── final_figure_report.md
```

## Conventions

* x-axis: time in **ns**, 0–25.
* SD: population standard deviation (ddof = 0), matching the thesis `resumen.txt` / `extract_stats.sh`.
* colour cycle fixed A1–A6 (A6 = black, bold, labelled "SimForge").
* vector output first (SVG + PDF), PNG for convenience.

## Known caveats to carry into the report

1. **HMG ligand RMSD** uses the fit-to-protein files (`HMG_CoA_R-A{n}_rmsd_ligand.xvg` for A1–A5;
   SimForge "LIG after lsq fit to Protein" for A6). The alternative `HMG-A4rmsd_ligand.xvg` /
   `HMG-A5rmsd_ligand.xvg` are fit-to-**ligand** and are deliberately not used.
2. **HMG target = 25 ns window** (`HMG-R-25ns-A6`), not the 200 ns run (`HMG-R-200ns-A6`).
3. **AG-A4** historical files are 501 points (50 ps stride) vs 2501 elsewhere — full 0–25 ns
   coverage, lower time resolution only.
4. **LP-A6 ligand RMSD** shows an excursion to ~7.7 nm while active-site contact is retained —
   flagged as a suspected periodic-image artefact in that one xvg (QC `WARN_ARTIFACT`).
