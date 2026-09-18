# xanthone_short (A1–A6) — final figure report

Date: 2026-09-08
Profile: `xanthone_short` / `1.0` (thesis-compatible short-study analysis)
Output root: `reports/xanthone_short_figures/2026-09-08/`
Builder: `scripts/build_short_figures.py` (re-runnable; reads only, writes only under this dir)

## 1. What was done

The complete presentation figure set for the short xanthone series (targets **AA, AG, HMG, LP**;
compounds **A1–A6**; 5 descriptor families; 0–25 ns) was assembled by **reading existing `.xvg`
outputs only**. No MD analysis was rerun, no source simulation file was touched, and no historical
result or figure was overwritten.

* A1–A5 time series come from the thesis corpus
  `Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM/`.
* A6 time series come from the SimForge short run
  `simforge/analysis_outputs/xanthone_short/{AA-A6, AG-A6, HMG-R-25ns-A6, LP-A6}/observables/`.
* All 120 target×compound×descriptor sources resolved and validated (`figure_inventory.json`).

## 2. Reuse vs. regeneration

| category | count | items |
| --- | --- | --- |
| **EXISTS_REUSABLE** (used as-is) | 0 | — |
| **EXISTS_BUT_INCOMPLETE** (historical A1–A5 figure existed; regenerated A1–A6) | 12 | protein RMSD ×4, ligand RMSD ×4, catalytic COM distance ×4 |
| **EXISTS_BUT_NEEDS_REGENERATION** | 0 | (subsumed into the above — historical figures are Grace/xmgrace exports, A1–A5, inconsistent axis units, so treated as incomplete) |
| **MISSING** (no comparative figure existed; generated new) | 12 | active-site min distance ×4, active-site contacts ×4, per-target heatmaps ×4 |

Historical figures that existed and whose **data** was reused (the figures themselves were *not*
edited or replaced; they remain in place):

```
25ns-DM/svg/AA-rmsd.svg  AG-rmsd.svg  HMG-rmsd.svg  LP-rmsd.svg            (protein RMSD, A1–A5)
25ns-DM/svg/AA-LIG-rmsd.svg  AG-rmsd-ligand.svg  HMG-rmsd-ligand.svg  LP-rmsd-ligand.svg  (ligand RMSD, A1–A5)
25ns-DM/svg/AA-dist.svg  AG-dist.svg  HMG-dist.svg  LP-dist.svg            (catalytic COM distance, A1–A5)
```

No historical figure was presentation-consistent with the required set (A1–A5 only, Grace styling,
mixed x-units), so every primary figure was regenerated with one unified style and the A6 trace
added. The historical `.svg` files are untouched.

## 3. Deliverables and exact paths

### A. Inventory / planning — `inventory/`
* `figure_inventory.json` — every source `.xvg`: path, type (historical / A6-simforge), n points,
  time coverage, detected x-unit, header subtitle, mean, SD, plausibility flag.
* `figure_gap_analysis.json` — 24 expected figures classified (12 incomplete + 12 missing) with
  the historical figure referenced and the action taken.
* `figure_generation_plan.md` — the plan, source mapping, filename-normalisation rules, caveats.

### B. Data tables — `tables/`
* `short_descriptor_timeseries_long.csv` — tidy long form, 290 120 rows:
  `target, compound, descriptor, time_ns, time_ps, value, source_type, source_path`
  (~46 MB uncompressed; ~2.6 MB gzipped).
* `short_descriptor_summary.csv` — 24 rows (4 targets × 6 compounds), columns:
  `target, compound, protein_rmsd_mean/sd, ligand_rmsd_mean/sd, active_mindist_mean/sd,
  active_contacts_mean/sd, catalytic_dist_mean/sd`.
* `short_descriptor_summary.md` — same table, human-readable.

SD convention: **population SD (ddof = 0)**, matching the thesis `resumen.txt` /
`extract_stats.sh`. Historical means reproduce `resumen.txt` exactly (spot-checked: AA-A1 protein
0.1533/0.0123, AA-A1 contacts 754.75/105.29, HMG_CoA_R-A1 ligand 0.6550/0.1699).

### C. Figures — `figures/`
`figures/primary/` — **20** comparative figures, each as `.svg` + `.pdf` + `.png`:

```
{AA,AG,HMG,LP}_protein_rmsd_A1-A6
{AA,AG,HMG,LP}_ligand_rmsd_A1-A6
{AA,AG,HMG,LP}_active_mindist_A1-A6
{AA,AG,HMG,LP}_active_contacts_A1-A6
{AA,AG,HMG,LP}_catalytic_dist_A1-A6
```
x = time (ns), 0–25; lines A1–A6 with a fixed colour cycle; A6 drawn bold/black and labelled
"SimForge".

`figures/summary/` — **4** per-target heatmaps + 1 combined bar chart, each `.svg` + `.pdf` + `.png`:

```
{AA,AG,HMG,LP}_heatmap_short_descriptors      rows A1–A6 × cols = 5 descriptors
                                              cell text = raw mean; colour = per-column z-score
                                              across A1–A6 within that target (clipped ±2.5σ)
ALL_targets_descriptor_bars_A1-A6             mean ± SD per target, one panel per descriptor
```

### D. Reports — root of this dir
* `qc/qc_report.json` — per-figure QC (see §5).
* `final_figure_report.md` — this file.

## 4. Descriptor / source mapping notes (filename robustness)

The resolver handles all observed naming inconsistencies:
`AA-A1rmsd_protein.xvg`, `AA-A1_rmsd_ligand.xvg`, `LP-A2_rmsd-ligand.xvg`, `AG-A4rmsd_ligand.xvg`,
`HMG_CoA_R-A1_rmsd_ligand.xvg`, hyphen/underscore variants, and x-axis values recorded in **ps**
(`25000`) or **ns** (`25.0`) — auto-detected per file and normalised to ns.

## 5. Quality control (`qc/qc_report.json`)

| check | result |
| --- | --- |
| all six compounds present in every primary figure | **20 / 20 PASS** |
| time coverage 0–25 ns for every trace | **PASS** (all traces tmin ≤ 0.1 ns, tmax = 25.0 ns) |
| source files readable | **PASS** (120 / 120) |
| no NaN-only traces | **PASS** |
| target / compound labels correct | **PASS** (validated against `.xvg` headers and working-dir strings) |
| HMG uses the short-analysis window | **PASS** — `HMG-R-25ns-A6` (0–25 ns), not the 200 ns run |
| point count = 2501 | 115 / 120 sources; the 5 exceptions are **AG-A4** (501 points, 50 ps stride) |
| physical-plausibility guard | **1 WARN_ARTIFACT**: `LP_ligand_rmsd_A1-A6` (LP-A6) |

Overall: **19 figures PASS, 1 WARN_ARTIFACT** (LP-A6 ligand RMSD, described below). No figure fails.

## 6. Caveats

### HMG-specific
1. **Ligand RMSD reference.** HMG ligand RMSD is taken from the fit-to-protein files
   `HMG_CoA_R-A{1..5}_rmsd_ligand.xvg` (subtitle *"LIG after lsq fit to Protein"*) for A1–A5, and
   from the SimForge A6 file with the same *"LIG after lsq fit to Protein"* fitting. The alternative
   thesis files `HMG-A4rmsd_ligand.xvg` / `HMG-A5rmsd_ligand.xvg` are *fit to ligand* and give very
   different values (e.g. A4: 0.18 nm vs 0.66 nm); they are **deliberately excluded** so all four
   targets and all six compounds use one consistent definition.
2. **Analysis window.** The HMG target for this set is the 25 ns system (`HMG-R-25ns-A6`). The
   200 ns system (`HMG-R-200ns-A6`) has no short-profile observables and is not used here.
3. HMG A1–A3 ligand files store time in ps, A4–A5 in ns (both 2501 points, 0–25 ns) — normalised.

### AG-specific
4. **AG-A4** historical `.xvg` files are 501 points (50 ps stride) for all five descriptors, versus
   2501 points (10 ps) elsewhere. Full 0–25 ns coverage; only the time resolution differs. The A6
   AG-A4… n/a (A6 is a single compound); AG-A6 is 2501 points.

### LP-specific
5. **LP-A6 ligand RMSD** rises to ~7.7 nm for most of the trajectory while the ligand–active-site
   **minimum distance stays ~0.2 nm and contacts stay ~1000** — i.e. the ligand is still bound.
   This inconsistency indicates a **periodic-image (PBC) artefact in that single `rmsd_ligand.xvg`**,
   not ligand egress. Per the "do not rerun analyses" constraint the file was not regenerated; the
   figure plots the trace with the off-scale portion masked and a QC caption, the heatmap prints
   the raw mean with a footnote, and `qc_report.json` records it as `WARN_ARTIFACT`. **Recommended
   follow-up (not done here): regenerate LP-A6 `rmsd_ligand.xvg` with `gmx trjconv -pbc mol
   -center` (or `-pbc nojump`) preprocessing, then re-plot `LP_ligand_rmsd_A1-A6`.**

## 7. Completeness statement

The short-study presentation figure set is **complete**:

* 20 / 20 primary comparative figures generated, each covering A1–A6 over 0–25 ns, in SVG + PDF + PNG.
* 4 / 4 per-target descriptor heatmaps generated.
* Combined summary table generated in CSV and Markdown, plus a tidy long-form timeseries CSV.
* 1 combined bar figure as an optional presentation aid.

One figure (`LP_ligand_rmsd_A1-A6`) carries a documented data-quality caveat (LP-A6 PBC artefact)
that is visualised honestly and flagged; it does not block presentation use of the set. Everything
else passes QC.

Mechanistic-profile figure generation was **not** started — out of scope for this phase.
