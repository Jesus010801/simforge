# Figure update report — LP-A6 ligand RMSD

**Date:** 2026-09-09
**Scope:** minimal — only the one thesis figure whose underlying data changed.

## What changed

| figure | status | why |
|---|---|---|
| `LP-LIG-rmsd-A1-A6.{svg,png,pdf}` (+ `svg/` copy) | **regenerated** | the A6 ligand-RMSD trace changed (raw `md.xtc` ~7.7 nm artefact → corrected `mdfit.xtc`, mean 0.24 nm) |
| the other **19** figures in `completadas_A6/` (`{AA,AG,HMG,LP}-{rmsd,dist,mindist,contacts}-A1-A6` and `{AA,AG,HMG}-LIG-rmsd-A1-A6`) | **NOT regenerated** | their underlying A6 data is identical to 4 dp on `md.xtc` vs `mdfit.xtc` (`short_old_vs_corrected.csv`); A1–A5 data unchanged |
| any aggregate / summary figure | **none exist** in `completadas_A6/` — all 20 are per-target × per-descriptor time series; nothing derives from ligand RMSD |

## Output location

New directory (originals preserved, never overwritten):

```
Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM/completadas_A6_corrected_2026-09-09/
├── LP-LIG-rmsd-A1-A6.svg / .png / .pdf
├── svg/LP-LIG-rmsd-A1-A6.svg
├── qc_report_LP-LIG-rmsd.json
└── regen_lp_lig_rmsd.py            (the generator — reuses the original build_thesis_style_A6.py verbatim)
```

`completadas_A6/SUPERSEDED_LP-LIG-rmsd.md` points to it. The original
`completadas_A6/LP-LIG-rmsd-A1-A6.*` files are **byte-identical** before and
after (sha256 verified).

## Preservation of A1–A5 traces

The regenerator resolves A1–A5 from the **same unchanged historical `.xvg`
files** (`LP-A1_rmsd_ligand.xvg` … `LP-A5_rmsd_ligand.xvg`) via the original
`resolve()` / `parse_xvg()` helpers, with the original colour map
(`A1 #FF00FF · A2 #FFFF00 · A3 #FFA500 · A4 #00FFFF · A5 #7221BC`), line width
(1.377 pt) and Grace style. Verified: A1–A5 source paths, point counts (2501)
and per-trace maxima are identical between the original and corrected figure
inventories. Only the A6 (`#00A000`) trace differs (vmax 7.7115 → 0.4310).

## Visible difference

| | original (`completadas_A6/`) | corrected (`completadas_A6_corrected_2026-09-09/`) |
|---|---|---|
| A6 trace | off-scale ~7.7 nm, **masked** above the clip line | normal, on-scale, ~0.24 nm |
| legend | `A6 *` (asterisk) | `A6` |
| caption | red "* A6: excursión a ~7.7 nm fuera de escala (artefacto PBC sospechado…)" | none |
| y-axis | clipped to ~1.3× the A1–A5 99.5th percentile | data-driven, no clip |
| QC | `WARN_ARTIFACT` | **`PASS`** |

The A1–A5 curves occupy different pixel positions **only** because the y-axis is
no longer clipped and the artefact caption no longer reserves bottom margin —
their data values are unchanged.

## QC / style reports

* `completadas_A6_corrected_2026-09-09/qc_report_LP-LIG-rmsd.json` — new: PASS,
  A6 vmax 0.431 nm, no artefacts, 0–25 ns coverage, 2501 pts all six compounds.
* The original `completadas_A6/qc_report.json` / `style_audit.md` /
  `figure_generation_report.md` are **left as historical record** (they
  correctly describe the 2026-09-08 build); `SUPERSEDED_LP-LIG-rmsd.md` is the
  pointer to the correction.
