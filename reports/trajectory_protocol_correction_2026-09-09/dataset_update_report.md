# Dataset update report — LP-A6 ligand RMSD backfill

**Date:** 2026-09-09
**Trigger:** the recompute checkpoint (`checkpoint_report.md`) established that
19/20 A6 short-study observables are numerically identical to 4 dp on `md.xtc`
vs `mdfit.xtc`; only **LP-A6 ligand RMSD** changed materially
(`INVALIDATED_WRONG_TRAJECTORY` → corrected on `mdfit.xtc`).
**Scope:** minimal, evidence-based — one cell per file.

## Files and cells changed

Repo: `/home/jesusxd/Escritorio/ML/xanthones_swissadme_datasets/xanthones_datasets/`

| file | row | column | before | after |
|---|---|---|---|---|
| `md_dataset_v2_A1-A6.csv` | `A6, pancreatic_lipase` (line 25) | `rmsd_ligand_nm` | *(blank)* | **0.24** |
| `master_multimodal_dataset_v2_A1-A6.csv` | `A6, pancreatic_lipase` (line 25) | `rmsd_ligand_nm` | *(blank)* | **0.24** |

`0.24` = mean of `analysis_outputs/xanthone_short/LP-A6/observables/ligand_rmsd/rmsd_ligand.xvg`
(corrected, `mdfit.xtc`, `xanthone_short/1.1`), 0–25 ns, 2501 points
(exact mean 0.24012 nm; median 0.24134; SD 0.05841; max 0.43099), rounded to the
column's 2-dp convention. Same 0–25 ns-mean rule as every other cell in the
column.

## Cells NOT changed (and why)

The other 19 A6 short-MD cells — `rmsd_protein_nm`, `rmsd_ligand_nm`,
`catalytic_distance_nm` for A6 × {alpha_amylase, alpha_glucosidase,
hmgcoa_reductase} and `rmsd_protein_nm` / `catalytic_distance_nm` for
A6/pancreatic_lipase — were left unchanged. Old (`md.xtc`) vs corrected
(`mdfit.xtc`) values are identical to 4 dp (see `short_old_vs_corrected.csv`):

| A6 row / column | dataset value | corrected `mdfit.xtc` mean | matches at 2 dp? |
|---|---|---|---|
| alpha_amylase / rmsd_protein_nm | 0.16 | 0.1553 | yes |
| alpha_amylase / rmsd_ligand_nm | 0.26 | 0.2634 | yes |
| alpha_amylase / catalytic_distance_nm | 1.41 | 1.4139 | yes |
| alpha_glucosidase / rmsd_protein_nm | 0.24 | 0.2404 | yes |
| alpha_glucosidase / rmsd_ligand_nm | 0.88 | 0.8782 | yes |
| alpha_glucosidase / catalytic_distance_nm | 0.91 | 0.9144 | yes |
| hmgcoa_reductase / rmsd_protein_nm | 0.19 | 0.1929 | yes |
| hmgcoa_reductase / rmsd_ligand_nm | 0.52 | 0.5153 | yes |
| hmgcoa_reductase / catalytic_distance_nm | 2.78 | 2.7763 | yes |
| pancreatic_lipase / rmsd_protein_nm | 0.23 | 0.2257 | yes |
| pancreatic_lipase / catalytic_distance_nm | 0.67 | 0.6687 | yes |

## A1–A5 integrity

* `md_dataset_v2_A1-A6.csv`: lines 2–21 (A1–A5 × 4 targets) **byte-identical**
  to the pre-edit file (verified by line diff). Shape 24×5 unchanged.
* `master_multimodal_dataset_v2_A1-A6.csv`: the edit `old_string` spanned
  `A6,pancreatic_lipase,…,Molecule 6,` and inserted `0.24` between two commas —
  the other 23 data rows and all 74 columns are untouched. Shape 24×74 unchanged.
* `master_multimodal_dataset_v1.csv` — not touched.

## QC interpretation updated (dataset repo `reports/`)

The "PBC artefact" interpretation was in report prose, not in the CSVs. Updated:

| file | change |
|---|---|
| `reports/targeted_A6_completion_2026-09-09/qc_and_unresolved_items.md` | Q3 retitled **"RESOLVED 2026-09-09 (`INVALIDATED_WRONG_TRAJECTORY`)"** with a SUPERSEDED banner; original prose preserved below it; blockers-summary row struck through and marked done |
| `reports/targeted_A6_completion_2026-09-09/master_v2_change_log.md` | "short MD" row: note updated (blank → 0.24, with `INVALIDATED_WRONG_TRAJECTORY` provenance); "Known consequences" item 3 struck through and marked RESOLVED |

Provenance preserved throughout: the old value came from **raw `md.xtc`**
(`INVALIDATED_WRONG_TRAJECTORY`); the corrected value comes from the canonical
PBC-corrected rot+trans fit **`mdfit.xtc`**.

## MM-PBSA — NOT touched this pass

`mmpbsa_dataset_v2_A1-A6.csv` and the MM-PBSA cells of
`master_multimodal_dataset_v2_A1-A6.csv` are **unchanged**. The audit conclusion
(`mmpbsa_equivalence_audit.md`) stands and is recorded for a separate dedicated
pass:

* scientific method / configuration is **equivalent** (same `mmpbsa.mdp`
  `954785fb1805…`, same fitted trajectory, same 20–25 ns interval, same PB/SA
  physics / dielectrics / T, same g_mmpbsa build; manual vs automated execution
  is not a comparability criterion);
* **A6 used 26 snapshots (`-dt 200`)** vs **51 snapshots (`-dt 100`) for A1–A5** —
  same window, coarser sampling only;
* the current `mmpbsa_dataset_v2` A1–A5 column contains **cross-batch values that
  do not reproduce from the currently available on-disk folders** (A5 rows,
  A1/lipase — earlier July-2025 charge batch).

MM-PBSA normalisation / recalculation → **separate dedicated pass** (recommended:
single-batch rerun of all 24 systems, then regenerate `mmpbsa_dataset` from
scratch).
