# Trajectory protocol correction — final report

**Date:** 2026-09-09
**Canonical analysis trajectory:** `mdfit.xtc`
(`md.xtc` → `gmx trjconv -pbc res -ur compact -center` → `mdcenter.xtc` →
`gmx trjconv -fit rot+trans` → `mdfit.xtc`). Raw `md.xtc` is the production
trajectory, never a scientific-observable input.

This report covers the `xanthone_short` A6 short-study correction end to end.
The mechanistic A1/A3/A6 suite is **planned, not executed** (separate authorised
phase). MM-PBSA: the canonical A1–A5 source has been **recovered** (thesis
heatmap) and shown **methodologically equivalent to A6** — **no rerun**
(§5, `mmpbsa_rerun_decision.md`).

---

## 1. Which short outputs were invalidated

**All 20** of the 2026-09-08 `xanthone_short` A6 outputs — 4 systems
(`AA-A6`, `AG-A6`, `HMG-R-25ns-A6`, `LP-A6`) × 5 observables (protein RMSD,
ligand RMSD, active-site minimum distance, active-site contacts, catalytic-site
COM distance). Every provenance record showed `"source_trajectory": ".../md.xtc"`.
**Verdict: `INVALIDATED_WRONG_TRAJECTORY` ×20.** Archived verbatim (nothing
deleted) to
`analysis_outputs/xanthone_short/_invalidated_md_xtc_2026-09-08/` + `INVALIDATED.md`.
`HMG-R-25ns-A6` was additionally pointed at the full 200 ns raw file.

## 2. Which were recomputed

All 20, on `mdfit.xtc`, protocol `xanthone_short/1.1`, via
`simforge study run Nuevos_sistemas --profile xanthone_short --force`:
**20 EXECUTED, 0 FAILED.** Each `mdfit.xtc` verified 2501 frames / 0–25 000 ps /
dt 10 ps; each output XVG 2501 points, 0–25 000 ps, all finite, 10 ps sampling,
valid semantics — **QC PASS ×20**. Window (`-b 0 -e 25000`), contact cutoff
(`-d 0.6`), catalytic COM semantics and name-based group resolution unchanged.

## 3. Corrected A6 statistics

Full table + CSV: `checkpoint_report.md` §5, `short_old_vs_corrected.csv`.

**19 of 20 observables are numerically identical to 4 decimal places** between
`md.xtc` and `mdfit.xtc` — `gmx rms` re-fits before computing RMSD, and
`gmx mindist` / `gmx distance` use the minimum image, so the raw PBC state does
not affect them. Corrected A6 means:

| system | protein RMSD | ligand RMSD | active-site min dist | active-site contacts | catalytic COM dist |
|---|---|---|---|---|---|
| AA-A6 | 0.155 | 0.263 | 0.223 | 678.1 | 1.414 |
| AG-A6 | 0.240 | 0.878 | 0.237 | 475.9 | 0.914 |
| HMG-R-25ns-A6 | 0.193 | 0.515 | 0.241 | 282.5 | 2.776 |
| LP-A6 | 0.226 | **0.240** | 0.204 | 1084.3 | 0.669 |

(nm, except contacts = count; 0–25 ns means, 2501 frames.)

## 4. LP-A6 ligand RMSD conclusion

| | old (`md.xtc`) | corrected (`mdfit.xtc`) |
|---|---|---|
| mean / median / max | 4.84 / 7.46 / 7.71 nm | **0.240 / 0.241 / 0.431 nm** |
| shape | bimodal, ~7.7 nm plateau | single well-behaved basin |

The ~7.7 nm excursion **disappears entirely**. In the same corrected run the
ligand is bound throughout (`mindist` 0.204 nm, `contacts` 1084). It was a
**periodic-image artefact of analysing the raw `md.xtc`** — the ligand had
jumped a box vector in the unwrapped production trajectory.

**Verdict:** reclassified from *"suspected PBC artefact"* to
**`INVALIDATED_WRONG_TRAJECTORY`**. The corrected value is a real, populatable
number.

## 5. MM-PBSA equivalence conclusion  *(revised 2026-09-09 after canonical-source recovery)*

Full reports: `mmpbsa_canonical_source_recovery.md`, `mmpbsa_rerun_decision.md`
(these supersede the "rerun 24" recommendation in `mmpbsa_equivalence_audit.md`).
**No MM-PBSA rerun; no dataset change.**

* **Canonical A1–A5 source recovered:** `/home/jesusxd/Escritorio/heatmap.py`
  (2026-05-12) → `Tesis-maestría/Heat.pdf` → thesis Fig. `heatmap_mmpbsa` →
  `mmpbsa_dataset.csv` / `_v2` — **all 20 A1–A5 values match to 1 dp**; controls
  match `Sustratos-inhibidores/*` (6/7).
* **Methodology documented and equivalent to A6:** thesis Methods —
  LigParGen/OPLS-AA `1.14*CM1A-LBCC`, SPC/E, 309.65 K, 25 ns, MM-PBSA over the
  **last 5 ns at 200 ps intervals** (~26 configs). A6 uses `-b 20000 -e 25000
  -dt 200` (26 configs). `mmpbsa.mdp` is physically identical everywhere (the two
  sha256 differ only by one comment character).
* **The disagreeing `/Escritorio/MMPBSA/` folders are a discarded preliminary
  `-dt 100` batch** (Jul-2025). 10/20 canonical values reproduce from a surviving
  folder; 7/20 (A5×4, HMG-A2, HMG-A4, LP-A1) were re-run for the thesis with the
  final folders not retained; 3/20 (A1/AG, A2/AG, A4/AG) have no folder anywhere.
  This is a **records-retention gap, not a methodological incompatibility.**
* **200 ns mechanistic MM-PBSA** (last 50 ns, ~250 configs) is a separate
  protocol; poorly-converged totals stay for per-residue interpretation only.
* **Decision:** retain A1–A5, keep A6, update provenance in docs only.
  **No 24-system rerun.** `mmpbsa_dataset_v2` / `master_v2` MM-PBSA cells already
  hold the canonical values and are unchanged.

## 6. Datasets updated

Repo `/home/jesusxd/Escritorio/ML/xanthones_swissadme_datasets/xanthones_datasets/`
— full detail in `dataset_update_report.md`.

| file | cell | before → after |
|---|---|---|
| `md_dataset_v2_A1-A6.csv` | `A6, pancreatic_lipase` · `rmsd_ligand_nm` | blank → **0.24** |
| `master_multimodal_dataset_v2_A1-A6.csv` | `A6, pancreatic_lipase` · `rmsd_ligand_nm` | blank → **0.24** |

Post-edit hashes:
`md_dataset` `564e4af6…`, `master` `9f2ff247…`. Shapes 24×5 / 24×74 unchanged.
A1–A5 rows byte-identical. `master_multimodal_dataset_v1.csv` untouched.
`mmpbsa_dataset_v2_A1-A6.csv` **unchanged.**

QC prose updated: `qc_and_unresolved_items.md` Q3 (→ RESOLVED,
`INVALIDATED_WRONG_TRAJECTORY`) and `master_v2_change_log.md` (short-MD note +
consequence 3). Provenance preserved: old value from raw `md.xtc`, corrected
value from canonical `mdfit.xtc`.

## 7. Figures updated

Full detail in `figure_update_report.md`.

* **Regenerated:** `LP-LIG-rmsd-A1-A6.{svg,png,pdf}` (+ `svg/` copy) into a new
  dir `Nuevos_sistemas/Resultados-Tesis-maestría/25ns-DM/completadas_A6_corrected_2026-09-09/`.
  QC `WARN_ARTIFACT` → **`PASS`**; A6 trace now on-scale (~0.24 nm); no mask, no
  caption, legend `A6` (was `A6 *`).
* **Not regenerated:** the other 19 `completadas_A6/` figures — their A6 data is
  identical to 4 dp; no aggregate/summary figure derives from ligand RMSD.
* **A1–A5 traces preserved exactly** — same historical `.xvg` sources, colours
  and style (verified: identical source paths, point counts, per-trace maxima).
* Originals in `completadas_A6/` **byte-identical before/after** (sha256
  verified); `completadas_A6/SUPERSEDED_LP-LIG-rmsd.md` points to the correction.

## 8. Standardized mechanistic protocol (planned)

`mechanistic_standardized_protocol.yaml`: one trajectory treatment
(`mdfit.xtc` for **all** systems incl. reprocessed `A1/APO/COA/COMP` references),
0–200 ns window / 100 ps stride, common atom selections, protein-backbone fit;
a **shared PCA basis** for A1/A3/A6 (single covariance on the concatenated,
commonly-aligned Cα trajectories; each system projected onto the shared
eigenvectors) and a **shared PC1/PC2 FEL** for direct comparison, plus
per-system historically-compatible FEL. Clustering and conformational-state
populations are marked **NEW — not historically present for A1**.

New command `simforge study mechanistic-plan <dir> [--reference <dir>] [--json]`
(planning only, never executes) + `tests/analysis/campaign/test_mechanistic_plan.py`.

## 9. What remains to execute for A3/A6

Dry-run: `mechanistic_gap_and_dryrun.md` + `mechanistic_dryrun.json`.

1. Reprocess `Mecanismo_inhibitorio/{A1,APO,COA,COMP}` to `mdfit.xtc`.
2. Recompute 12 descriptors on `mdfit.xtc` for A1/APO/COA/COMP
   (`SUPERSEDED_PENDING` → done): protein RMSD, ligand RMSD, RMSF, Rg, SASA,
   H-bonds, active-site min distance, active-site contacts, catalytic COM
   distance, PCA, FEL, DCCM.
3. Compute the same 12 for the new systems (`A3-HMG-R`, `HMG-R-200ns-A6`,
   `system_A3_COA`, `system_A6_COA`) — all `PLANNED`.
4. MM-PBSA + per-residue: already `AVAILABLE` on `mdfit.xtc` for A3-HMG-R /
   HMG-R-200ns-A6; equalise stride vs the A1 reference if a quantitative
   comparison is needed. 200 ns totals for A3/A6 are poorly converged
   (SD > |mean|) → per-residue interpretation only.
5. Build the shared PCA basis + shared PC1/PC2 FEL for A1/A3/A6.
6. NEW extensions: `gmx cluster` + conformational-state populations.

**No mechanistic GROMACS analysis executed in this work.**

## 10. Test results

| suite | result |
|---|---|
| focused (`test_xanthone_short` + `test_legacy_commands` + `test_mechanistic_plan` + `test_legacy`) | 56 passed |
| campaign (`tests/analysis/campaign/`) | 326 passed (was 315; +11 new) |
| full repository (`pytest -q`) | **2549 passed, 29 skipped, 9 xfailed, 0 failed** (re-run after the dataset + figure edits) |

New regression tests enforce: short analyses require `mdfit.xtc`; mechanistic
analyses require `mdfit.xtc`; no fallback to raw `md.xtc`; missing `mdfit.xtc`
blocks execution; the analysis window is independent of production length;
semantic index roles stay dynamic (by name); corrected results supersede
invalidated ones; source files remain unchanged; MM-PBSA manual provenance can
be scientifically equivalent when parameters match; manual-vs-automated is not
itself a comparability criterion.

## 11. Source-integrity verification

* `run_report.json` → `source_integrity.verdict` =
  `{ok: true, added: [], removed: [], modified: [], index_changed: []}`.
* Independent `sha256` of every `md.xtc` / `md.tpr` / `index.ndx` for the 4
  systems — **all match** the fingerprints in the archived 2026-09-08 provenance
  (HMG-R-25ns-A6/md.xtc > 2 GiB: size + mtime match).
* No raw trajectory, `.tpr` or source `index.ndx` was modified anywhere.
* Writes: `simforge/analysis/` + `tests/` + `cli.py` (code) ·
  `simforge/analysis_outputs/xanthone_short/` (results + archive) ·
  `simforge/reports/trajectory_protocol_correction_2026-09-09/` ·
  the two ML dataset CSVs (one cell each) + their QC prose + a new dated report ·
  `.../25ns-DM/completadas_A6_corrected_2026-09-09/` (one figure) +
  `completadas_A6/SUPERSEDED_LP-LIG-rmsd.md`.
