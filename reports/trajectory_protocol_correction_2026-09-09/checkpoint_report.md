# Trajectory protocol correction — recompute checkpoint

**Date:** 2026-09-09
**Pass scope:** provenance audit · SimForge patch enforcing `mdfit.xtc` ·
regression tests · recompute the 20 A6 short-study observables on `mdfit.xtc` ·
archive/invalidate the old `md.xtc` outputs · validate corrected XVGs ·
old-vs-corrected statistics · LP-A6 ligand-RMSD verdict · source-integrity
verification · MM-PBSA provenance/config audit (report only) · mechanistic
protocol + dry-run (prepared, not executed).

**STOP after this checkpoint.** No ML dataset, figure, MM-PBSA dataset or
mechanistic analysis was touched.

---

## 1. Provenance audit — all 20 previous outputs used raw `md.xtc`

Every one of the 20 A6 short-study outputs from 2026-09-08
(`AA-A6`, `AG-A6`, `HMG-R-25ns-A6`, `LP-A6` × {protein RMSD, ligand RMSD,
active-site min distance, active-site contacts, catalytic-site COM distance})
records `"source_trajectory": ".../md.xtc"`.

**Verdict: all 20 = `INVALIDATED_WRONG_TRAJECTORY`.** HMG-R-25ns-A6 additionally
had the analysis pointed at the full 200 ns raw file
(`trajectory_measured_end_ps = 200000`).

Full table: `short_analysis_provenance_audit.md`.

## 2. SimForge patch (protocol `xanthone_short/1.0` → `1.1`)

| file | change |
|---|---|
| `analysis/campaign/xanthone_short.py` | new `_fit_trajectory()` resolves the analysis trajectory to `mdfit.xtc` only (discovered DERIVED artifact named `mdfit.xtc`, else `<dir>/mdfit.xtc`); **never** `md.xtc` / `mdcenter.xtc` / `*.trr`. `run_system()` binds it; if absent → every analysis `REVIEW_REQUIRED` ("mdfit.xtc not found … raw md.xtc is not an acceptable analysis input"). `_coverage_ok()` measures **`mdfit.xtc`** span directly via `gmx check` (falls back to discovery inspection only when gmx is unavailable / dry-run). Provenance now carries `source_trajectory = mdfit.xtc`, `analysis_trajectory_role`, `raw_production_trajectory` (recorded, unused), `preprocessing_chain`. `assess_eligibility()` reports a missing `mdfit.xtc`. |
| `analysis/campaign/legacy_commands.py` | `build_xanthone_commands` binds the five analysis rows to `mdfit.xtc`; missing → `REVIEW_REQUIRED` naming `mdfit.xtc`; the `-pbc`/`-fit` preprocessing rows still consume the raw `md.xtc`. |
| `analysis/campaign/mechanistic.py` (new) | `plan_mechanistic()` — dry-run planner for the standardized mechanistic profile; never executes GROMACS. |
| `analysis/campaign/cli.py`, `cli.py` | new `simforge study mechanistic-plan <dir> [--reference <dir>] [--json]`. |

**Not changed:** the trajectory stage classifier (correct — `mdfit.xtc` *is* a
derived trajectory), discovery, the legacy XVG `study analyze` path, the analysis
window (`-b 0 -e 25000`), contact cutoff (`-d 0.6`), catalytic COM semantics, and
name-based semantic-group resolution.

### Dry-run against the real dataset

`simforge study run Nuevos_sistemas --profile xanthone_short --dry-run` →
all 20 analyses for the 4 A6 systems `PLANNED`, each argv `-f` = `mdfit.xtc`;
the other 4 systems stay `REVIEW_REQUIRED` (identity/cofactor, unchanged);
source integrity `UNCHANGED`.

## 3. Regression tests

New tests (`tests/analysis/campaign/`):

* `test_xanthone_short.py`: analysis binds `mdfit.xtc` not raw `md.xtc`;
  missing `mdfit.xtc` → all `REVIEW_REQUIRED` and **no** gmx analysis call is
  made; provenance records fit + raw separately with `protocol_version 1.1`;
  window stays `-b 0 -e 25000` on a 200 ns fitted trajectory; semantic index
  roles still resolve by name.
* `test_legacy_commands.py`: analysis rows require `mdfit.xtc`, never fall back
  to `md.xtc`; adding `mdfit.xtc` flips them to `PLANNED` bound to `mdfit.xtc`;
  preprocessing still consumes raw `md.xtc`.
* `test_mechanistic_plan.py` (new file): planner never executes; uses
  `mdfit.xtc`; missing `mdfit.xtc` → `REVIEW_REQUIRED`; clustering /
  state-populations are `NEW_EXTENSION`; historical raw output → `SUPERSEDED_PENDING`;
  profile independent of production length.

Results:

* focused (`test_xanthone_short` + `test_legacy_commands` + `test_mechanistic_plan` + `test_legacy`): **56 passed**
* campaign suite (`tests/analysis/campaign/`): **326 passed** (was 315; +11 new)
* full repository suite (`pytest -q`): **2549 passed, 29 skipped, 9 xfailed, 0 failed** (6:06)

## 4. Recompute — 20 observables on `mdfit.xtc`

Old outputs copied verbatim to
`analysis_outputs/xanthone_short/_invalidated_md_xtc_2026-09-08/` + `INVALIDATED.md`
(nothing deleted). `simforge study run Nuevos_sistemas --profile xanthone_short
--force` → **20 EXECUTED, 0 FAILED**. New provenance: `xanthone_short/1.1`,
`source_trajectory` ends `mdfit.xtc`, every XVG 2501 pts, 0–25000 ps, all finite,
10 ps sampling.

All four `mdfit.xtc` inputs verified 2501 frames / 0–25000 ps / dt 10 ps
(`gmx check`). HMG-R-25ns-A6's `mdfit.xtc` is already the 0–25 ns window (not the
200 ns raw file).

## 5. Old vs corrected statistics (all 20)

CSV: `short_old_vs_corrected.csv`. **19 of 20 observables are numerically
identical to 4 dp** between `md.xtc` and `mdfit.xtc` — `gmx rms` re-fits and
`gmx mindist` / `gmx distance` use the minimum image, so the raw PBC state does
not affect them. **Only LP-A6 ligand RMSD changed materially.**

| system | analysis | old mean | corr mean | old SD | corr SD | old median | corr median | old max | corr max | Δmean | %Δ | pts | QC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AA-A6 | protein RMSD | 0.1553 | 0.1553 | 0.0106 | 0.0106 | 0.1563 | 0.1563 | 0.1823 | 0.1823 | 0.0000 | ~0 | 2501 | PASS |
| AA-A6 | ligand RMSD | 0.2635 | 0.2634 | 0.0516 | 0.0516 | 0.2585 | 0.2584 | 0.5215 | 0.5216 | −0.0000 | ~0 | 2501 | PASS |
| AA-A6 | active-site min dist | 0.2229 | 0.2229 | 0.0137 | 0.0137 | 0.2219 | 0.2220 | 0.2870 | 0.2877 | 0.0000 | ~0 | 2501 | PASS |
| AA-A6 | active-site contacts | 678.09 | 678.10 | 113.40 | 113.39 | 701 | 700 | 925 | 925 | +0.01 | ~0 | 2501 | PASS |
| AA-A6 | catalytic COM dist | 1.4139 | 1.4139 | 0.0906 | 0.0906 | 1.407 | 1.407 | 1.745 | 1.745 | 0.0000 | ~0 | 2501 | PASS |
| AG-A6 | protein RMSD | 0.2404 | 0.2404 | 0.0419 | 0.0419 | 0.2431 | 0.2431 | 0.3360 | 0.3360 | 0.0000 | ~0 | 2501 | PASS |
| AG-A6 | ligand RMSD | 0.8781 | 0.8782 | 0.3078 | 0.3078 | 0.9958 | 0.9960 | 1.5035 | 1.5036 | +0.0001 | ~0 | 2501 | PASS |
| AG-A6 | active-site min dist | 0.2369 | 0.2369 | 0.0256 | 0.0256 | 0.2361 | 0.2361 | 0.6355 | 0.6349 | −0.0000 | ~0 | 2501 | PASS |
| AG-A6 | active-site contacts | 475.83 | 475.85 | 121.25 | 121.25 | 510 | 510 | 718 | 719 | +0.02 | ~0 | 2501 | PASS |
| AG-A6 | catalytic COM dist | 0.9144 | 0.9144 | 0.1319 | 0.1319 | 0.924 | 0.924 | 1.624 | 1.624 | −0.0000 | ~0 | 2501 | PASS |
| HMG-R-25ns-A6 | protein RMSD | 0.1929 | 0.1929 | 0.0234 | 0.0234 | 0.1969 | 0.1969 | 0.2334 | 0.2334 | 0.0000 | ~0 | 2501 | PASS |
| HMG-R-25ns-A6 | ligand RMSD | 0.5153 | 0.5153 | 0.0783 | 0.0783 | 0.5277 | 0.5277 | 0.6637 | 0.6637 | −0.0000 | ~0 | 2501 | PASS |
| HMG-R-25ns-A6 | active-site min dist | 0.2413 | 0.2413 | 0.0174 | 0.0174 | 0.2411 | 0.2411 | 0.3426 | 0.3426 | −0.0000 | ~0 | 2501 | PASS |
| HMG-R-25ns-A6 | active-site contacts | 282.52 | 282.53 | 35.86 | 35.88 | 284 | 284 | 428 | 428 | +0.01 | ~0 | 2501 | PASS |
| HMG-R-25ns-A6 | catalytic COM dist | 2.7763 | 2.7763 | 0.0442 | 0.0442 | 2.783 | 2.783 | 2.887 | 2.887 | 0.0000 | ~0 | 2501 | PASS |
| LP-A6 | protein RMSD | 0.2257 | 0.2257 | 0.0252 | 0.0252 | 0.2252 | 0.2252 | 0.3178 | 0.3178 | −0.0000 | ~0 | 2501 | PASS |
| **LP-A6** | **ligand RMSD** | **4.8364** | **0.2401** | **3.5226** | **0.0584** | **7.4578** | **0.2413** | **7.7115** | **0.4310** | **−4.5963** | **−95.0%** | 2501 | PASS |
| LP-A6 | active-site min dist | 0.2040 | 0.2040 | 0.0166 | 0.0166 | 0.2035 | 0.2036 | 0.2595 | 0.2597 | −0.0000 | ~0 | 2501 | PASS |
| LP-A6 | active-site contacts | 1084.33 | 1084.26 | 94.21 | 94.23 | 1078 | 1078 | 1379 | 1379 | −0.06 | ~0 | 2501 | PASS |
| LP-A6 | catalytic COM dist | 0.6687 | 0.6687 | 0.0411 | 0.0411 | 0.665 | 0.665 | 0.855 | 0.855 | 0.0000 | ~0 | 2501 | PASS |

All corrected time ranges 0–25000 ps; all point counts 2501; all values finite;
all XVG semantics valid. QC = PASS for all 20.

## 6. LP-A6 ligand RMSD — scientific verdict

| | old (`md.xtc`) | corrected (`mdfit.xtc`) |
|---|---|---|
| mean | 4.84 nm | **0.240 nm** |
| SD | 3.52 nm | 0.058 nm |
| median | 7.46 nm | **0.241 nm** |
| max | 7.71 nm | **0.431 nm** |
| shape | bimodal, ~7.7 nm plateau | single well-behaved basin |

The ~7.7 nm excursion **disappears completely** when computed from `mdfit.xtc`.
In the same corrected run the ligand is demonstrably bound throughout
(`mindist_lig_active` mean 0.204 nm, `contacts_lig_active` mean 1084). The
excursion was a **periodic-image artefact of analysing the raw `md.xtc`** — the
ligand had jumped a box vector in the unwrapped production trajectory.

**Verdict:** the previous LP-A6 ligand RMSD is reclassified from *"suspected PBC
artefact"* to **`INVALIDATED_WRONG_TRAJECTORY`**. The corrected value
(mean 0.24 nm, median 0.24 nm) is a real, populatable number — recorded here
only; **not written to any dataset in this pass**.

## 7. Source-integrity verification

* `run_report.json` → `source_integrity.verdict` = `{ok: true, added: [], removed: [], modified: [], index_changed: []}` (before/after `snapshot_tree` of `Nuevos_sistemas`).
* Independent `sha256` of every `md.xtc` / `md.tpr` / `index.ndx` for the 4
  systems — **all match** the fingerprints recorded in the archived 2026-09-08
  provenance (HMG-R-25ns-A6/md.xtc > 2 GiB: size + mtime match).
* Writes this pass are confined to: `analysis/` + `tests/` + `cli.py` (code),
  `analysis_outputs/xanthone_short/` (results + `_invalidated_md_xtc_2026-09-08/`
  archive), and `reports/trajectory_protocol_correction_2026-09-09/`.
* **No raw trajectory, `.tpr`, or source `index.ndx` was modified.**

## 8. MM-PBSA — provenance/config audit only (no rerun, no dataset change)

> **UPDATE 2026-09-09:** the "dataset A1–A5 not internally reproducible → 24-system
> rerun" conclusion below is **superseded**. The canonical source was found —
> `/home/jesusxd/Escritorio/heatmap.py` → thesis `Heat.pdf` → `mmpbsa_dataset.csv`,
> all 20 A1–A5 values matching; methodology documented in the thesis (last 5 ns,
> 200 ps) and **equivalent to A6**. The `/Escritorio/MMPBSA/` folders that disagree
> are a discarded preliminary `-dt 100` batch. **Decision: retain A1–A5, no rerun,
> no dataset change.** See `mmpbsa_canonical_source_recovery.md` and
> `mmpbsa_rerun_decision.md`.

Full report: `mmpbsa_equivalence_audit.md`. Summary:

* **Configuration is method-equivalent** between the A6 runs and the on-disk
  A1–A5 runs sharing `mmpbsa.mdp` `954785fb1805…` (17 systems incl. all 4 A6):
  same fitted `mdfit.xtc`, same 20–25 ns interval, same `Protein`/`LIG`
  selection, same PB (pdie 2 / sdie 78.3 / lpbe), same SASA-only apolar model,
  same T = 309.65 K, same g_mmpbsa build. **Manual vs automated launch is not a
  comparability criterion** — every run (A1–A6) was a manual `g_mmpbsa`
  invocation.
* **One real parameter difference:** A6 used `-dt 200` (26 snapshots), A1–A5
  used `-dt 100` (51 snapshots) over the identical window. Sampling density
  only; equalise in any future rerun.
* **Separate batch:** HMGR-A2..A5 used a different `mmpbsa.mdp` and 815 vs 1615
  residues — do not pool with HMGR-A1 / HMG-R-25ns-A6.
* **Dataset-level problem (reported, not fixed):** `mmpbsa_dataset_v2` A1–A5 is a
  cross-batch mix — the A5 rows and A1/lipase do **not** reproduce from any
  on-disk folder (20–80 kJ/mol more positive; earlier July-2025 charge batch,
  not on disk). The fix is a **single-batch 24-system rerun**, reported here and
  **not performed**. `mmpbsa_dataset_v2` and `master_v2` MM-PBSA cells
  **unchanged**.

## 9. Mechanistic protocol + dry-run (prepared, not executed)

* `mechanistic_standardized_protocol.yaml` — one trajectory treatment
  (`mdfit.xtc` for **all** systems incl. reprocessed references), one 0–200 ns
  window / 100 ps stride, common atom selections, a **shared PCA basis** for
  A1/A3/A6 (single covariance on concatenated commonly-aligned Cα trajectories;
  project each system onto the shared eigenvectors) and a shared PC1/PC2 FEL for
  direct comparison, plus per-system historically-compatible FEL. Clustering and
  conformational-state populations are marked **NEW — not historically present
  for A1**.
* `simforge study mechanistic-plan` (new command) + `test_mechanistic_plan.py`.
* `mechanistic_gap_and_dryrun.md` + `mechanistic_dryrun.json`: reference systems
  A1/APO/COA/COMP → 12 descriptors each `SUPERSEDED_PENDING` (raw `md.xtc`, no
  `mdfit.xtc` yet); new A3/A6 HMG-R systems → 12 × `PLANNED` + MM-PBSA/decomp
  `AVAILABLE` (already on `mdfit.xtc`); 24 × `NEW_EXTENSION`.
* **No mechanistic GROMACS analysis executed.**

---

## STOP — awaiting review

Nothing downstream has been propagated. Before the next pass, please review:

1. the corrected A6 statistics and the **LP-A6 ligand RMSD = 0.24 nm**
   (`INVALIDATED_WRONG_TRAJECTORY`, not PBC artefact) conclusion;
2. the **MM-PBSA** conclusion (method-equivalent by config; dataset A1–A5 column
   is a non-reproducible cross-batch mix — recommend a 24-system rerun);
3. the standardized **mechanistic protocol** and dry-run.

On approval, the deferred pass will: update the 3 MD columns (+ LP-A6 backfill)
in `md_dataset_v2_A1-A6.csv` and `master_multimodal_dataset_v2_A1-A6.csv`
(A1–A5 untouched); regenerate the `completadas_A6/` thesis figures; and write
`dataset_update_report.md`, `figure_update_report.md`, `final_report.md`.
Mechanistic execution and any MM-PBSA rerun remain separate authorised phases.
