# A6 200 ns mechanistic correction — `md_final.xtc` (2026-09-10)

> **Note (later same-day update):** the "n/a (no 200 ns run)" A1 cell in the
> table below is superseded — A1's 200 ns MM-PBSA was recovered later on
> 2026-09-10 from `HMG_CoA_R-A1-R2` (−236.5 ± 33.7 kJ/mol). This document is
> kept as the historical record of the A6 correction specifically; see
> `a1_source_reconciliation_2026-09-10/` and the updated `final_mechanistic_report.md`
> §8 for the current three-way A1/A3/A6 200 ns MM-PBSA comparison.

**Instruction (verbatim intent):** use `Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc`
directly, with `md.tpr` + `index.ndx`, as the canonical analysis trajectory for
**all** HMG-R-200ns-A6 200 ns mechanistic analyses — not the invalidated
`mdfit.xtc` / `mdcenter.xtc`, and without inferring or re-deriving a fit. This
closes the single open gap flagged in the 2026-09-09 pass
(`final_mechanistic_report.md` §2/§10: *"HMG-R-200ns-A6 — obtain / re-extract a
PBC-clean 200 ns single-ligand trajectory... The one true gap."*).

## 1. What was wrong before

`analysis_outputs/xanthone_mechanistic/HMG-R-200ns-A6/traj/mdfit.xtc` (built
2026-09-09 with `-pbc nojump → -pbc mol -center → -fit rot+trans` from the only
full-length trajectory then on disk, `mdcenter.xtc`, itself `-pbc res`-built)
still split the tetramer: Backbone RMSD reached 5–7 nm and Rg doubled
(3.3 → 6.5 nm) — an imaging artefact, not real dynamics. No raw `md.xtc` existed
to rebuild from. The system was marked **REVIEW_REQUIRED**; 200 ns whole-protein
observables, PCA/FEL/DCCM, clustering and the shared A1/A3/A6 basis were **not
computed** for A6.

## 2. The corrected trajectory

| | |
|---|---|
| path | `/home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc` |
| sha256 | `9309606e79aece25fe44eed546e787718e8c49b2a572b96787578daf5d4a721d` |
| size | 991 004 796 bytes |
| frames | 2001 (native **100 ps** stride, 0–200 000 ps) |
| atoms | 132 474 (matches `md.tpr`, matches the other HMG-R systems) |
| used with | `Nuevos_sistemas/HMG-R-200ns-A6/md.tpr`, `index.ndx` (unmodified, read-only) |

**PBC sanity check** (`gmx rms`/`gmx gyrate` on the full trajectory, Backbone):
protein RMSD mean **0.196 ± 0.025 nm**, max **0.255 nm**; Rg **3.318 ± 0.008 nm**
— clean, matches the other 7 clean systems in the study (A1: RMSD 0.173, Rg
3.324; A3: RMSD 0.222, Rg 3.330). No tetramer-splitting artefact.

## 3. How it was integrated

`md_final.xtc` is used **directly** — it is *not* renamed to `mdfit.xtc`, *not*
re-fitted, and no new derivative trajectory is written for it. A thin adapter
(`scripts/a6_common.py`) monkey-patches the trajectory-resolution function used
by the existing phase-2..6 modules (`core_descriptors.py`,
`essential_dynamics.py`, `clustering.py`, `cluster_states.py`) so that, for
`HMG-R-200ns-A6` only, they resolve their analysis trajectory to `md_final.xtc`
with correctly-labelled provenance (`source_trajectory`,
`source_trajectory_sha256`, `source_trajectory_designation`). Every other
system is completely unaffected. GROMACS tools that need a least-squares fit
(`gmx rms`, `gmx rmsf`, `gmx covar`) perform it **internally, on the fly** — the
same way they do for every other system in this pipeline; no separate fitted
trajectory file was produced.

The previously-contaminated managed outputs were archived, not deleted:
`HMG-R-200ns-A6/_pbc_contaminated_2026-09-09/` (already existed) and
`HMG-R-200ns-A6/_superseded_2026-09-10_pre_md_final/` (the last `mdfit.xtc`
attempt). Raw sources
(`Nuevos_sistemas/HMG-R-200ns-A6/{md_final.xtc,md.tpr,index.ndx}`) were only
read.

## 4. What was (re)computed

All from `md_final.xtc`, full 200 ns, 2001 frames:

- **Core descriptors** (9 xvg): protein RMSD, ligand RMSD, RMSF (Cα),
  Rg, SASA, protein H-bonds, ligand–active-site H-bonds, ligand–active-site
  min-distance + contacts (< 0.6 nm), ligand–catalytic-COM distance.
- **PCA** (`gmx covar`/`gmx anaeig`, per-system, Cα): PC1 29.3 %, PC2 9.3 %.
- **FEL** (per-system PC1/PC2, 2D histogram, ΔG = −kT ln P): 12 basins < 3 kJ/mol.
- **DCCM** (Cα–Cα normalised cross-correlation, MDAnalysis): mean |C| 0.095.
- **Clustering** (`gmx cluster` gromos, Cα/Backbone RMSD, 0.15 nm cutoff):
  5 clusters, top states 75.2 % / 16.7 % / 7.8 % / 0.2 %.
- **Macro-state scan** (adaptive 0.15/0.25/0.35 nm): converges at 0.15 nm
  (same 5 states — consistent with every other system in the study; 0.25 nm
  collapses to 1 cluster everywhere).
- **Windowed analysis** (0–50/50–100/100–150/150–200 ns): stationarity + OLS
  slope for every observable — folded into `windowed_analysis_A1_A3_A6.csv`.
- **Shared A1/A3/A6 Cα PCA basis + shared FEL**: rebuilt to include A6
  (previously A1-vs-A3 only). A6's Cα trajectory (2001 frames) is subsampled
  1-in-2 (in memory, no new trajectory file) to ~1001 frames so it carries the
  same statistical weight as A1/A3 in the **pooled basis fit only**; A6's own
  per-system PCA/FEL/DCCM above use the full native 2001-frame series.
- **DCCM comparisons**: ΔDCCM(A6−A1) and ΔDCCM(A6−A3) added alongside the
  existing ΔDCCM(A3−A1).
- **Single-ligand vs ternary comparison**: A6-solo vs `system_A6_COA`
  (previously used the 25 ns short-study as a proxy for A6-solo; now uses the
  real 200 ns A6-solo data).
- **200 ns MM-PBSA + per-residue decomposition**: **incorporated, not rerun** —
  `Nuevos_sistemas/HMG-R-200ns-A6/{energy_summary.csv,
  residues_energy_summary.csv}` were already completed by the user on
  `md_final.xtc` (`g_mmpbsa run -f md_final.xtc -s md.tpr -b 150000 -e 200000
  -dt 500 ...`, 101 configurations, 150–200 ns window). See
  `mmpbsa_200ns_A1_A3_A6.md`.

126 → 153 provenance JSON files in total across the study (27 new, all for A6,
`source_trajectory` ending in `md_final.xtc` by design — see `qc_report.md`).

## 5. Headline results now available for A6 (200 ns, real data — not a proxy)

| observable | A1 | A3 | **A6 (200 ns, md_final.xtc)** |
|---|---|---|---|
| protein RMSD (nm) | 0.173 ± 0.017 | 0.222 ± 0.026 | **0.196 ± 0.025** |
| ligand RMSD (nm) | 0.66 ± 0.20 | 0.64 ± 0.09 | **0.47 ± 0.06** |
| RMSF Cα mean (nm) | 0.089 | 0.093 | **0.087** |
| Rg (nm) | 3.324 ± 0.006 | 3.330 ± 0.015 | **3.318 ± 0.008** |
| SASA (nm²) | 527 ± 5 | 531 ± 8 | **528 ± 5** |
| protein H-bonds | 1212 ± 17 | 1209 ± 21 | **1236 ± 17** |
| active-site min-dist (nm) | 0.239 ± 0.032 | 0.257 ± 0.027 | **0.240 ± 0.018** |
| active-site contacts | 342 ± 214 | 267 ± 92 | **282 ± 41** |
| catalytic-COM distance (nm) | 2.603 ± 0.099 | 3.027 ± 0.127 | **2.764 ± 0.041** |
| clustering (0.15 nm) top states | 3 cl., 93.4/3.9 % | 5 cl., 76.2/15.6 % | **5 cl., 75.2/16.7 %** |
| per-system PCA PC1/PC2 | 21.2/11.6 % | 26.4/9.8 % | **29.3/9.3 %** |
| MM-PBSA 200 ns (150–200 ns, kJ/mol) | n/a (no 200 ns run) | −127.8 ± 110.9 (poorly converged) | **−232.2 ± 28.1 (well converged)** |

A6's own 200 ns numbers are close to the 25 ns short-study proxy that stood in
for it previously (25 ns: ligand RMSD 0.52, active-site contacts 283, catalytic
COM 2.78 nm) — the 25 ns picture was a reasonable stand-in, now confirmed
directly.

**New, previously unavailable:**

- **Shared A1/A3/A6 PCA basis** (`shared_pca_summary.md`): PC1 27.9 %, PC2
  13.3 %, cum(PC1,2) 41.2 %. **A6 occupies a third, distinct region of the
  shared conformational plane** — centroid (0.01, 2.64) nm, essentially on the
  A1↔A3 PC1 midpoint (A1 −3.30, A3 +3.29) but displaced along **PC2**, which
  A1 and A3 barely populate (−1.26, −1.38). Pairwise Bhattacharyya overlap is
  low for all three pairs (A1–A3 0.009, A1–A6 0.009, A3–A6 **0.001**, the
  lowest) — A6's conformational ensemble is the most distinct of the three.
- **ΔDCCM(A6−A1)**: mean |Δ| 0.114, 4.8 % of Cα pairs shift > 0.3 — comparable
  in magnitude to ΔDCCM(A3−A1) (0.121, 5.7 %). **ΔDCCM(A6−A3)**: mean |Δ| 0.106,
  3.6 % — A6 and A3 are somewhat closer to each other dynamically than either
  is to A1.
- **A6-solo vs A6+CoA (real data)**: catalytic-COM distance solo 2.764 nm →
  ternary 3.342 nm (**Δ +0.578 nm, the largest of the three** — vs A1 +0.321,
  A3 +0.296); active-site contacts solo 282 → ternary 161; **conformational
  confinement solo 75.2 % → ternary 71.2 %** (essentially unchanged / slightly
  *less* confined) — confirms directly, no longer via proxy, that unlike A1
  (93.4 → 95.7 %) and A3 (76.2 → **97.4 %**), **CoA binding does not lock the
  enzyme when A6 is present**; A6 is also pushed furthest from the catalytic
  centre by CoA of the three.
- **200 ns MM-PBSA, well converged** (SD/|mean| 0.12 vs A3's 0.87): **A6
  −232.2 ± 28.1 kJ/mol**, more favourable than A1's 25 ns canonical value
  (−190.9) and far more favourable than A3's 200 ns value (−127.8, poorly
  converged). Per-residue decomposition (`mmpbsa_200ns_A1_A3_A6.md`) shows A6's
  binding is anchored by a hydrophobic cluster (LEU-853 −36.3, VAL-683 −29.8,
  HIS-752 −26.1, SER-684 −26.1, LEU-857 −23.2 kJ/mol) **plus** real, if modest,
  catalytic-region electrostatic engagement (LYS-691 −3.7, LYS-692 −3.0,
  **ARG-702 −2.2**, ARG-646 −1.2) — a refinement of the earlier "diffuse,
  non-catalytically-focused" read: A6 does engage ARG-702, just far more weakly
  than A1 (−9.9 kJ/mol at 25 ns) and mostly through a large hydrophobic anchor
  elsewhere in the pocket.

## 6. What did NOT change

- A1, APO, COA, COMP, A3-HMG-R, system_A3_COA, system_A6_COA — none of their
  analyses, trajectories or outputs were touched.
- No raw simulation input was modified (`md.xtc`, `md.tpr`, `index.ndx`,
  `md_final.xtc` — read-only; `source_integrity.json` / `qc_report.md`
  confirm `Mecanismo_inhibitorio` unchanged and the only `Nuevos_sistemas`
  deltas are the user's own new `HMG-R-200ns-A6` files (`md_final.xtc`, the
  refreshed MM-PBSA outputs) plus unrelated report/figure writes already
  nested under that tree from other passes).
- The 25 ns `xanthone_short` short-MD study and its figures are unaffected.

## 7. Outputs touched by this correction

- `analysis_outputs/xanthone_mechanistic/HMG-R-200ns-A6/{observables,essential_dynamics,clustering}/**`
  (new, 27 provenance JSON + xvg/npy/png/pdb)
- `analysis_outputs/xanthone_mechanistic/shared_pca/**` (rebuilt; A1-vs-A3-only
  backed up to `_A1_A3_only_2026-09-09/`)
- `analysis_outputs/xanthone_mechanistic/mmpbsa_200ns/**` (new)
- `reports/xanthone_mechanistic_2026-09-10/{mechanistic_descriptor_summary.csv,
  clustering_state_populations.csv, shared_pca_summary.md, fel_state_summary.md,
  dccm_summary.md, windowed_analysis_A1_A3_A6.csv, execution_summary.md,
  qc_report.md, final_mechanistic_report.md}` (regenerated / updated)
- new: `mmpbsa_200ns_A1_A3_A6.md`, `solo_vs_ternary_A1_A3_A6.md`,
  `mechanistic_A6_md_final_correction.md` (this file)
- `analysis_outputs/xanthone_mechanistic/HMG-R-200ns-A6/REVIEW_REQUIRED.md` and
  `analysis_outputs/xanthone_mechanistic/_NEW_SYSTEM_PBC_ISSUE.md` — resolution
  banners added, original findings kept.
- `system_compatibility.csv` — HMG-R-200ns-A6 row updated to OK.

Scripts: `reports/xanthone_mechanistic_A6_2026-09-10/scripts/` (`a6_common.py`,
`run_a6_core.py`, `shared_pca_a6.py`, `run_a6_shared.py`, `finalize_a6.py`,
`mmpbsa_200ns_incorporate.py`, `solo_vs_ternary_a6.py`).
