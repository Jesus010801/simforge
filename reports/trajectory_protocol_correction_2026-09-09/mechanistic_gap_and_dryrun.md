# Mechanistic profile — gap analysis and dry-run

**Date:** 2026-09-09
**Command:** `simforge study mechanistic-plan /home/jesusxd/Escritorio/Nuevos_sistemas
--reference /home/jesusxd/Escritorio/Mecanismo_inhibitorio --json`
**Machine-readable output:** `mechanistic_dryrun.json` (this directory).
**Nothing was executed.** The planner never invokes GROMACS.

The standardized protocol this plan is measured against:
`mechanistic_standardized_protocol.yaml` (this directory).

## Status vocabulary

| status | meaning |
|---|---|
| `AVAILABLE` | a valid output already exists **and** was produced on `mdfit.xtc` |
| `PLANNED` | `mdfit.xtc` present; the command can be built + run when authorised |
| `REVIEW_REQUIRED` | `mdfit.xtc` / TPR / index missing, or no ligand for a ligand observable |
| `SUPERSEDED_PENDING` | a historical output exists but was computed on raw `md.xtc` — recompute on `mdfit.xtc` |
| `NEW_EXTENSION` | analysis with no historical A1 precedent (clustering / state populations) |

## Totals

| status | count |
|---|---|
| PLANNED | 100 |
| SUPERSEDED_PENDING | 43 |
| NEW_EXTENSION | 24 |
| REVIEW_REQUIRED | 13 |
| AVAILABLE | 12 |

## Per system

### Historical reference study (`Mecanismo_inhibitorio/`)

| system | ligand | `mdfit.xtc`? | descriptor outputs | verdict |
|---|---|---|---|---|
| A1 | A1 (xanthone) | **no** | rmsd / rmsf / rg / sasa / hb / mindist / contacts / dist / PCA (eigenval, pc1, pc2) / FEL / DCCM — all present, all from raw `md.xtc` (200 ns, 100 ps) | **12 × SUPERSEDED_PENDING**; mmpbsa + per-residue REVIEW_REQUIRED (MM-PBSA for A1 lives in `/Escritorio/MMPBSA/HMGR-A1/`, not this dir) |
| APO | none | no | rmsd / rmsf / rg / sasa / hb / PCA / FEL / DCCM | 8 × SUPERSEDED_PENDING; 6 × REVIEW_REQUIRED (4 need a ligand, 2 = mmpbsa/decomp) |
| COA | CoA | no | full descriptor set | 12 × SUPERSEDED_PENDING; 2 × REVIEW_REQUIRED (mmpbsa/decomp) |
| COMP | A1 + CoA | no | full descriptor set (ligand + CoA variants) | 11 × SUPERSEDED_PENDING; 3 × REVIEW_REQUIRED |

**Action for the reference systems:** generate `mdcenter.xtc → mdfit.xtc` with the
standard chain, then recompute every descriptor on `mdfit.xtc`. Per the user
decision, the references are reprocessed so that A1/A3/A6 share one trajectory
treatment.

### New mechanistic systems (`Nuevos_sistemas/`)

| system | role | `mdfit.xtc`? | verdict |
|---|---|---|---|
| A3-HMG-R | A3 / HMG-CoA reductase, 200 ns | yes | 12 × PLANNED (rmsd…dccm), 2 × AVAILABLE (MM-PBSA + per-residue — `-f mdfit.xtc`, `-b 150000 -e 200000 -dt 200`), 2 × NEW_EXTENSION |
| HMG-R-200ns-A6 | A6 / HMG-CoA reductase, 200 ns | yes (raw `md.xtc` absent from folder) | 12 × PLANNED, 2 × AVAILABLE, 2 × NEW_EXTENSION |
| system_A3_COA | A3 + CoA ternary, ~200 ns | yes | 14 × PLANNED, 2 × NEW_EXTENSION — MM-PBSA exists but in `A3-MMPBSA/` sub-dirs (planner scans the top dir only; treat as AVAILABLE after confirming `-f mdfit.xtc` in those headers) |
| system_A6_COA | A6 + CoA ternary, ~200 ns | yes | 14 × PLANNED, 2 × NEW_EXTENSION — same MM-PBSA-in-subdir note |
| AA-A6 / AG-A6 / HMG-R-25ns-A6 / LP-A6 | 25 ns short systems (in the same tree) | yes | 12 × PLANNED, 2 × AVAILABLE (MM-PBSA on `mdfit.xtc`), 2 × NEW_EXTENSION — **but these are `xanthone_short` systems, not mechanistic**; a full 200 ns mechanistic suite is not in scope for them |

## What remains to execute for the A1/A3/A6 mechanistic comparison

1. **Reprocess the 4 reference systems** to `mdfit.xtc`
   (`Mecanismo_inhibitorio/{A1,APO,COA,COMP}`).
2. **Recompute 12 descriptors** on `mdfit.xtc` for A1/APO/COA/COMP
   (`SUPERSEDED_PENDING` → done): protein RMSD, ligand RMSD, RMSF, Rg, SASA,
   H-bonds, active-site min distance, active-site contacts, catalytic COM
   distance, PCA, FEL, DCCM.
3. **Compute the same 12 descriptors** for the new systems
   (`A3-HMG-R`, `HMG-R-200ns-A6`, and the two COA ternaries) — all `PLANNED`.
4. **MM-PBSA + per-residue decomposition:** already `AVAILABLE` on `mdfit.xtc`
   for A3-HMG-R and HMG-R-200ns-A6; for A1 use `/Escritorio/MMPBSA/HMGR-A1/`
   (already `mdfit.xtc`, but `-dt 100` vs the 200 ns runs' `-dt 200` — equalise
   in a rerun if quantitative comparison is needed). 200 ns totals for A3/A6 are
   poorly converged (SD > |mean|) → per-residue interpretation only.
5. **Shared PCA basis** for A1/A3/A6 (same receptor): one covariance on the
   concatenated, commonly-aligned Cα trajectories; project each system onto the
   shared eigenvectors; build a shared PC1/PC2 FEL for direct comparison
   (see `mechanistic_standardized_protocol.yaml` §`shared_pca_basis` / `shared_fel`).
   Keep per-system FEL as the historically-compatible output.
6. **NEW extensions** (24 × `NEW_EXTENSION`): `gmx cluster` (gromos) +
   conformational-state populations — flagged explicitly as **not historically
   present for A1**.

## Not done in this pass

No mechanistic GROMACS analysis was executed. No reference system was
reprocessed. The `simforge study mechanistic-plan` command and the protocol YAML
are delivered for review; execution is a separate, explicitly-authorised phase.
