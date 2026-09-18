# MM-PBSA — rerun decision

**Date:** 2026-09-09
**Decision: DO NOT rerun the 24 MM-PBSA systems.**

## Basis (task criteria 5 & 6)

> 5. If a canonical final A1–A5 source is recovered and its methodology is
>    equivalent to A6 → retain A1–A5, add A6, update provenance, do NOT rerun.
> 6. Only recommend a 24-system rerun if the canonical source cannot be
>    recovered **or** there is a genuine methodological incompatibility.

| test | outcome |
|---|---|
| canonical final A1–A5 source recovered? | **YES** — `/home/jesusxd/Escritorio/heatmap.py` (2026-05-12) → `Tesis-maestría/Heat.pdf` / thesis Fig. `heatmap_mmpbsa` → `Tesis-JPA.pdf`. All 20 A1–A5 values match `mmpbsa_dataset.csv` / `_v2` to 1 dp. Controls recovered too (`Sustratos-inhibidores/*`). Details: `mmpbsa_canonical_source_recovery.md`. |
| methodology documented? | **YES** — thesis §"Dinámica molecular de los complejos enzima–ligando": LigParGen/OPLS-AA `1.14*CM1A-LBCC`, SPC/E, 309.65 K, 25 ns, MM-PBSA over **the last 5 ns at 200 ps intervals** (~26 configs); controls under the same parameters. |
| methodology equivalent to A6? | **YES** — A6 25 ns MM-PBSA uses `-b 20000 -e 25000 -dt 200` (26 configs) on `mdfit.xtc`, physically-identical `mmpbsa.mdp` (the two mdp sha256 differ only by one comment character). Same window, same stride, same config count, same PB/SA physics, same T. |
| genuine methodological incompatibility found? | **NO** — one ligand-parameter protocol, one PB/SA setup, one 20–25 ns window across A1–A6. The differences among on-disk folders are (a) stride `-dt 100` in a **discarded preliminary batch** and (b) per-system re-runs whose output folders were not retained. |

**Criterion 5 is met; criterion 6 is not.** → retain A1–A5, keep A6, no rerun.

## Residual caveat (does not change the decision)

For **7 of 20** A1–A5 cells (A5 × 4, HMG-A2, HMG-A4, LP-A1) and **3 of 20**
(A1/AG, A2/AG, A4/AG) there is **no surviving g_mmpbsa output directory** — the
canonical value is documented in `heatmap.py` and the thesis but cannot be
recomputed from retained files. The remaining **10 of 20** reproduce exactly from
a surviving preliminary folder.

This is a **records-retention gap**, not a methodological problem. It does **not**
warrant re-running all 24 systems. If bit-for-bit reproducibility of those 10
specific cells is later required for publication, the minimal action is to
**re-derive only those 10** on their existing `mdfit.xtc` at `-b 20000 -e 25000
-dt 200` — a targeted 10-run job, not 24 — and confirm they land within ~2 kJ/mol
of the thesis values. That is optional and out of scope here.

## Dataset action

**None.** `mmpbsa_dataset_v2_A1-A6.csv` and the `mmpbsa_dg_kjmol` column of
`master_multimodal_dataset_v2_A1-A6.csv` already hold the canonical thesis
values (A1–A5) plus the A6 values. No cell is added, changed or removed.
`mmpbsa_dataset.csv` (base) is likewise already correct.

Provenance is updated in documentation only:

* `mmpbsa_canonical_source_recovery.md`, `mmpbsa_value_provenance_A1-A6.csv`,
  `mmpbsa_heatmap_vs_dataset_comparison.csv` (this directory) — new.
* `mmpbsa_equivalence_audit.md`, `checkpoint_report.md` §8, `final_report.md` §5
  (this directory) — annotated: the earlier "cross-batch, non-reproducible,
  recommend 24-system rerun" conclusion is **superseded**.
* `ML/…/reports/targeted_A6_completion_2026-09-09/qc_and_unresolved_items.md` Q4
  and `master_v2_change_log.md` — annotated: A1–A6 MM-PBSA **is comparable**;
  the "different batch = incompatible" framing is withdrawn.

## A1–A6 MM-PBSA comparability — final statement

A1–A6 MM-PBSA (`mmpbsa_dg_kjmol`) is **comparable within the 25 ns dataset**:
one ligand-parameter protocol, one PB/SA configuration, one 20–25 ns / 200 ps
analysis window, `mdfit.xtc` input throughout. The 200 ns mechanistic MM-PBSA
(A3-HMG-R, HMG-R-200ns-A6, the COA ternaries) is a **separate** protocol
(last 50 ns, 200 ps, ~250 configs) and its poorly-converged totals
(SD > |mean|) remain for per-residue interpretation only.
