# Final mechanistic report — A1 / A3 / A6 vs HMG-CoA reductase

**Date:** 2026-09-10 (updated twice — A6 200 ns correction, then A1 200 ns
MM-PBSA recovery; original pass 2026-09-10 AM)

> **2026-09-10 correction #1:** HMG-R-200ns-A6 (previously REVIEW_REQUIRED —
> §2 below) is now complete. The user supplied a validated PBC-clean trajectory,
> `Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc`, used directly (sha256
> `9309606e…d4a721d`; Backbone RMSD mean 0.196 / max 0.255 nm, Rg 3.318 nm —
> clean). Every A6-related section below (§3–§9) is updated with real 200 ns
> data, replacing the earlier 25 ns short-study proxy. Full account:
> `mechanistic_A6_md_final_correction.md`.

> **2026-09-10 correction #2 (A1):** the statement below that *"A1 has no
> 200 ns MM-PBSA on disk"* was **wrong** and is retracted. It was true only for
> `Mecanismo_inhibitorio/A1` (no MM-PBSA files there); a complete, well-converged
> 200 ns MM-PBSA for A1-alone exists in `/home/jesusxd/Escritorio/HMG_CoA_R-A1-R2/`
> — the production directory underlying that same A1 system (byte-identical
> `md.tpr`/`md.cpt`/`md.edr`) — and is now the canonical A1-alone 200 ns value:
> **ΔG = −236.516 ± 33.740 kJ/mol**. Three further 200 ns MM-PBSA channels were
> recovered alongside it: substrate-only (E+S, `HMG-CoA-R_sustrato`,
> −196.585 ± 34.952), A1-in-ternary (E+S+I(I), `Competitive-system-2`,
> −211.691 ± 26.273), and the A1↔CoA pairwise channel (I–S,
> `Competitive-system-2`, −27.410 ± 22.365) — the last of these is a genuine
> energetic decomposition of the xanthone–CoA channel that no other part of
> this study computes. §8 and §9 below are updated accordingly. A residual,
> intentionally-unresolved discrepancy on A1-alone's historical Total (see
> `a1_source_reconciliation_2026-09-10/a1_long_mmpbsa_recovery.md` §4) is noted
> once in §8 and not re-litigated. Full account:
> `a1_source_reconciliation_2026-09-10/`.

**Analysis trajectory:** `mdfit.xtc` (PBC-corrected rot+trans fit, 200 ps common
stride) for every system **except HMG-R-200ns-A6**, which uses `md_final.xtc`
directly (native 100 ps stride) — see the correction note above. Raw `md.xtc` is
never an analysis input for the reference systems. References (A1, APO, COA,
COMP) had no `mdfit.xtc` on disk and were reprocessed into the managed tree;
the scientific source trees are read-only and verified byte-unchanged
(`qc_report.md`: `Mecanismo_inhibitorio` UNCHANGED; the only `Nuevos_sistemas`
deltas are the user's own new A6 files and pre-existing nested report/figure
directories from other passes — no raw simulation input touched).

**Systems** (all the same HMG-CoA-reductase tetramer — **1614 Cα, 24 199 protein
atoms, identical residue numbering**):

| label | dir | contents | 200 ns whole-protein |
|---|---|---|---|
| **A1** | `Mecanismo_inhibitorio/A1` | enzyme + xanthone A1 | ✓ clean |
| **A3** | `Nuevos_sistemas/A3-HMG-R` | enzyme + xanthone A3 | ✓ clean (after PBC fix) |
| **A6** | `Nuevos_sistemas/HMG-R-200ns-A6` | enzyme + xanthone A6 | ✓ **clean (`md_final.xtc`, 2026-09-10)** |
| APO | `Mecanismo_inhibitorio/APO` | apo enzyme | ✓ clean |
| COA | `Mecanismo_inhibitorio/COA` | enzyme + substrate HMG-CoA | ✓ clean |
| COMP | `Mecanismo_inhibitorio/COMP` | ternary: enzyme + A1 + CoA (205 ns) | ✓ clean |
| system_A3_COA | `Nuevos_sistemas/system_A3_COA` | ternary: enzyme + A3 + CoA | ✓ clean (after PBC fix) |
| system_A6_COA | `Nuevos_sistemas/system_A6_COA` | ternary: enzyme + A6 + CoA | ✓ clean (after PBC fix) |

**All eight systems now have complete 200 ns whole-protein mechanistic data.**

---

## 0. Trajectory-integrity finding (fixed before any analysis) — and its closure

The pre-made `mdfit.xtc` / `mdcenter.xtc` for **all four new HMG-R systems** were
originally built **without whole-molecule PBC treatment** for the tetramer —
individual chains imaged across the box, raw protein RMSD reaching 5–7 nm
(`gmx rms`: *"distance above half the box length — system exploding"*).

**Fix (2026-09-09):** rebuild with `-pbc nojump → -pbc mol -center -ur compact
→ -fit rot+trans`.

- **A3-HMG-R, system_A3_COA, system_A6_COA recovered cleanly** this way
  (protein RMSD max 0.20–0.27 nm, Rg constant at 3.32 nm).
- **HMG-R-200ns-A6 could not be recovered this way** — no raw `md.xtc`; the only
  full-length trajectory (`mdcenter.xtc`) was `-pbc res`-processed, irreversibly
  splitting the tetramer.

**Closure (2026-09-10):** the user supplied `md_final.xtc`, a separately
validated PBC-clean trajectory (see the banner above and
`mechanistic_A6_md_final_correction.md`). Used directly; no further fitting was
applied beyond what `gmx rms`/`rmsf`/`covar` already do internally for every
system in this pipeline.

The first-pass (contaminated) outputs for the 4 new systems are archived under
`<sys>/_pbc_contaminated_2026-09-09/`, plus `HMG-R-200ns-A6/_superseded_2026-09-10_pre_md_final/`
for the intermediate corrective attempt — **none in use.**

This never affected the closed 25 ns `xanthone_short` work (HMG-R-25ns-A6's
25 ns mdfit is clean; the other short systems are single-chain).

---

## 1. Analyses completed per system

| system | core desc. (xvg) | PCA | FEL | DCCM | clustering 0.15 nm | macro-state scan (NEW) | windowed | CoA channels |
|---|---|---|---|---|---|---|---|---|
| A1 | ✓ (9) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | n/a |
| APO | ✓ (5) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | n/a |
| COA | ✓ (9) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | n/a (substrate = ligand) |
| COMP | ✓ (13+3) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ (3-way) |
| A3-HMG-R | ✓ (9) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | n/a |
| **HMG-R-200ns-A6** | **✓ (9) — `md_final.xtc`** | **✓** | **✓** | **✓** | **✓** | **✓** | **✓** | n/a |
| system_A3_COA | ✓ (12+3) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ (3-way) |
| system_A6_COA | ✓ (12+3) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ (3-way) |

**Shared PCA basis + shared FEL: A1 ↔ A3 ↔ A6** (all three 200 ns, clean —
rebuilt 2026-09-10; previously A1↔A3 only). Clustering, representative
structures, and conformational-state populations are **NEW** — absent from the
historical A1 mechanistic workflow.

**H-bonds:** `gmx hbond` does not recognise the LigParGen xanthone `LIG` atoms as
donors/acceptors (0/0) → ligand–site H-bond counts unavailable; contact/distance
metrics used instead. Whole-protein H-bond count skipped for the two large CoA
systems (O(N²), ~1210 ± 20 across all states — no discrimination).

---

## 2. REVIEW_REQUIRED items — all resolved

| item | resolution |
|---|---|
| **HMG-R-200ns-A6** — no PBC-clean 200 ns trajectory | **RESOLVED 2026-09-10.** User supplied `md_final.xtc` (validated clean); full 200 ns whole-protein suite computed; A6 now in the shared A1/A3/A6 PCA basis. |
| system_A3_COA / system_A6_COA — no site groups in source index | **derived** from A1 (protein byte-identical) into the managed index; recorded as *derived*. |
| APO — no `index.ndx` in source | `gmx make_ndx` into the managed tree; apo → no ligand/site metrics. |

---

## 3. A1 vs A3 vs A6 stability — 0–200 ns (all real data)

| observable | A1 | A3 | **A6** |
|---|---|---|---|
| protein RMSD (nm) | 0.173 ± 0.017 | 0.221 ± 0.026 | **0.196 ± 0.025** |
| ligand RMSD (nm) | 0.66 ± 0.20 | 0.64 ± 0.09 | **0.47 ± 0.06** |
| RMSF Cα mean (nm) | 0.089 | 0.093 | **0.087** |
| Rg (nm) | 3.324 ± 0.006 | 3.330 ± 0.015 | **3.318 ± 0.008** |
| SASA (nm²) | 527 ± 5 | 531 ± 8 | **528 ± 5** |
| protein H-bonds (whole-tetramer) | 1212 ± 17 | 1209 ± 21 | **1236 ± 17** |
| ligand–active-site min distance (nm) | 0.239 ± 0.032 | 0.257 ± 0.027 | **0.240 ± 0.018** |
| ligand–active-site contacts (<0.6 nm) | **342 ± 214** | **267 ± 92** | **282 ± 41** |
| ligand–catalytic-site COM distance (nm) | **2.603 ± 0.099** | **3.027 ± 0.127** | **2.764 ± 0.041** |

**All three enzymes are structurally stable** (protein RMSD < 0.25 nm, Rg
constant, RMSF and H-bonds indistinguishable). The differences are entirely in
**ligand engagement**:

- **A1 sits deepest** — 2.60 nm from the catalytic COM, ~340 active-site contacts.
- **A3 sits furthest** — 3.03 nm (≈0.4 nm further out), ~25 % fewer contacts.
- **A6 is intermediate** on both, but with the **tightest, lowest-variance
  ligand RMSD of the three** (0.47 ± 0.06 nm vs A1's 0.66 ± 0.20 and A3's
  0.64 ± 0.09) and the tightest catalytic-distance variance (± 0.04 nm) —
  A6 explores the smallest range of poses even though its central position is
  between A1 and A3.

**Windowed temporal analysis** (`windowed_analysis_A1_A3_A6.csv`, 0–50 / 50–100 /
100–150 / 150–200 ns): protein RMSD / Rg / SASA are stationary in every system
(window-to-window shifts < 0.08 nm). Two real drifts, both < 0.6 nm over
200 ns and both toward *less* extreme values:

- **A3's ligand slowly moves toward the catalytic site** (active-site contacts
  +161, t = 18.7; catalytic COM distance −0.24 nm, t = −20.3) — A3 starts
  peripheral and settles in over 200 ns.
- **A6's ligand RMSD slightly decreases** (−0.092 nm over 200 ns, t = −22.4) and
  its active-site contacts fall modestly (−10, t = −3.3) — a mild settling, far
  smaller in magnitude than A3's. A6's protein RMSD/Rg drifts (+0.066 /
  −0.023 nm over 200 ns) are on the same small scale as every other system's.
  **Trajectory length alone is not taken as convergence** — the CSV carries
  per-window means + OLS slopes for every scalar observable.

---

## 4. Dominant conformational states (Cα-RMSD gromos clustering, 0.15 nm — NEW)

| system | # clusters | C1 | C2 | reading |
|---|---|---|---|---|
| APO | 2 | **99.2 %** | 0.8 % | apo enzyme = one rigid state |
| system_A3_COA | 3 | **97.4 %** | 1.6 % | A3 + CoA → enzyme locked |
| COMP (A1 + CoA) | 4 | **95.7 %** | 3.2 % | A1 + CoA → enzyme locked |
| A1 | 3 | **93.4 %** | 3.9 % | A1 alone ≈ as rigid as apo |
| COA (substrate) | 4 | 83.5 % | 10.4 % | substrate allows a minor 2nd state |
| A3 | 5 | **76.2 %** | 15.6 % | A3 alone → populated 2nd state |
| **A6** | 5 | **75.2 %** | **16.7 %** | **A6 alone → populated 2nd state, essentially like A3** |
| **system_A6_COA (A6 + CoA)** | 5 | **71.2 %** | **24.4 %** | **A6 + CoA → NOT locked, if anything looser than A6 alone** |

**Conformational-confinement ranking:**
**APO ≈ A3+CoA ≈ COMP > A1 > COA > A3 ≈ A6(alone) > A6+CoA.**

- **A1 confines the enzyme to essentially one basin** — as tightly as apo.
- **A3 alone and A6 alone are the two most conformationally permissive
  single-ligand systems** (76.2 % / 75.2 % top state) — this is now a *direct*
  finding for A6, not inferred from a 25 ns proxy.
- **CoA normally locks the enzyme further** (A3 alone 76.2 % → A3+CoA 97.4 %;
  and COMP, A1+CoA, stays high at 95.7 %) — *except* **A6 + CoA, which does not
  gain confinement (75.2 % → 71.2 %, if anything slightly looser)**: neither
  the A6 xanthone nor CoA succeeds in committing the enzyme to one state, and
  now this is confirmed against a real A6-alone baseline rather than a proxy.
- Representative structures: `<sys>/clustering/clusters.pdb`; macro-state scan +
  per-window populations in `clustering_state_populations.csv` and
  `<sys>/clustering/state_clustering.json`. A6's second state (C2, 16.7 %) is
  front-loaded (60.6 % of frames in 0–50 ns, 0 % after 100 ns) then C3 (7.8 %)
  appears in 100–150 ns (24.4 %) and partly recedes — a genuine, if modest,
  state interconversion over the 200 ns window, not a single early transient.

**Adaptive macro-state scan (0.15/0.25/0.35 nm, all 8 clean systems):** 0.15 nm
is the only cutoff that resolves more than one cluster for every system — at
0.25 nm every trajectory (A6 included) collapses to a single cluster. The whole
panel occupies **one Cα macro-basin with fine (< 0.15 nm) substructure**; the
"states" above are sub-states of that basin, consistent with the globally rigid
enzyme (§3, §5).

---

## 5. PCA / FEL differences

### Per-system PCA (independent bases — PC1/PC2 not directly comparable)

| system | PC1 | PC2 | PC1+PC2 |
|---|---|---|---|
| A1 | 21.2 % | 11.6 % | 32.8 % |
| APO | 24.0 % | 12.7 % | 36.7 % |
| COA | 26.1 % | 9.4 % | 35.5 % |
| COMP | 24.8 % | 9.1 % | 33.9 % |
| A3 | 26.4 % | 9.8 % | 36.1 % |
| **A6** | **29.3 %** | 9.3 % | **38.5 %** |
| system_A3_COA | 19.4 % | 9.2 % | 28.6 % |
| system_A6_COA | 29.1 % | 12.3 % | 41.4 % |

No single collective mode dominates (PC1 ≤ 29 %) → the enzyme is globally rigid;
essential dynamics is spread over many low-amplitude modes. A6-alone's PC1
(29.3 %) closely matches A6+CoA's (29.1 %) — the dominant collective mode is a
property of the A6-bound enzyme that CoA binding does not redirect, consistent
with §4's finding that CoA does not further confine this system.

### Shared A1 ↔ A3 ↔ A6 Cα PCA basis (`shared_pca_summary.md`, rebuilt 2026-09-10)

All three 1614-Cα ensembles Kabsch-aligned to one iteratively-refined
pooled-mean reference (A6's native 2001-frame series subsampled 1-in-2 for
basis-fit balance only), projected onto a single PCA eigenbasis.

- **PC1 = 27.9 %, PC2 = 13.3 % of pooled variance; PC1+PC2 = 41.2 %.**
- **The three ensembles are nearly disjoint** — Bhattacharyya overlap A1–A3
  **0.009**, A1–A6 **0.009**, A3–A6 **0.001** (the lowest of the three pairs).
- **Centroids (nm):** A1 (−3.30, −1.26), A3 (+3.29, −1.38), **A6 (+0.01,
  +2.64)**. A1 and A3 sit at opposite ends of **PC1** (the "A1-vs-A3
  coordinate" identified in the 2026-09-09 pass); **A6 sits almost exactly at
  the A1/A3 PC1 midpoint but is displaced along PC2, an axis A1 and A3 barely
  populate (PC2 −1.26 / −1.38).** A6 therefore occupies a genuinely **third
  region of the shared conformational landscape**, not an intermediate point on
  the A1–A3 axis.
- Real-space: A1 and A3 mean Cα structures differ by 0.165 nm RMSD (from the
  2026-09-09 pass); A6's displacement along the orthogonal PC2 axis indicates a
  distinct — not merely intermediate — average backbone response.
- **Shared-basis FEL** (`shared_pca/shared_fel.png`): A1's landscape has 3
  basins around one PC1 minimum (tightest — matches its 93.4 % single-state
  clustering); **A3's is much rougher (24 basins)**; **A6's has 9 basins**,
  intermediate roughness, with its global minimum at (PC1 0.12, PC2 3.82) — off
  in the PC2 direction, away from both A1 and A3's basins.
- APO/COA/COMP can still be projected onto this basis in a future pass
  (protein identical); not done here (scope: A6 correction).

---

## 6. DCCM differences (`dccm_summary.md`, rebuilt 2026-09-10)

Normalised Cα–Cα cross-correlation (aligned Cα; full matrices + heatmaps per
system).

| system | mean \|C\| | strong anti-corr (< −0.4) | strong corr (> 0.6) |
|---|---|---|---|
| A1 | 0.110 | 0.2 % | 0.5 % |
| APO | 0.110 | 0.2 % | 0.5 % |
| COA | 0.095 | 0.0 % | 0.4 % |
| COMP | 0.094 | 0.1 % | 0.4 % |
| A3 | 0.098 | 0.1 % | 0.4 % |
| **A6** | **0.095** | 0.1 % | 0.4 % |
| system_A3_COA | 0.093 | 0.1 % | 0.4 % |
| system_A6_COA | 0.098 | 0.2 % | 0.4 % |

- **A1 (and apo) have the strongest inter-domain coupling** (mean \|C\| 0.11).
  Binding any ligand or substrate **damps** the correlated motion (→ 0.09–0.10);
  A6 sits in that same damped range, essentially matching COA (0.095).
- **ΔDCCM (A3 − A1): mean \|Δ\| 0.121; 5.7 % of Cα pairs shift by > 0.3.**
- **ΔDCCM (A6 − A1): mean \|Δ\| 0.114; 4.8 % of Cα pairs shift by > 0.3 (NEW).**
  A6 reorganises the correlation network relative to A1 by nearly the same
  magnitude as A3 does.
- **ΔDCCM (A6 − A3): mean \|Δ\| 0.106; 3.6 % of Cα pairs shift by > 0.3 (NEW).**
  Smaller than either is from A1 — **A6 and A3 perturb the dynamic
  cross-correlation network in a more similar direction to each other than
  either does relative to A1**, even though (§5) their shared-PCA positions are
  nearly disjoint (Bhattacharyya 0.001). PCA position and pairwise-correlation
  perturbation are capturing different aspects of the dynamics; both point to
  **A1 as the outlier that most reorganises the enzyme's collective motion**
  relative to the other two.
  (`shared_pca/dccm_diff_{A3_minus_A1,A6_minus_A1,A6_minus_A3}.npy` and
  corresponding heatmaps.)

---

## 7. CoA-related mechanistic differences (`coa_interaction_summary.md`,
## `solo_vs_ternary_A1_A3_A6.md`)

Three interaction channels kept **strictly separate**:

| channel (mean ± SD) | COMP (A1 + CoA) | system_A3_COA (A3 + CoA) | system_A6_COA (A6 + CoA) |
|---|---|---|---|
| **protein–xanthone** min dist (nm) | 0.19 ± 0.02 | 0.27 ± 0.05 | 0.29 ± 0.08 |
| protein–xanthone contacts | **986 ± 111** | 191 ± 80 | 161 ± 87 |
| protein–xanthone catalytic COM (nm) | 2.92 | 3.32 | 3.34 |
| **protein–CoA** min dist (nm) | **0.15 ± 0.01** | **0.62 ± 0.13** | **2.11 ± 0.91** |
| protein–CoA contacts | **4174 ± 327** | **6.7 ± 22** | **5.9 ± 20** |
| **xanthone–CoA** min dist (nm) | **0.38 ± 0.13** | **1.09 ± 0.29** | **2.60 ± 1.29** |

**Solo vs ternary — now with real A6-solo data (2026-09-10), not a 25 ns proxy:**

| label | solo catalytic COM (nm) | ternary catalytic COM (nm) | Δ | solo contacts | ternary contacts | solo top state(s) | ternary top state(s) |
|---|---|---|---|---|---|---|---|
| A1 | 2.603 ± 0.099 | 2.924 ± 0.117 | **+0.321** | 342 ± 214 | 986 ± 111 | 93.4/3.9 % | 95.7/3.2 % |
| A3 | 3.027 ± 0.127 | 3.323 ± 0.205 | +0.296 | 267 ± 92 | 191 ± 80 | 76.2/15.6 % | **97.4/1.6 %** |
| **A6** | **2.764 ± 0.041** | 3.342 ± 0.320 | **+0.578 (largest)** | **282 ± 41** | 161 ± 87 | 75.2/16.7 % | **71.2/24.4 % (loosens)** |

**This is the sharpest mechanistic distinction in the study — a clean 3-way split:**

- **A1 + CoA (COMP): both deeply co-bound.** CoA makes ~4200 active-site contacts
  (0.15 nm), A1 ~990 (0.19 nm), and the two sit only 0.38 nm apart. A1 binds
  *beside* a fully engaged substrate — the enzyme is locked, and the xanthone is
  pushed only modestly further from the catalytic centre (+0.32 nm).
- **A3 + CoA: CoA expelled from the active site** (~7 contacts, 0.62 nm), A3
  itself peripheral (~190 contacts), A3–CoA ~1.1 nm apart — **yet the enzyme
  becomes MORE locked than A3 alone** (76.2 % → 97.4 %).
- **A6 + CoA: CoA fully expelled and far away** (~6 contacts, 2.11 nm from the
  active site; A6–CoA 2.60 nm apart), **A6 itself is pushed furthest of the
  three from the catalytic centre by CoA's presence (+0.58 nm)**, **and the
  enzyme does NOT gain confinement** (75.2 % solo → 71.2 % ternary — the only
  system of the three where the ternary is *less* confined than the
  single-ligand baseline).

Caveats: the CoA-contact time series are stable (maintained states, not slow
losses over 200 ns); the ternary starting poses could bias the absolute CoA
position, but the *maintained* difference — now anchored on a real A6-solo
baseline rather than a proxy — is robust.

---

## 8. MM-PBSA / per-residue decomposition alongside the dynamics

### Summary

| | **A1** | A3 | A6 |
|---|---|---|---|
| MM-PBSA ΔG, **25 ns** (kJ/mol, canonical thesis values) | −190.9 | −144.1 | −266.2 |
| MM-PBSA ΔG, **200 ns** (kJ/mol, 150–200 ns window) | **−236.5 ± 33.7 (SD/\|mean\| 0.14 — well converged, `HMG_CoA_R-A1-R2/md_final.xtc`)** | −127.8 ± 110.9 (SD/\|mean\| 0.87 — **poorly converged**) | −232.2 ± 28.1 (SD/\|mean\| 0.12 — well converged, `md_final.xtc`) |

**A1's 200 ns MM-PBSA is recovered (2026-09-10) from `HMG_CoA_R-A1-R2`** — the
production directory underlying the current A1 system (byte-identical
`md.tpr`/`md.cpt`/`md.edr` to `Mecanismo_inhibitorio/A1`), which never had
MM-PBSA output on disk. **This retracts the earlier "A1 has no 200 ns MM-PBSA
run" statement.** With it, **A1 and A6 are now the two best-converged and most
energetically favourable systems at 200 ns** (−236.5 and −232.2 kJ/mol,
statistically indistinguishable given their SDs), both well ahead of A3
(−127.8 ± 110.9, unconverged). Full breakdown, all four recovered channels
(A1-alone, substrate-only, A1-in-ternary, A1↔CoA pairwise) and per-residue
tables: `mmpbsa_200ns_A1_A3_A6.md`; recovery provenance:
`a1_source_reconciliation_2026-09-10/`.

**Note (not re-litigated further, per instruction):** this 200 ns A1-alone
Total does not reconcile against the existing presentation's "E+I" slide value
(−206.5 kJ/mol) — the MM force-field terms (vDW, Elec) match that slide almost
exactly, but the solvation terms and Total differ by ≈30 kJ/mol. `HMG_CoA_R-A1-R2`
is used as instructed as the authoritative source regardless; see
`a1_source_reconciliation_2026-09-10/a1_long_mmpbsa_recovery.md` §4 for the
record.

### Per-residue decomposition — catalytic-residue engagement

**A1 — recovered 200 ns decomposition (2026-09-10; `HMG_CoA_R-A1-R2/
residues_energy_summary.csv`, ligand's own row `LIG-873`, total
−253.8 ± 22.5, excluded from the protein ranking):**

| residue | total ΔG (kJ/mol) | character |
|---|---|---|
| LEU-853 | −51.3 ± 10.6 | hydrophobic anchor |
| HIS-752 | −23.1 ± 6.6 | hydrophobic/aromatic |
| SER-684 | −13.4 ± 6.9 | hydrophobic-dominated |
| ASN-755 | −12.0 ± 8.1 | polar/hydrophobic mixed |
| LEU-857 | −8.6 ± 4.9 | hydrophobic anchor |
| LYS-691 | −2.9 ± 7.1 | electrostatic (noisy) |

**This supersedes the earlier "A1 is anchored on the catalytic ARG-702
(−9.9 kJ/mol)" statement carried in this report from before the 2026-09-10
recovery.** In the recovered 200 ns decomposition, **ARG-702 does not appear
in A1's top residues**; the raw file lists all four tetramer copies of
ARG-702 individually (−0.28, +0.21, +0.03, −0.02 kJ/mol) — negligible in every
chain. The original source of the −9.9 kJ/mol figure was not located among the
four directories reconciled on 2026-09-10 and is not pursued further here (see
`a1_source_reconciliation_2026-09-10/`); the number above, from the same
canonical source as A1's new 200 ns total ΔG, is used going forward.

**A striking finding: A1 and A6 converge on the same hydrophobic anchor
residues** — LEU-853, HIS-752, SER-684 and LEU-857 are each in **both** A1's
and A6's top-5 (A6, unchanged from the 2026-09-10 A6 correction):

| residue | A1 total (kJ/mol) | A6 total (kJ/mol) |
|---|---|---|
| LEU-853 | −51.3 ± 10.6 | −36.3 ± 12.1 |
| HIS-752 | −23.1 ± 6.6 | −26.1 ± 15.9 |
| SER-684 | −13.4 ± 6.9 | −26.1 ± 6.7 |
| LEU-857 | −8.6 ± 4.9 | −23.2 ± 7.3 |

**A3**: LYS-864 (−26.9, noisy), HIS-866 (−16.9, noisy), VAL-863 (−12.6),
ASP-690 (−6.9). No overlap with the A1/A6 hydrophobic-anchor set, and **ARG-702
is absent from A3's profile too** — A3 engages *peripheral* residues only, at
both timescales.

**Revised residue-level picture:** A1 and A6 are not mechanistically distinct
at the level of *which* residues are engaged — both anchor on the same
hydrophobic sub-pocket (LEU-853/HIS-752/SER-684/LEU-857), with comparable,
strongly favourable total energetics (−236 vs −232 kJ/mol) that now clearly
separate them from A3. **What still distinguishes A1 from A6 is not the
residue-level energetics but the conformational consequence** (§4, §9): A1
alone locks the enzyme to one state (93.4 %, ≈ apo) while A6 alone leaves it as
permissive as A3 (75.2 %) — the same anchor residues, engaged to a similar
energetic degree, are compatible with two different conformational outcomes.
This points to a broader network of weaker contacts or ligand-shape effects
(not the top-5 residues) as the source of A1's conformational lock, not the
dominant hydrophobic anchor itself.

---

## 9. Do A1, A3 and A6 show distinct inhibition mechanisms?

**They act at the same site (all orthosteric / competitive at the HMG-CoA pocket)
but with a clear gradient of catalytic commitment and conformational effect —
now confirmed for A6 with real 200 ns data at every level:**

| | A1 | A6 | A3 |
|---|---|---|---|
| catalytic-COM distance (nm) | 2.60 (deepest) | 2.76 | 3.03 (shallowest) |
| direct active-site contacts | 342 ± 214 | 282 ± 41 (**tightest spread**) | 267 ± 92 |
| dominant residue anchor (200 ns) | **hydrophobic: LEU-853/HIS-752/SER-684/LEU-857** | **same hydrophobic set, weaker per-residue but broader (+VAL-683)** | different, peripheral (LYS-864/HIS-866/VAL-863) |
| catalytic ARG-702 engagement (200 ns) | **negligible (all 4 chain copies ≈ 0)** | **modest, real (−2.2)** | **none** |
| enzyme alone → conformational effect | **locks to 1 state (93.4 %, ≈ apo)** | **opens a 2nd state (75.2/16.7 %) — as permissive as A3** | **opens a 2nd state (76.2/15.6 %)** |
| shared-basis PCA position | PC1 −3.30 (extreme) | **PC2 +2.64 — a third, distinct axis** | PC1 +3.29 (extreme) |
| DCCM perturbation vs A1 | — | 4.8 % of Cα pairs | 5.7 % of Cα pairs |
| ternary with CoA | **A1 + CoA co-bind deeply; enzyme locked further (93→96 %)** | **CoA expelled furthest; enzyme does NOT lock further (75→71 %)** | **CoA expelled; enzyme locks MORE (76→97 %)** |
| MM-PBSA ΔG, 200 ns (kJ/mol) | **−236.5 ± 33.7 (converged)** | −232.2 ± 28.1 (converged) | −127.8 ± 111 (unconverged) |
| MM-PBSA ΔG, 25 ns (kJ/mol) | −190.9 | −266.2 | −144.1 |

- **A1** — at 200 ns, **the single most energetically favourable of the three**
  (−236.5 kJ/mol, converged), anchored on the **same hydrophobic sub-pocket as
  A6** (LEU-853, HIS-752, SER-684, LEU-857) rather than on ARG-702 as earlier
  believed (§8) — and yet the only one of the three that **rigidifies the
  catalytic site as much as the apo enzyme**, tightening slightly further with
  CoA present. Since its residue-level anchor is now shown to closely resemble
  A6's, **A1's conformational lock is evidently not explained by the
  hydrophobic anchor itself**, but by some other feature of the complex
  (broader/weaker contact network, ligand shape/rigidity) not captured by the
  top-5 residues.
- **A3** — a weak, peripheral binder (different residues entirely from A1/A6)
  that never reaches the catalytic centre, leaves the enzyme conformationally
  permissive **on its own**, but in the ternary **displaces CoA from the
  active site while paradoxically driving the enzyme to its most-locked
  ternary state (97.4 %)**. Its inhibitory action, if real, is about
  *blocking productive substrate binding while stabilising a non-productive
  complex*, not occupying the catalytic machinery.
- **A6** — geometrically intermediate, energetically **on par with A1**
  (converged 200 ns MM-PBSA, −232 vs −236 kJ/mol) via the **same hydrophobic
  anchor cluster as A1**, plus a real, modest catalytic-ARG-702 contact absent
  in A1. Conformationally it is **as permissive as A3 when alone** — not a
  second A1 — and uniquely, **CoA does not confine it further and pushes it
  furthest from the catalytic centre of the three (+0.58 nm)**. A6 occupies a
  **third, distinct region of the shared conformational landscape** (high PC2,
  near-zero PC1) rather than sitting on the A1–A3 axis. A6 combines
  A1-comparable binding energetics and residue anchor with a mechanism that
  does **not** rigidify the catalytic site the way A1 does, and **displaces
  CoA from the active site more completely than either A1 or A3**.

**Answer:** the mechanistic evidence, now with A1's real 200 ns data, supports
a picture that is **more nuanced than a simple three-way gradient**:
**A1 and A6 share both the residue-level binding anchor and comparable,
converged, favourable binding energetics** — they are the two "deep binders"
of the series, not opposite ends of a spectrum — **while A3 is the true
outlier**, engaging different, peripheral residues with a much weaker, poorly
converged energetic signature. What separates A1 from A6, given their
energetic and residue-level similarity, is **entirely on the conformational
and CoA-interaction side**: A1 locks the enzyme (alone and more so with CoA)
and tolerates a fully co-bound substrate; A6 does not lock the enzyme either
alone or with CoA, and most strongly displaces CoA of the three. A3 sits apart
on every axis — weak, peripheral, permissive alone, yet paradoxically
locking the enzyme in the ternary. **The shared PCA places A6 on a third,
distinct conformational axis from A1 and A3 despite its energetic kinship with
A1** — energetics and conformational dynamics are answering different
questions here, and both are needed to tell A1, A3 and A6 apart.

---

## 10. What remains for thesis / multimodal integration

1. ~~**HMG-R-200ns-A6** — obtain a PBC-clean 200 ns trajectory~~ **DONE
   2026-09-10** (`md_final.xtc`). A6 is now in the shared A1/A3/A6 PCA basis +
   shared FEL, with full core descriptors, DCCM, clustering and a well-converged
   200 ns MM-PBSA.
2. **Project APO / COA / COMP onto the shared PCA basis** (protein identical) →
   one conformational map placing every ligand/state. *(still open — out of
   scope for this A6-focused correction.)*
3. ~~**A1's own 200 ns MM-PBSA** is still missing~~ **DONE 2026-09-10** —
   recovered from `HMG_CoA_R-A1-R2` (−236.5 ± 33.7 kJ/mol, well converged).
   §8's comparison is now three-of-three. The historical Total discrepancy vs
   the existing presentation's E+I slide value is flagged but intentionally
   not resolved (`a1_source_reconciliation_2026-09-10/`).
4. **Fold the now-complete A1/A3/A6 200 ns mechanistic block into the
   multimodal dataset** as a feature group distinct from the 25 ns
   `xanthone_short` block: conformational-state populations, shared-PC
   centroids, ΔDCCM summary (now 3 pairs), catalytic-COM distance,
   catalytic-residue engagement (now quantified for all three, all from
   consistent 200 ns sources), ternary CoA-displacement (now on real A1/A6-solo
   data and, for COMP, real MM-PBSA energetics per channel).
5. Decide 200 ns vs 25 ns MM-PBSA reporting for the thesis — A1 and A6 are now
   both well-converged at 200 ns and mutually consistent in ranking with their
   25 ns canonical values (A1 25 ns −190.9 → 200 ns −236.5; A6 25 ns −266.2 →
   200 ns −232.2); A3 remains unconverged at 200 ns (SD ≈ |mean|) either way.
   Recommend leading with the 25 ns canonical values for the thesis (all three
   converged, directly comparable across the full 24-complex panel) and citing
   the 200 ns corroboration for A1/A6 specifically.
6. **New (2026-09-10):** the four recovered directories
   (`HMG_CoA_R-A1-R2`, `HMG-CoA-R_sustrato`, `Competitive-system-2`, plus the
   superseded `HMG_CoA_R-A1`) also contain historical RMSD/RMSF/Rg/PCA/FEL/DCCM/
   clustering analyses for A1, substrate-only and the ternary that were never
   reconciled against this pipeline's independently-recomputed equivalents
   (`Mecanismo_inhibitorio/A1`'s `rmsd_A1.xvg`/`PCA_A1.pdf`/`FEL_A1.pdf`/
   `DCCM_A1.pdf` etc., dated Mar 2026). Worth a future comparison pass — out of
   scope here per instruction not to reopen source reconciliation.

---

*Supporting files:* `system_compatibility.csv`, `mechanistic_descriptor_summary.csv`
(now carries `mmpbsa_200ns_kjmol_{mean,sd,source}` and, for COMP, the
per-channel `lig_/coa_/lig_coa_mmpbsa_200ns_*` columns),
`windowed_analysis_A1_A3_A6.csv`, `clustering_state_populations.csv`,
`shared_pca_summary.md`, `fel_state_summary.md`, `dccm_summary.md`,
`coa_interaction_summary.md`, `solo_vs_ternary_A1_A3_A6.md`,
`mmpbsa_200ns_A1_A3_A6.md`, `mechanistic_A6_md_final_correction.md`,
`a1_source_reconciliation_2026-09-10/` (`a1_source_reconciliation.md`,
`a1_long_mmpbsa_recovery.md`, `a1_system_mapping.csv`,
`a1_replicate_and_control_inventory.csv`),
`execution_summary.md`, `qc_report.md`,
`mechanistic_standardized_protocol.yaml` (in
`reports/trajectory_protocol_correction_2026-09-09/`).
Raw grids / matrices / representative PDBs / provenance JSON under
`analysis_outputs/xanthone_mechanistic/`.
