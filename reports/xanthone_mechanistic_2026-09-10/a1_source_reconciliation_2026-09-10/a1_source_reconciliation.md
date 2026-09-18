# A1 source reconciliation (2026-09-10)

**Scope (as instructed):** inspect only
`/home/jesusxd/Escritorio/{HMG_CoA_R-A1, HMG_CoA_R-A1-R2, HMG-CoA-R_sustrato,
Competitive-system-2}` and reconcile them against the current
`Mecanismo_inhibitorio/{A1,APO,COA,COMP}` mechanistic records. No broader
filesystem audit was performed. **Read-only** — nothing in any of these six
directories was modified, and no current mechanistic result
(`final_mechanistic_report.md`, `mechanistic_descriptor_summary.csv`, etc.) was
overwritten. This document, `a1_long_mmpbsa_recovery.md`, `a1_system_mapping.csv`
and `a1_replicate_and_control_inventory.csv` are the only new artifacts.

---

## 1. Identity of each directory

| directory | biological system |
|---|---|
| `HMG_CoA_R-A1` | **A1 alone** — HMG-CoA reductase + xanthone A1, first/preliminary production (100 ns, smaller box) |
| `HMG_CoA_R-A1-R2` | **A1 alone, replicate/refined production** — same biological system, 200 ns, full box; this is the run underlying the current `Mecanismo_inhibitorio/A1` |
| `HMG-CoA-R_sustrato` | **substrate-only** — HMG-CoA reductase + HMG-CoA/CoA substrate, no inhibitor (E+S) |
| `Competitive-system-2` | **competitive ternary** — HMG-CoA reductase + xanthone A1 + CoA substrate, simultaneously bound (E+S+I) |

Full per-directory inventory: `a1_system_mapping.csv`,
`a1_replicate_and_control_inventory.csv`.

---

## 2. Reconciliation against `Mecanismo_inhibitorio/{A1,APO,COA,COMP}`

| current dir | original production folder | status |
|---|---|---|
| `Mecanismo_inhibitorio/A1` | `HMG_CoA_R-A1-R2` | **processed/partial copy** of the R2 production: `md.tpr` (12 142 572 B), `md.cpt` (6 832 564 B) and `md.edr` (13 441 484 B) are byte-size-identical between the two directories, and both report **284 597 atoms**. `Mecanismo_inhibitorio/A1/md.xtc` (2.15 GB) is a much smaller **stripped/reduced trajectory** than `HMG_CoA_R-A1-R2/md.xtc` (21.37 GB raw) — almost certainly a water-stripped subset of the same production, not a separate run. `Mecanismo_inhibitorio/A1` also already carries a full historical RMSD/RMSF/Rg/SASA/H-bond/PCA/FEL/DCCM/site-descriptor analysis (`rmsd_A1.xvg`, `rmsf_A1.xvg`, `rg_A1.xvg`, `sasa_A1.xvg`, `hb_A1.xvg`, `pc1/2-A1.xvg`, `PCA_A1.pdf`, `FEL_A1.pdf`, `DCCM_A1.pdf`, `mindist_lig_active.xvg`, `contacts_lig_active.xvg`, `dist_lig_catalytic.xvg`, dated Mar 2026) — **but no MM-PBSA files of any kind.** The MM-PBSA that exists for this exact production sits only in `HMG_CoA_R-A1-R2`, never copied over. |
| `Mecanismo_inhibitorio/APO` | *(none of the 4 directories)* | no candidate directory for the apo system was among the 4 inspected; `Mecanismo_inhibitorio/APO` remains the only known source — unaffected by this reconciliation. |
| `Mecanismo_inhibitorio/COA` | `HMG-CoA-R_sustrato` | `md.tpr` is byte-size-identical (8 507 044 B) between the two. `md_200ns.xtc` (4.95 GB, `HMG-CoA-R_sustrato`) and `md.xtc` (4.89 GB, `Mecanismo_inhibitorio/COA`) are close but **not** byte-identical in size — same production, likely a slightly different concatenation/frame-count of the same continuation chain (`md.cpt`/`md_continue.*` present in the source), not two independent runs. `HMG-CoA-R_sustrato` carries a complete 200 ns MM-PBSA that `Mecanismo_inhibitorio/COA` does not have on disk. |
| `Mecanismo_inhibitorio/COMP` | `Competitive-system-2` | `md_final.xtc` (`Competitive-system-2`) and `md.xtc` (`Mecanismo_inhibitorio/COMP`) are **byte-identical in size** (8 428 580 692 B) — for all practical purposes the same file / a direct copy. `Competitive-system-2` carries three complete 200 ns MM-PBSA decompositions (A1-in-ternary, CoA-in-ternary, A1↔CoA pairwise) that `Mecanismo_inhibitorio/COMP` does not have on disk. |

**Pattern across all three matched systems:** `Mecanismo_inhibitorio/{A1,COA,COMP}`
hold the *dynamical* analysis (or, for A1, a stripped trajectory plus a full
historical dynamical analysis) but **none of the three carry any MM-PBSA
output** — the MM-PBSA for all three was computed directly in the original
production directories (`HMG_CoA_R-A1-R2`, `HMG-CoA-R_sustrato`,
`Competitive-system-2`) and never copied into `Mecanismo_inhibitorio/`. This is
why the 2026-09-10 mechanistic pipeline (which only ever read from
`Mecanismo_inhibitorio/` and `Nuevos_sistemas/`) reported "A1 has no 200 ns
MM-PBSA run" — true for the directory it looked in, false for the underlying
production.

**`HMG_CoA_R-A1` (no "-R2" suffix)** does not map to any current
`Mecanismo_inhibitorio/` directory — it is a **preliminary, superseded, 100 ns,
smaller-box** A1 production with no MM-PBSA, extensive but exploratory
FEL/clustering work, and a distinct (non-matching) `md.tpr`. Nothing in
`Mecanismo_inhibitorio/` derives from it.

---

## 3. Is "A1 has no 200 ns MM-PBSA run" false?

**Yes, false.** A complete, well-converged 200 ns MM-PBSA for A1-alone exists
in `HMG_CoA_R-A1-R2/` (ΔG = −236.516 ± 33.740 kJ/mol, SD/\|mean\| = 0.14).
Two more directly A1-relevant 200 ns MM-PBSA results (A1-in-ternary and
A1↔CoA-pairwise) exist in `Competitive-system-2/`. Full recovery, including
the reconciliation of these numbers against the current presentation's
existing MM-PBSA decomposition slide (4 of 5 rows match exactly; the A1-alone
row does not — flagged, not resolved): see `a1_long_mmpbsa_recovery.md`.

---

## 4. What current report statements must be corrected

These are **not applied yet** (task 6: "do not overwrite current mechanistic
results") — listed here for review and later action:

1. **`final_mechanistic_report.md` §8, §10.3** currently states *"A1: no 200 ns
   MM-PBSA on disk; the 25 ns canonical value stands"* and lists this as an
   open item ("A1's own 200 ns MM-PBSA is still missing... consider running it
   for full three-way parity"). **This is incorrect as written** — a 200 ns
   A1-alone MM-PBSA exists (`HMG_CoA_R-A1-R2`); it does not need to be *run*,
   it needs to be *recovered, reconciled and, if trusted, incorporated* the
   same way A6's was.
2. **`mmpbsa_200ns_A1_A3_A6.md`** (`SRC` dict in `mmpbsa_200ns_incorporate.py`)
   only defines sources for A3 and A6; A1 should be added once the E+I
   discrepancy in §4 of `a1_long_mmpbsa_recovery.md` is resolved (i.e., once it
   is determined whether −236.5 or the slide's −206.5 kJ/mol is the value to
   carry forward, or whether a still-different A1-alone run is the true
   source).
3. **`mechanistic_descriptor_summary.csv` / the mechanistic pipeline's `SYS`
   dict** (`mech_common.py`) treats `Mecanismo_inhibitorio/A1` as the sole A1
   source and reprocesses its `md.xtc` from scratch. Given
   `Mecanismo_inhibitorio/A1` already contains a complete historical
   RMSD/RMSF/Rg/SASA/HB/PCA/FEL/DCCM/site-descriptor analysis (dated Mar 2026,
   apparently produced directly from `HMG_CoA_R-A1-R2`), the pipeline's
   from-scratch recomputation is not wrong, but it is **duplicated effort**
   against an existing historical analysis that was never reconciled against
   it — worth a future comparison pass (out of scope here).
4. **The existing presentation's MM-PBSA decomposition slide** (E+I / E+S /
   E+S+I(I) / E+S+I(S) / I–S) is now traced to exact source files for 4 of its
   5 rows (`HMG-CoA-R_sustrato`, `Competitive-system-2`) — worth citing that
   provenance explicitly in the thesis/slide notes. The 5th row (E+I) has a
   located candidate source (`HMG_CoA_R-A1-R2`) whose MM force-field terms
   match but whose solvation terms and total do not — **do not cite
   `HMG_CoA_R-A1-R2` as the confirmed source of the slide's −206.5 kJ/mol
   value** until that is resolved.
5. **A new, previously-undocumented MM-PBSA channel exists**: A1↔CoA pairwise
   interaction energy (`Competitive-system-2/A1-COA*`, ΔG = −27.410 ±
   22.365 kJ/mol) — an actual energetic decomposition of the xanthone–CoA
   channel that the current mechanistic pipeline (`coa_decomposition.py`) only
   characterises geometrically (mindist/contacts/COM-distance), never
   energetically. Worth incorporating into `coa_interaction_summary.md` in a
   future pass.

---

## 5. Final answers (task 8)

- **Does A1 200 ns MM-PBSA exist?** **Yes** — in `HMG_CoA_R-A1-R2/`
  (ΔG = −236.516 ± 33.740 kJ/mol, 150–200 ns, dt 500 ps, 101 configs,
  well-converged, canonical `mmpbsa.mdp`). The current claim that it does not
  exist is **false** and refers only to the absence of MM-PBSA files inside
  `Mecanismo_inhibitorio/A1`.
- **Should A1-R2 be treated as a replicate or excluded?** **Neither, exactly**
  — `HMG_CoA_R-A1-R2` is not an independent statistical replicate to be
  averaged with a separate A1 run; its `md.tpr`/`md.cpt`/`md.edr` are
  byte-size-identical to `Mecanismo_inhibitorio/A1`, so it should be treated
  as **the same underlying production as the current A1 record** — the
  authoritative location that additionally holds the MM-PBSA and clustering
  outputs `Mecanismo_inhibitorio/A1` lacks. It should **not** be excluded; it
  should be **reconciled in** (once the E+I discrepancy above is settled).
  `HMG_CoA_R-A1` (no suffix, 100 ns, smaller box, no MM-PBSA) **should be
  excluded** from the main mechanistic comparison — it is a superseded,
  preliminary production.
- **Canonical source for A1-alone:** `HMG_CoA_R-A1-R2` (200 ns, matches
  `Mecanismo_inhibitorio/A1`'s topology exactly, carries the only known
  200 ns A1-alone MM-PBSA).
- **Canonical source for substrate-only:** `HMG-CoA-R_sustrato` (matches
  `Mecanismo_inhibitorio/COA`'s topology exactly; its MM-PBSA reproduces the
  existing presentation's E+S row exactly).
- **Canonical source for competitive A1+substrate:** `Competitive-system-2`
  (its `md_final.xtc` is byte-identical in size to
  `Mecanismo_inhibitorio/COMP/md.xtc`; its three MM-PBSA channels reproduce
  the existing presentation's E+S+I(I), E+S+I(S) and I–S rows exactly).
- **What must be corrected:** see §4 above — principally, the
  `final_mechanistic_report.md` statement that A1 has no 200 ns MM-PBSA, and
  the omission of A1 from `mmpbsa_200ns_A1_A3_A6.md`, both pending resolution
  of the E+I total discrepancy (§4 of `a1_long_mmpbsa_recovery.md`) before any
  number is written into the accepted mechanistic record.

**No current mechanistic result was changed by this pass.** This is a
provenance/reconciliation layer only, per instruction.
