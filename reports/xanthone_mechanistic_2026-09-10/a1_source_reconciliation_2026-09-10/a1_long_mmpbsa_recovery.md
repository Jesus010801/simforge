# A1 long-MD MM-PBSA recovery (2026-09-10)

**Scope:** the 4 directories named in this pass only —
`HMG_CoA_R-A1`, `HMG_CoA_R-A1-R2`, `HMG-CoA-R_sustrato`, `Competitive-system-2` —
reconciled against `Mecanismo_inhibitorio/{A1,APO,COA,COMP}`. No broader
filesystem audit performed. **Read-only**: nothing below was rerun; nothing in
any of these directories or in `Mecanismo_inhibitorio/` was modified.

## Headline answer

**The statement "A1 has no 200 ns MM-PBSA run" is FALSE.** It is true only for
`Mecanismo_inhibitorio/A1` specifically (that directory has no MM-PBSA output
files at all). A complete, 200 ns, well-formed MM-PBSA run for A1-alone exists
in `/home/jesusxd/Escritorio/HMG_CoA_R-A1-R2/` — a directory the previous
mechanistic pass never inspected (it only reprocessed
`Mecanismo_inhibitorio/A1/md.xtc`).

Three further complete 200 ns MM-PBSA results, all directly relevant to A1,
were also recovered from `HMG-CoA-R_sustrato/` and `Competitive-system-2/` —
these correspond to the **substrate-only (E+S)** and **ternary (E+S+I)**
systems and reproduce, almost to the last decimal, the MM-PBSA decomposition
table already printed on the existing presentation slide *"Comparación
estructural de A1 y HMG-CoA en el sitio activo"*.

## 1. A1-alone, 200 ns — `HMG_CoA_R-A1-R2/`

| | |
|---|---|
| source trajectory | `md_final.xtc` (per `energy_MM.xvg` command line) |
| topology | `md.tpr`, 284 597 atoms — **byte-size-identical** to `Mecanismo_inhibitorio/A1/md.tpr` (12 142 572 B both) |
| index | `index.ndx` (also `index_full.ndx`, `index_old.ndx`) |
| g_mmpbsa command | `g_mmpbsa run -f md_final.xtc -s md.tpr -b 150000 -e 200000 -dt 500 -i mmpbsa.mdp -pbsa -mm energy_MM.xvg -pol polar.xvg -apol apolar.xvg -decomp -mmcon contrib_MM.dat -pcon contrib_pol.dat -apcon contrib_apol.dat -incl_14 -n index.ndx -unit1 'group "Protein"' -unit2 'group "LIG"'` |
| window / stride / n | **150–200 ns**, dt = 500 ps, **101 configurations** |
| mmpbsa.mdp | sha256 `a08eb840aa52c8aa5bef668b2ce2f3526a7c51e3d662fa96be6758fd184ad0a7` — the *"ISCB+/A2 / Sustratos-inhibidores"* canonical variant (comment-character-only difference from the `954785fb1805…` variant used for A3/A6/`HMG-R-200ns-A6`; PB dielectrics pdie 2/sdie 78.3/vdie 1, chgm spl4, srfm smol, temp 309.65 K, gamma 0.0226778, sasaconst 3.84928 — physically identical) |

### Total ΔG (`energy_summary.csv`)

| term | mean (kJ/mol) | SD |
|---|---|---|
| vDW | −157.166 | 11.466 |
| Electrostatic | −44.399 | 16.394 |
| Polar-solvation | +102.289 | 15.533 |
| Non-polar-solvation | −137.239 | 22.242 |
| **Total** | **−236.516** | **33.740** |

SD/\|mean\| = 0.143 — **well converged** (comparable to A6's 200 ns run, 0.12).

### Per-residue decomposition

`residues_energy_summary.csv`, 1615 residues (1616 lines incl. header). Not
re-tabulated here (out of this pass's scope — file is intact and readable at
the path above); available for a future top-N pull the same way it was done
for A3/A6.

### Duration / integrity checks

- `md.mdp`: `integrator = md`, `nsteps = 200000000`, `dt = 0.001` → 200 000 ps =
  **200 ns**, standard classical GROMACS 2025.2 MD (comments in the mdp file
  itself are stale/mislabelled — "1000 ps (1 ns)" / "2 fs" — but do not match
  the actual computed duration or the observed data span; the real duration is
  confirmed independently by three consistent sources: nsteps×dt, the
  `binding_energy.xvg` time axis (ends 200000.000 ps), and `rmsd-A1.xvg`
  (ends 200.0 ns)).
- A QM/MM side-attempt exists in this directory (`qmmm.inp`, `qm_atoms.txt`,
  `QM_MM/`) but **`qmmm.out` shows an ABORTED CP2K run** ("invalid value for
  enumeration: CHARMm" / "GROMACS") — it never produced results and does not
  affect the classical MD or the MM-PBSA above, both computed with the
  standard `md` integrator.
- `gmx cluster` was already run in this directory: `cluster.log` reports
  **"Found 4 clusters"** (own clustering, not reconciled against the current
  pipeline's independently-computed A1 clustering from
  `Mecanismo_inhibitorio/A1` — see the reconciliation note below).

## 2. A1-in-ternary, CoA-in-ternary, A1↔CoA pairwise — `Competitive-system-2/`

Same `md_final.xtc` (8 428 580 692 B — **byte-size-identical** to
`Mecanismo_inhibitorio/COMP/md.xtc`) and `md_200ns.tpr`, three independent
`g_mmpbsa` invocations:

| channel | files | unit1 / unit2 | window | dt | n cfg | Total ΔG (kJ/mol) |
|---|---|---|---|---|---|---|
| A1 in ternary | `energy_summary-A1.csv`, `A1energy_MM.xvg ` (note trailing space in filename) | Protein / LIG | 150–200 ns | 100 | 501 | **−211.691 ± 26.273** |
| CoA in ternary | `COA-energy_summary.csv`, `COA-energy_MM.xvg` | Protein / COA | 150–200 ns | 500 | 101 | **−133.832 ± 58.683** |
| A1 ↔ CoA pairwise | `A1-COAenergy_summary.csv`, `A1-COAenergy_MM.xvg` | COA / LIG | 150–200 ns | 500 | 101 | **−27.410 ± 22.365** |

The pairwise channel's `A1-COAresidues_energy_summary.csv` (2 rows: `LIG-873`
−15.882, `COA-874` −11.528) sums to −27.410 — internally consistent with the
channel total.

**Note (flagged, not resolved in this pass):** the A1-in-ternary channel used
`-dt 100` (501 configurations) while the other two channels in the same
directory used `-dt 500` (101 configurations) — same 150–200 ns window, finer
sampling for A1 specifically. Not an error per se, but an inconsistency in
sampling density across channels worth flagging for anyone averaging or
comparing convergence across the three.

## 3. Substrate-only (E+S), 200 ns — `HMG-CoA-R_sustrato/`

| | |
|---|---|
| source trajectory | `md_final.xtc` |
| topology | `md.tpr`, 8 507 044 B — **byte-size-identical** to `Mecanismo_inhibitorio/COA/md.tpr` |
| g_mmpbsa command | `-b 150000 -e 200000 -dt 500`, unit1 Protein / unit2 COA |
| window / stride / n | 150–200 ns, 500 ps, 101 configurations |

| term | mean (kJ/mol) | SD |
|---|---|---|
| vDW | −115.888 | 24.827 |
| Electrostatic | −123.753 | 20.407 |
| Polar-solvation | +103.517 | 16.003 |
| Non-polar-solvation | −60.462 | 15.272 |
| **Total** | **−196.585** | **34.952** |

## 4. Reconciliation against the current presentation's MM-PBSA table

The existing presentation slide *"Comparación estructural de A1 y HMG-CoA en
el sitio activo"* carries a 5-row MM-PBSA decomposition table
(ΔEvdW / ΔEelec / ΔGtotal, kJ/mol) for E+I, E+S, E+S+I(I), E+S+I(S), I–S(E+S+I).
Comparing it against what was recovered here:

| row | slide ΔEvdW / ΔEelec / ΔGtotal | recovered (this pass) | match |
|---|---|---|---|
| E+I (A1 alone) | −157.2 / −44.4 / **−206.5** | `HMG_CoA_R-A1-R2`: −157.166 / −44.399 / **−236.516** | **vDW & Elec match to 0.1 kJ/mol; Total does NOT match (Δ ≈ 30 kJ/mol)** |
| E+S (substrate alone) | −115.9 / −123.8 / −196.6 | `HMG-CoA-R_sustrato`: −115.888 / −123.753 / −196.585 | **exact match** |
| E+S+I (I) (A1 in ternary) | −179.4 / −21.8 / −211.7 | `Competitive-system-2` A1-in-ternary: −179.427 / −21.768 / −211.691 | **exact match** |
| E+S+I (S) (CoA in ternary) | −100.3 / −102.5 / −133.8 | `Competitive-system-2` CoA-in-ternary: −100.293 / −102.501 / −133.832 | **exact match** |
| I–S (E+S+I) (A1↔CoA pairwise) | −27.0 / +1.7 / −27.4 | `Competitive-system-2` A1↔CoA: −27.011 / +1.734 / −27.410 | **exact match** |

**Four of the five rows match to 2–3 decimal places** — this confirms
`HMG-CoA-R_sustrato` and `Competitive-system-2` are indeed the authoritative
sources already used to build the existing presentation table, and their raw
files (not just the printed slide numbers) are now recovered and provenance-
tracked.

**The fifth row (E+I, A1 alone) does not reconcile on Total**, despite the MM
force-field terms (vDW, Elec — which depend only on the trajectory and force
field, not on the implicit-solvent model) matching almost exactly. Since vDW
and Elec are identical but Polar+Apolar solvation are not (implied: slide's
Polar+Apolar ≈ −4.9 kJ/mol combined vs `HMG_CoA_R-A1-R2`'s actual −34.95
kJ/mol combined), the discrepancy sits entirely in the **implicit-solvent
(APBS) step**, not in the classical MD. Two explanations are consistent with
the evidence and were NOT adjudicated in this pass (would require inspecting
further files/timestamps beyond the 4 directories in scope):

1. `HMG_CoA_R-A1-R2`'s MM-PBSA (dated **1 Mar 2026**) is a **later recomputation**
   superseding whatever produced the slide's −206.5 value, using the same
   trajectory/frames but a different APBS/solvation configuration or software
   state than whatever generated the original slide number; or
2. the slide's E+I row was sourced from a **different A1-alone run entirely**
   (not yet located within the 4 directories inspected here).

**This discrepancy is flagged, not resolved — no current report value should
be silently replaced by −236.516 without further review.**

## 5. What was NOT done in this pass

- No MM-PBSA was rerun.
- No file in `Mecanismo_inhibitorio/`, `HMG_CoA_R-A1`, `HMG_CoA_R-A1-R2`,
  `HMG-CoA-R_sustrato`, or `Competitive-system-2` was modified.
- `mechanistic_descriptor_summary.csv`, `final_mechanistic_report.md`,
  `mmpbsa_200ns_A1_A3_A6.md` and every other current mechanistic report file
  were **left untouched** — see `a1_source_reconciliation.md` §"Corrections
  required" for exactly what should change and where, once reviewed and
  approved.
