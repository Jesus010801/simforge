# CoA-system interaction decomposition

Three channels kept separate for COMP / system_A3_COA / system_A6_COA.

**2026-09-10 addition:** COMP now also has a real **MM-PBSA energetic**
decomposition for all three channels (not just geometric mindist/contacts),
recovered from `Competitive-system-2` (see
`a1_source_reconciliation_2026-09-10/`, `mmpbsa_200ns_A1_A3_A6.md`). This is
incorporated below per channel; it is **not** produced by `coa_decomposition.py`
(which computes geometry only) — re-running `build_reports.py`/`coa_decomposition.py`
will regenerate the geometric lines above this note but will not remove or
regenerate the MM-PBSA lines, which must be re-added from
`mmpbsa_200ns_incorporate.py` if this file is rebuilt from scratch.
**system_A3_COA / system_A6_COA have no equivalent energetic decomposition** —
only the geometric channels below, unchanged.

## COMP
- protein–xanthone: active-site mindist 0.188 ± 0.016 nm, contacts 985.510 ± 111.343, catalytic COM 2.924 ± 0.117 nm
- protein–CoA: active-site mindist 0.149 ± 0.005 nm, contacts 4173.646 ± 326.519, catalytic COM 3.977 ± 0.128 nm
- xanthone–CoA: mindist 0.379 ± 0.130 nm, contacts(<0.6 nm) 87.223 ± 88.742, COM distance 1.347 ± 0.155 nm
- **MM-PBSA (200 ns, 150–200 ns window, `Competitive-system-2`):**
  protein–xanthone(A1) ΔG **−211.7 ± 26.3 kJ/mol** (dt 100 ps, 501 cfg);
  protein–CoA ΔG **−133.8 ± 58.7 kJ/mol** (dt 500 ps, 101 cfg);
  xanthone(A1)–CoA pairwise ΔG **−27.4 ± 22.4 kJ/mol** (dt 500 ps, 101 cfg;
  additive per-ligand share: A1 −15.9 ± 8.0, CoA −11.5 ± 21.1).
  Energetically, A1 binds the enzyme roughly as strongly as CoA binds the
  enzyme when both are present together (−211.7 vs −133.8 kJ/mol) while the
  two ligands' direct pairwise interaction is comparatively weak (−27.4
  kJ/mol) — consistent with the geometric picture (both close to the protein,
  0.38 nm from each other) but quantifies for the first time that A1's own
  enzyme engagement dominates over the A1–CoA contact itself.

## system_A3_COA
- protein–xanthone: active-site mindist 0.268 ± 0.053 nm, contacts 191.292 ± 80.014, catalytic COM 3.323 ± 0.205 nm
- protein–CoA: active-site mindist 0.616 ± 0.133 nm, contacts 6.671 ± 22.227, catalytic COM 3.806 ± 0.207 nm
- xanthone–CoA: mindist 1.086 ± 0.286 nm, contacts(<0.6 nm) 5.416 ± 22.633, COM distance 2.001 ± 0.356 nm

## system_A6_COA
- protein–xanthone: active-site mindist 0.291 ± 0.076 nm, contacts 160.564 ± 87.073, catalytic COM 3.342 ± 0.320 nm
- protein–CoA: active-site mindist 2.111 ± 0.910 nm, contacts 5.881 ± 20.457, catalytic COM 4.412 ± 0.887 nm
- xanthone–CoA: mindist 2.600 ± 1.286 nm, contacts(<0.6 nm) 13.725 ± 38.980, COM distance 3.464 ± 1.290 nm
