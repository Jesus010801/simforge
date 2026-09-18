# Execution summary — mechanistic A1/A3/A6 pass (2026-09-10)

Canonical trajectory: **mdfit.xtc** (PBC-corrected rot+trans fit), resampled to a common **200 ps** stride (0–200 ns), for every system **except HMG-R-200ns-A6**. Raw `md.xtc` never analysed for the reference systems. References (A1/APO/COA/COMP) had no `mdfit.xtc` — reprocessed into the managed tree (`md.xtc → trjconv -pbc mol -center -ur compact → trjconv -fit rot+trans`); source trees read-only.

**2026-09-10 correction — HMG-R-200ns-A6:** the earlier `mdfit.xtc`/`mdcenter.xtc` for this system were `-pbc res`-built and split the tetramer (see the superseded REVIEW_REQUIRED item below). The user supplied a validated, PBC-clean replacement, `Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc` (2001 frames, native 100 ps stride, 0–200 ns; sha256 `9309606e79aece25fe44eed546e787718e8c49b2a572b96787578daf5d4a721d`). Sanity: Backbone RMSD mean 0.196 / max 0.255 nm, Rg 3.318 ± 0.008 nm — clean. This trajectory is used **directly** for every HMG-R-200ns-A6 analysis below (no new fitted trajectory written; `gmx rms`/`rmsf`/`covar` perform their usual on-the-fly least-squares fit, as for every other system). `HMG-R-200ns-A6` is now **fully complete** (was REVIEW_REQUIRED) and included in the shared A1/A3/A6 PCA basis, shared FEL, DCCM comparisons and clustering. See `mechanistic_A6_md_final_correction.md` for the full account.

## Analyses completed per system

| system | label | role | core descriptors | PCA | FEL | DCCM | clustering | macro-states | windowed | CoA channels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1 | A1 | reference | 9 xvg | yes | yes | yes | yes | yes (NEW) | yes | — |
| APO |  | reference | 5 xvg | yes | yes | yes | yes | yes (NEW) | yes | — |
| COA |  | reference | 9 xvg | yes | yes | yes | yes | yes (NEW) | yes | n/a (CoA = ligand) |
| COMP |  | reference | 16 xvg | yes | yes | yes | yes | yes (NEW) | yes | yes (3-way) |
| A3-HMG-R | A3 | new | 9 xvg | yes | yes | yes | yes | yes (NEW) | yes | — |
| HMG-R-200ns-A6 | A6 | new | 9 xvg | yes | yes | yes | yes | yes (NEW) | yes | — |
| system_A3_COA |  | new | 15 xvg | yes | yes | yes | yes | yes (NEW) | yes | yes (3-way) |
| system_A6_COA |  | new | 15 xvg | yes | yes | yes | yes | yes (NEW) | yes | yes (3-way) |

## REVIEW_REQUIRED items

- **system_A3_COA / system_A6_COA**: source `index.ndx` lacked `ActiveSite_HMG` / `Catalytic_HMG`. Resolved by deriving those groups from A1 (protein is byte-identical: 1614 Cα, 24199 atoms, same numbering) into the managed index. Recorded as *derived*.
- **APO**: no `index.ndx` in source → `gmx make_ndx` from `md.tpr` into the managed tree (default groups; apo enzyme → no ligand/site observables).
- **HMG-R-200ns-A6 — RESOLVED 2026-09-10**: the original pre-made trajectory derivatives (`mdfit.xtc`, `mdcenter.xtc`) were built with `-pbc res` and irreversibly split the tetramer (protein RMSD 5-7 nm, Rg doubles — an imaging artefact, not real dynamics). The user supplied `md_final.xtc`, a validated PBC-clean replacement; the full 200 ns whole-protein suite (core descriptors, PCA/FEL/DCCM, clustering, macro-states, windowed) is now computed from it, and A6 is included in the shared A1/A3/A6 PCA basis. See `HMG-R-200ns-A6/REVIEW_REQUIRED.md` for the original finding and `mechanistic_A6_md_final_correction.md` for the resolution. The closed short-study `HMG-R-25ns-A6` mdfit.xtc was never affected (clean).
- **ligand–active-site H-bonds**: `gmx hbond` reports the LigParGen xanthone `LIG` as 0 donors / 0 acceptors (non-standard atom names) → ligand–site H-bond count is not available from `gmx hbond`; ligand–site contact/distance metrics used instead.

See `system_compatibility.csv` for the full pre-flight.