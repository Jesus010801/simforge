# Short-study (A6) provenance audit — analysis trajectory

**Date:** 2026-09-09
**Scope:** the 20 `xanthone_short` A6 short-study outputs generated 2026-09-08
(protocol `xanthone_short/1.0`, commit `2d2501d`).
**Question:** was each output computed from raw `md.xtc` or from the canonical
PBC-corrected rot+trans fit `mdfit.xtc`?

## Method

Each output directory carries a machine-written provenance record
(`analysis_outputs/xanthone_short/<system>/observables/<task>/provenance.<analysis>.json`).
The `source_trajectory` field and the `command_argv` `-f` argument were read
directly from all 20 records (archived verbatim at
`analysis_outputs/xanthone_short/_invalidated_md_xtc_2026-09-08/`).

## Result — all 20 outputs used raw `md.xtc`

| # | system | analysis | `source_trajectory` | `command_argv -f` | window | groups (by name) | verdict |
|---|---|---|---|---|---|---|---|
| 1 | AA-A6 | protein RMSD | `md.xtc` | `.../AA-A6/md.xtc` | -b 0 -e 25000 | Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 2 | AA-A6 | ligand RMSD | `md.xtc` | `.../AA-A6/md.xtc` | -b 0 -e 25000 | LIG vs Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 3 | AA-A6 | active-site min distance | `md.xtc` | `.../AA-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_AA | `INVALIDATED_WRONG_TRAJECTORY` |
| 4 | AA-A6 | active-site contacts (<0.6 nm) | `md.xtc` | `.../AA-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_AA | `INVALIDATED_WRONG_TRAJECTORY` |
| 5 | AA-A6 | catalytic-site COM distance | `md.xtc` | `.../AA-A6/md.xtc` | -b 0 -e 25000 | COM(LIG)–COM(Catalytic_AA) | `INVALIDATED_WRONG_TRAJECTORY` |
| 6 | AG-A6 | protein RMSD | `md.xtc` | `.../AG-A6/md.xtc` | -b 0 -e 25000 | Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 7 | AG-A6 | ligand RMSD | `md.xtc` | `.../AG-A6/md.xtc` | -b 0 -e 25000 | LIG vs Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 8 | AG-A6 | active-site min distance | `md.xtc` | `.../AG-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_AG | `INVALIDATED_WRONG_TRAJECTORY` |
| 9 | AG-A6 | active-site contacts (<0.6 nm) | `md.xtc` | `.../AG-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_AG | `INVALIDATED_WRONG_TRAJECTORY` |
| 10 | AG-A6 | catalytic-site COM distance | `md.xtc` | `.../AG-A6/md.xtc` | -b 0 -e 25000 | COM(LIG)–COM(Catalytic_AG) | `INVALIDATED_WRONG_TRAJECTORY` |
| 11 | HMG-R-25ns-A6 | protein RMSD | `md.xtc` | `.../HMG-R-25ns-A6/md.xtc` | -b 0 -e 25000 | Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 12 | HMG-R-25ns-A6 | ligand RMSD | `md.xtc` | `.../HMG-R-25ns-A6/md.xtc` | -b 0 -e 25000 | LIG vs Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 13 | HMG-R-25ns-A6 | active-site min distance | `md.xtc` | `.../HMG-R-25ns-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_HMG | `INVALIDATED_WRONG_TRAJECTORY` |
| 14 | HMG-R-25ns-A6 | active-site contacts (<0.6 nm) | `md.xtc` | `.../HMG-R-25ns-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_HMG | `INVALIDATED_WRONG_TRAJECTORY` |
| 15 | HMG-R-25ns-A6 | catalytic-site COM distance | `md.xtc` | `.../HMG-R-25ns-A6/md.xtc` | -b 0 -e 25000 | COM(LIG)–COM(Catalytic_HMG) | `INVALIDATED_WRONG_TRAJECTORY` |
| 16 | LP-A6 | protein RMSD | `md.xtc` | `.../LP-A6/md.xtc` | -b 0 -e 25000 | Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 17 | LP-A6 | ligand RMSD | `md.xtc` | `.../LP-A6/md.xtc` | -b 0 -e 25000 | LIG vs Protein | `INVALIDATED_WRONG_TRAJECTORY` |
| 18 | LP-A6 | active-site min distance | `md.xtc` | `.../LP-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_LIP | `INVALIDATED_WRONG_TRAJECTORY` |
| 19 | LP-A6 | active-site contacts (<0.6 nm) | `md.xtc` | `.../LP-A6/md.xtc` | -b 0 -e 25000 | LIG–ActiveSite_LIP | `INVALIDATED_WRONG_TRAJECTORY` |
| 20 | LP-A6 | catalytic-site COM distance | `md.xtc` | `.../LP-A6/md.xtc` | -b 0 -e 25000 | COM(LIG)–COM(Catalytic_LP) | `INVALIDATED_WRONG_TRAJECTORY` |

**0 of 20** outputs used `mdfit.xtc`. **All 20** are classified
`INVALIDATED_WRONG_TRAJECTORY` and must be recomputed on `mdfit.xtc`.

## Secondary finding — HMG-R-25ns-A6 measured span

The 5 HMG-R-25ns-A6 records also carry `trajectory_measured_end_ps = 200000`:
the analysis was pointed at the **full 200 ns** raw production trajectory
(`md.xtc`, 9.9 GB, 200 000 ps). The `-e 25000` argument still bounded the
computation to 0–25 ns, so the *window* was correct, but the *trajectory* was
both raw and the 200 ns file. The system's `mdfit.xtc` is already the 0–25 ns
fitted window (2501 frames, 0–25000 ps, dt 10 ps), which the corrected run uses.

## Semantic groups (unchanged by the correction)

Group resolution was and remains by **name** from each system's `index.ndx`
(numeric ids recorded for provenance only, never fed to GROMACS):

| system | ligand | protein | active site | catalytic site | LIG id | site id | cat id |
|---|---|---|---|---|---|---|---|
| AA-A6 | LIG | Protein | ActiveSite_AA | Catalytic_AA | 13 | 21 | 22 |
| AG-A6 | LIG | Protein | ActiveSite_AG | Catalytic_AG | 13 | 21 | 22 |
| HMG-R-25ns-A6 | LIG | Protein | ActiveSite_HMG | Catalytic_HMG | 13 | 21 | 22 |
| LP-A6 | LIG | Protein | ActiveSite_LIP | Catalytic_LP | 13 | 21 | 22 |

## Source fingerprints recorded in the old provenance (for integrity checks)

| system | file | bytes | sha256 |
|---|---|---|---|
| AA-A6 | md.xtc | 516 690 612 | `69c040d2d260a1f046f0428647883477fa1876c513324721941f888f7961f496` |
| AA-A6 | md.tpr | 3 041 296 | `42a573952c64dddd2277015ec9d808a8d5f94f9490375766ea4bb1b2f4bcf17f` |
| AA-A6 | index.ndx | 1 686 006 | `0dbdfd727b2d8a51aadb1eefa9a259b052f7af08580a001016b1669063a985b2` |
| AG-A6 | md.xtc | 1 148 822 852 | `6749c8562f5e4f0d5a81c2112c24b3d54a2f83ae51f2ca1a1915348632e402e2` |
| AG-A6 | md.tpr | 6 247 368 | `155ca5e568a66bfb2da96798fdc01f74d68b2ee1646604aa553d597a807caec5` |
| AG-A6 | index.ndx | 3 887 924 | `5f3891e5cd308a5f92d8237fdfcd83f3368d4d5eff8910bccefe66fe4f245b1c` |
| HMG-R-25ns-A6 | md.xtc | 9 899 197 968 | *(>2 GiB — size+mtime fingerprint only)* |
| HMG-R-25ns-A6 | md.tpr | 8 491 596 | `dbfcb22d6c21ba6c37f319cead9a0c416cf9a6e94e2ce90ce9ad3ca8ee7948f8` |
| HMG-R-25ns-A6 | index.ndx | 4 308 079 | `a1b2e8d51cb06629bab9551d963949931acefaecb6308fde86708492afd92164` |
| LP-A6 | md.xtc | 521 771 640 | `534c9917724352dc7f9bf26f73f887047f767bc299eb0a5159e2e7dfddb84a47` |
| LP-A6 | md.tpr | 2 881 992 | `179134773db1d36d894ded4ff701ea04c14ffb7965c12cc099afe047b866528f` |
| LP-A6 | index.ndx | 1 696 880 | `89c54716349d5ff90ff9120afa6e81b6220bbbcf75be3983d4e875ae6f388048` |

These are cross-checked against the live files in Part 6 (source-integrity
verification) — the raw trajectories, TPRs and index files must be byte-identical
before and after the recompute.
