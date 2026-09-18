# Single-ligand vs ternary (E+I+CoA) — A1 / A3 / A6 (2026-09-10 update)

A6-solo now has full 200 ns data (`HMG-R-200ns-A6`, `md_final.xtc`), replacing the 25 ns proxy used previously. Same clustering method (gmx cluster gromos, Cα/Backbone RMSD, 0.15 nm cutoff) for solo and ternary.

| label | solo: catalytic COM (nm) | ternary: xanthone catalytic COM (nm) | Δ (ternary−solo) | solo: active-site contacts | ternary: active-site contacts | solo: clusters (top state %) | ternary: clusters (top state %) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A1 | 2.603±0.099 | 2.924±0.117 | +0.321 | 342±214 | 986±111 | 3 (93.4/3.9%) | 4 (95.7/3.2%) |
| A3 | 3.027±0.127 | 3.323±0.205 | +0.296 | 267±92 | 191±80 | 5 (76.2/15.6%) | 3 (97.4/1.6%) |
| A6 | 2.764±0.041 | 3.342±0.320 | +0.578 | 282±41 | 161±87 | 5 (75.2/16.7%) | 5 (71.2/24.4%) |