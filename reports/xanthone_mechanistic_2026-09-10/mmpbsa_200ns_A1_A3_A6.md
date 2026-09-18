# MM-PBSA, 200 ns window — A1 / A3 / A6 and HMG-CoA-reductase mechanistic controls

Incorporated from already-completed `g_mmpbsa` runs; **not rerun** by this pass.
**A1-alone now has a canonical 200 ns value** (source: `HMG_CoA_R-A1-R2`, recovered 2026-09-10 — see `a1_source_reconciliation_2026-09-10/`). The earlier statement "A1 has no 200 ns MM-PBSA run" is **retracted**: it was true only for `Mecanismo_inhibitorio/A1`, which has no MM-PBSA output on disk; the run exists in `HMG_CoA_R-A1-R2`, the production directory underlying that same A1 system (byte-identical `md.tpr`/`md.cpt`/`md.edr`).

## A1 / A3 / A6 — single-ligand, 200 ns

| system | ΔG total (kJ/mol) | vDW | Elec | Polar-solv | Non-polar-solv | window | stride | n cfg | source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1 | **-236.5 ± 33.7** | -157.2 | -44.4 | 102.3 | -137.2 | 150-200 ns | 500 ps | 101 | `/home/jesusxd/Escritorio/HMG_CoA_R-A1-R2` (md_final.xtc (HMG_CoA_R-A1-R2 — canonical A1-alone source per 2026-09-10 reconciliation; supersedes the earlier 'no 200 ns MM-PBSA on disk' statement, which was true only for Mecanismo_inhibitorio/A1)) |
| A3 | **-127.8 ± 110.9** | -81.4 | -16.8 | 48.4 | -78.0 | 150-200 ns | 200 ps | 251 | `/home/jesusxd/Escritorio/Nuevos_sistemas/A3-HMG-R` (mdfit.xtc (A3's own, PBC-clean)) |
| A6 | **-232.2 ± 28.1** | -163.7 | -37.0 | 116.4 | -147.9 | 150-200 ns | 500 ps | 101 | `/home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6` (md_final.xtc) |

## Convergence

- **A1**: SD/|mean| = 0.14 — well converged.
- **A3**: SD/|mean| = 0.87 — poorly converged (SD approaching |mean|).
- **A6**: SD/|mean| = 0.12 — well converged.

A1 (0.14) and A6 (0.12) are both well converged and now directly comparable at 200 ns; A3 (0.87) remains poorly converged — its 200 ns value is reported for completeness but should not be over-interpreted on its own.

## HMG-CoA-reductase mechanistic controls — substrate-only and competitive ternary (200 ns; recovered 2026-09-10, see `a1_source_reconciliation_2026-09-10/`)

| channel | ΔG total (kJ/mol) | vDW | Elec | Polar-solv | Non-polar-solv | window | stride | n cfg | source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| substrate-only (E+S) | **-196.6 ± 35.0** | -115.9 | -123.8 | 103.5 | -60.5 | 150-200 ns | 500 ps | 101 | `/home/jesusxd/Escritorio/HMG-CoA-R_sustrato` |
| A1 in ternary (E+S+I, I channel) | **-211.7 ± 26.3** | -179.4 | -21.8 | 86.2 | -96.7 | 150-200 ns | 100 ps | 501 | `/home/jesusxd/Escritorio/Competitive-system-2` |
| CoA in ternary (E+S+I, S channel) | **-133.8 ± 58.7** | -100.3 | -102.5 | 148.1 | -79.2 | 150-200 ns | 500 ps | 101 | `/home/jesusxd/Escritorio/Competitive-system-2` |
| A1<->CoA pairwise (I-S, E+S+I) | **-27.4 ± 22.4** | -27.0 | 1.7 | 22.3 | -24.4 | 150-200 ns | 500 ps | 101 | `/home/jesusxd/Escritorio/Competitive-system-2` |

These reproduce (4 of 4 checked rows, to 2-3 decimals) the existing presentation's MM-PBSA decomposition table for E+S, E+S+I(I), E+S+I(S) and I-S — see `a1_source_reconciliation_2026-09-10/a1_long_mmpbsa_recovery.md` §4 for the row-by-row check.

## A1 — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)

*(ligand's own decomposition row, LIG-873: total -253.8 ± 22.5 — excluded from the ranking below, which is protein residues only)*

| residue | total | vdW | elec | polar | apolar |
| --- | --- | --- | --- | --- | --- |
| LEU-853 | -51.32 ± 10.62 | -10.63 | -1.72 | 0.74 | -39.71 |
| HIS-752 | -23.12 ± 6.57 | -5.34 | -10.60 | 10.92 | -18.11 |
| SER-684 | -13.37 ± 6.88 | -2.04 | 0.11 | 2.56 | -14.00 |
| ASN-755 | -12.05 ± 8.08 | -2.14 | -1.54 | 2.82 | -11.19 |
| LEU-857 | -8.62 ± 4.90 | -3.30 | -0.41 | 0.59 | -5.49 |
| LYS-691 | -2.87 ± 7.14 | -0.97 | 6.09 | 0.58 | -8.56 |
| ALA-751 | -1.90 ± 2.07 | -0.77 | 0.17 | 0.39 | -1.70 |
| ARG-568 | -1.16 ± 1.26 | -0.30 | -0.54 | -0.27 | -0.05 |
| THR-566 | -0.97 ± 0.31 | -0.45 | -0.72 | 0.20 | 0.00 |
| GLY-849 | -0.95 ± 1.25 | -0.85 | -0.40 | 0.75 | -0.44 |

## A3 — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)

*(ligand's own decomposition row, LIG-873: total -126.6 ± 115.5 — excluded from the ranking below, which is protein residues only)*

| residue | total | vdW | elec | polar | apolar |
| --- | --- | --- | --- | --- | --- |
| LYS-864 | -26.89 ± 51.42 | -4.51 | -5.41 | 2.26 | -19.23 |
| HIS-866 | -16.86 ± 51.70 | -3.92 | -1.93 | 4.30 | -15.31 |
| VAL-863 | -12.62 ± 23.55 | -2.02 | -0.47 | 0.39 | -10.52 |
| ASP-690 | -6.85 ± 14.88 | -1.00 | -5.29 | 7.38 | -7.95 |
| MET-867 | -4.36 ± 15.40 | -0.91 | 0.14 | 0.26 | -3.85 |
| LYS-691 | -2.99 ± 4.89 | -0.26 | 0.64 | -2.55 | -0.83 |
| ARG-590 | -1.58 ± 26.61 | -2.18 | 7.65 | 2.80 | -9.85 |
| LYS-662 | -1.52 ± 2.17 | -0.26 | -0.20 | -0.93 | -0.13 |
| GLU-850 | -1.06 ± 1.36 | -0.01 | -0.91 | -0.16 | 0.02 |
| PRO-543 | -1.04 ± 6.87 | -0.00 | 0.00 | -0.00 | -1.04 |

## A6 — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)

*(ligand's own decomposition row, LIG-873: total -269.8 ± 23.8 — excluded from the ranking below, which is protein residues only)*

| residue | total | vdW | elec | polar | apolar |
| --- | --- | --- | --- | --- | --- |
| LEU-853 | -36.25 ± 12.10 | -5.86 | -0.12 | 0.67 | -30.93 |
| VAL-683 | -29.77 ± 8.03 | -5.93 | -1.31 | 4.11 | -26.63 |
| HIS-752 | -26.14 ± 15.87 | -4.25 | 1.02 | 3.47 | -26.38 |
| SER-684 | -26.09 ± 6.72 | -3.67 | -0.41 | 4.97 | -26.98 |
| LEU-857 | -23.23 ± 7.33 | -4.38 | -0.92 | 1.22 | -19.15 |
| LYS-691 | -3.73 ± 1.36 | -0.23 | -3.37 | -0.05 | -0.08 |
| LYS-692 | -3.03 ± 1.31 | -0.10 | -6.49 | 3.58 | -0.02 |
| ARG-702 | -2.16 ± 0.55 | -0.03 | -3.47 | 1.35 | -0.01 |
| ARG-646 | -1.22 ± 0.34 | -0.01 | -2.00 | 0.78 | 0.01 |
| LEU-562 | -1.13 ± 2.43 | -1.76 | -0.56 | 0.09 | 1.10 |

## substrate-only (E+S) — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)

*(ligand's own decomposition row, COA-873: total -347.2 ± 33.7 — excluded from the ranking below, which is protein residues only)*

| residue | total | vdW | elec | polar | apolar |
| --- | --- | --- | --- | --- | --- |
| LYS-691 | -141.15 ± 12.80 | 4.87 | -185.57 | 74.76 | -35.21 |
| LYS-864 | -68.57 ± 21.99 | 5.08 | -153.44 | 86.58 | -6.78 |
| GLU-559 | -61.82 ± 15.33 | 6.16 | -165.95 | 137.47 | -39.50 |
| MET-655 | -38.09 ± 11.03 | -16.06 | -4.64 | 15.00 | -32.39 |
| LYS-692 | -32.58 ± 5.05 | -0.16 | -43.91 | 11.54 | -0.06 |
| GLN-766 | -26.15 ± 6.56 | -8.05 | -16.72 | 11.59 | -12.96 |
| GLY-656 | -22.80 ± 6.76 | -8.11 | 0.25 | 1.56 | -16.50 |
| LYS-735 | -21.73 ± 2.38 | -0.03 | -29.13 | 7.46 | -0.03 |
| ARG-702 | -21.36 ± 2.84 | -0.32 | -25.01 | 4.15 | -0.18 |
| HIS-866 | -18.26 ± 8.13 | -3.35 | -34.83 | 12.63 | 7.29 |

## CoA in ternary (E+S+I, S channel) — top stabilising PROTEIN residues (200 ns, total ΔG, kJ/mol)

*(ligand's own decomposition row, COA-874: total 51.2 ± 63.5 — excluded from the ranking below, which is protein residues only)*

| residue | total | vdW | elec | polar | apolar |
| --- | --- | --- | --- | --- | --- |
| GLY-524 | -37.23 ± 14.31 | -8.11 | -1.49 | 5.76 | -33.38 |
| MET-523 | -30.72 ± 18.11 | -8.93 | -3.55 | 7.24 | -25.48 |
| GLU-559 | -29.25 ± 4.96 | -0.27 | -33.45 | 4.61 | -0.14 |
| VAL-863 | -19.86 ± 23.39 | -4.06 | -1.30 | 2.50 | -17.01 |
| ALA-525 | -17.81 ± 9.06 | -3.94 | -0.14 | -0.03 | -13.70 |
| HIS-866 | -17.34 ± 14.57 | -3.32 | -3.66 | 1.80 | -12.17 |
| GLU-528 | -16.82 ± 5.80 | -0.38 | -25.74 | 10.08 | -0.78 |
| HIS-861 | -15.62 ± 18.44 | -1.69 | -1.66 | 4.39 | -16.67 |
| ASP-690 | -15.60 ± 3.20 | -0.02 | -15.10 | -0.21 | -0.27 |
| ASP-767 | -15.17 ± 2.02 | -0.03 | -17.39 | 2.30 | -0.04 |

## A1<->CoA pairwise per-ligand share

`A1-COAresidues_energy_summary.csv` decomposes the pairwise ΔG additively between the two ligand "residues" themselves:

| entry | total ΔG (kJ/mol) |
| --- | --- |
| LIG-873 | -15.88 ± 8.00 |
| COA-874 | -11.53 ± 21.10 |
