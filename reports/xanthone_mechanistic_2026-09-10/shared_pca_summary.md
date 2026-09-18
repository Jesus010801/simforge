# Shared PCA basis — A1 / A3 / A6

shared Cα PCA: iterative pooled-mean reference, Kabsch alignment, sklearn PCA on the concatenated A1/A3/A6 aligned Cα ensemble, each system projected onto the shared eigenvectors (Å->nm). A6 subsampled 1-in-2 (2001->~1001 frames) for basis-fit/projection balance only; A6's own per-system PCA/FEL/DCCM use the full native series.

- Cα per system: **1614** (identical topology / residue numbering)
- frames: {'A1': 1001, 'A3': 1001, 'A6': 1001}
- kT = 2.5746 kJ/mol (309.65 K)

## Variance explained (shared eigenbasis)

| PC | fraction | cumulative |
| --- | --- | --- |
| PC1 | 0.279 | 0.279 |
| PC2 | 0.133 | 0.412 |
| PC3 | 0.070 | 0.482 |
| PC4 | 0.042 | 0.524 |
| PC5 | 0.038 | 0.563 |
| PC6 | 0.029 | 0.591 |

PC1+PC2 = **41.2%**, PC1–PC5 = 56.3%

## PC1/PC2 projections (nm)

| system | PC1 centroid | PC2 centroid | PC1 range | PC2 range |
| --- | --- | --- | --- | --- |
| A1 | -3.3 | -1.263 | [-5.193, -0.336] | [-2.002, -0.117] |
| A3 | 3.291 | -1.379 | [-2.48, 5.795] | [-2.326, -0.36] |
| A6 | 0.008 | 2.642 | [-2.241, 1.329] | [-0.742, 4.683] |

## Overlap / separation (shared plane)

| pair | Bhattacharyya overlap | centroid separation (nm) |
| --- | --- | --- |
| A1-A3 | 0.009 | 6.592 |
| A1-A6 | 0.009 | 5.118 |
| A3-A6 | 0.001 | 5.192 |

## Dominant conformational regions (shared-basis FEL basins, ΔG < 2.5 kJ/mol)

- **A1**: 3 basin(s) — (PC1 -3.6195177365938824, PC2 -1.3017871016263964) ΔG 0.0; (PC1 -3.4495595955848697, PC2 -1.515349794427554) ΔG 0.0; (PC1 -3.279601454575857, PC2 -1.1950057552258175) ΔG 0.621
- **A3**: 24 basin(s) — (PC1 2.1590590577125544, PC2 -1.3017871016263964) ΔG 0.0; (PC1 4.028598608811696, PC2 -0.8746617160240812) ΔG 0.43; (PC1 4.368514890829721, PC2 -1.7289124872287116) ΔG 0.947; (PC1 4.538473031838734, PC2 -1.8356938336292905) ΔG 0.947
- **A6**: 9 basin(s) — (PC1 0.11956136560440056, PC2 3.8237175256013862) ΔG 0.0; (PC1 0.45947764762242604, PC2 2.435560022393862) ΔG 1.121; (PC1 0.45947764762242604, PC2 2.5423413687944407) ΔG 1.121; (PC1 0.6294357886314388, PC2 3.7169361792008075) ΔG 1.121

Figures: `analysis_outputs/xanthone_mechanistic/shared_pca/shared_pc12.png`, `shared_fel.png`.
