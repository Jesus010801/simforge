#!/usr/bin/env python3
"""Phase 4 — shared PCA basis + shared FEL for A1 / A3 / A6.

The three HMG-CoA-reductase complexes have byte-identical protein topology
(1614 Cα, 24199 protein atoms, same residue numbering).  We:

  1. extract the Cα trajectory of each (already done in essential_dynamics);
  2. align every frame of every system to ONE shared reference (iteratively
     refined mean Cα structure of the pooled ensemble);
  3. pool all aligned Cα coordinate vectors and run PCA (sklearn) -> the SHARED
     eigenbasis + variance explained;
  4. project each system separately onto the shared PC1/PC2;
  5. build per-system FELs on COMMON bin edges (shared axes/limits) + an overlay;
  6. report overlap / separation between A1, A3, A6 in the shared PC1-PC2 plane.

Run with system python3 (MDAnalysis + sklearn).  Writes to
analysis_outputs/xanthone_mechanistic/shared_pca/.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")
SHARED = OUT / "shared_pca"
SHARED.mkdir(parents=True, exist_ok=True)
# HMG-R-200ns-A6 excluded: no PBC-clean 200 ns trajectory on disk (see
# analysis_outputs/xanthone_mechanistic/HMG-R-200ns-A6/REVIEW_REQUIRED.md).
# Shared basis is A1 vs A3 (both 200 ns, clean); A6 can be projected later.
TRIO = {"A1": "A1", "A3-HMG-R": "A3"}
KT = 0.0083144626 * 309.65
GMX = "/usr/local/gromacs/bin/gmx"


def load_ca(system):
    import MDAnalysis as mda
    d = OUT / system / "essential_dynamics"
    u = mda.Universe(str(d / "ca.gro"), str(d / "ca.xtc"))
    ca = u.select_atoms("name CA")
    X = np.array([ca.positions.copy() for _ in u.trajectory])  # (T,N,3)
    t = np.array([ts.time for ts in u.trajectory]) if False else \
        np.arange(X.shape[0]) * 200.0
    return X, t


def kabsch(P, Q):
    """rotate P onto Q (both centred), return rotated P."""
    H = P.T @ Q
    V, S, Wt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1, 1, dsign])
    R = V @ D @ Wt
    return P @ R


def align_to(X, ref):
    ref_c = ref - ref.mean(0)
    out = np.empty_like(X)
    for i in range(X.shape[0]):
        c = X[i] - X[i].mean(0)
        out[i] = kabsch(c, ref_c)
    return out


def main():
    raw = {}
    for s in TRIO:
        d = OUT / s / "essential_dynamics"
        if not (d / "ca.xtc").is_file():
            print(f"MISSING ca.xtc for {s}; run essential_dynamics first"); return
        X, t = load_ca(s)
        raw[s] = (X, t)
        print(f"{s}: {X.shape[0]} frames, {X.shape[1]} Cα")

    N = next(iter(raw.values()))[0].shape[1]
    assert all(v[0].shape[1] == N for v in raw.values()), "Cα count mismatch"

    # iterative reference: start from A1 mean, refine on pooled aligned ensemble
    ref = raw["A1"][0].mean(0)
    for it in range(3):
        pooled = []
        for s in TRIO:
            pooled.append(align_to(raw[s][0], ref))
        pooled_all = np.concatenate(pooled, axis=0)
        newref = pooled_all.mean(0)
        shift = np.sqrt(((newref - ref) ** 2).sum(1).mean())
        ref = newref
        print(f"  ref iter {it}: mean Cα shift {shift:.4f} Å")
    aligned = {s: align_to(raw[s][0], ref) for s in TRIO}

    Xall = np.concatenate([aligned[s].reshape(aligned[s].shape[0], -1) for s in TRIO], axis=0)
    Xall -= Xall.mean(0)
    from sklearn.decomposition import PCA
    p = PCA(n_components=10)
    p.fit(Xall)
    var = p.explained_variance_ratio_

    proj = {}
    for s in TRIO:
        M = aligned[s].reshape(aligned[s].shape[0], -1) - Xall.mean(0) * 0  # already centred by fit mean below
        Mc = aligned[s].reshape(aligned[s].shape[0], -1) - (np.concatenate(
            [aligned[k].reshape(aligned[k].shape[0], -1) for k in TRIO], axis=0)).mean(0)
        pc = p.transform(Mc)[:, :2] / 10.0   # Å -> nm for readability
        proj[s] = pc
        np.save(SHARED / f"proj_{TRIO[s]}.npy", pc)

    # shared FEL: common bin edges from pooled projections
    allpc = np.concatenate([proj[s] for s in TRIO], axis=0)
    x0, x1 = np.percentile(allpc[:, 0], [0.5, 99.5])
    y0, y1 = np.percentile(allpc[:, 1], [0.5, 99.5])
    xe = np.linspace(x0, x1, 61); ye = np.linspace(y0, y1, 61)
    fels = {}
    basins = {}
    for s in TRIO:
        H, _, _ = np.histogram2d(proj[s][:, 0], proj[s][:, 1], bins=[xe, ye])
        P = H / H.sum()
        with np.errstate(divide="ignore"):
            G = -KT * np.log(P)
        G[~np.isfinite(G)] = np.nan
        G -= np.nanmin(G)
        fels[TRIO[s]] = G
        gg = np.where(np.isfinite(G), G, 1e9)
        bs = []
        for i in range(1, gg.shape[0] - 1):
            for j in range(1, gg.shape[1] - 1):
                if gg[i, j] == gg[i-1:i+2, j-1:j+2].min() and gg[i, j] < 2.5:
                    bs.append(dict(pc1=float((xe[i]+xe[i+1])/2), pc2=float((ye[j]+ye[j+1])/2),
                                   G=round(float(gg[i, j]), 3)))
        basins[TRIO[s]] = sorted(bs, key=lambda b: b["G"])
    np.savez(SHARED / "shared_fel.npz", xedges=xe, yedges=ye, kT=KT,
             **{f"G_{k}": v for k, v in fels.items()})

    # overlap: Bhattacharyya-ish on the 2D histograms; centroid separation
    def hist(pc):
        h, _, _ = np.histogram2d(pc[:, 0], pc[:, 1], bins=[xe, ye])
        return h / h.sum()
    hh = {TRIO[s]: hist(proj[s]) for s in TRIO}
    overlap = {}
    labs = list(hh)
    for i in range(len(labs)):
        for j in range(i + 1, len(labs)):
            bc = float(np.sum(np.sqrt(hh[labs[i]] * hh[labs[j]])))
            cen_i = proj[[k for k in TRIO if TRIO[k] == labs[i]][0]].mean(0)
            cen_j = proj[[k for k in TRIO if TRIO[k] == labs[j]][0]].mean(0)
            sep = float(np.linalg.norm(cen_i - cen_j))
            overlap[f"{labs[i]}-{labs[j]}"] = dict(bhattacharyya_overlap=round(bc, 3),
                                                   centroid_separation_nm=round(sep, 3))

    summary = dict(
        method="shared Cα PCA: iterative pooled-mean reference, Kabsch alignment, "
               "sklearn PCA on the concatenated A1/A3/A6 aligned Cα ensemble, each "
               "system projected onto the shared eigenvectors (Å->nm).",
        n_calpha=int(N), kT_kjmol=round(KT, 4),
        frames={TRIO[s]: int(raw[s][0].shape[0]) for s in TRIO},
        shared_variance_explained=[round(float(v), 4) for v in var],
        cum_variance_pc1_2=round(float(var[:2].sum()), 4),
        cum_variance_pc1_5=round(float(var[:5].sum()), 4),
        pc_ranges_nm={TRIO[s]: dict(pc1=[round(float(proj[s][:, 0].min()), 3), round(float(proj[s][:, 0].max()), 3)],
                                    pc2=[round(float(proj[s][:, 1].min()), 3), round(float(proj[s][:, 1].max()), 3)])
                      for s in TRIO},
        pc_centroids_nm={TRIO[s]: [round(float(proj[s][:, 0].mean()), 3), round(float(proj[s][:, 1].mean()), 3)]
                         for s in TRIO},
        pairwise_overlap=overlap,
        shared_fel_basins=basins,
        outputs=["proj_{A1,A3,A6}.npy", "shared_fel.npz", "shared_pca_summary.json",
                 "shared_pc12.png", "shared_fel.png"],
    )
    (SHARED / "shared_pca_summary.json").write_text(json.dumps(summary, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        col = {"A1": "#1f77b4", "A3": "#ff7f0e", "A6": "#2ca02c"}
        plt.figure(figsize=(5.6, 4.8))
        for s in TRIO:
            l = TRIO[s]
            plt.scatter(proj[s][:, 0], proj[s][:, 1], s=6, alpha=0.45, label=l, color=col[l])
            c = proj[s].mean(0)
            plt.scatter(*c, s=120, marker="X", color=col[l], edgecolor="k", zorder=5)
        plt.xlabel("shared PC1 (nm)"); plt.ylabel("shared PC2 (nm)")
        plt.title(f"Shared PCA — A1/A3/A6  (PC1 {var[0]*100:.0f}%, PC2 {var[1]*100:.0f}%)")
        plt.legend(); plt.tight_layout(); plt.savefig(SHARED / "shared_pc12.png", dpi=150); plt.close()

        labs2 = list(TRIO.values())
        fig, ax = plt.subplots(1, len(labs2), figsize=(4.3 * len(labs2), 4),
                               sharex=True, sharey=True, squeeze=False)
        ax = ax[0]
        for k, a in zip(labs2, ax):
            cs = a.contourf(xe[:-1], ye[:-1], np.ma.masked_invalid(fels[k]).T, levels=25, cmap="viridis_r",
                            vmin=0, vmax=np.nanmax([np.nanmax(fels[q]) for q in fels]))
            a.set_title(f"FEL {k} (shared PC basis)"); a.set_xlabel("PC1 (nm)")
        ax[0].set_ylabel("PC2 (nm)")
        fig.colorbar(cs, ax=ax, label="ΔG (kJ/mol)", shrink=0.85)
        plt.savefig(SHARED / "shared_fel.png", dpi=150); plt.close()
    except Exception as e:  # noqa
        (SHARED / "plot_error.txt").write_text(str(e))

    print("shared PCA:", summary["shared_variance_explained"][:3],
          "cum(PC1,2)=", summary["cum_variance_pc1_2"])
    print("overlap:", overlap)


if __name__ == "__main__":
    main()
