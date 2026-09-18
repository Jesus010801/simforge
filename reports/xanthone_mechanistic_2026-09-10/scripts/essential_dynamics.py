#!/usr/bin/env python3
"""Phase 3 — per-system essential dynamics: Cα extraction, PCA, FEL, DCCM.

PCA  : gmx covar + gmx anaeig on C-alpha  -> eigenvalues, PC1/PC2 projection,
       variance explained.
FEL  : 2D histogram of (PC1,PC2) -> G = -kT ln P, G -= min  (kT at 309.65 K).
DCCM : normalised Cα-Cα cross-correlation matrix (MDAnalysis).

All read the managed mdfit.xtc (200 ps).  Writes raw grids/matrices (npy/xvg) +
PDF + provenance.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, gmx, mdfit, tpr, ndx, read_xvg, provenance, sha, run, GMX

KT = 0.0083144626 * 309.65   # kJ/mol


def ed_dir(system):
    d = OUT / system / "essential_dynamics"
    d.mkdir(parents=True, exist_ok=True)
    return d


def extract_ca(system):
    d = ed_dir(system)
    ca_xtc, ca_gro = d / "ca.xtc", d / "ca.gro"
    if not ca_xtc.is_file():
        gmx(["trjconv", "-s", tpr(system), "-f", mdfit(system), "-n", ndx(system),
             "-o", ca_xtc], stdin="C-alpha\n", log=d / "ca_xtc.log")
    if not ca_gro.is_file():
        gmx(["trjconv", "-s", tpr(system), "-f", mdfit(system), "-n", ndx(system),
             "-dump", "0", "-o", ca_gro], stdin="C-alpha\n", log=d / "ca_gro.log")
    return ca_xtc, ca_gro


def pca(system):
    d = ed_dir(system)
    ca_xtc, ca_gro = extract_ca(system)
    eig = d / "eigenval.xvg"; vec = d / "eigenvec.trr"; avg = d / "average.pdb"
    proj2d = d / "pc1_pc2.xvg"; proj = d / "proj.xvg"
    if not eig.is_file():
        gmx(["covar", "-s", ca_gro, "-f", ca_xtc, "-o", eig, "-v", vec, "-av", avg,
             "-l", d / "covar.log", "-xpma", d / "covar.xpm"],
            stdin="System\nSystem\n", log=d / "covar_run.log")
    if not proj2d.is_file():
        gmx(["anaeig", "-v", vec, "-f", ca_xtc, "-s", ca_gro, "-first", "1", "-last", "2",
             "-proj", proj, "-2d", proj2d],
            stdin="System\nSystem\n", log=d / "anaeig.log")
    ev = read_xvg(eig)[1][:, 0] if eig.is_file() else np.array([])
    ev = ev[ev > 0]
    var = (ev / ev.sum()).tolist() if ev.size else []
    # gmx anaeig -2d writes "PC1  PC2" per line (no time column): x=PC1, y[:,0]=PC2
    x, y = read_xvg(proj2d) if proj2d.is_file() else (np.array([]), np.zeros((0, 1)))
    pc1 = np.asarray(x, float)
    pc2 = y[:, 0] if (y.ndim == 2 and y.shape[1] >= 1) else np.array([])
    if pc2.size != pc1.size:
        m = min(pc1.size, pc2.size)
        pc1, pc2 = pc1[:m], pc2[:m]
    np.save(d / "pc12.npy", np.column_stack([pc1, pc2]) if pc1.size else np.zeros((0, 2)))
    prov = provenance(system, "pca",
                      argv=["gmx covar", "gmx anaeig", "-first 1 -last 2 (C-alpha)"],
                      stdin="System", outputs=[eig, proj2d, avg],
                      extra=dict(n_eigenvalues=int(ev.size),
                                 variance_explained_pc1_5=[round(v, 4) for v in var[:5]],
                                 cum_variance_pc1_2=round(sum(var[:2]), 4) if len(var) >= 2 else None,
                                 pc1_range=[float(pc1.min()), float(pc1.max())] if pc1.size else None,
                                 pc2_range=[float(pc2.min()), float(pc2.max())] if pc2.size else None,
                                 n_frames=int(pc1.size)))
    (d / "provenance.pca.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return pc1, pc2, var


def fel(system, pc1, pc2, bins=40):
    d = ed_dir(system)
    if pc1.size == 0:
        return
    H, xe, ye = np.histogram2d(pc1, pc2, bins=bins)
    P = H / H.sum()
    with np.errstate(divide="ignore"):
        G = -KT * np.log(P)
    G[~np.isfinite(G)] = np.nan
    G -= np.nanmin(G)
    np.savez(d / "fel.npz", G=G, xedges=xe, yedges=ye, kT=KT)
    # basins: strict local minima in a 5x5 window, populated, below 1.5 kT,
    # separated by >= 2 bins (merge nearby).
    gg = np.where(np.isfinite(G), G, 1e9)
    thr = 1.5 * KT
    raw = []
    R = 2
    for i in range(R, gg.shape[0] - R):
        for j in range(R, gg.shape[1] - R):
            if H[i, j] < 3:
                continue
            w = gg[i-R:i+R+1, j-R:j+R+1]
            if gg[i, j] == w.min() and gg[i, j] <= thr:
                raw.append((gg[i, j], (xe[i]+xe[i+1])/2, (ye[j]+ye[j+1])/2))
    raw.sort()
    mins = []
    for g, x, y in raw:
        if all((x - m["pc1"])**2 + (y - m["pc2"])**2 > (2*(xe[1]-xe[0]))**2 for m in mins):
            mins.append(dict(pc1=float(x), pc2=float(y), G=float(g)))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5.2, 4.4))
        cs = plt.contourf(xe[:-1], ye[:-1], np.ma.masked_invalid(G).T, levels=25, cmap="viridis_r")
        plt.colorbar(cs, label="ΔG (kJ/mol)")
        plt.xlabel("PC1 (nm)"); plt.ylabel("PC2 (nm)")
        plt.title(f"FEL {system} (per-system PCA)")
        plt.tight_layout(); plt.savefig(d / "fel.png", dpi=150); plt.close()
    except Exception as e:  # noqa
        (d / "fel_plot_error.txt").write_text(str(e))
    (d / "fel_basins.json").write_text(json.dumps(
        dict(system=system, kT_kjmol=KT, n_basins=len(mins), basins=sorted(mins, key=lambda m: m["G"])),
        indent=2))


def dccm(system):
    """Cα-Cα normalised cross-correlation via MDAnalysis (system python3)."""
    d = ed_dir(system)
    ca_xtc, ca_gro = extract_ca(system)
    script = d / "_dccm_run.py"
    script.write_text(f'''
import numpy as np, MDAnalysis as mda
from MDAnalysis.analysis import align
u = mda.Universe("{ca_gro}", "{ca_xtc}")
ref = mda.Universe("{ca_gro}")
align.AlignTraj(u, ref, select="name CA", in_memory=True).run()
ca = u.select_atoms("name CA")
coords = np.array([ca.positions.copy() for _ in u.trajectory])   # (T, N, 3)
mean = coords.mean(axis=0)
d = coords - mean
N = d.shape[1]
num = np.einsum("tik,tjk->ij", d, d)
var = np.einsum("tik,tik->i", d, d)
den = np.sqrt(np.outer(var, var))
C = num / den
np.save("{d}/dccm.npy", C)
print("dccm", C.shape, "range", float(C.min()), float(C.max()))
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(5,4.2)); plt.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
    plt.colorbar(label="Cα-Cα correlation"); plt.xlabel("residue"); plt.ylabel("residue")
    plt.title("DCCM {system}"); plt.tight_layout(); plt.savefig("{d}/dccm.png", dpi=150); plt.close()
except Exception as e:
    open("{d}/dccm_plot_error.txt","w").write(str(e))
''')
    r = run(["python3", str(script)], log=d / "dccm.log", timeout=7200)
    C = np.load(d / "dccm.npy") if (d / "dccm.npy").is_file() else None
    prov = provenance(system, "dccm", argv=["python3", "MDAnalysis Cα cross-correlation"],
                      stdin="", outputs=[d / "dccm.npy", d / "dccm.png"],
                      extra=dict(matrix_shape=list(C.shape) if C is not None else None,
                                 method="normalised covariance C_ij = <Δr_i·Δr_j>/sqrt(<Δr_i²><Δr_j²>) on aligned Cα",
                                 returncode=r.returncode))
    (d / "provenance.dccm.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")


def run_system(system):
    if not Path(mdfit(system)).is_file():
        print(f"  {system}: mdfit not ready, SKIP"); return
    print(f"[{system}] essential dynamics")
    pc1, pc2, var = pca(system)
    fel(system, pc1, pc2)
    dccm(system)
    print(f"  {system}: pca var(PC1,PC2)={[round(v,3) for v in var[:2]]}  frames={pc1.size}")


if __name__ == "__main__":
    for s in (sys.argv[1:] or list(SYS)):
        run_system(s)
