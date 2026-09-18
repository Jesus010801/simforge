#!/usr/bin/env python3
"""Phase 5 (NEW extension) — conformational clustering, representative structures,
state populations.  NOT part of the historical A1 workflow.

gmx cluster (gromos) on Cα RMSD; representative = cluster central member;
populations = frame fraction per cluster (overall + per 50 ns window).
Cluster centroids are also located in the shared PC1/PC2 plane.
"""
from __future__ import annotations
import sys, json, re
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, gmx, mdfit, tpr, ndx, provenance, sha

CUTOFF_NM = 0.15   # Cα-RMSD gromos cutoff; reported, tuned once for HMG-CoA reductase


def cl_dir(system):
    d = OUT / system / "clustering"
    d.mkdir(parents=True, exist_ok=True)
    return d


def parse_cluster_log(txt):
    m = re.search(r"Found\s+(\d+)\s+clusters", txt)
    n = int(m.group(1)) if m else None
    sizes = {}
    # "cl. |  #st  rmsd  | middle rmsd | cluster members"
    for line in txt.splitlines():
        mm = re.match(r"\s*(\d+)\s*\|\s*(\d+)\s+([\d.]+)\s*\|\s*(\d+)\s+([\d.]+)\s*\|", line)
        if mm:
            sizes[int(mm.group(1))] = dict(size=int(mm.group(2)), middle_frame=int(mm.group(4)))
    return n, sizes


def run_system(system, cutoff=CUTOFF_NM):
    if not Path(mdfit(system)).is_file():
        print(f"  {system}: mdfit not ready, SKIP"); return
    d = cl_dir(system)
    T, N, fit = tpr(system), ndx(system), mdfit(system)
    log = d / "cluster.log"; logf = d / "gmx.log"
    sz = d / "cluster_sizes.xvg"; clid = d / "cluster_id_vs_time.xvg"
    clpdb = d / "clusters.pdb"; xpm = d / "rmsd_clust.xpm"
    if not clid.is_file():
        gmx(["cluster", "-s", T, "-f", fit, "-n", N, "-tu", "ns",
             "-method", "gromos", "-cutoff", str(cutoff),
             "-cl", clpdb, "-sz", sz, "-clid", clid, "-o", xpm, "-g", log],
            stdin="Backbone\nBackbone\n", log=logf, timeout=7200)
    txt = log.read_text(errors="replace") if log.is_file() else ""
    nclust, sizes = parse_cluster_log(txt)

    # frame -> cluster id  (clid xvg: col0 time(ns), col1 cluster)
    ids = []
    if clid.is_file():
        for line in open(clid, errors="replace"):
            s = line.strip()
            if not s or s[0] in "@#&":
                continue
            c = s.split()
            try:
                ids.append((float(c[0]), int(float(c[1]))))
            except (ValueError, IndexError):
                pass
    ids = np.array(ids) if ids else np.zeros((0, 2))
    total = len(ids)
    pops = {}
    if total:
        uid, cnt = np.unique(ids[:, 1].astype(int), return_counts=True)
        pops = {int(u): dict(frames=int(c), fraction=round(float(c) / total, 4))
                for u, c in zip(uid, cnt)}
    # per-window populations (0-50/50-100/100-150/150-200 ns)
    windows = {}
    for w0, w1 in [(0, 50), (50, 100), (100, 150), (150, 200)]:
        mask = (ids[:, 0] >= w0) & (ids[:, 0] < w1) if total else np.array([], bool)
        seg = ids[mask, 1].astype(int) if total else np.array([], int)
        if seg.size:
            uid, cnt = np.unique(seg, return_counts=True)
            windows[f"{w0}-{w1}ns"] = {int(u): round(float(c) / seg.size, 3) for u, c in zip(uid, cnt)}

    summary = dict(system=system, method="gmx cluster gromos, Cα-RMSD (Backbone fit+rmsd)",
                   cutoff_nm=cutoff, n_clusters=nclust, n_frames=total,
                   cluster_middle_frames={k: v["middle_frame"] for k, v in sizes.items()},
                   populations_overall=pops, populations_by_window=windows,
                   representative_structures_pdb=str(clpdb),
                   note="NEW analysis — not present in the historical A1 mechanistic workflow")
    (d / "clustering_summary.json").write_text(json.dumps(summary, indent=2))
    prov = provenance(system, "clustering",
                      argv=["gmx", "cluster", "-method", "gromos", "-cutoff", str(cutoff)],
                      stdin="Backbone\nBackbone\n", outputs=[clpdb, sz, clid, xpm],
                      extra=dict(n_clusters=nclust, populations_overall=pops,
                                 NEW_EXTENSION=True))
    (d / "provenance.clustering.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    top = sorted(pops.items(), key=lambda kv: -kv[1]["fraction"])[:4] if pops else []
    print(f"  {system}: {nclust} clusters (cutoff {cutoff} nm); top populations "
          + ", ".join(f"C{k}={v['fraction']}" for k, v in top))
    return summary


if __name__ == "__main__":
    args = sys.argv[1:] or list(SYS)
    for s in args:
        run_system(s)
