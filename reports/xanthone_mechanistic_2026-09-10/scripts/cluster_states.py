#!/usr/bin/env python3
"""Adaptive 'state' clustering — a coarser cutoff giving an interpretable number
of macro-states (target 2-10) for representative structures + state populations,
in addition to the fixed-cutoff heterogeneity count from clustering.py.

Writes clustering/state_clustering.json (+ state_clusters.pdb, state_clid.xvg).
"""
from __future__ import annotations
import sys, json, re
from pathlib import Path
import numpy as np
from mech_common import SYS, OUT, gmx, mdfit, tpr, ndx, provenance

CUTOFFS = [0.15, 0.25, 0.35]


def n_from_log(txt):
    m = re.search(r"Found\s+(\d+)\s+clusters", txt)
    return int(m.group(1)) if m else None


def run_system(system):
    if not Path(mdfit(system)).is_file():
        return
    d = OUT / system / "clustering"
    d.mkdir(parents=True, exist_ok=True)
    T, N, fit = tpr(system), ndx(system), mdfit(system)
    chosen = None
    trace = []
    for cut in CUTOFFS:
        log = d / f"state_probe_{cut}.log"
        clid = d / f"state_probe_{cut}.xvg"
        clpdb = d / f"state_probe_{cut}.pdb"
        if not log.is_file():
            gmx(["cluster", "-s", T, "-f", fit, "-n", N, "-tu", "ns", "-method", "gromos",
                 "-cutoff", str(cut), "-clid", clid, "-g", log,
                 "-o", d / f"state_probe_{cut}.xpm"],
                stdin="Backbone\nBackbone\n", log=d / f"gmx_state_{cut}.log", timeout=3600)
        nc = n_from_log(log.read_text(errors="replace")) if log.is_file() else None
        trace.append(dict(cutoff=cut, n_clusters=nc))
        if nc is not None and 2 <= nc <= 10 and chosen is None:
            chosen = (cut, clid, clpdb, nc)
    if chosen is None:
        # fall back to the coarsest tried
        for cut in reversed(CUTOFFS):
            log = d / f"state_probe_{cut}.log"
            nc = n_from_log(log.read_text(errors="replace")) if log.is_file() else None
            if nc:
                chosen = (cut, d / f"state_probe_{cut}.xvg", d / f"state_probe_{cut}.pdb", nc)
                break
    if chosen is None:
        print(f"  {system}: no clustering converged"); return
    cut, clid, clpdb, nc = chosen
    # representative structures for the chosen cutoff only
    if not Path(clpdb).is_file():
        gmx(["cluster", "-s", T, "-f", fit, "-n", N, "-tu", "ns", "-method", "gromos",
             "-cutoff", str(cut), "-cl", clpdb, "-clid", d / f"state_final_{cut}.xvg",
             "-g", d / f"state_final_{cut}.log", "-o", d / f"state_final_{cut}.xpm"],
            stdin="Backbone\nBackbone\n", log=d / f"gmx_state_final.log", timeout=3600)
        clid = d / f"state_final_{cut}.xvg"
    ids = []
    for line in open(clid, errors="replace"):
        s = line.strip()
        if not s or s[0] in "@#&":
            continue
        c = s.split()
        try:
            ids.append((float(c[0]), int(float(c[1]))))
        except (ValueError, IndexError):
            pass
    ids = np.array(ids)
    total = len(ids)
    uid, cnt = np.unique(ids[:, 1].astype(int), return_counts=True)
    pops = {int(u): dict(frames=int(c), fraction=round(float(c) / total, 4)) for u, c in zip(uid, cnt)}
    windows = {}
    for w0, w1 in [(0, 50), (50, 100), (100, 150), (150, 200)]:
        m = (ids[:, 0] >= w0) & (ids[:, 0] < w1)
        seg = ids[m, 1].astype(int)
        if seg.size:
            u2, c2 = np.unique(seg, return_counts=True)
            windows[f"{w0}-{w1}ns"] = {int(u): round(float(x) / seg.size, 3) for u, x in zip(u2, c2)}
    out = dict(system=system, state_cutoff_nm=cut, n_states=nc,
              cutoff_scan=trace, populations_overall=pops, populations_by_window=windows,
              representative_structures_pdb=str(clpdb),
              note="adaptive coarse cutoff for interpretable macro-states; the fixed-0.15 nm "
                   "count in clustering_summary.json is the conformational-heterogeneity metric")
    (d / "state_clustering.json").write_text(json.dumps(out, indent=2))
    prov = provenance(system, "state_clustering",
                      argv=["gmx", "cluster", "gromos", f"-cutoff {cut}"],
                      stdin="Backbone\nBackbone\n", outputs=[clpdb, clid],
                      extra=dict(NEW_EXTENSION=True, n_states=nc, state_cutoff_nm=cut))
    (d / "provenance.state_clustering.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    top = sorted(pops.items(), key=lambda kv: -kv[1]["fraction"])[:4]
    print(f"  {system}: {nc} macro-states @ {cut} nm; "
          + ", ".join(f"S{k}={v['fraction']}" for k, v in top))


if __name__ == "__main__":
    for s in (sys.argv[1:] or list(SYS)):
        run_system(s)
