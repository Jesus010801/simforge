#!/usr/bin/env python3
"""CORRECTIVE reprocessing — the pre-made mdfit/mdcenter for the 4 new HMG-R
systems were built WITHOUT whole-molecule PBC treatment for the tetramer:
protein RMSD reaches 5-7 nm and Rg doubles (3.3 -> 6.5 nm) purely from chains
imaged across the box.  Rebuild from the rawest available trajectory with
`-pbc mol -center -ur compact` (the treatment that gives clean references).

A3-HMG-R / system_A3_COA / system_A6_COA : from md.xtc.
HMG-R-200ns-A6 : no md.xtc -> from mdcenter.xtc with `-pbc nojump` then
                 `-pbc mol -center` (mdcenter is residue-wrapped, not mol-whole).
"""
from __future__ import annotations
import subprocess, time, json, hashlib
from pathlib import Path

GMX = "/usr/local/gromacs/bin/gmx"
OUT = Path("/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_mechanistic")
NUEV = "/home/jesusxd/Escritorio/Nuevos_sistemas"
STRIDE = 200

JOBS = {
 "A3-HMG-R":       dict(raw="md.xtc",       kind="raw"),
 "system_A3_COA":  dict(raw="md.xtc",       kind="raw"),
 "system_A6_COA":  dict(raw="md.xtc",       kind="raw"),
 "HMG-R-200ns-A6": dict(raw="mdcenter.xtc", kind="mdcenter"),
}


def sh(argv, stdin, log):
    t0 = time.time()
    p = subprocess.run(argv, input=stdin, text=True, capture_output=True, timeout=28800)
    Path(log).write_text(f"$ {' '.join(argv)}\nstdin {stdin!r}\nrc={p.returncode} {time.time()-t0:.0f}s\n\n"
                         f"STDERR:\n{p.stderr[-10000:]}\n")
    return p


def sha(p, lim=2 * 1024**3):
    p = Path(p)
    if not p.is_file() or p.stat().st_size > lim:
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def check(traj):
    import re
    p = subprocess.run([GMX, "check", "-f", str(traj)], capture_output=True, text=True, timeout=14400)
    t = p.stderr + p.stdout
    n = re.search(r"^(?:Coords|Step)\s+(\d+)", t, re.M)
    lf = re.search(r"Last frame\s+\d+\s+time\s+([\d.eE+]+)", t)
    at = re.search(r"# Atoms\s+(\d+)", t)
    return dict(n_frames=int(n.group(1)) if n else None,
               last_ps=float(lf.group(1)) if lf else None,
               atoms=int(at.group(1)) if at else None)


def run(name, cfg):
    src = Path(f"{NUEV}/{name}")
    tj = OUT / name / "traj"
    tj.mkdir(parents=True, exist_ok=True)
    logs = OUT / name / "prep_logs"
    logs.mkdir(exist_ok=True)
    tpr, ndx = src / "md.tpr", src / "index.ndx"
    final = tj / "mdfit.xtc"
    # archive the contaminated version
    if final.is_file():
        bad = tj / "mdfit_PBC_CONTAMINATED.xtc"
        if not bad.is_file():
            final.rename(bad)
    for stale in ("mdcenter_dt200.xtc", "nojump.xtc", "molcen.xtc", "ref.gro"):
        (tj / stale).unlink(missing_ok=True)

    # Tetramer in a tight rectangular box: `-pbc mol` alone leaves chains imaged
    # apart.  Recipe: -pbc nojump (keep the complex together as at frame 0)
    # -> -pbc mol -center -ur compact -> -fit rot+trans.
    nj = tj / "nojump_dt200.xtc"
    sh([GMX, "trjconv", "-s", str(tpr), "-f", str(src / cfg["raw"]), "-n", str(ndx),
        "-pbc", "nojump", "-dt", str(STRIDE), "-o", str(nj)],
       "System\n", logs / "R1_nojump.log")
    mc = tj / "molcen_dt200.xtc"
    sh([GMX, "trjconv", "-s", str(tpr), "-f", str(nj), "-n", str(ndx),
        "-pbc", "mol", "-center", "-ur", "compact", "-o", str(mc)],
       "Protein\nSystem\n", logs / "R2_molcenter.log")
    sh([GMX, "trjconv", "-s", str(tpr), "-f", str(mc), "-n", str(ndx),
        "-fit", "rot+trans", "-o", str(final)],
       "Backbone\nSystem\n", logs / "R3_fit.log")
    nj.unlink(missing_ok=True)

    ref = tj / "ref.gro"
    sh([GMX, "trjconv", "-s", str(tpr), "-f", str(final), "-dump", "0", "-o", str(ref)],
       "System\n", logs / "ref_gro.log")
    (tj / "molcen_dt200.xtc").unlink(missing_ok=True)

    info = check(final)
    # quick RMSD sanity
    rlog = logs / "sanity_rmsd.log"
    sh([GMX, "rms", "-s", str(tpr), "-f", str(final), "-n", str(ndx), "-tu", "ns",
        "-o", str(logs / "sanity_rmsd.xvg")], "Backbone\nBackbone\n", rlog)
    rv = []
    if (logs / "sanity_rmsd.xvg").is_file():
        for l in open(logs / "sanity_rmsd.xvg", errors="replace"):
            if l.strip() and l[0] not in "@#&":
                try:
                    rv.append(float(l.split()[1]))
                except (ValueError, IndexError):
                    pass
    import numpy as np
    rmsd_max = float(np.max(rv)) if rv else None
    prov = dict(system=name, corrective=True, source_dir=str(src),
                reason="pre-made mdfit/mdcenter lacked -pbc mol (tetramer chains imaged apart)",
                rebuilt_from=cfg["raw"], kind=cfg["kind"], stride_ps=STRIDE,
                analysis_trajectory=str(final), sha256=sha(final), check=info,
                sanity_protein_rmsd_max_nm=round(rmsd_max, 3) if rmsd_max else None,
                pbc_fix_ok=(rmsd_max is not None and rmsd_max < 0.6))
    (OUT / name / "prep_corrective.json").write_text(json.dumps(prov, indent=2))
    print(f"{name}: frames={info['n_frames']} last={info['last_ps']} atoms={info['atoms']} "
          f"RMSD_max={prov['sanity_protein_rmsd_max_nm']} nm  "
          f"{'OK' if prov['pbc_fix_ok'] else '*** STILL BAD ***'}")
    return prov


if __name__ == "__main__":
    import sys
    for n in (sys.argv[1:] or list(JOBS)):
        run(n, JOBS[n])
