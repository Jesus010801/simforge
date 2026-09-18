#!/usr/bin/env python3
"""Phase 1 — build the managed canonical analysis trajectory for every system.

Output: analysis_outputs/xanthone_mechanistic/<sys>/traj/mdfit.xtc
        (+ ref.gro, index.ndx for APO, prep.json provenance)

* references (A1/APO/COA/COMP): md.xtc -> trjconv -pbc mol -center -ur compact
  -> trjconv -fit rot+trans   (Backbone lsq fit), resampled to 200 ps.
* A3-HMG-R / system_A?_COA: source mdfit.xtc resampled to 200 ps.
* HMG-R-200ns-A6: source mdfit.xtc is 1000 ps (201 frames) -> regenerate from
  mdcenter.xtc (PBC-corrected) with -fit rot+trans at 200 ps.

The source tree is READ ONLY. All writes go to the managed tree.
"""
from __future__ import annotations
import json, subprocess, sys, time, hashlib, os
from pathlib import Path

GMX = "/usr/local/gromacs/bin/gmx"
ROOT = Path("/home/jesusxd/Escritorio/simforge")
OUT = ROOT / "analysis_outputs/xanthone_mechanistic"
STRIDE = 200
MECH = "/home/jesusxd/Escritorio/Mecanismo_inhibitorio"
NUEV = "/home/jesusxd/Escritorio/Nuevos_sistemas"

PLAN = {
 "A1":            dict(src=f"{MECH}/A1",   mode="reprocess",  raw="md.xtc"),
 "APO":           dict(src=f"{MECH}/APO",  mode="reprocess",  raw="md.xtc", make_ndx=True),
 "COA":           dict(src=f"{MECH}/COA",  mode="reprocess",  raw="md.xtc"),
 "COMP":          dict(src=f"{MECH}/COMP", mode="reprocess",  raw="md.xtc"),
 "A3-HMG-R":      dict(src=f"{NUEV}/A3-HMG-R",       mode="resample", fit="mdfit.xtc"),
 "HMG-R-200ns-A6":dict(src=f"{NUEV}/HMG-R-200ns-A6", mode="refit",    fit="mdcenter.xtc"),
 "system_A3_COA": dict(src=f"{NUEV}/system_A3_COA",  mode="resample", fit="mdfit.xtc"),
 "system_A6_COA": dict(src=f"{NUEV}/system_A6_COA",  mode="resample", fit="mdfit.xtc"),
}


def sh(argv, stdin=None, log=None, timeout=21600):
    t0 = time.time()
    p = subprocess.run(argv, input=stdin, text=True, capture_output=True, timeout=timeout)
    if log:
        log.write_text(f"$ {' '.join(argv)}\nstdin: {stdin!r}\nrc={p.returncode} dt={time.time()-t0:.0f}s\n\n"
                       f"STDERR:\n{p.stderr[-8000:]}\n\nSTDOUT:\n{p.stdout[-4000:]}\n")
    return p


def sha(path, limit=2 * 1024**3):
    if not Path(path).is_file() or Path(path).stat().st_size > limit:
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def check(traj):
    p = subprocess.run([GMX, "check", "-f", str(traj)], capture_output=True, text=True, timeout=14400)
    txt = p.stderr + p.stdout
    import re
    n = re.search(r"^(?:Coords|Step)\s+(\d+)\s+([\d.]+)?", txt, re.M)
    lf = re.search(r"Last frame\s+\d+\s+time\s+([\d.eE+]+)", txt)
    at = re.search(r"# Atoms\s+(\d+)", txt)
    return dict(n_frames=int(n.group(1)) if n else None,
               dt_ps=float(n.group(2)) if n and n.group(2) else None,
               last_ps=float(lf.group(1)) if lf else None,
               atoms=int(at.group(1)) if at else None)


def prep(name, cfg):
    src = Path(cfg["src"])
    d = OUT / name
    tj = d / "traj"
    tj.mkdir(parents=True, exist_ok=True)
    logs = d / "prep_logs"
    logs.mkdir(exist_ok=True)
    tpr = src / "md.tpr"
    ndx = src / "index.ndx"
    steps = []

    if cfg.get("make_ndx"):
        ndx = d / "index.ndx"
        if not ndx.is_file():
            r = sh([GMX, "make_ndx", "-f", str(tpr), "-o", str(ndx)], stdin="q\n",
                   log=logs / "make_ndx.log")
            steps.append(("make_ndx", r.returncode))

    final = tj / "mdfit.xtc"
    if final.is_file():
        print(f"  {name}: mdfit.xtc already present, skipping build")
    elif cfg["mode"] == "reprocess":
        raw = src / cfg["raw"]
        cen = tj / "mdcenter_dt200.xtc"
        r1 = sh([GMX, "trjconv", "-s", str(tpr), "-f", str(raw), "-n", str(ndx),
                 "-pbc", "mol", "-center", "-ur", "compact", "-dt", str(STRIDE),
                 "-o", str(cen)], stdin="Protein\nSystem\n", log=logs / "1_center.log")
        steps.append(("center", r1.returncode))
        r2 = sh([GMX, "trjconv", "-s", str(tpr), "-f", str(cen), "-n", str(ndx),
                 "-fit", "rot+trans", "-o", str(final)],
                stdin="Backbone\nSystem\n", log=logs / "2_fit.log")
        steps.append(("fit", r2.returncode))
    elif cfg["mode"] == "resample":
        r = sh([GMX, "trjconv", "-s", str(tpr), "-f", str(src / cfg["fit"]), "-n", str(ndx),
                "-dt", str(STRIDE), "-o", str(final)], stdin="System\n", log=logs / "resample.log")
        steps.append(("resample", r.returncode))
    elif cfg["mode"] == "refit":
        r = sh([GMX, "trjconv", "-s", str(tpr), "-f", str(src / cfg["fit"]), "-n", str(ndx),
                "-fit", "rot+trans", "-dt", str(STRIDE), "-o", str(final)],
               stdin="Backbone\nSystem\n", log=logs / "refit.log")
        steps.append(("refit", r.returncode))

    # reference structure (frame 0) matching mdfit
    ref = tj / "ref.gro"
    if final.is_file() and not ref.is_file():
        sh([GMX, "trjconv", "-s", str(tpr), "-f", str(final), "-dump", "0", "-o", str(ref)],
           stdin="System\n", log=logs / "ref_gro.log")

    info = check(final) if final.is_file() else {}
    prov = dict(system=name, mode=cfg["mode"], source_dir=str(src),
                source_read_only=True, managed_tpr=str(tpr), managed_index=str(ndx),
                analysis_trajectory=str(final), analysis_trajectory_sha256=sha(final),
                stride_ps=STRIDE, steps=steps, check=info,
                raw_production_trajectory=str(src / "md.xtc") if (src / "md.xtc").is_file() else None,
                preprocessing_chain=(
                    ["md.xtc", "trjconv -pbc mol -center -ur compact -dt 200 -> mdcenter_dt200.xtc",
                     "trjconv -fit rot+trans (Backbone) -> mdfit.xtc"] if cfg["mode"] == "reprocess"
                    else ["source mdfit.xtc", f"trjconv -dt {STRIDE} -> mdfit.xtc"] if cfg["mode"] == "resample"
                    else ["mdcenter.xtc (PBC-corrected)", f"trjconv -fit rot+trans -dt {STRIDE} -> mdfit.xtc"]))
    (d / "prep.json").write_text(json.dumps(prov, indent=2))
    print(f"  {name}: frames={info.get('n_frames')} dt={info.get('dt_ps')} "
          f"last={info.get('last_ps')} atoms={info.get('atoms')} sha={str(sha(final))[:12]}")
    return prov


if __name__ == "__main__":
    only = sys.argv[1:] or list(PLAN)
    for name in only:
        print(f"[{name}] {PLAN[name]['mode']}")
        prep(name, PLAN[name])
    print("done")
