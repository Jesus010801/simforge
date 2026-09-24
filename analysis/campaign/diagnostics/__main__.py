"""Developer entry point (read-only):

    python -m analysis.campaign.diagnostics TRAJ [-s md.tpr] [-n semantic.ndx]
        [--group receptor=Receptor ...] [--auto-groups] [--cache-dir D] [-o OUT_DIR]
"""
from __future__ import annotations

import argparse
import json
import sys

from analysis.campaign.diagnostics import (
    DiagnosticContext, groups_from_semantic_index, run_diagnostics, write_report,
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Read-only trajectory diagnostics.")
    ap.add_argument("trajectory")
    ap.add_argument("-s", "--structure", default=None)
    ap.add_argument("-n", "--index", default=None)
    ap.add_argument("--group", action="append", default=[], help="role=GroupName")
    ap.add_argument("--auto-groups", action="store_true",
                    help="map standard semantic group names found in the index")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("-o", "--out-dir", default=None)
    ns = ap.parse_args(argv)
    groups = groups_from_semantic_index(ns.index) if (ns.auto_groups and ns.index) else {}
    for g in ns.group:
        role, _, name = g.partition("=")
        groups[role] = name
    ctx = DiagnosticContext(trajectory_path=ns.trajectory, structure_path=ns.structure,
                            index_path=ns.index, groups=groups, gmx=ns.gmx,
                            cache_dir=ns.cache_dir)
    report = run_diagnostics(ctx)
    if ns.out_dir:
        print(write_report(report, ns.out_dir), file=sys.stderr)
    print(json.dumps(report.to_dict(), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
