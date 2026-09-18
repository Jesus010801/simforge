#!/usr/bin/env python3
"""Re-run every analysis for the given new systems on the PBC-corrected mdfit.xtc."""
import sys
NEW = ["A3-HMG-R", "HMG-R-200ns-A6", "system_A3_COA", "system_A6_COA"]


def main(systems):
    import core_descriptors, essential_dynamics, clustering, cluster_states
    for s in systems:
        print(f"\n########## {s} ##########", flush=True)
        core_descriptors.run_system(s)
        essential_dynamics.run_system(s)
        clustering.run_system(s)
        cluster_states.run_system(s)
    print("\nRUN_ALL_NEW COMPLETE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or NEW)
