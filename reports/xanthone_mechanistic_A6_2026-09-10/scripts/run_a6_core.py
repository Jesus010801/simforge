#!/usr/bin/env python3
"""Run core_descriptors + essential_dynamics (PCA/FEL/DCCM) + clustering +
cluster_states for HMG-R-200ns-A6 on md_final.xtc (via a6_common patch)."""
import time
import a6_common  # noqa: F401  (applies the mdfit/provenance patch)
import core_descriptors, essential_dynamics, clustering, cluster_states

A6 = "HMG-R-200ns-A6"

t0 = time.time()
print(f"\n===== core_descriptors [{A6}] =====", flush=True)
core_descriptors.run_system(A6)
print(f"  elapsed {time.time()-t0:.0f}s", flush=True)

t1 = time.time()
print(f"\n===== essential_dynamics (PCA/FEL/DCCM) [{A6}] =====", flush=True)
essential_dynamics.run_system(A6)
print(f"  elapsed {time.time()-t1:.0f}s", flush=True)

t2 = time.time()
print(f"\n===== clustering (fixed 0.15 nm) [{A6}] =====", flush=True)
clustering.run_system(A6)
print(f"  elapsed {time.time()-t2:.0f}s", flush=True)

t3 = time.time()
print(f"\n===== cluster_states (adaptive macro-states) [{A6}] =====", flush=True)
cluster_states.run_system(A6)
print(f"  elapsed {time.time()-t3:.0f}s", flush=True)

print(f"\nRUN_A6_CORE COMPLETE  total {time.time()-t0:.0f}s", flush=True)
