# Phase 10 — Root-Cause Architectural Audit: Membrane Packing Failure

**Date:** 2026-06-25  
**Scope:** Why SimForge membrane embedding fails to achieve correct protein–lipid packing despite seven layers of corrections.  
**Status:** Analysis only. No fixes implemented.

---

## Executive Summary

The current pipeline detects the problem perfectly and fails to fix it for a fundamental reason: **every actuator it possesses either moves lipids vertically, removes lipids, or globally rescales coordinates**. Not one of them can add lipids to, or laterally redistribute lipids toward, an under-covered TM surface. The closed-loop optimizer (Phase 9B) searches a parameter space that is structurally incapable of containing a valid solution to annular gap under-coverage. The sensor (Phase 9A) is correct. The builder is searching the wrong manifold.

---

## 1. Sensors vs. Actuators

### Current sensors

| Sensor | What it measures | Works? |
|--------|-----------------|--------|
| TM burial score | Fraction of TM residues whose CA is between the headgroup planes | Yes |
| Void fraction (local) | Fraction of the TM footprint volume that contains no lipid atoms | Yes |
| Trapped lipid diagnostic | Lipid residues whose centroid lies inside the protein cavity | Yes |
| Surface interference gate | Lipid atoms overlapping TM surface within a cutoff | Yes |
| Protein–membrane interface coverage (Phase 9A) | `fraction_covered`, `fraction_exposed_gap`, `p90_nearest_lipid_distance` for TM-surface atoms | Yes — detects the exact problem |

**All sensors work.** The GLP-1R run confirms: `fraction_covered ≈ 0.41`, `fraction_exposed_gap ≈ 0.41`, `p90 ≈ 1.21 nm`. The Phase 9A evaluator has a correct read on what the optimizer must fix.

---

### Current actuators

| Actuator | Physical effect | Can it increase TM-surface lipid coverage? |
|----------|----------------|-------------------------------------------|
| **Bilayer Z-shift** (optimizer) | Translates all lipid atoms uniformly along Z | Only if the gap is caused by Z-misalignment. Irrelevant for a lateral (XY) annular gap. |
| **TM-exclusion mask padding** (optimizer + tm_aware backend) | Deletes lipid residues whose any atom falls within `base_cutoff + pad` of any TM atom | Makes coverage **worse**, not better. Widens the annular gap. |
| **Lipid exclusion cutoff change** | Same as mask padding — controls deletion threshold | Same problem. |
| **Energy minimization** (steepest descent, per-iteration in shrink loop) | Relaxes van der Waals clashes locally; displaces atoms by ~Å | Cannot move a lipid headgroup across a 1.2 nm gap. Lipids already far from the surface have no gradient pulling them inward. |
| **InflateGRO shrink loop** | Globally scales all lipid XY coordinates by factor 0.95 per iteration until APL converges | No. Explained in detail in §5. |
| **TM-aware backend** (tm_aware) | Classifies and deletes lipids inside TM cavity, then runs GROMACS topology rebuild | Reduces invasion. Does not fill the annular gap. |

**Conclusion:** Not one actuator can increase lateral lipid coverage around an under-covered TM surface. The correction mechanism set has no element that performs the required physical operation.

---

## 2. The Optimizer Searches a Space That Cannot Contain the Solution

The Phase 9B optimizer generates candidates by varying:

- `z_shift_offset_nm ∈ {-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6}` — vertical translation of lipids
- `mask_padding_nm ∈ {0.0, 0.1, 0.2, 0.3}` — additional radial exclusion around TM

The annular gap in a real GPCR embedding is a **lateral (XY) packing problem**. The TM bundle occupies a finite XY footprint (~2.5 × 2.5 nm for GLP-1R). After the bilayer is merged and the exclusion filter removes overlapping lipids, a ring-shaped void forms in the XY plane around the bundle. Lipids on the outer edge of this ring are ~0.9–1.2 nm away laterally.

Z-shift changes which Z-slice of the protein is at the membrane midplane. This is a useful correction for misaligned TM center vs. hydrophobic slab — but it does not change the lateral distribution of lipids. A lipid at XY position (3.0, 5.0) stays at (3.0, 5.0) regardless of what Z-shift is applied to the bilayer.

Mask padding systematically removes lipids that are too close to the TM surface. The optimizer's penalty term for this (`lip_del_penalty`) tries to discourage excessive removal, but the fundamental effect of any `mask_padding > 0` candidate is to **widen the annular gap**. The optimizer correctly discounts these candidates — but it has no candidate that can narrow the gap.

**The optimizer is topologically correct and computationally sound. It is searching a feasible region that does not intersect the set of physically valid solutions.**

---

## 3. Deletion Can Fix Invasion but Cannot Fix Under-Coverage

It is important to distinguish four failure modes, because they require opposite corrective actions:

| Failure mode | What it looks like | Correct actuator |
|---|---|---|
| **Lipid invasion / trapping** | Lipid residues whose centroid or most atoms are inside the TM helix bundle cavity (z-channel, lumen) | **Deletion** |
| **Lipid overcrowding** | Too many lipids per nm² near the protein; APL well below target; headgroups forced into the protein | **Deletion + inflation** |
| **Lipid under-coverage (annular gap)** | TM surface atoms have no lipid neighbors within 0.7 nm; large contiguous void patch | **Lateral placement / diffusion** |
| **Lateral packing failure** | The overall bilayer APL is correct but the lipid distribution is non-uniform; some patches are dense, some sparse | **Lateral redistribution** |

The trapped-lipid detector and the TM-exclusion mask are both deletion instruments. They are designed for modes 1 and 2. When applied to a system already in mode 3 or 4, they actively worsen coverage: every lipid removed from near the TM is one less lipid contributing to `fraction_covered`.

The GLP-1R result (`fraction_covered ≈ 0.41`, no trapped lipids) is unambiguously mode 3. All deletion-based actuators are contra-indicated.

---

## 4. Energy Minimization Cannot Close a Large Annular Gap

Steepest-descent and L-BFGS minimization follow the gradient of the potential energy landscape. The restoring force experienced by a lipid far from the protein surface is:
- zero from electrostatics (lipid tails are non-polar; the TM surface contributes minimal charge)
- attractive van der Waals from protein atoms — but this force decays as ~r⁻⁶ and is effectively zero beyond 0.5 nm

A lipid headgroup sitting 1.2 nm from the nearest TM atom (as reported in the GLP-1R run) is beyond the force cutoff for practical energy minimization. There is no gradient to follow. The lipid sits in a local minimum on the flat part of the potential surface.

Even if there were a gradient, steepest descent moves atoms by ~0.01–0.1 Å per step. Closing a 1.2 nm gap would require 100–1000 minimization steps per lipid, but the force decays toward zero as the lipid approaches the surface — the gradient never builds to drive the lipid across the gap.

The per-iteration minimization in the shrink loop is there to relax VDW clashes after each deflation step, not to diffuse lipids. It does exactly what it is designed for and nothing more.

---

## 5. InflateGRO Shrink Cannot Guarantee TM Surface Coverage

### How inflategro actually works

InflateGRO.pl inflates the XY coordinates of all lipid atoms by multiplying their X and Y components by the inflation factor (default 4.0) relative to some expansion origin (typically the system XY centroid). It then iteratively deflates (multiplies by 0.95) and minimizes until the APL converges to target.

### What this means geometrically

When lipids are inflated by 4×, they spread radially outward from the box centroid. As they deflate, they converge radially **toward the box centroid** — not toward any molecular surface.

For a protein whose TM bundle is at or near the box XY centroid (after box matching), the deflation should in principle bring lipids toward the protein. In practice, two things prevent this:

1. **The exclusion cutoff (Å parameter passed to inflategro) removes overlapping lipids** at each inflation step. Any lipid that would have ended up filling the annular gap is deleted the moment it inflates into the cutoff zone. The gap is re-created at every iteration.

2. **The convergence criterion is APL, not surface coverage.** The loop terminates when the mean area per lipid across the full bilayer matches the target (e.g. ≤ 64 Å²). This global criterion is satisfied even if the lipid distribution is deeply inhomogeneous — a dense outer ring with a vacant annular zone can have the same mean APL as a uniform bilayer.

The shrink loop optimizes a scalar (mean APL) that does not capture the spatial distribution of lipids around the protein surface. It can produce a system with a perfect global APL and a 1.2 nm protein–lipid gap at every TM surface atom.

---

## 6. The Missing Operation

The pipeline has seven layers of measurement and correction, none of which performs the following operation:

> **Identify specific patches of TM surface with no nearby lipids, and place or move lipid molecules into those patches.**

This operation has three required components:

1. **Surface patch identification** — find contiguous regions of TM surface atoms with `min_nearest_lipid_distance > threshold` (already available from Phase 9A gap clusters).

2. **Lipid introduction into the gap region** — either by:
   a. Selecting lipid molecules from the existing bilayer that are in over-packed regions and translating them into the gap; or
   b. Inserting new lipid molecules copied from a template library; or
   c. Running MD with enough thermal energy to diffuse lipids into the gap.

3. **Clash resolution and topology consistency** — removing introduced lipids that overlap with protein atoms, updating residue counts in the topology.

No current component performs step 2. The pipeline is a complete sensor system attached to an incomplete actuator system: it can see the problem in full detail, but it has no arm that reaches into the gap.

---

## 7. A General Solution Architecture

The solution must work for:
- GPCRs (helical TM bundles, ~2.5 nm diameter)
- Ion channels (often homo-oligomers with larger footprints, ~4–6 nm)
- Transporters (TM barrels and bundles of varying shape)
- Oligomeric proteins (irregular XY footprints; non-circular gaps)
- Proteins with large soluble domains (ECD/ICD that extend beyond the bilayer; XY centroid displaced from TM center)

The general algorithm must operate on the **TM surface geometry**, not on any protein-specific residue numbering or shape. The following interface is sufficient:

- Input: a `.gro` file with protein + bilayer, a set of TM residue numbers, a lipid residue name
- Derived automatically: gap patch coordinates (from Phase 9A gap_clusters), the per-residue `mean_distance_by_residue` field
- Output: a modified `.gro` with improved lateral coverage, consistent atom/residue count and box

---

## 8. Proposed Architecture: Protein–Membrane Interface Builder (Phase 10)

### Components

#### 10A. Hydrophobic Surface Map
From the TM residue set and their `_Atom` records (already available), compute the set of "TM surface atoms": those with fewer than 18 protein-neighbor atoms within 0.5 nm (the current Phase 9A definition). Group them into contiguous gap clusters.

#### 10B. Gap Patch Target Generator
For each gap cluster, compute:
- the cluster's 3D centroid
- an inward surface normal from the cluster center toward the bilayer interior
- a target XY position for a lipid headgroup: centroid displaced inward by `target_distance_nm` (0.45 nm) along the surface normal projected to XY

This yields a set of `N_gap_patches` target XY positions where lipids should be placed.

#### 10C. Lipid Placement Module (one of three strategies — see §10)

#### 10D. Clash Filter
After placement, evaluate each introduced lipid atom against all protein atoms with a hard-sphere distance check (cutoff = 0.12 nm). Remove any introduced lipid that violates this check. Prefer removing the whole lipid residue atomically rather than individual atoms.

#### 10E. Topology Update
Count the final number of each lipid species. Update the `[molecules]` section of `topol.top` accordingly. This is already implemented in `update_topology_lipid_count`.

#### 10F. Short Restrained Relaxation
Run a brief position-restrained minimization + 500 ps NVT MD with:
- protein heavy atoms restrained (1000 kJ/mol/nm²)
- lipid and water atoms free
- small time step (1 fs) and moderate temperature (298 K)

This allows lipids to settle into the gap without enough thermal energy to escape the membrane midplane.

#### 10G. Interface Re-evaluation
Run Phase 9A `evaluate_protein_membrane_interface` on the relaxed structure. If `quality_passed` is True, proceed. If not, report the remaining gap and warn (policy = "warn") or block (policy = "strict").

---

## 9. Design Sketch

```
embed_in_bilayer  →  [existing merge + initial exclusion]
        │
        ▼
membrane_embedding  →  [existing shrink loop → converged.gro]
        │
        ▼
  ┌─────────────────────────────────────────────────────────┐
  │  Phase 10: Protein–Membrane Interface Builder (NEW)     │
  │                                                         │
  │  10A: hydrophobic surface map + gap cluster detection   │
  │  10B: target XY positions from gap centroids           │
  │  10C: lipid placement (Strategy A, B, or C)            │
  │  10D: clash filter                                      │
  │  10E: topology count update                             │
  │  10F: short restrained NVT relaxation                   │
  │  10G: Phase 9A re-evaluation + quality gate            │
  │                                                         │
  │  Output: interface_packed.gro                           │
  └─────────────────────────────────────────────────────────┘
        │
        ▼
  solvate_membrane  →  [existing]
```

The Interface Builder runs **after** the shrink loop converges, so it operates on a bilayer that already has the correct global APL. Its job is to fix local under-coverage without disturbing the global packing.

---

## 10. Three Corrective Strategies Evaluated

### Strategy A — Lateral Redistribution of Existing Lipids

**Mechanism:** Identify overpacked XY regions (APL locally below target) and under-packed regions (the gap clusters). Translate selected lipid molecules laterally from the overpacked source region into the gap target region.

| Criterion | Assessment |
|-----------|-----------|
| Physical realism | Moderate. Lipid diffusion is the natural process, but moving a whole residue as a rigid body is not. The resulting configuration will have strained bonds and require relaxation. |
| Implementation difficulty | Moderate. Requires local APL mapping, a donor selection algorithm, rigid-body translation, and clash filtering. ~500 lines of Python. |
| Determinism | High. Given the same gap map and the same donor selection rule, the output is reproducible. |
| Topology implications | None. Same number of lipids, same residues. No `.itp` changes. |
| Risk of artifacts | Moderate. The moved lipid may be sterically strained in its new position even after minimization. The donor region becomes transiently under-packed. |
| Suitability for automation | High. No external tool calls. No new parameters. Runs inside a Python script. |

**Best when:** The bilayer has locally dense regions (common with flat bilayer + irregular protein shape). May not help if the bilayer is already uniformly packed at the correct APL — there are no donor lipids.

---

### Strategy B — Template Insertion from a Lipid Library

**Mechanism:** Copy a pre-equilibrated lipid molecule from a template library (e.g., a single DPPC from `dppc512_whole.gro`). Place its headgroup at each gap-target XY position, at the correct Z (bilayer midplane ± half-thickness for upper/lower leaflet). Assign the correct residue number. Apply a random rotation around Z.

| Criterion | Assessment |
|-----------|-----------|
| Physical realism | High. The lipid conformation is taken from a real bilayer equilibrium structure. The inserted lipid is physically valid. |
| Implementation difficulty | Moderate. Requires a lipid template parser, a leaflet assignment rule (Z > midplane → upper leaflet), a clash filter, and topology update. ~400 lines of Python plus a template GRO. |
| Determinism | Moderate. The Z-rotation introduces randomness. Determinism can be recovered by fixing the random seed, but makes testing easier. |
| Topology implications | **New lipid residues are added.** The topology count increases. The `topol.top` `[molecules]` section must be updated. If the `.itp` file uses a fixed molecule count, it must also be updated. This is manageable but requires care. |
| Risk of artifacts | Low if clash filtering is strict. An improperly placed lipid (overlapping with protein) will be caught by the clash filter and removed. The remaining lipids are drawn from a real bilayer. |
| Suitability for automation | High. The bilayer template file is already staged in `membrane_assets/` by the workspace builder. |

**Best when:** The annular gap is large (many missing lipids), the bilayer is already at target APL (no donor regions available), and the topology system can tolerate additional lipid residues.

---

### Strategy C — Short Restrained Lateral-Relaxation MD

**Mechanism:** Run a short NVT MD simulation (500 ps – 2 ns) with protein heavy atoms fully restrained (1000–5000 kJ/mol/nm²). Lipid and water atoms are free. Thermal motion at physiological temperature (298–310 K) provides enough kinetic energy for lipids to diffuse into the annular gap via the natural Brownian-motion-driven lateral diffusion of membrane lipids.

| Criterion | Assessment |
|-----------|-----------|
| Physical realism | **Highest.** This is exactly what happens in real membranes after an insertion event. No artificial moves, no rigid-body translations, no library interpolation. |
| Implementation difficulty | High. Requires GROMACS, a position-restraint `.itp` for the protein, an NVT MDP file tuned for lipid relaxation, and topology integrity. This is an entire GROMACS run, not a Python script. |
| Determinism | Low without velocity seeds. MD trajectories diverge from initial conditions within picoseconds. The final structure depends on the random seed. This makes automated testing difficult. |
| Topology implications | None. No molecules added or removed. |
| Risk of artifacts | Low with proper restraints. The main risks are: (a) protein conformational drift if restraint force constants are too low; (b) lipid flip-flop (headgroup passing through the hydrophobic core) at elevated temperatures; (c) water entry into membrane gaps if the relaxation is too long. All manageable. |
| Suitability for automation | Moderate. GROMACS is already required for the main workflow. However, the relaxation MD adds ~10–60 minutes of wall time per run (GPU dependent). The output structure must be post-processed to remove water that entered the membrane. |

**Best when:** The existing bilayer is uniform but the gap is too large for Strategy A. Time is available. The force field is already parametrized for the lipid type (which it is — DPPC + OPLS-AA is the current target).

---

## 11. Recommendation for Phase 10

### Immediate priority: Strategy B (template insertion) as the primary actuator

**Rationale:**

1. Strategy A requires donor regions. A bilayer at correct APL after a full shrink loop may not have sufficient local excess density to donate lipids to the gap. Strategy B does not depend on the existing lipid distribution.

2. Strategy C is physically ideal but adds significant runtime and introduces non-determinism that is incompatible with automated testing and reproducible workspace generation.

3. Strategy B is deterministic (with a fixed seed), compatible with the current workspace structure (the template bilayer is already staged), and requires no external tool invocations beyond GROMACS minimization.

4. The topology update is handled by `update_topology_lipid_count`, which is already implemented.

5. Strategy C should be **provided as an optional backend** behind a `interface_builder_relaxation_md: true` flag for users who want maximum physical realism and can tolerate the runtime cost.

### Implementation sequence

1. **Phase 10A/B** — Gap target generator: reads Phase 9A `gap_clusters` from `protein_membrane_interface_report.json`, computes per-cluster centroid and target XY positions. Pure Python, no GROMACS. Unit-testable.

2. **Phase 10C (Strategy B)** — Lipid template inserter: reads `membrane_assets/dppc512_whole.gro`, extracts one lipid residue, places N copies at target positions (upper/lower leaflet from Z), generates output GRO. Pure Python.

3. **Phase 10D** — Clash filter: atoms within 0.12 nm of any protein atom → remove whole residue. Pure Python.

4. **Phase 10E** — Topology count update: already implemented, wire it in.

5. **Phase 10F** — Restrained minimization: reuse the existing `minim_shrink.mdp` + one GROMACS call. The existing `_write_minim_mdp` method is sufficient.

6. **Phase 10G** — Re-evaluation: call `evaluate_protein_membrane_interface` on the output; assert `quality_passed` or issue policy-controlled warning.

### What Phase 10 must NOT do

- It must not call inflategro again (would re-create the gap via the exclusion filter).
- It must not delete lipids as a corrective action for under-coverage.
- It must not depend on GLP-1R-specific residue numbers — all operations derive from gap_cluster coordinates.
- It must not require a separate force field or lipid template beyond what is already in `membrane_assets/`.

---

## 12. Why Did All Previous Attempts Fail?

Every phase from Phase 1 through Phase 9B addressed a real problem, but none of them addressed annular gap under-coverage:

| Phase | What it fixed | Did it touch lateral packing? |
|-------|--------------|-------------------------------|
| 1 — TM-center alignment | Z-misalignment of hydrophobic core | No |
| 2 — TM-aware exclusion | Lipid invasion into TM cavity | No — made lateral gap worse |
| 3 — Surface interference gate | Lipid clashes with EC/IC domains | No |
| 4 — Box match / lateral margin | System size mismatch | No |
| 5 — Trapped lipid diagnostic | Reporting lipids in protein cavity | No |
| 6 — Embedding quality (APL, void, burial) | Global APL convergence and TM burial | No — APL is a global metric |
| 7 — Clean water gate | Water in membrane interior | No |
| 8 — TM exclusion mask | Lipid overlaps with TM surface | No — also worsens lateral gap |
| 9A — Interface evaluator | **Detects** the annular gap correctly | No — sensor only |
| 9B — Embedding optimizer | Searches Z-shift + mask-padding grid | No — search space cannot contain solution |

The root cause is not a bug in any individual phase. It is an architectural gap: **the pipeline was built incrementally to correct specific failure modes observed during development (invasion, Z-misalignment, clashes), and the under-coverage failure mode was not addressed because the diagnostic for it (Phase 9A) was added only after all the other corrections were already in place.**

Now that Phase 9A makes the problem observable with quantitative metrics, Phase 10 can be designed to fix it precisely.

---

## 13. What the General Algorithm Should Be

```
INPUT:  converged.gro  (post-shrink-loop system with correct global APL)
        protein_membrane_interface_report.json  (gap_clusters, per-residue distances)
        membrane_assets/dppc512_whole.gro  (template bilayer for lipid insertion)
        tm_residues  (set of residue integers)

ALGORITHM:
  1. Read gap_clusters from Phase 9A report.
     If len(gap_clusters) == 0 or fraction_exposed_gap <= threshold:
         exit 0  # no action needed

  2. For each gap cluster c in gap_clusters:
     a. Compute cluster centroid (x_c, y_c, z_c)
     b. Assign leaflet: z_c > bilayer_midplane_z → upper, else lower
     c. Compute target headgroup position:
            x_t = x_c + 0.45 * outward_normal_x
            y_t = y_c + 0.45 * outward_normal_y
            z_t = leaflet_headgroup_z  (from current bilayer headgroup mean Z)
     d. Extract one lipid molecule from template bilayer
     e. Translate lipid so its headgroup reference atom is at (x_t, y_t, z_t)
     f. Apply random Z-rotation (seed fixed for reproducibility)
     g. Check clash with all protein atoms (cutoff 0.12 nm)
        → if clash: try up to 4 rotations, then discard if all clash
     h. If accepted: append to atom list, assign new residue number

  3. Write output GRO with original atoms + accepted new lipids

  4. Update topol.top lipid count += n_accepted

  5. Run steepest-descent minimization (minim_shrink.mdp, maxcycles=2000)

  6. Run Phase 9A re-evaluation on minimized structure

  7. Write interface_packed.gro

OUTPUT: interface_packed.gro
        interface_packing_report.json
```

This algorithm is:
- Protein-shape-agnostic (all geometry comes from gap cluster coordinates)
- Topology-consistent (count update is explicit)
- Testable (fixed seed, deterministic output)
- Composable (wraps around the existing Phase 9A evaluator for the quality gate)
- Independent of inflategro (does not re-run the shrink loop)
