# Membrane Protein MD — Architectural Design

> Analysis of `Prot-Memb_FILES/` to guide SimForge membrane pipeline implementation.
> **Status:** Knowledge capture — no code written yet. Implementation targets: `membrane_knowledge.py` → `adapters/` → validators → `MembraneWorkflowOPLSAA`.

---

## 1. Scientific reasoning behind each workflow step

The 24-step tutorial is not a checklist — it encodes physical constraints:

| Step | Scientific Why |
|---|---|
| Orientation (editconf -princ) | Align protein inertia axes; TM helices must be Z-parallel for bilayer embedding |
| Rotation (270° around Y) | Conclusion from structure, not a parameter — correct orientation means cytosolic face at −Z, extracellular at +Z |
| Box size match | Bilayer pre-built at fixed XY (12.84 × 12.89 nm); protein box must match exactly or PBC fails |
| Concatenate + atom count fix | GROMACS .gro format: atom count on line 2 must be exact; any mismatch aborts grompp |
| MoveMemb (Z-shift) | Place bilayer below protein without overlap; amount is protein-height-dependent |
| Strong position restraints | Prevent protein migration during shrink loop; 100000 kJ/mol/nm² — far stiffer than standard posre |
| InflateGRO (scale 4.0) | Expand lipid XY to create space around protein; factor ~4 for 512-lipid DPPC |
| Shrink loop (scale 0.95) | Iterative compaction; each iteration = minimize + deflate + measure APL; converges when APL ≤ target |
| APL convergence criterion | Physical observable: DPPC at 298K → ~62 Å²/lipid; above = too much space, below = over-compressed |
| Solvate | Add TIP3P/SPC water; `spc216.gro` used for SPC water placement |
| water_deletor | Waters placed inside bilayer core by solvate must be removed — bilayer interior is hydrophobic |
| Ion neutralization | Standard protocol; physiological 0.154 M NaCl |
| Minimization (with POSRES) | Relax clashes from concatenation; STRONG_POSRES keeps protein/lipid geometry intact |
| NVT 100ps | Thermalize system at 298K while keeping structure restrained |
| NPT 1ns | Allow box + bilayer area to equilibrate; semiisotropic coupling is non-negotiable here |
| Production MD (500ns) | Release restraints; switch to Nosé-Hoover + Parrinello-Rahman for NVT ensemble accuracy |

---

## 2. Implicit scientific decisions (not visible as params)

These are the decisions the workflow crystallizes but never states explicitly:

**Rotation angle is a conclusion, not a free parameter.**
The `-rotate 0 270 0` is the result of visualizing the structure after `-princ`. A different protein may need 90°, 0°, or no rotation. SimForge cannot hardcode this — it must be detected (DeepTMHMM → helix axis → box axis alignment) or asked.

**APL target is lipid+temperature dependent.**
- DPPC at 298K → 62 Å²
- DPPC at 323K → ~68 Å²
- POPC at 298K → ~65 Å²
- POPE at 298K → ~56 Å²
These are experimental values, not computed. They must live in `membrane_knowledge.py`.

**Inflation factor scales with system size and protein footprint.**
Factor 4.0 is calibrated for 512-lipid DPPC with a single-pass TM helix. Larger proteins (GPCRs, channels) may need factor 5–6. The protocol states "increase for large proteins" — heuristic, not formula.

**Z-displacement for MoveMemb is manual.**
The Fortran program is interactive — it asks for nm to shift. The correct value is `z_center_protein - z_center_bilayer`. SimForge must compute this from the .gro coordinates or provide a validator that checks for overlap before embedding.

**Semiisotropic pressure coupling is non-negotiable for membranes.**
Isotropic coupling would compress XY and Z together, collapsing the bilayer. XY pressure ref = 0.5 bar (membrane tension), Z = 0.5 bar — different physical meaning from soluble protein NPT. This is not a tunable knob.

**dt = 0.001 fs in production (not 0.002).**
Lipid force fields (OPLS-AA membrane) require smaller timestep than protein-only. The NVT/NPT use 0.002, production drops to 0.001. This is a stability constraint from the OPLS-AA lipid parameterization.

**tc-grps = system (not Protein + Non-Protein).**
During equilibration, the membrane system uses a single thermostat group. Splitting into Protein + Lipid + Solvent risks "hot solvent / cold protein" artifacts during the shrink loop when lipid degrees of freedom are suppressed.

---

## 3. Failure signals the workflow is designed to prevent

Each failure mode maps to a specific protocol step:

| Failure | Signal | Prevented by |
|---|---|---|
| TM helix not in bilayer | Helix Z-coordinates outside bilayer Z-range | Rotation + orientation validation |
| Lipid clash with protein | GROMACS atom-atom overlap crash at grompp | InflateGRO expansion + shrink loop |
| APL non-convergence | APL oscillates above target, never converges | Correct inflation factor; may need manual restart with higher factor |
| Water in bilayer core | Water oxygen between headgroup and tail Z-coords | water_deletor.pl |
| Topology atom count mismatch | grompp fatal: "atoms in .top != atoms in .gro" | Manual atom count update in concat step |
| LINCS warnings at equilibration | Too many bond constraint failures | Strong POSRES during NVT; reduce dt |
| Lipid flip-flop at NPT | Lipid migrates from one leaflet to other | Adequate equilibration; check if APL converged |
| Incorrect box Z dimension | Protein extends outside box | Box Z must be > protein height + bilayer thickness + 2× water layer |

---

## 4. Remediations the human performs when things fail

These are the expert judgment calls that the shrink loop script cannot make:

- **APL stalls above target**: Increase inflation factor (4.0 → 5.0), restart from pre-inflation geometry.
- **APL drops below target and bounces**: Deflation factor too aggressive (0.95 → 0.97).
- **MoveMemb displacement wrong**: Re-inspect concatenated .gro in VMD, estimate visual overlap, re-run with different nm value.
- **Protein tilts during shrink**: Increase STRONG_POSRES force constant (100000 → 500000) or add backbone restraints.
- **Water count wrong in topology**: Open topol.top, find `SOL` line, update integer manually.
- **Rotation wrong**: Re-visualize after `-princ`, adjust rotation angles, re-run from step 2.
- **LINCS failure during NVT**: Reduce dt from 0.002 to 0.001 in nvt.mdp; inspect .log for worst-offending bonds.

---

## 5. SimForge module mapping

Where each piece of the workflow belongs in the SimForge architecture:

```
DeepTMHMM analysis     → core/scientific_planner.py (PROTOCOL_SELECTION question)
                          adapters/deepTMHMM_adapter.py (web service call)

Protein orientation    → builders/step_builders/membrane_orient_builder.py
  - editconf -princ
  - rotation detection/question

Box sizing             → builders/step_builders/membrane_box_builder.py
  - reads bilayer .gro dimensions
  - sets matching box

Concatenate + fix      → builders/step_builders/membrane_embed_builder.py
  - cat protein + bilayer .gro
  - Python atom count fix

MoveMemb (Z-shift)     → adapters/movememb_adapter.py
  - computes Z-overlap from coordinates
  - runs Fortran or equivalent Python reimplementation

Strong POSRES          → builders/step_builders/membrane_embed_builder.py
  - gmx genrestr with fc 100000 100000 100000

InflateGRO             → adapters/inflategro_adapter.py
  - wraps inflategro-Jorge.pl
  - inputs: gro, scale, lipid_name, cutoff, output, gridsize, area_dat

Shrink loop            → builders/step_builders/shrink_loop_builder.py (meta-step)
  - generates shrink_loop.sh that mimics ScriptCamilo-Jorge.sh
  - APL convergence criterion from membrane_knowledge.py
  - internal: grompp → mdrun → inflategro → AperR → check → repeat

water_deletor          → adapters/water_deletor_adapter.py
  - wraps water_deletor.pl
  - inputs: in, out, ref atom (O33), middle atom (C50), nwater

Ion neutralization     → builders/step_builders/solvation_builder.py (shared with protein pipeline)

MDP files              → builders/step_builders/membrane_equilibration_builder.py
  - generates nvt-memb.mdp (semiisotropic, tc-grps=system)
  - generates npt-memb.mdp (Berendsen, semiisotropic, ref_p=0.5)
  - generates md-memb.mdp (Nosé-Hoover + Parrinello-Rahman, dt=0.001)

APL monitoring         → validators/membrane_validators.py
  - apl_converged(area_dat_path, target_apl, tolerance=2.0) → bool
  - protein_in_bilayer(gro, bilayer_z_range) → bool
  - no_water_in_bilayer(gro, headgroup_atom, tail_atom) → bool
```

---

## 6. OPLS-AA specific vs membrane-specific vs generalizable

| Feature | OPLS-AA specific | Membrane specific | Generalizable |
|---|---|---|---|
| `oplsaa_membrane.ff` directory | ✓ | | |
| DPPC lipid residue name "DPP" | | ✓ (DPPC=DPP, POPC=POP) | |
| APL target 62 Å² | ✓ (OPLS-AA DPPC 298K) | | |
| Headgroup atom O33 | ✓ OPLS-AA DPPC numbering | | |
| Tail atom C50 | ✓ OPLS-AA DPPC numbering | | |
| rcoulomb/rvdw = 1.2 nm | | ✓ membrane standard | |
| DispCorr = EnerPres | | ✓ membrane standard | |
| dt = 0.001 in production | ✓ OPLS-AA lipid constraint | | |
| Semiisotropic coupling | | ✓ | |
| ref_p = 0.5 bar (XY) | | ✓ DPPC surface tension | |
| tc-grps = system | | ✓ shrink loop phase | |
| Nosé-Hoover + PR in production | | | ✓ any long MD |
| 512-lipid bilayer geometry | | ✓ | |
| InflateGRO algorithm | | ✓ | |

---

## 7. The DAG loop problem

The shrink loop is fundamentally not a DAG step. It violates DAG acyclicity:

```
minimize → deflate → measure_APL → [if APL > target: go back to minimize]
```

**Chosen approach: meta-step.**

Treat the entire InflateGRO shrink loop as a single opaque DAG node with:
- `step_type: "shrink_loop"`  
- Internal convergence managed by the generated shell script
- Inputs: inflated system, APL target, deflation factor, max iterations
- Outputs: converged system GRO + final APL value
- The builder generates a self-contained `shrink_loop.sh` analogous to `ScriptCamilo-Jorge.sh`

This is the correct v1 choice. Dynamic DAG rewriting (hot workflow mutation) is explicitly deferred to v3.

The convergence criterion lives in `membrane_knowledge.py` so the builder can query it:
```python
membrane_knowledge.apl_target(lipid="DPPC", forcefield="opls-aa", temperature_K=298) → 62.0
```

---

## 8. Proposed DAG for `MembraneWorkflowOPLSAA`

```
predict_tm_helices          [optional, DeepTMHMM adapter]
       ↓
orient_protein              [editconf -princ + rotation]
       ↓
size_box_to_bilayer         [editconf -box matching bilayer XY]
       ↓
embed_in_bilayer            [cat + atom count fix + MoveMemb + strong POSRES]
       ↓
generate_topology           [pdb2gmx with oplsaa_membrane.ff]
       ↓
shrink_loop (meta-step)     [inflate → minimize loop until APL converged]
       ↓
solvate                     [gmx solvate]
       ↓
clean_water                 [water_deletor.pl]
       ↓
add_ions                    [gmx genion]
       ↓
minimize_full               [em.mdp with -DPOSRES -DSTRONG_POSRES]
       ↓
equilibrate_nvt             [nvt-memb.mdp, 100ps, V-rescale, single group]
       ↓
equilibrate_npt             [npt-memb.mdp, 1ns, Berendsen, semiisotropic]
       ↓
production_md               [md-memb.mdp, 500ns, NH+PR, dt=0.001]
       ↓
analysis                    [APL over time, RMSD, bilayer thickness, order params]
```

---

## 9. New layers required (not yet in SimForge)

### `core/membrane_knowledge.py`
Single source of truth for lipid physical constants. No builder or validator should hardcode these.

Key contents:
- `apl_target(lipid, forcefield, temperature_K) → float`
- `lipid_residue_name(lipid, forcefield) → str`  (e.g. DPPC+opls-aa → "DPP")
- `headgroup_atom(lipid, forcefield) → str`  (e.g. "O33")
- `tail_atom(lipid, forcefield) → str`  (e.g. "C50")
- `n_lipids_for_box(box_xy_nm) → int`  (512 for 12.84×12.89)
- Inflation factors, convergence thresholds, deflation rates

### `adapters/` layer
New abstraction: `ExternalToolAdapter` base class for tools that are not GROMACS.

Contracts (before implementation):
- `InflateGROAdapter.run(gro_in, scale, lipid_name, cutoff, gro_out, gridsize, area_dat)`
- `WaterDeletorAdapter.run(gro_in, gro_out, ref_atom, middle_atom, nwater)`
- `MoveMembAdapter.compute_z_shift(protein_gro, bilayer_gro) → float`; `run(gro_in, z_shift_nm, gro_out)`
- `DeepTMHMMAdapter.predict(fasta_path) → List[TMHelix]`  (web service, async)

### `validators/membrane_validators.py`
Checks before automation, not after:
- `validate_orientation(gro, tm_helices) → ValidationResult`
- `validate_apl_convergence(area_dat_path, target, tolerance) → ValidationResult`
- `validate_no_water_in_bilayer(gro, ref_atom, middle_atom) → ValidationResult`
- `validate_protein_in_bilayer(gro, helix_residues, bilayer_z_range) → ValidationResult`

---

## 10. Implementation priority (explicit ordering)

1. **`core/membrane_knowledge.py`** — knowledge first, before any builder consumes it
2. **`adapters/` contracts** — base class + interface definitions; stub implementations
3. **`validators/membrane_validators.py`** — build validators against knowledge layer
4. **`pipelines/membrane_pipeline.py`** — `MembraneWorkflowOPLSAA` specifically, no premature generalization
5. **Builders** — one builder per DAG step, consuming knowledge + adapters
6. **`core/scientific_planner.py` extensions** — membrane-specific questions (lipid, orientation, protocol)

**Explicitly deferred:**
- Multi-lipid / multi-bilayer generalization
- Adaptive DAG (hot workflow mutation)
- DeepTMHMM web adapter (implement stub, manual for now)
- `MembranePipeline` generalized superclass

---

## 11. Key MDP parameters reference

### Shrink-loop minimization (`minim.mdp`)
```
define      = -DPOSRES -DSTRONG_POSRES
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
```

### NVT equilibration (100 ps)
```
nsteps      = 25000    ; dt=0.002 → 100ps
tcoupl      = V-rescale
tc-grps     = system
tau_t       = 0.1
ref_t       = 298
pcoupl      = no
gen_vel     = yes
```

### NPT equilibration (1 ns)
```
nsteps      = 150000   ; dt=0.002 → 1ns  (note: script says 500000 but nsteps is 150000)
tcoupl      = V-rescale
pcoupl      = Berendsen
pcoupltype  = semiisotropic
ref_p       = 0.5  0.5
compressibility = 4.5e-5  4.5e-5
```

### Production MD (500 ns)
```
nsteps      = 50000000 ; dt=0.001 → 500ns
tcoupl      = Nose-Hoover
tc-grps     = System
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
ref_p       = 1.0  1.0
dt          = 0.001    ; critical: OPLS-AA lipid stability constraint
constraints = h-bonds  ; not all-bonds in production
```

---

*Last updated: 2026-05-22. Source: `Prot-Memb_FILES/` analysis.*
