# Design Document: Phase 7 — TM-Aware Membrane Packing Backend

## 1. Context & Motivation

Our scientific audit of the G-protein coupled receptor (GPCR) GLP-1R membrane embedding identified two fundamental limitations in the classic `InflateGRO` method:
1. **Soluble Domain Interference**: Large, bulky extracellular (ECD) and intracellular (ICD) soluble domains extend far beyond the transmembrane (TM) bundle in the XY plane. Uniform radial box scaling causes these domains to contact and clash with lipids prematurely, halting deflation before the bilayer midplane lipids can pack against the narrower TM bundle. This leaves large persistent vacuum voids around the TM helices.
2. **Simplified Overlap Criteria**: Checking only headgroup Phosphorus (P) to protein C-alpha (CA) distances misses steric clashes of lipid hydrocarbon tails, side chains, and cavity-trapped lipids.

To resolve these issues, we design a new **TM-Aware Membrane Packing** backend. This backend isolates membrane-plane packing from soluble domain footprints and uses 3D voxel-based excluded volume reconstruction to ensure correct lipid packing and cavity-aware lipid exclusion.

---

## 2. Configuration & Legacy Support

The classic `InflateGRO` backend remains fully available to maintain backwards compatibility. We introduce a new configuration parameter:

```yaml
membrane:
  embedding:
    backend: "tm_aware"      # Options: "inflategro" (classic) | "tm_aware" (default in Phase 7)
    cutoff_nm: 0.14          # Cutoff distance for local steric checks
    safety_limit: 50         # Max lipids allowed to be excluded before triggering safety gate
```

---

## 3. Scientific Model

The new backend operates on a Z-segmented physical model:
- **Transmembrane (TM) Excluded Volume**: Only the protein residues annotated as part of the transmembrane segments (`tm_residues`) define the lateral membrane exclusion zone.
- **Soluble Domain Exemption**: Intracellular and extracellular domains that lie outside the bilayer core ($Z < z_{\text{core\_bot}}$ or $Z > z_{\text{core\_top}}$) are excluded from the lateral membrane footprint and do not prevent lipids from packing in XY.
- **Lipid Packing Quality**: Packing quality is no longer evaluated using a global Area Per Lipid (APL). Instead, it is validated using localized density profiles, 2D grid-based void analysis at the midplane, TM burial fraction, and angular-gap cavity-trapping checks.

```
       Extracellular Domain (ECD) - Wide Footprint
             [=========]         <-- Excluded from lateral packing footprint
     ~~~~~~~~~~~~~~~~~~~~~~~~~~  <-- Bilayer Upper Leaflet (Packs under ECD)
     |  Lipids pack   |  Lipids| 
     |  against TM    |  pack  |  <-- Bilayer Midplane (No Voids!)
     |  helices       |  here  | 
     ~~~~~~~~~~~~~~~~~~~~~~~~~~  <-- Bilayer Lower Leaflet (Packs over ICD)
             [=========]
       Intracellular Domain (ICD) - Wide Footprint
```

---

## 4. Algorithm v1: Voxel-Masking & Local Relaxation

The `tm_aware` backend executes the following algorithmic steps:

### Step 1: TM Domain Isolation & Z-Core Definition
1. Identify the Z-boundaries of the TM domain based on the C-alpha coordinates of `tm_residues`:
   $$z_{\text{min, TM}} = \min(z_{\text{CA}}), \quad z_{\text{max, TM}} = \max(z_{\text{CA}})$$
2. Align the bilayer coordinates so its midplane $z_{\text{mid}}$ matches the centroid of the TM domain:
   $$z_{\text{mid}} = \frac{z_{\text{min, TM}} + z_{\text{max, TM}}}{2}$$

### Step 2: 3D Voxel Mask Construction
1. Define a 3D grid of voxels of spacing $0.1$ nm over the simulation box.
2. Filter the protein atoms to keep only those belonging to `tm_residues` or lying inside $Z \in [z_{\text{min, TM}}, z_{\text{max, TM}}]$.
3. For each Z-slice of the grid:
   - Identify voxels within the van der Waals radius of any TM atom (steric exclusion).
   - Compute the 2D alpha-shape (concave hull) of the TM atoms in the slice.
   - Mark all voxels *inside* the concave hull as part of the **TM Excluded Volume** (this automatically encloses internal cavities, channels, and pores).

### Step 3: Pre-Shrink Lipid Exclusion
1. For each lipid molecule in the initial bilayer:
   - Check if any of its atoms lie inside the 3D **TM Excluded Volume**.
   - If an overlap is detected, flag the entire lipid residue for deletion.
2. Remove the flagged lipid residues from the coordinates (`system.gro`).
3. Renumber residues to ensure sequential numbering starting from 1.

### Step 4: GROMACS Topology Generation
- Run `gmx pdb2gmx` directly on the cleaned coordinate file. The resulting GROMACS topology (`topol.top`) will automatically match the correct number of lipids, avoiding manual topology hacks.

### Step 5: Lateral Pressure Relaxation
1. Run GROMACS energy minimization with strong position restraints on all protein atoms to resolve local atomic clashes.
2. Run a short (e.g. 200 ps) restrained NPT equilibration using a semi-isotropic barostat (keeping $P_z$ decoupled or fixed, while adjusting $P_{xy}$ to 1 bar). 
   - Under lateral pressure, the lipids flow dynamically to pack tightly against the irregular outer surface of the TM domain, naturally closing any local gaps without requiring an artificial geometric shrink loop.

---

## 5. Safety Gates

If any of the following gates fail, the step exits with code 1 and blocks downstream steps:
1. **Lipid Deletion Limit**: The number of excluded lipids exceeds the `safety_limit` (typically 50). This prevents bulk bilayer deletion due to incorrect orientation or wrong annotations.
2. **Cavity Trapping**: Any lipid residue COM is detected inside the TM bundle (angular gap of surrounding TM C-alpha atoms $< 185^\circ$) after relaxation.
3. **High Local Void Fraction**: The total unoccupied void area at the bilayer midplane exceeds $1.0\text{ nm}^2$.
4. **TM Burial Failure**: The fraction of TM C-alpha atoms vertically buried inside the hydrophobic core ($Z \in [z_{\text{core\_bot}}, z_{\text{core\_top}}]$) is $< 80\%$.
5. **Soluble Domain Penetration**: More than 10 non-TM protein atoms lie vertically inside the bilayer core.

---

## 6. Generated Reports

The backend will produce three JSON reports in the step directory:

### A. `tm_mask_report.json`
Details on the generated TM mask (voxel size, Z-slices, area per slice, cavity volume):
```json
{
  "grid_spacing_nm": 0.1,
  "tm_z_min": 5.152,
  "tm_z_max": 9.240,
  "total_excluded_voxels": 45120,
  "estimated_tm_volume_nm3": 45.12,
  "cavity_volume_nm3": 8.45
}
```

### B. `lipid_exclusion_report.json`
Details on lipids deleted during pre-shrink exclusion:
```json
{
  "n_lipids_checked": 512,
  "n_lipids_removed": 24,
  "removed_residue_ids": [453, 454, 455, 460, 514, 527, 565, 572, 577, 580, 584, 593, 656, 674, 688, 699, 703, 713, 717, 718, 825, 842, 843, 904],
  "removal_reasons": {
    "direct_clash_with_tm_helices": 12,
    "geometrically_trapped_in_cavity": 12
  }
}
```

### C. `embedding_quality_report.json`
Final quality metrics of the relaxed system:
```json
{
  "expanded_state_ok": true,
  "shrink_convergence_profile": [
    {
      "step": "initial_embed",
      "void_area_nm2": 4.12,
      "trapped_lipid_count": 0
    },
    {
      "step": "post_minimization",
      "void_area_nm2": 1.25,
      "trapped_lipid_count": 0
    },
    {
      "step": "final_relaxation",
      "void_area_nm2": 0.18,
      "trapped_lipid_count": 0
    }
  ],
  "final_void_score": 0.18,
  "final_trapped_lipid_count": 0,
  "final_tm_burial_score": 0.94,
  "soluble_domain_core_penetration": 2,
  "warnings": []
}
```

---

## 7. Verification Test Suite

We will add the following synthetic test cases in `tests/test_embedding_quality.py`:

1. **GPCR with Large ECD**:
   - Construct a synthetic protein consisting of a narrow cylindrical TM region (4 helices) and a wide extracellular dome (width 3x TM region).
   - Place a uniform lipid bilayer with overlapping and cavity-trapped lipids.
   - Verify that the `tm_aware` backend deletes the overlapping and cavity-trapped lipids but does NOT delete lipids adjacent to the wide extracellular dome, resulting in 0 voids and 0 trapped lipids after relaxation.
2. **Channel-Like Cylindrical Protein**:
   - Construct a hollow cylinder representing an ion channel pore.
   - Verify that the voxel-mask correctly identifies the internal channel volume and excludes lipids from the pore.
3. **Oligomeric Protein**:
   - Construct a multimeric assembly with asymmetric interfaces.
   - Verify that the concave hull detects inter-subunit cavities and properly excludes lipids from interfaces.
4. **Soluble-Domain False Exclusion Case**:
   - Construct a soluble protein with no TM annotation.
   - Verify that the safety gate triggers or defaults to the classic legacy path without performing bulk deletions.
