# SimForge — Declarative System Build

`simforge build` turns a single YAML *system specification* into a
GROMACS-ready system (topology + coordinates + provenance + validation),
replacing a 15-line `simforge ligand integrate` invocation.

It is a thin, declarative front end to the already-validated assembly engine
(`ligand.integrate.assemble_system_multi`): protein + **N** independently
parameterized non-protein components, with automatic atomtype namespace
isolation. No assembly logic is duplicated.

```
simforge build system.yaml
simforge build system.yaml --out competitive_system/
simforge build system.yaml --no-grompp        # skip the gmx grompp check
simforge build system.yaml --verbose           # print full atomtype rename map
```

Related commands:

| Command | Purpose |
|---|---|
| `simforge doctor` | Inspect the runtime environment (GROMACS, RDKit, force fields, conda env). |
| `simforge validate-system SYSTEM_DIR` | Structured validation of a built system, without modifying it. |
| `simforge inspect-run RUN_DIR` | Summarize a built system's `provenance.json`. |

---

## The specification

The document root may contain the keys directly, or nest them under a
`system:` key. Every relative path is resolved **relative to the YAML file's
own directory**, never the current working directory.

```yaml
system:
  name: hmgcoa_competitive          # optional — used to name the output dir

  protein:
    structure: protein_only.gro     # .pdb / .ent  -> raw mode (runs pdb2gmx)
                                    # .gro          -> preparameterized mode
    topology:                       # preparameterized mode only
      - topol_Protein_chain_A.itp
      - topol_Protein_chain_B.itp
      - topol_Protein_chain_C.itp
      - topol_Protein_chain_D.itp
    restraints:                     # optional, preparameterized mode only
      - posre_Protein_chain_A.itp
      - posre_Protein_chain_B.itp
      - posre_Protein_chain_C.itp
      - posre_Protein_chain_D.itp

  components:
    - id: A1
      topology: A1.itp
      coordinates: A1.gro
      role: competitive_ligand      # metadata only

    - id: COA
      topology: COA_original.itp
      coordinates: COA.gro
      role: cofactor

  forcefield: oplsaa
  water_model: spce
  box:                              # optional, currently metadata only
    distance: 1.2
    type: triclinic
```

### Raw vs. pre-parameterized protein

* **Raw mode** — `structure:` is a `.pdb`/`.ent`. SimForge runs `gmx pdb2gmx`
  (with the adaptive `-ignh` fallback). Do **not** give `topology:` /
  `restraints:` — pdb2gmx generates them.
* **Pre-parameterized mode** — `structure:` is a `.gro` that has *already* been
  through pdb2gmx. You supply the per-chain `topol_Protein_chain_*.itp` files
  directly. pdb2gmx is **never** invoked. This is a first-class input, not a
  workaround: re-deriving a raw PDB from a processed `.gro` does not reliably
  round-trip through pdb2gmx (post-processing renames atoms the force field's
  `.rtp` templates no longer recognize).

### Auto-discovery

```yaml
protein:
  structure: protein_only.gro
  topology: auto
  restraints: auto
```

`auto` discovers sibling files next to `structure`:

* `topol_Protein_chain_*.itp` (preferred), else a single `topol_Protein.itp`
* `posre_Protein_chain_*.itp` (or `posre_Protein.itp`)

Auto-discovery is **only** used when it is deterministic. If chain files and a
generic `topol_Protein.itp` both exist, the build fails and lists the
candidates — it never guesses. Every auto-discovery decision is recorded in
`provenance.json` under `discovery_notes`.

### parameterization topology vs. system-specific coordinates

A component's `.itp` is its **parameterization** — bonded/nonbonded terms and
its own `[ atomtypes ]` block, produced once by LigParGen. Its `.gro` is the
**system-specific pose** for *this* system. SimForge copies each component's
coordinates verbatim (pose resolution is a separate, upstream concern — see
`ligand.pose_rewriter`) and merges topologies with atomtype namespace
isolation.

> **LigParGen components must keep their original `[ atomtypes ]` block.**
> When two independently-run LigParGen submissions both define, say,
> `opls_800` with different parameters, SimForge namespaces the conflicting
> type per component (`A1_opls_800`, `COA_opls_800`). A component that ships
> *no* `[ atomtypes ]` section but references local `opls_*` names cannot be
> namespaced — the build fails in **preflight**, before assembly:
>
> ```
> [A1] A1.itp references atomtype(s) ['opls_800', ...] but contains NO
> [ atomtypes ] section. Re-export the component from LigParGen with its
> [ atomtypes ] definitions included.
> ```

---

## Invalid / valid combinations

| | |
|---|---|
| ❌ | protein `.pdb` **and** `topology:` set |
| ❌ | protein `.gro` **without** `topology:` (and not `auto`) |
| ❌ | component with `topology:` but no `coordinates:` |
| ❌ | a referenced file that does not exist |
| ❌ | `topology: auto` with ambiguous candidates |
| ✅ | raw protein `.pdb` + one component |
| ✅ | pre-parameterized protein `.gro` + ligand + cofactor |
| ✅ | pre-parameterized 4-chain protein + N components |

---

## Outputs

```
<out>/
  topol.top                  master topology (correct include order)
  complex.gro                merged protein + component coordinates
  protein.gro                protein coordinates (copied verbatim in preparam mode)
  component_atomtypes.itp    merged / deduplicated / namespaced atomtypes
  <component>.itp            each component, cleaned + renamed
  topol_Protein_chain_*.itp  copied verbatim (preparam mode)
  assembly_report.yaml       assembly-engine report
  provenance.json            full reproducibility record (see below)
  validation_report.json     structured validation checks
  system_build_report.json   human-oriented build summary
```

### `provenance.json`

```
schema, generated_utc
environment          simforge/python/platform/conda_env versions,
                     git commit + branch + dirty flag,
                     gromacs / rdkit / openbabel versions, timestamp
command              the exact argv
spec_path            absolute path to the YAML
system_spec          normalized specification
forcefield, water_model, protein_mode, pdb2gmx_used, pdb2gmx_ignh_fallback
discovered_protein_topology, discovery_notes
protein_molecules, component_molecule_names
atomtype_renames     {component: {old_type: new_type}}
total_atom_count, warnings
inputs / outputs     {label: {path, sha256, bytes}}
validation           readiness + every check
```

---

## Validation (`simforge validate-system`)

Structured, read-only checks over a built system directory:

| check | severity |
|---|---|
| `topology_present` | blocking |
| `coordinate_present` | blocking |
| `includes_resolve` — every local `#include` on disk | blocking |
| `molecule_table` — `[ molecules ]` names all have a `[ moleculetype ]` | blocking |
| `atomtype_resolution` — every referenced atomtype defined locally or by the force field | blocking (advisory if GROMACS force-field data is not inspectable) |
| `atom_counts` — Σ(mol atoms × count) == coordinate atom count | blocking |
| `coordinate_validity` — no NaN/inf, sane box, exact-overlap report | blocking / advisory |
| `grompp_dry_run` — `gmx grompp -maxwarn 0` (when `gmx` is in PATH) | blocking |

Exit code is non-zero when any **blocking** check fails.
`validation_report.json` is written next to the system.
