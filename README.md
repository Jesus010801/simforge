# SimForge

**Reproducible construction, validation and provenance for GROMACS
molecular-simulation systems.**

> **Status: research preview / alpha (`0.1.0a1`).**
> The protein–ligand and multi-component system-construction layers are
> functional and covered by real-`gmx grompp` tests. The membrane workflow is
> **beta** and requires manual QC. Interfaces may change. This is **not yet a
> stable scientific release** — see [Known limitations](#known-limitations).

SimForge turns a single declarative YAML *system specification* into a
GROMACS-ready system — merged topology, merged coordinates, a structured
validation report, and a `provenance.json` that records every input hash, tool
version, and automated decision. It is a thin front end over an assembly engine
that handles the parts most likely to go **silently wrong** when you build a
multi-component system by hand.

---

## What SimForge does

- **`simforge build system.yaml`** — compile a `SystemSpec` into
  `topol.top` + `complex.gro` + `provenance.json` + `validation_report.json`.
- **`simforge validate-system DIR`** — 8 structured, read-only checks on a
  built system, ending with a blocking `gmx grompp -maxwarn 0`. Non-zero exit
  on any blocking failure.
- **`simforge inspect-run DIR`** — render the provenance / reproducibility
  manifest for a built system.
- **`simforge doctor`** — check the runtime environment (GROMACS, force fields,
  optional RDKit / Open Babel, Python, conda env).
- **`simforge campaign-build campaign.yaml`** — build many systems that share
  one parameterized protein.
- **`simforge ligand …`** — prepare a ligand from a docked complex, assess
  protonation, export for LigParGen, assemble.
- **`simforge fel …`** — 2D free-energy-landscape from two GROMACS XVG
  observables (user-chosen collective variables).
- **`simforge study …`** / **`simforge analyze …`** — trajectory-first
  comparative MD analysis over a directory of finished runs.

The original workflow compiler (`simforge compile` / `run` / `status`, YAML →
DAG → GROMACS scripts) is still present but is not the focus of this preview.

## Key differentiators (vs. a hand-written shell script)

| | |
|---|---|
| **Atomtype namespace isolation** | Two independently parameterized components (e.g. two LigParGen exports) that both define `opls_800` for *different* atoms are namespaced per component (`A_opls_800` / `B_opls_800`); the rename map is recorded. A naive `cat A.itp B.itp` silently corrupts one molecule. |
| **Deterministic preflight** | A component that references an unresolvable atomtype fails **before** assembly, with the component, the missing types, and the fix. |
| **Pre-parameterized protein as a first-class input** | An already-`pdb2gmx`'d protein is reused verbatim — SimForge never round-trips a processed `.gro` back through `pdb2gmx`. |
| **`provenance.json`** | Schema-versioned: SHA256 of every input/output, GROMACS / RDKit / Open Babel versions, git commit + dirty flag, the exact argv, the atomtype rename map, the embedded validation report. |
| **`validate-system` as a standalone gate** | Re-runnable, exit-code-meaningful; drop it in a Makefile or CI. |
| **`doctor`** | One command turns "why doesn't this work on the cluster" into a table. |

## What SimForge does **not** do

- It does **not** parameterize ligands. LigParGen (or your own parameters) is
  an external, manual step; SimForge prepares the input and validates the
  output.
- It does **not** solvate or add ions. `build` stops at the assembled,
  pre-solvation system; you run `gmx solvate` / `genion` and your own
  equilibration protocol afterwards (the build output prints the next command).
- It does **not** choose collective variables for the FEL — you do.
- It does **not** yet produce production-ready membrane systems unattended.

## Scientific scope and assumptions

- Validated primarily against **GROMACS 2025.2**. The conda recipe tracks
  recent releases; other versions are untested.
- **OPLS-AA + SPC/E** is the primary tested force-field / water-model scope
  (defaults throughout).
- Component `.itp` files must keep their own `[ atomtypes ]` block (this is how
  collisions are namespaced) — the preflight enforces it.
- The membrane builder follows the Berger-lipid / InflateGRO lineage and is
  **beta**: visual + metric QC (APL, thickness, tilt, pore hydration) is
  **mandatory** before production MD.
- FEL surfaces are descriptive (Boltzmann inversion over user-chosen CVs), not
  predictive.

---

## Installation

Requires **Python ≥ 3.11** and a **GROMACS** install with `gmx` on `PATH`
(plus the `oplsaa.ff` force field bundled with GROMACS).

```bash
git clone https://github.com/Jesus010801/simforge
cd simforge

# conda (recommended — also gets GROMACS, numpy, matplotlib, optional RDKit)
conda env create -f environment.yml
conda activate simforge
pip install -e .

# or pip only (you provide GROMACS yourself)
pip install -e ".[dev]"
```

Optional extras: `pip install -e ".[chem]"` (RDKit), `".[analysis]"`
(matplotlib for FEL plots). SimForge runs without RDKit on a documented
degraded path.

### Verify the environment

```bash
simforge doctor
```

```
   Component     Level      Version   Status
   SimForge      REQUIRED   0.1.0a1   ✓ PASS
   Python        REQUIRED   3.13.9    ✓ PASS
   GROMACS       REQUIRED   2025.2    ✓ PASS
   OPLSAA        REQUIRED   —         ✓ PASS
   RDKit         OPTIONAL   —         ✗ MISSING   (features degrade to heuristics)
   Open Babel    OPTIONAL   3.1.1     ✓ AVAILABLE

  Overall status: READY
```

---

## Quick start

```bash
simforge doctor

simforge build examples/multicomponent/system.yaml --out multi_system

simforge validate-system multi_system

simforge inspect-run  multi_system
```

`examples/multicomponent/` ships every input (a 6-residue peptide + two
independently parameterized small molecules with a deliberate `opls_800`
collision). The build assembles the topology, namespaces the collision, and
passes `gmx grompp -maxwarn 0`; `inspect-run` shows the rename map in the
provenance. See [`examples/`](examples/) for details — these are software
demonstrations, not scientific benchmarks.

### Example specification

```yaml
system:
  name: multicomponent_demo

  protein:
    structure: protein.gro       # .gro -> preparameterized (no pdb2gmx)
    topology: auto               # discover sibling topol_Protein_chain_*.itp
    restraints: auto

  components:
    - id: MTH
      topology: MTH.itp          # keeps its own [ atomtypes ]
      coordinates: MTH.gro       # system-specific pose
      role: ligand
    - id: EOL
      topology: EOL.itp
      coordinates: EOL.gro
      role: cofactor

  forcefield: oplsaa
  water_model: spce
```

A raw `.pdb` under `structure:` instead switches to raw mode (SimForge runs
`gmx pdb2gmx`); then you do **not** supply `topology:` / `restraints:`.

### Outputs

```
multi_system/
  topol.top                  master topology (correct include order)
  complex.gro                merged protein + component coordinates
  component_atomtypes.itp    merged / deduplicated / namespaced atomtypes
  provenance.json            full reproducibility record
  validation_report.json     the 8 structured checks
  system_build_report.json   human-oriented build summary
```

---

## Documentation

- [`docs/system_build.md`](docs/system_build.md) — the `SystemSpec`, raw vs.
  pre-parameterized mode, auto-discovery, the validation checks, the provenance
  schema.
- [`docs/study_campaign.md`](docs/study_campaign.md) — the MD study / campaign
  analysis layer.
- [`examples/`](examples/) — runnable systems.
- [`CHANGELOG.md`](CHANGELOG.md) · [`CONTRIBUTING.md`](CONTRIBUTING.md)

## Known limitations

- **Membrane generation is beta.** 9 pre-shrink lipid-exclusion tests are
  marked `xfail` / `membrane_wip` (synthetic-fixture contract issues, not
  regressions — real-data regression tests pass). Manual QC is required before
  production MD.
- Only **GROMACS 2025.2** is actively validated.
- The full end-to-end scientific milestone (a 4-chain protein + inhibitor +
  cofactor taken through assembly → grompp → solvation → ions → NVT without
  manual topology repair) has been achieved **manually**; the repository
  protects the *topology mechanism* with a small real-`grompp` regression
  (`tests/test_examples.py`), not the full 24k-atom system.
- `cli.py` is a large single module (known refactor debt).
- Some `pydantic` v1-style config warnings remain.

## Testing

```bash
pytest                       # ~2400 tests, a few minutes
pytest -m "not gmx"          # skip tests needing a real GROMACS binary
```

CI (GitHub Actions) runs the non-GROMACS suite on Python 3.11–3.13 plus an
advisory real-GROMACS smoke job.

## Citing

If you use SimForge, please cite it — see [`CITATION.cff`](CITATION.cff).
A Zenodo DOI will be minted for the first tagged release.

## License

MIT — see [`LICENSE`](LICENSE). Some bundled membrane-setup helpers under
`docs/Prot-Memb_FILES/` are third-party community tools with their own terms;
see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
