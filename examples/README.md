# SimForge examples

Small, fully self-contained systems that ship every input needed to run them
from a fresh clone. They are **software demonstrations**, not curated
scientific benchmarks — the "protein" is a 6-residue peptide fragment and the
"ligands" are methane / methanol. The point is to exercise the
build → validate → provenance path end to end.

| Example | What it shows |
|---|---|
| [`protein_ligand_min/`](protein_ligand_min/) | Minimal preparameterized protein + one component. `build`, `validate-system`, `inspect-run`. |
| [`multicomponent/`](multicomponent/) | Same protein + **two** independently parameterized components whose ITPs both declare `opls_800` for different elements. SimForge namespaces the collision (`MTH_opls_800` / `EOL_opls_800`) and records the map in `provenance.json`. |

## Run it

```bash
simforge doctor                                   # check GROMACS / force field

simforge build examples/multicomponent/system.yaml --out multi_system
simforge validate-system multi_system             # 8 checks incl. gmx grompp -maxwarn 0
simforge inspect-run  multi_system                # provenance / reproducibility manifest
```

`multi_system/` will contain `topol.top`, `complex.gro`, `component_atomtypes.itp`
(merged + namespaced), `provenance.json`, `validation_report.json` and
`system_build_report.json`.

## Provenance of the inputs

- `protein.gro`, `topol_Protein_chain_A.itp`, `posre_Protein_chain_A.itp` —
  produced by `gmx pdb2gmx` (OPLS-AA / SPC-E) from a 6-residue fragment
  (SER-GLY-PHE-LEU-ARG-ASP). No external structure is redistributed.
- `MTH.itp` / `EOL.itp` — hand-written minimal ITPs using published OPLS-AA
  non-bonded values and OPLS/CHARMM22 bonded values; each carries its own
  `[ atomtypes ]` block, deliberately shaped like an independent LigParGen
  export. See the header comment in each file.
