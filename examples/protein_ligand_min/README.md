# protein_ligand_min

Minimal preparameterized protein + one component. Software demonstration.

```bash
simforge build examples/protein_ligand_min/system.yaml --out example_system
simforge validate-system example_system
simforge inspect-run  example_system
```

## Inputs (all shipped)

| File | What it is |
|---|---|
| `protein.gro` | 6-residue peptide (SER-GLY-PHE-LEU-ARG-ASP), already through `gmx pdb2gmx` (OPLS-AA / SPC-E), 96 atoms |
| `topol_Protein_chain_A.itp` | the chain topology pdb2gmx wrote (discovered via `topology: auto`) |
| `posre_Protein_chain_A.itp` | position restraints (discovered via `restraints: auto`) |
| `MTH.itp` | methane, as a stand-in ligand; self-contained `[ atomtypes ]` |
| `MTH.gro` | a pose for MTH in this system |
| `system.yaml` | the declarative spec |

## Expected result

`example_system/` with `topol.top`, `complex.gro` (101 atoms), the merged
`component_atomtypes.itp`, `provenance.json`, and `validation_report.json`
showing all 8 checks `PASS`, including `gmx grompp -maxwarn 0`.

This is **not solvated** — it is a preprocessed system ready for you to run
`gmx editconf` / `solvate` / `genion` and your own equilibration protocol.
