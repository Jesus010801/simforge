# Contributing to SimForge

SimForge is a research-preview tool for declarative construction, validation, and
provenance of GROMACS molecular-simulation systems. Contributions are welcome —
please read this short guide first.

## Development setup

Requires Python >= 3.11. Either:

```bash
# Option A: conda (recommended if you want the optional chem/analysis stack)
conda env create -f environment.yml
conda activate simforge

# Option B: pip, editable install with all extras
pip install -e ".[dev,chem,analysis]"
```

RDKit is an **optional** dependency. SimForge has a documented degraded path that
runs without it; install the `chem` extra (or a conda `rdkit` build) to exercise
the full ligand layer.

## Running tests

```bash
pytest                      # full suite (~2400 tests, a few minutes)
pytest -m "not gmx"         # skip tests that need a real GROMACS binary
pytest -m "not rdkit"       # skip tests that need RDKit
```

Marker reference:

- `gmx` — needs a real `gmx` (GROMACS) binary on `PATH`.
- `rdkit` — needs RDKit importable.
- `membrane_wip` — known-incomplete membrane packing; expected to fail/xfail
  until the membrane pipeline leaves beta.

CI runs `pytest -m "not gmx and not rdkit"` on Python 3.11–3.13, plus an advisory
real-GROMACS smoke job.

## Coding expectations

- Match the style of the surrounding code. Use type hints on new functions.
- No new hard runtime dependencies without discussion first (open an issue).
  Optional features belong behind an extra and a graceful import guard.
- Keep changes focused; add tests for new behaviour.

## Scientific-regression requirement

Any change that touches **topology assembly, atomtype handling, or coordinate
merging** MUST add or update a test that runs a real
`gmx grompp -maxwarn 0` on a small fixture and asserts it succeeds (and, where
relevant, checks atom counts / total charge / box). Changes that alter what
SimForge hands to GROMACS are not accepted without this coverage.

## Reporting a scientific-correctness concern

If you think SimForge produced a scientifically wrong topology, system, or
analysis, open a **"Scientific concern"** issue rather than a plain bug report.
Attach the spec YAML, `provenance.json`, `validation_report.json`, the full
`gmx grompp` output, and a clear expected-vs-got description.
