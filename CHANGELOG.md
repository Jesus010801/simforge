# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0a1] - 2026-09-03

First public research-preview release.

### Added

- Declarative `SystemSpec` model and `simforge build` to compile a single
  specification into GROMACS topologies and coordinates.
- `simforge doctor` environment diagnostics (Python, GROMACS, optional
  dependencies, force-field availability).
- `simforge validate-system` structured validation with a blocking
  `gmx grompp -maxwarn 0` gate and a machine-readable `validation_report.json`.
- `simforge inspect-run` provenance summary for a completed run.
- `simforge campaign-build` for multi-system campaigns that share a single
  parameterized protein.
- `provenance.json` output with SHA256 hashing of every input and output file
  plus environment, version, and git capture.
- Support for N independently-parameterized non-protein components with
  deterministic atomtype namespace isolation.
- Deterministic atomtype preflight that surfaces conflicts before assembly.
- Pre-parameterized protein support: an already-parameterized protein is reused
  without re-running `pdb2gmx`.
- Adaptive `pdb2gmx -ignh` fallback for structures with missing/extra hydrogens.
- Ligand `prepare` / `integrate` / batch tooling.
- LigParGen export/import helpers for external ligand parameterization.
- FEL (free-energy-landscape) analysis built from XVG observables.
- Trajectory-first MD study and campaign analysis.

### Known limitations

- Membrane system generation is **beta** and requires manual QC; 9
  known-incomplete pre-shrink-exclusion tests are marked `xfail` /
  `membrane_wip`.
- Validated primarily against GROMACS 2025.2.
- OPLS-AA with SPC/E water is the primary tested force-field / water-model
  scope.
- LigParGen parameterization is an external, manual step outside the pipeline.

[Unreleased]: https://github.com/Jesus010801/simforge/compare/v0.1.0a1...HEAD
[0.1.0a1]: https://github.com/Jesus010801/simforge/releases/tag/v0.1.0a1
