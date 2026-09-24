 Phase 1 — Pre-Implementation Regression Baseline

**Purpose:** establish the repository test state before any functional SSKI Phase 1 code is introduced.

## Focused regression baseline

The Phase 0 audit-defined focused regression suite completed with:

770 passed  
2 skipped  
5 warnings  
0 failed

This focused baseline covers existing parsing, decision, inference, structural annotation,
MD knowledge, compiler/workflow, system contracts, runtime/recovery, ligand identity,
provenance, topology, CLI, and campaign behavior.

## Repository-wide non-GROMACS / non-RDKit baseline

Command:

`pytest -m "not gmx and not rdkit" -q -p no:cacheprovider`

Result before Phase 1 implementation:

4 failed  
3062 passed  
6 skipped  
4 deselected  
11 xfailed  
1 xpassed  
68 warnings

## Pre-existing failures

All four failures are confined to `ligand/test_integrate.py`:

1. `TestPrepareCLI::test_prepare_with_ambiguous_h_ratio_blocks_without_rdkit`
2. `TestHydrogenationCLIOption::test_hydrogenation_none_unknown_status_exits_nonzero`
3. `TestHydrationStatus::test_status_unknown_for_plausible_ratio`
4. `TestHydrationStatus::test_complex_with_h_status_is_unknown`

Observed behavior:

- expected hydrogenation status: `unknown`
- actual hydrogenation status: `complete`

The associated CLI tests therefore exit successfully where the current tests expect a
non-zero result.

These failures existed before Phase 1 SSKI implementation and are outside the SSKI scope.

## Phase 1 regression rule

Phase 1 MUST NOT introduce additional regressions.

A successful Phase 1 validation therefore requires:

- all new SSKI Phase 1 tests to pass;
- the focused existing regression suite to remain at least as healthy as this baseline;
- the repository-wide gate to introduce no new failures;
- the four failures listed above to be reported separately as pre-existing unless fixed independently outside the SSKI scope.

No attempt to repair these ligand failures is authorized as part of Phase 1.
