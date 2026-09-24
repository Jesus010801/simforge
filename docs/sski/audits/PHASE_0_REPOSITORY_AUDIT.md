PHASE 0 — SSKI REPOSITORY ARCHITECTURE AUDIT

Sí mi lord.

**Verdict:** `simforge/knowledge/` is the correct location. Phase 1 can proceed to domain design without restructuring SimForge. **No unresolved BLOCKER was identified for that limited phase.** Existing execution, artifact, identity, and evidence models must not become SSKI domain dependencies.

All requested normative documents, all ten accepted ADRs, and the Phase 0 prompt were read. This audit used file inspection and searches only. No files were modified, no tests or project code were executed, and no commits were created. Initial and final working-tree status listings matched. Existing uncommitted work was treated as user-owned context, not SSKI work.

## A. Repository topology

### A1 — MINOR: SimForge is a multi-package distribution

[pyproject.toml](/home/jesusxd/Escritorio/simforge/pyproject.toml) defines:

| Property | Actual configuration |
|---|---|
| Distribution | `simforge`, version `0.1.0a1` |
| Build backend | `setuptools.build_meta` |
| Python | `>=3.11` |
| Console command | `simforge = "simforge.cli:app"` |
| Standalone module | Root `cli.py`, explicitly included through `py-modules` |
| Package discovery | Repository root, with explicit package-name patterns |
| Required dependencies | Typer, Rich, Pydantic 2+, PyYAML, NumPy |
| Optional dependencies | RDKit; Matplotlib; SciPy; development/testing dependencies |

These directories are actual top-level Python packages with `__init__.py`:

```text
simforge/
core/
builders/
executors/
runtime/
pipelines/
descriptors/
validators/
workflows/
adapters/
ligand/
utils/
analysis/
```

All are included by the packaging configuration. `benchmarks/` is also a Python package, but is not included by those discovery patterns.

Relevant nested packages include:

```text
simforge/knowledge/
core/md_knowledge/
builders/step_builders/
validators/ligand_parsers/
validators/membrane_water/
validators/membrane_water_v2/
analysis/md/
analysis/fel/
analysis/campaign/
analysis/campaign/observables/
analysis/campaign/orchestration/
analysis/campaign/structure/
analysis/campaign/trajectory/
```

`tests/knowledge/` currently contains a README, not implemented SSKI tests. `simforge/knowledge/` contains its instructions and a bootstrap `__init__.py`, not production domain models.

### A2 — MINOR: The apparent duplication is an entry-point arrangement

There are not two equivalent implementations of the application:

- `core/` contains established application models, parsing, inference, planning, and compilation.
- `simforge/__init__.py` currently defines the package version.
- Root `cli.py` implements the CLI.
- `simforge/cli.py` imports `cli.cli` as `app`.
- `simforge/__main__.py` invokes that application.

The entry-point chain is:

```text
simforge console command / python -m simforge
    → simforge.cli
    → root cli.py
    → existing top-level application packages
```

Imports such as `from core.models import SystemState` are absolute imports of independently packaged top-level modules. They are not shorthand for `simforge.core`.

Source-tree execution resolves them when the repository root is on Python’s import path. Installed distributions include these packages and the root CLI module explicitly. This audit verified the configuration, not a built wheel.

Some generated worker scripts additionally insert the repository root into `sys.path`, notably in assembly, embedding, and match-box builders. That is existing execution behavior, not a packaging pattern SSKI should copy.

**Consequence:** no package migration, CLI relocation, or import rewrite is necessary for Phase 1.

## B. Existing concepts relevant to SSKI

### B1 — MAJOR: Existing system models describe simulation inputs, not canonical scientific entities

[core/models.py](/home/jesusxd/Escritorio/simforge/core/models.py:214) defines `ComponentModel` with a local string ID, role, file, biological context, validation, descriptors, and reasoning. `SystemState` aggregates simulation configuration and accumulated inference.

[core/system_spec.py](/home/jesusxd/Escritorio/simforge/core/system_spec.py:81) defines `ProteinSpec` as coordinate/topology/restraint inputs. `ComponentSpec` and `SystemSpec` describe assembly inputs, force fields, and water models.

These are valuable downstream contracts. They are not replacements for:

- Provider-independent `EntityId`.
- Biological `Protein`.
- Exact `ProteinSequence`.
- Versioned external mappings.
- Separate chemical identity and species.

**Required treatment:** preserve existing models; later translate resolved SSKI context into them on the SimForge side.

### B2 — MAJOR: Existing decision-making is simulation planning

[core/decision_engine.py](/home/jesusxd/Escritorio/simforge/core/decision_engine.py) builds `SimulationPlan`, workflow policy, steps, and engine parameters.

[core/scientific_planner.py](/home/jesusxd/Escritorio/simforge/core/scientific_planner.py) produces interactive planning questions and applies answers as patches to `SystemState`.

[core/semantic_inference.py](/home/jesusxd/Escritorio/simforge/core/semantic_inference.py) normalizes objectives, expands presets, and infers membrane-related simulation hints.

These functions must remain owned by SimForge. SSKI policy resolution must not duplicate them or promote their inferred hints into externally established scientific facts.

### B3 — MAJOR: Existing evidence lacks SSKI’s full provenance semantics

[core/md_knowledge/evidence.py](/home/jesusxd/Escritorio/simforge/core/md_knowledge/evidence.py:12) defines `Evidence` as an observable value, temporal pattern, simulation-state vote, confidence, and explanatory message. `EvidenceBundle` accumulates votes.

[core/structural_annotation.py](/home/jesusxd/Escritorio/simforge/core/structural_annotation.py:139) defines `OrientationEvidence` with source, confidence, reference, and notes. Its residue annotations use local residue-range representations.

These do not provide the complete SSKI chain of source release, raw artifact hash, importer version, claim references, and decision policy version. Local residue ranges also do not establish mappings among sequence, author numbering, prepared coordinates, and topology indices.

Campaign `ComponentEvidence` similarly supports component classification with weighted observations. It is not a canonical claim/evidence graph.

### B4 — MAJOR: Existing artifact tracking is deliberately mutable and workspace-oriented

[runtime/artifacts.py](/home/jesusxd/Escritorio/simforge/runtime/artifacts.py:25) defines a mutable `ArtifactRef` containing:

```text
path, checksum, semantic_role, step_id, created_at, size_bytes
```

Its registry indexes lineages by absolute path and persists workspace metadata. `checksum()` returns an empty string for a missing file. That behavior is explicitly covered by an existing regression test.

[core/semantic_artifacts.py](/home/jesusxd/Escritorio/simforge/core/semantic_artifacts.py) includes topology mutation tracking and engine-format-specific artifacts.

These are appropriate execution concepts but cannot directly satisfy immutable, hash-validated bundle artifact semantics.

**Do not “fix” the existing runtime artifact contract as part of SSKI.** Introduce a separate qualified SSKI `ArtifactRef`.

### B5 — MAJOR: Existing provenance and manifests have different purposes

Existing implementations include:

- `core/provenance.py`: assembly/toolchain provenance, including GROMACS and optional chemistry tool versions.
- `core/workspace_fingerprint.py`: builder/build/template freshness.
- `builders/workspace_builder.py`: execution manifest with workflow steps.
- `core/variant_compiler.py`: comparative variant/workspace manifest.
- `runtime/trajectory_ingestor.py`: discovered trajectory files.
- `ligand/campaign.py`: campaign preparation manifests.
- `analysis/campaign/manifest.py`: editable study discovery and assignment manifest.
- `analysis/md/provenance.py`, `analysis/fel/provenance.py`, `analysis/campaign/provenance.py`: analysis provenance.
- `runtime/journal.py`: execution events.

None is already a `KnowledgeBundleManifest`, `KnowledgeSnapshot`, or complete `RunDerivationLedger`.

Existing runtime records can supply information to a future derivation integration. They must not be renamed or reused as though their guarantees already match SSKI.

### B6 — MAJOR: Chemical identity already exists, but with narrower semantics

[ligand/campaign.py](/home/jesusxd/Escritorio/simforge/ligand/campaign.py:309) defines a frozen `ChemicalIdentity` derived from a reference file, with SMILES, InChIKey, heavy-atom count, and elements. Its computation uses RDKit.

`ligand.normalization.LigandIdentity` instead manages component IDs and GROMACS-safe residue/molecule names.

Neither establishes the SSKI separation:

```text
ChemicalIdentity → ChemicalSpecies → SimulationMicrostate
```

A future adapter must preserve the method and provenance of existing chemical observations. It must not equate a residue label, structural signature, or available InChIKey with an automatically confirmed SSKI identity.

## C. Reuse / Adapt / Isolate / New matrix

“ADAPT” means a future explicit boundary conversion, **not a Phase 1 implementation or import permission**.

| Concept / implementation | Classification | Decision |
|---|---|---|
| Existing package discovery and entry points | REUSE | Already accommodate `simforge.knowledge.domain`. |
| Pydantic 2 model mechanism | REUSE | Available convention; use independently of existing application models. |
| Standard-library SHA256 | REUSE | Use `hashlib`; do not import permissive runtime hashing helpers. |
| `ComponentModel`, `SystemState`, `SystemSpec` | ADAPT | Future downstream consumption of resolved scientific inputs. |
| `ProteinSpec`, `ComponentSpec` | ADAPT | Assembly-facing representations, not entity definitions. |
| Structural annotations / `OrientationEvidence` | ADAPT | Preserve provenance and explicit coordinate mappings during conversion. |
| `ligand.campaign.ChemicalIdentity` | ADAPT | Reference-derived chemical observation, requiring validated mapping. |
| `LigandIdentity` | ISOLATE | Simulation naming and parameterization contract. |
| Campaign canonical system IDs | ISOLATE | Study/system identity, not global scientific entity identity. |
| Decision engine and `WorkflowPolicy` | ISOLATE | Simulation planning and execution parameters. |
| Semantic inference and objective presets | ISOLATE | Application intent normalization. |
| Scientific planner / planning session | ISOLATE | Interactive simulation decisions; do not clone into SSKI. |
| MD `Evidence`, heuristics, state voting | ISOLATE | Simulation interpretation; later explicit ingestion could produce evidence. |
| Runtime artifact registry / semantic artifacts | ADAPT | Future derived-output linkage; never the immutable bundle registry. |
| Build and analysis provenance | ADAPT | Potential derivation inputs with explicit schema conversion. |
| Existing execution/study/campaign manifests | ISOLATE | Retain their existing lifecycles and meanings. |
| Workflow graph, compiler, executors, builders | ISOLATE | Existing SimForge execution architecture. |
| SSKI identities and scientific entity models | NEW | Provider-independent scientific vocabulary. |
| SSKI `Claim` and `Evidence` | NEW | Separate assertions and supporting records. |
| SSKI `ArtifactRef` | NEW | Strict digest-bearing scientific reference. |
| `PolicyRef`, `Decision`, `DecisionContext` | NEW | Version-pinned, four-state scientific resolution records. |
| `KnowledgeBundleManifest` | NEW | Immutable domain representation in Phase 1. |
| Actual bundle, snapshot machinery, derivation ledger | NEW, later | Do not implement their infrastructure in Phase 1. |

## D. Naming and semantic collisions

Searches covered class declarations and typed-definition patterns, plus related concepts in the inspected Python packages.

### D1 — MAJOR: Exact-name collisions require qualified imports

| Proposed name | Existing name | Semantic difference |
|---|---|---|
| `ArtifactRef` | `runtime.artifacts.ArtifactRef` | Mutable workspace/path reference versus strict scientific artifact reference. |
| `Evidence` | `core.md_knowledge.evidence.Evidence` | Observable/state-voting evidence versus source-provenanced claim support. |
| `ChemicalIdentity` | `ligand.campaign.ChemicalIdentity` | File-derived chemistry summary versus canonical scientific entity. |

Keep both meanings in their own modules. Future bridge code should use explicit aliases such as `KnowledgeArtifactRef` and `RuntimeArtifactRef`. Do not introduce global compatibility aliases that make them interchangeable.

### D2 — MINOR: Related names must not be conflated

| Requested concept | Existing related concepts |
|---|---|
| Entity / EntityId | Component IDs, campaign IDs, ligand internal IDs |
| ExternalIdentifier | Source/reference strings, accessions in annotations |
| Artifact | `BaseArtifact`, `AnalysisArtifact`, `TrajectoryArtifact`, `TopologyState` |
| Claim | Annotation/inference assertions, without a dedicated SSKI claim record |
| Decision | `CleanupDecision`, planning answers, execution reasoning |
| DecisionContext | `SystemContext`, `AnalysisContext`, `SystemState` |
| Manifest | `VariantManifest`, `CampaignManifest`, `StudyManifest`, `TrajectoryManifest`, execution-manifest dictionaries |
| Snapshot | `runtime.metrics.SystemSnapshot`, source/build fingerprints |
| Provenance | `analysis.md.provenance.Provenance`, other analysis/build provenance records |
| Policy | `WorkflowPolicy`, membrane cleanup/contact policies |
| Protein | `ProteinSpec`, protein validation results, protein-role components |
| ChemicalSpecies | Local species/water classification concepts |

No exact declarations of `Entity`, `EntityId`, `ExternalIdentifier`, `Claim`, `DecisionContext`, `Protein`, or `ChemicalSpecies` were found in the inspected application code.

**Duplicate-domain-model risk:** manageable if SSKI owns scientific identity and provenance while existing models retain ownership of preparation, execution, and analysis. It becomes serious if developers attempt implicit conversion based on matching names or string IDs.

## E. Architectural boundary analysis

### E1 — MAJOR: The hard boundary is feasible but not implemented

The existing import topology permits:

```text
simforge.knowledge.domain
    ↓ future immutable contracts
KnowledgeBundle + DecisionContext
    ↓ future SimForge-side reader/conversion
existing preparation, planning, compilation and execution
```

`simforge/__init__.py` does not eagerly import the application execution stack, so importing the proposed domain package need not load it.

Phase 1 should introduce **no imports from existing SimForge application packages**. Later conversion belongs on the consuming side, not inside scientific domain models.

### E2 — MAJOR: Current `SimulationPlan` is not engine-independent

[core/execution_models.py](/home/jesusxd/Escritorio/simforge/core/execution_models.py:104) requires `SimulationStep.engine` and exposes an explicitly engine-specific `params` dictionary.

The decision engine emits:

- `gromacs:pdb2gmx`, `gromacs:solvate`, and `gromacs:genion`.
- GROMACS MD steps.
- Parameters such as `integrator`, `nsteps`, `constraints`, and output frequencies.

This is a real gap against ADR-010’s target architecture. It does **not** prevent pure domain models from starting safely, because Phase 1 neither consumes nor changes `SimulationPlan`.

Do not claim the existing plan already satisfies ADR-010. Resolve that integration gap in a separately approved later scope.

### E3 — MAJOR: Convenient reuse would create inappropriate coupling

Particularly unsuitable domain imports are:

- `core.semantic_artifacts` → `runtime.artifacts`.
- `core.compiler` → pipelines and workflow compilation.
- `core.provenance` → tool probing, subprocesses, and execution environment.
- Ligand identity computation → RDKit.
- `core.execution_models` → execution status and engine parameter contracts.
- CLI modules → broad application command registration.

The targeted network search found documentation URLs rather than active HTTP clients in the inspected core/execution paths. That is **not proof** of worker isolation: subprocesses and generated scripts require separate later runtime enforcement.

### E4 — MAJOR: Ordinary pytest execution cannot prove domain import isolation

Root [conftest.py](/home/jesusxd/Escritorio/simforge/conftest.py) eagerly imports the parser, decision engine, compiler, workspace builder, and shell executor.

Consequently, a test asserting that the current pytest process contains no execution modules would fail for unrelated reasons—or provide misleading evidence if inadequately written.

Use static dependency checks and a fresh child interpreter to test the SSKI import boundary. Do not change root `conftest.py` in Phase 1.

## F. SSKI package-location decision

### F1 — MINOR: Approve `simforge/knowledge/domain/`

This path:

- Matches both repository instruction files and the implementation plan.
- Falls under the existing `simforge*` discovery rule.
- Preserves established `core.*` imports.
- Separates source knowledge from `core/md_knowledge/`, which interprets simulation observables.
- Requires no packaging or CLI change.

Reject `core/knowledge/`, `core/md_knowledge/`, root `knowledge/`, and a new `src/` layout for Phase 1. They would introduce unnecessary divergence or restructuring.

## G. Required regression tests

### G1 — MINOR: Use the existing regression suites, without editing them

The focused Phase 1 regression set should include these exact existing paths:

| Coverage | Required paths |
|---|---|
| Parsing, inference, decisions | `core/test_parser.py`, `core/test_decision_engine.py`, `core/test_semantic_inference.py`, `core/test_structural_annotation.py`, `core/test_md_knowledge.py` |
| Compilation and workflow contracts | `core/test_compiler.py`, `workflows/test_workflow_graph.py`, `builders/test_workspace_builder.py`, `builders/test_workspace_selfcontained.py` |
| System input/build contracts | `core/test_system_spec.py`, `core/test_system_build.py`, `core/test_system_validate.py` |
| Runtime artifacts and recovery | `runtime/test_runtime.py`, `runtime/test_checkpoint_recovery.py`, `executors/test_executor.py`, `executors/test_dag_blocking.py`, `tests/test_workflow_resume.py` |
| Chemical-identity overlap | `ligand/test_campaign.py`, `ligand/test_normalization.py`, `ligand/test_chemical_perception.py` |
| Provenance and topology | `tests/test_pdb2gmx_provenance.py`, `tests/test_topology_assembly.py` |
| Application imports and CLI | `test_cli_smoke.py`, `test_cli_platform.py`, `test_cli_ligand.py`, `test_cli_integrate.py` |
| Campaign overlap | `tests/analysis/campaign/test_manifest.py`, `tests/analysis/campaign/test_canonical_ids.py`, `tests/analysis/campaign/test_fingerprint.py` |

The existing repository-wide CI regression command is:

```bash
pytest -m "not gmx and not rdkit" -q -p no:cacheprovider
```

Run that gate after the focused tests. It also exercises user-owned campaign and membrane work without authorizing edits to those areas. Compare failures with a pre-implementation baseline; do not attribute every existing failure to SSKI.

The existing GROMACS smoke job remains:

```bash
pytest -m gmx -q -p no:cacheprovider -k "build or validate or integrate or spec"
```

Preserve that CI job. Report unavailable optional environments rather than claiming those checks passed.

No dedicated scientific-planner behavioral suite was found; the smoke suite checks its importability. Adding unrelated planner tests is outside Phase 1.

### G2 — MAJOR: Required new architecture and invariant tests

| Requirement | Proposed verification |
|---|---|
| No engine dependency | AST/import-graph checks and clean-process imports reject application execution packages and engine bindings. |
| No network dependency | Reject HTTP libraries and standard-library network imports; guard network access in the clean process. |
| No persistence dependency | Reject SQLAlchemy, database drivers, SQLite, DuckDB, and persistence operations. |
| Internal identity | Reject external accessions as entity IDs; verify external identifiers link to separately supplied internal IDs. |
| Claim ≠ Evidence | Separate types, distinct serialization, explicit claim references and evidence provenance. |
| Chemical distinctions | Identity and species remain separate; distinct species may reference the same parent identity. |
| Experimental structure semantics | Predicted records cannot silently validate as experimental structures. |
| SHA256 semantics | Reject missing, empty, wrong-length, and nonhexadecimal digests; define case normalization explicitly. |
| Deterministic serialization | Equivalent supported values produce identical canonical output, including across fresh processes. |
| Deep immutability | Reject field assignment and nested mutation; changing caller-owned collections cannot alter stored contracts. |
| Decision states | Preserve all four states; `DEFERRED` carries explicit operation requirements rather than meaning “unknown.” |
| Namespace coexistence | Existing and SSKI collision types coexist without aliasing or import side effects. |

**Phase limitation:** Phase 1 can prove immutable manifest/context value semantics. It cannot prove an actual filesystem `KnowledgeBundle` is immutable, complete, or hash-verified before execution. Those tests belong to bundle/compiler/runtime phases. A frozen manifest or `immutable: true` field alone is not that proof.

## H. Specification conflicts

| Severity | Finding | Required disposition |
|---|---|---|
| MAJOR | ADR-010 expects an engine-independent plan; the current plan carries engines and backend parameters. | Record the gap; defer integration redesign. Keep domain models independent. |
| MAJOR | Existing runtime artifacts permit mutation and missing-file empty hashes, unlike bundle requirements. | Use separate SSKI contracts. Preserve existing runtime behavior. |
| MAJOR | Existing evidence and annotations do not satisfy SSKI provenance and coordinate-mapping requirements. | Require explicit future conversion; no direct substitution. |
| MAJOR | Complete bundle immutability is requested for testing, but the Phase 1 object list includes only `KnowledgeBundleManifest`. | Test domain immutability now; reserve physical bundle guarantees for later phases. |
| MINOR | `tests/knowledge/README.md` anticipates `architecture/`, while Phase 1 prompts allow `tests/knowledge/domain/`. | Place Phase 1 architecture tests inside `domain/`; no scope expansion is necessary. |
| MINOR | The architecture describes databases, federation, and bundle files not present today. | These are future components, not missing Phase 1 dependencies. |
| MINOR | Per-phase process calls for plan updates, while the implementation prompt narrowly allows domain files and minimal exports. | Explicitly include a documentation-only ledger update in the future approved scope. No update during this audit. |

**BLOCKER: none for isolated Phase 1 domain design.** The MAJOR findings constrain that design and later integration; they do not authorize refactoring existing modules.

## PROPOSED PHASE 1 PATCH SCOPE

This is a proposed allowlist, not implementation authorization. The Phase 1 domain-design prompt requires a reviewed design followed by an implementation instruction.

### Files allowed to create

```text
simforge/knowledge/domain/__init__.py
simforge/knowledge/domain/identity.py
simforge/knowledge/domain/entities.py
simforge/knowledge/domain/proteins.py
simforge/knowledge/domain/chemistry.py
simforge/knowledge/domain/structures.py
simforge/knowledge/domain/artifacts.py
simforge/knowledge/domain/claims.py
simforge/knowledge/domain/evidence.py
simforge/knowledge/domain/decisions.py
simforge/knowledge/domain/bundle.py
simforge/knowledge/domain/serialization.py

tests/knowledge/domain/test_identity.py
tests/knowledge/domain/test_entities.py
tests/knowledge/domain/test_proteins.py
tests/knowledge/domain/test_chemistry.py
tests/knowledge/domain/test_structures.py
tests/knowledge/domain/test_artifacts.py
tests/knowledge/domain/test_claims_evidence.py
tests/knowledge/domain/test_decisions.py
tests/knowledge/domain/test_bundle_manifest.py
tests/knowledge/domain/test_serialization.py
tests/knowledge/domain/test_architecture.py
```

These files cover only the sixteen objects enumerated in Phase 1, plus necessary private validation/serialization helpers and supporting value types.

No registry, source adapter, importer, resolver, policy evaluator, bundle compiler, filesystem bundle, or worker reader.

### Existing files allowed to edit

- `simforge/knowledge/__init__.py`: minimal domain exports only, if needed.
- `docs/sski/09_IMPLEMENTATION_PLAN.md`: Phase 1 status and validation record only, after completion and review.

No other existing-file edits.

### Allowed import paths and dependencies

Production imports:

- Modules within `simforge.knowledge.domain`.
- Pure standard-library facilities such as `typing`, `dataclasses`, `enum`, `collections.abc`, `datetime`, `uuid`, `re`, `json`, and `hashlib`.
- Existing `pydantic>=2`, if selected as the consistent model mechanism.

**Recommendation:** reuse Pydantic 2’s validation convention, but define independent SSKI models with recursively immutable values and explicit deterministic serialization. Frozen models alone do not freeze nested lists or dictionaries.

Tests may additionally use pytest and standard-library inspection/subprocess facilities. Legacy-type coexistence checks may import existing models in a separate test process; production domain code may not.

No dependency additions or version changes.

### Forbidden paths

Everything outside the explicit allowlist, including:

```text
core/
builders/
executors/
runtime/
pipelines/
descriptors/
validators/
workflows/
adapters/
ligand/
utils/
analysis/
reports/
scripts/
docs/design/
configs/
cli.py
simforge/cli.py
simforge/__main__.py
simforge/__init__.py
pyproject.toml
conftest.py
.github/
```

Also forbidden:

- Existing tests outside the new `tests/knowledge/domain/` files.
- Other SSKI implementation directories or future-phase scaffolding.
- Normative specifications and ADRs, except the narrowly allowed implementation-ledger update.
- All user-owned membrane, campaign, validator, report, and script changes.
- Moving, renaming, staging, reverting, or cleaning unrelated work.

### Forbidden dependencies

- Requests, HTTPX, source API clients, and standard-library network access.
- SQLAlchemy, PostgreSQL drivers, SQLite, DuckDB, and other persistence backends.
- RDKit and chemistry execution tooling.
- GROMACS, OpenMM, AMBER, and engine-specific packages.
- Existing application/runtime/build/analysis packages.
- NumPy, SciPy, PyYAML, Typer, and Rich in the pure domain package; their presence elsewhere does not establish a Phase 1 need.
- Subprocess execution or filesystem persistence in production domain models.

### Mandatory tests

- All new domain and architecture tests listed above.
- The focused existing regression set in section G.
- The existing non-GROMACS/non-RDKit CI gate.
- Preservation of existing optional CI checks, with explicit reporting of skips and unavailable environments.
- Baseline comparison for failures involving pre-existing user work.

### Phase 1 exit criteria

1. Only the approved domain objects and necessary supporting value types exist.
2. Internal IDs are independent of external accessions.
3. Claims, evidence, chemical identities, species, and structural provenance remain semantically distinct.
4. Artifact references enforce explicit SHA256 syntax.
5. Manifest/context values are deeply immutable and deterministically serialized.
6. All four decision states are represented without embedding policy evaluation or execution logic.
7. Dependency tests prove domain isolation in a clean interpreter.
8. Existing regression behavior is preserved; baseline failures are separately documented.
9. No existing application model, CLI, dependency configuration, or user-owned work is changed.
10. No actual bundle compiler, storage, federation, registry, or runtime integration is introduced.
11. Architecture review has no unresolved Phase 1 BLOCKER or MAJOR finding.
12. Changed files, test results, assumptions, and deferred physical-bundle guarantees are recorded precisely.
