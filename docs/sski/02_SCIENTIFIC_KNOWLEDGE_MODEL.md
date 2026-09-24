# Scientific Knowledge Model v0

This document defines the initial domain vocabulary. It is intentionally smaller than the future ontology.

## Entity

A uniquely identified SSKI object. Required conceptual fields: `sf_id`, `entity_type`, and lifecycle/version metadata where applicable.

An Entity is not identified by an external database accession.

## Protein domain

### Protein
A biological protein concept. It is not an exact amino-acid sequence, PDB chain, topology molecule or single structure.

### ProteinSequence
An exact ordered amino-acid sequence with deterministic content fingerprinting.

### ProteinIsoform
A biologically meaningful sequence form related to a Protein.

### ProteinVariant
A sequence differing from a defined reference under explicit variant semantics.

## Chemical domain

### ChemicalIdentity
Connectivity/stereochemical identity independent of a specific simulation microstate.

### ChemicalSpecies
Chemically explicit species with formal charge/protonation/tautomer state where known.

### SimulationMicrostate
Microstate actually chosen for simulation under a defined context such as pH, ionic strength, temperature or binding environment.

ChemicalIdentity and ChemicalSpecies MUST NOT be collapsed. SimulationMicrostate may be derived locally by SimForge and belongs in the RunDerivationLedger when calculated there.

## Structural domain

### ExperimentalStructure
Structure derived from an experimental structural biology source. Preserve source accession, revision/version metadata, method, relevant quality metadata and artifact hash.

### PredictedStructure
Computational structural model. Must remain epistemically distinguishable from ExperimentalStructure.

### BiologicalAssembly
Source-defined or explicitly inferred biological assembly.

### Chain
Concrete polymer instance within a structural artifact.

### ResidueMapping
Explicit mapping among residue coordinate/sequence spaces:

```text
UniProt residue
↔ source polymer entity
↔ source chain residue
↔ prepared structure residue
↔ topology residue
```

Mappings preserve gaps, missing residues and uncertainty.

## Membrane domain

### MembraneContext
Requested scientific membrane environment: taxonomic context, tissue/cell, organelle, temperature and known lipid information where available.

### MembraneTemplate
Practical fallback assembly model. It is not biological truth. It must include scope, composition, evidence level, intended use and limitations.

## Epistemic domain

### Claim
Machine-readable scientific assertion.

### Evidence
Support/context/observation associated with Claims. It retains source, release/revision, evidence class and import provenance.

### ConfidenceAssessment
Contextual assessment attached to an assertion or decision. No universal fixed weight is assigned solely from evidence class.

## Artifact domain

### Artifact
Physical/logical object used by a workflow. Every bundle artifact must have a cryptographic content hash.

### DerivedScientificArtifact
Artifact generated from scientific inputs, e.g. protonated structure, molecular topology, rebuilt loop or assembled membrane system. Derivation provenance is mandatory.

## Decision domain

### Policy
Versioned deterministic rules.

### Decision
Allowed states: `RESOLVED`, `RESOLVED_WITH_WARNING`, `DEFERRED`, `ABORT`.

### DecisionContext
Workflow-scoped record of resolved entities, selected artifacts, warnings, deferred operations and policy outputs.

### KnowledgeBundle
Immutable compiled scientific context consumed by SimForge.

### KnowledgeSnapshot
Reproducibility record spanning INPUT, KNOWLEDGE, POLICY, DECISION and ARTIFACTS.

### RunDerivationLedger
Append-only record of transformations produced after the SSKI/SimForge boundary.

## Explicit non-equivalences

```text
Protein != ProteinSequence
ProteinSequence != PDB chain
PDB chain != prepared structure
ChemicalIdentity != ChemicalSpecies
ChemicalSpecies != SimulationMicrostate
Claim != Evidence
ExternalIdentifier != SSKI primary key
MembraneTemplate != experimentally established membrane composition
PredictedStructure != ExperimentalStructure
KnowledgeBundle != RunDerivationLedger
```
