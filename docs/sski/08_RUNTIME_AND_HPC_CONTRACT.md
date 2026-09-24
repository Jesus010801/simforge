# Runtime and HPC Contract

## RUNTIME_INVARIANT_001

A SimForge compute worker is strictly forbidden from:

1. calling external scientific APIs
2. querying the SSKI Control Plane
3. querying the Global Knowledge Plane
4. mutating its KnowledgeBundle
5. dynamically federating missing scientific knowledge

Mandatory missing information must fail fast with an explicit error such as:

```text
KnowledgeBundleIncompleteError
```

before expensive simulation work begins.

## Allowed behavior

Workers MAY:

- read the local bundle
- validate hashes
- execute an authorized PreparationContract
- calculate a deferred microstate
- generate derived topology
- assemble a simulation system
- write derived artifacts
- append to RunDerivationLedger
- execute a SimulationPlan

These operations are not scientific federation.

## Data staging

```text
Control Plane
    ↓
validated KnowledgeBundle
    ↓
shared project/object storage
    ↓
job staging
    ↓
node-local scratch
    ↓
worker validation
    ↓
preparation/execution
```

## Runtime storage roles

### SQLite
Compact mappings, lookup, bundle metadata and provenance indices.

### DuckDB / Parquet
Large local analytical tables and candidate collections.

### Filesystem
mmCIF/PDB, SDF/MOL, FASTA, transforms, parameterization inputs and physical artifacts.

## RunDerivationLedger

Append-only. Each operation records parent bundle/hash, operation, software/version, input hashes, parameters, output hashes, warnings, decision result and timestamps.

## Job submission gate

Submit compute jobs only after:

```text
bundle compiled
AND bundle validated
AND mandatory decisions != ABORT
AND required artifacts present
```
