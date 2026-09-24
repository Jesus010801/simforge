# KnowledgeBundle Specification

## Purpose

KnowledgeBundle is the immutable contract crossing the SSKI/SimForge boundary. It contains all workflow-required scientific knowledge, decisions, policies and physical artifacts needed to begin without scientific federation.

## Canonical layout

```text
knowledge_bundle/
├── manifest.json
├── runtime.sqlite
│
├── knowledge/
│   ├── entities.parquet
│   ├── structures.parquet
│   ├── compounds.parquet
│   ├── claims.parquet
│   └── provenance.jsonl
│
├── contracts/
│   ├── decision_context.json
│   ├── policies.yaml
│   └── preparation_contract.yaml
│
└── artifacts/
    ├── structures/
    ├── compounds/
    ├── sequences/
    └── membrane/
```

Not every bundle requires every optional table/artifact class.

## Invariants

A compiled bundle MUST:

- have a unique identity
- be immutable after compilation
- pin the KnowledgeSnapshot
- pin policy versions
- hash all physical artifacts
- contain all required dependencies
- validate without Internet access
- fail compilation if mandatory data are unavailable

## Five reproducibility pillars

1. **INPUT** — exact normalized request
2. **KNOWLEDGE** — entities and claims consulted/selected
3. **POLICY** — exact policy IDs/versions
4. **DECISION** — resolution outputs and warnings
5. **ARTIFACTS** — exact hashes, source versions and bundle paths

## Manifest concept

```json
{
  "bundle_id": "sfkb:...",
  "schema_version": "0.1.0",
  "created_at": "...",
  "knowledge_snapshot_id": "...",
  "immutable": true,
  "artifacts": [
    {
      "artifact_id": "sf:artifact:...",
      "sha256": "...",
      "role": "receptor_structure",
      "path": "artifacts/structures/..."
    }
  ]
}
```

## Knowledge vs contracts

`knowledge/` answers: **what did SSKI know/include?**

`contracts/` answers: **what did SSKI resolve or instruct SimForge to do?**

## Validation gate

Before job submission verify:

- manifest/schema compatibility
- artifact presence
- SHA256 integrity
- required mappings
- required policy data
- preparation contract completeness
- no mandatory ABORT state

## No mutation

Workers must not append generated topologies/protonated structures into the bundle. Derived outputs belong to the run workspace and RunDerivationLedger.
