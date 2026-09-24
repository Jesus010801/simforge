# Evidence and Provenance Specification

## Principle

External database content is evidence-bearing source material, not automatic truth.

SSKI must be able to answer:

- Who asserted this?
- From which source/release?
- Was it experimental, curated, predicted or inferred?
- Which importer/version normalized it?
- Which claims support or contradict it?
- Which claim was used by a decision?

## Claim

Conceptual representation:

```text
claim_id
subject
predicate
object/value
qualifiers
scope
provenance links
```

Claims should be granular enough that disagreement can be represented explicitly.

## Evidence classes

Initial classes:

- experimental
- curated
- computational_prediction
- computational_simulation
- inferred
- user_supplied

These classes do not carry universal fixed numeric weights.

## Evidence record

Conceptual fields:

```text
evidence_id
claim_id(s)
source_id
source_record_id
source_release/revision
evidence_class
method
publication/reference
extraction_method
importer_version
imported_at
```

## Provenance chain

```text
source record
    ↓
raw payload/artifact hash
    ↓
KnowledgeImporter version
    ↓
normalized Entity / Claim / Evidence
    ↓
Policy version
    ↓
Decision
```

## Contradictions

Contradictory claims are preserved. The resolver may select under policy, preserve multiple candidates, mark uncertainty, defer resolution or abort.

Conflict resolution is workflow-context dependent.

## Derived SimForge results

Simulation outputs are not inserted as external truth. They may later be explicitly ingested as `computational_simulation` evidence after a defined curation/import step.

## Confidence

Confidence is contextual and should preserve target assertion/decision, method, inputs, limitations and policy/version. Avoid unexplained universal scalars.
