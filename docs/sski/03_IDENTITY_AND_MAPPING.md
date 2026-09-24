# Identity and Mapping Specification

## Internal identity

SSKI uses internal identifiers independent of external providers.

```text
sf:protein:<opaque-id>
sf:sequence:<opaque-id>
sf:chemical:<opaque-id>
sf:structure:<opaque-id>
sf:artifact:<opaque-id>
sf:claim:<opaque-id>
```

External accessions are mappings, never SSKI primary keys.

## ExternalIdentifier

Conceptual fields:

```text
entity_id
namespace
accession
version/revision
valid_from
valid_until
source metadata
```

Namespaces may include UniProt, PDB, PDB CCD, PubChem, ChEMBL, ChEBI, Rhea and NCBI taxonomy.

## Sequence identity

Exact amino-acid sequences SHOULD receive deterministic content fingerprints. Canonicalization rules must be explicit and versioned.

Fuzzy similarity may generate `CANDIDATE_MAPPING`, never silent identity.

## Chemical identity

Chemical resolution must preserve where relevant:

- connectivity
- stereochemistry
- isotopes
- parent identity
- formal species distinction

Different protonation states must not be merged solely because they share a parent compound.

```text
external record
    ↓
ChemicalIdentity
    ↓
ChemicalSpecies
    ↓
SimulationMicrostate
```

## Structure mappings

```text
Protein / ProteinSequence
        ↓
source polymer entity
        ↓
source chain instance
        ↓
prepared structure
        ↓
topology molecule
```

Residue mappings must support source sequence coordinates, label IDs, author numbering, insertion codes, missing residues, engineered mutations, prepared-system renumbering and topology indices.

No implicit `residue N == residue N` assumption is allowed across coordinate spaces.

## Alias and supersession

Registry semantics must support synonyms, deprecated accessions, merged/split records, source supersession and revised records.

Historical workflows must continue to resolve the exact pinned version/revision recorded in their KnowledgeSnapshot.

## Candidate vs confirmed

Use at least:

```text
CANDIDATE_MAPPING
CONFIRMED_MAPPING
```

ML, embeddings and similarity may generate candidates. Promotion requires an explicit deterministic rule or sufficient recorded evidence.
