# Federation and Scientific Sources

## Federation architecture

### SourceAdapter — transport layer

Owns:

- network/file transport
- authentication if required
- rate-limit handling
- cache access
- integrity checks
- source/release metadata
- raw payload acquisition

It MUST NOT make scientific decisions.

### KnowledgeImporter — semantic layer

Owns:

- parsing
- unit normalization
- entity emission
- external identifier emission
- claim emission
- evidence/provenance emission

It MUST NOT perform network access.

## Initial MVP sources

### UniProt
Protein identity, canonical nomenclature and sequence mapping.

### RCSB PDB / wwPDB data
Experimental structure candidates, assemblies, chains/polymer entities, quality metadata and revision-aware artifacts.

### PDB Chemical Component Dictionary
Chemical component identity/base connectivity for structure-associated non-polymers.

### PubChem
Fallback/enrichment for chemical identity; not mandatory if CCD is sufficient.

### Membrane sources (conditional/later)
OPM, MemProtMD, LIPID MAPS and other curated membrane/lipid resources.

## Later source families

Potential future integration:

- ChEMBL
- ChEBI
- Rhea
- Reactome
- InterPro/Pfam
- AlphaFold DB/model repositories
- BindingDB
- Guide to Pharmacology
- ClinVar/dbSNP/gnomAD
- interaction resources
- taxonomy resources
- specialized lipidomics resources

Licensing/redistribution must be reviewed before packaging redistributed content.

## Refresh strategies

TTL is source/resource specific.

```yaml
rcsb_entry:
  strategy: revision_aware
uniprot_record:
  strategy: release_based
predicted_model:
  strategy: version_check
bulk_dataset:
  strategy: release_manifest
```

Do not hardcode universal `max_age_days`.

## Cache behavior

Prefer verified cached content for repeat retrieval. Raw payloads/artifacts should be content-hashed. HTTP 429/5xx uses bounded retry/backoff; infinite retry is forbidden.

Federation failures occur before compute execution.

## Bulk vs on-demand

Use bulk/release synchronization for large stable datasets when appropriate, on-demand enrichment for workflow-specific detail, and incremental updates when supported.
