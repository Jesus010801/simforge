# SimForge Scientific Knowledge Infrastructure Protocol (SSKI)

**Documento:** Protocolo de arquitectura para acceso universal, actualización continua y razonamiento trazable sobre conocimiento molecular  
**Proyecto:** SimForge  
**Versión:** 1.0-draft  
**Fecha de diseño:** 2026-09-23  
**Estado:** Arquitectura propuesta; no implica implementación existente  

---

## 0. Resumen ejecutivo

SimForge no debe intentar convertirse en una copia física de UniProt, PDB, PubChem, ChEMBL, AlphaFold DB, Gene Ontology, Reactome, LIPID MAPS y el resto del ecosistema biomolecular. Tampoco debe depender de un modelo de lenguaje que “memorice” esos datos.

La arquitectura óptima es una **infraestructura científica federada y versionada** que mantenga localmente una identidad universal mínima, relaciones normalizadas, evidencia, procedencia, fingerprints y snapshots de decisiones; recupere el detalle desde fuentes autoritativas cuando haga falta; y conserve artefactos grandes únicamente cuando sean necesarios para reproducibilidad o trabajo activo.

El sistema propuesto se denomina en este documento:

> **SimForge Scientific Knowledge Infrastructure (SSKI)**

SSKI debe funcionar como cinco cosas simultáneamente:

1. **Universal Entity Registry** — sabe qué entidad científica es cada cosa y cómo se relacionan sus identificadores externos.
2. **Federation Engine** — sabe dónde consultar la información más reciente y por qué vía.
3. **Scientific Memory** — recuerda conocimiento ya resuelto, evidencias, resultados propios y decisiones previas.
4. **Evidence & Provenance Engine** — nunca convierte una fuente externa en verdad automática; conserva quién afirma qué, cuándo y con qué evidencia.
5. **Knowledge Resolver** — transforma información potencialmente contradictoria en contexto estructurado para el Decision Engine de SimForge.

Principio fundacional:

> **SimForge no almacena el conocimiento universal completo. SimForge sabe identificarlo, localizarlo, verificarlo, relacionarlo, versionarlo, recordar lo necesario y reconstruir exactamente de dónde provino.**

---

# 1. Objetivos

## 1.1 Objetivo principal

Dar a SimForge acceso práctico al mayor volumen posible de conocimiento científico relevante para simulación molecular, biofísica computacional, química computacional, biología estructural y diseño molecular sin requerir una réplica completa de las bases de datos externas.

## 1.2 Dominios que SSKI debe cubrir

SSKI debe estar diseñado desde el inicio para representar, como mínimo:

- genes;
- proteínas;
- isoformas;
- secuencias;
- variantes y mutaciones;
- complejos macromoleculares;
- estructuras experimentales;
- modelos estructurales predichos;
- ensamblajes biológicos;
- cadenas y dominios;
- sitios activos y sitios de unión;
- péptidos;
- ligandos;
- metabolitos;
- fármacos;
- iones;
- cofactores;
- solventes;
- lípidos;
- membranas;
- composiciones de membrana;
- proteínas transmembranales;
- transportadores;
- canales;
- receptores;
- reacciones;
- actividades enzimáticas;
- vías metabólicas y de señalización;
- interacciones proteína-proteína;
- interacciones proteína-ligando;
- afinidades y bioactividades;
- modificaciones postraduccionales;
- ontologías funcionales;
- taxonomía;
- enfermedades y fenotipos;
- evidencia bibliográfica;
- parámetros físico-químicos;
- estados de protonación y tautomería;
- topologías de simulación;
- parámetros de force field;
- sistemas de simulación;
- resultados de docking, MD, free-energy, QM y QM/MM;
- conocimiento producido localmente por SimForge.

## 1.3 Requisitos no negociables

El sistema debe ser:

- **actualizable**;
- **reproducible**;
- **agnóstico a una sola base de datos**;
- **tolerante a contradicciones**;
- **trazable**;
- **versionado**;
- **cacheable**;
- **capaz de trabajar offline parcialmente**;
- **compatible con millones o cientos de millones de entidades**;
- **extensible a nuevas fuentes sin modificar el núcleo**;
- **capaz de distinguir hechos, anotaciones, predicciones e inferencias**;
- **capaz de preservar conocimiento local sin mezclarlo con conocimiento público**;
- **capaz de servir al Decision Engine sin entregarle “verdades” no justificadas**.

---

# 2. Lo que NO se debe construir

## 2.1 No construir un espejo universal monolítico

No descargar y mantener permanentemente todos los archivos de todas las fuentes.

Problemas:

- almacenamiento enorme;
- duplicación de información;
- sincronización compleja;
- licencias heterogéneas;
- versiones incompatibles;
- alto coste de actualización;
- información que nunca será utilizada.

## 2.2 No usar un LLM como base de datos

Un modelo puede ayudar a interpretar, clasificar o generar hipótesis, pero no debe ser la memoria factual principal.

Un modelo no garantiza:

- exactitud factual;
- fecha de actualización;
- procedencia;
- versión;
- reproducibilidad;
- ausencia de alucinaciones.

## 2.3 No convertir identificadores externos en claves internas

Incorrecto:

```text
protein.id = P00533
```

Correcto:

```text
sf_entity_id = sf:protein:01...
external_identifier = UniProtKB:P00533
```

Esto permite sobrevivir a merges, obsolescencias, cambios de accession y discrepancias entre bases.

## 2.4 No permitir que el Decision Engine consulte directamente fuentes externas

Prohibido arquitectónicamente:

```text
DecisionEngine -> requests.get(UniProt)
```

Ruta obligatoria:

```text
DecisionEngine
      ↑
DecisionContext
      ↑
KnowledgeResolver
      ↑
Scientific Memory / Federation Engine
```

---

# 3. Arquitectura general

```text
                     EXTERNAL SCIENTIFIC ECOSYSTEM

  UniProt   PDB   AlphaFold   InterPro   GO   PubChem   ChEMBL
     |       |        |          |       |       |         |
  ChEBI   UniChem   Rhea   Reactome   BindingDB   LIPID MAPS
     |       |        |       |           |            |
  ClinVar  NCBI   MemProtMD  ModelArchive  ... future sources
     \       |       /          /          /
      \      |      /          /          /
       +-----+-----+----------+----------+
                       |
                 SOURCE ADAPTERS
                       |
                 SOURCE REGISTRY
                       |
               FEDERATION ENGINE
                       |
                 NORMALIZATION
                       |
              ENTITY RESOLUTION
                       |
        +--------------+---------------+
        |                              |
 UNIVERSAL ENTITY REGISTRY      EVIDENCE / CLAIM STORE
        |                              |
        +--------------+---------------+
                       |
                KNOWLEDGE GRAPH
                       |
              SCIENTIFIC MEMORY
                       |
              KNOWLEDGE RESOLVER
                       |
                DECISION CONTEXT
                       |
                 DECISION ENGINE
                       |
                    SIMFORGE
                       |
          LOCAL SCIENTIFIC RESULTS
                       |
             EVIDENCE / CLAIM STORE
```

---

# 4. Las cinco capas de SSKI

## 4.1 Capa A — Universal Entity Registry (UER)

Mantiene identidad mínima y mappings.

Debe responder:

- ¿qué entidad es esto?;
- ¿qué identificadores externos tiene?;
- ¿qué tipo de entidad es?;
- ¿qué fingerprints permiten reconocerla?;
- ¿ha sido sustituida, dividida o fusionada?;
- ¿qué fuentes conocen esta entidad?;

No debe almacenar necesariamente el contenido completo de las fuentes.

### Registro mínimo

```yaml
sf_entity_id: sf:protein:01J...
entity_type: protein
canonical_label: EGFR
status: active

external_ids:
  - namespace: UniProtKB
    accession: P00533
  - namespace: NCBI_Gene
    accession: "1956"

identity_fingerprints:
  sequence_sha256: ...

source_presence:
  - UniProtKB
  - PDB
  - AlphaFoldDB
  - InterPro

created_at: ...
updated_at: ...
```

## 4.2 Capa B — Federation Engine

Decide qué fuente consultar, en qué orden, bajo qué contrato y con qué fallback.

Ejemplo:

```text
Need: experimental protein structure

1. RCSB/wwPDB
2. PDBe/PDB mappings if needed
3. ModelArchive only if request accepts computational models
4. AlphaFold DB only if predicted models are acceptable
```

El motor debe comprender que fuentes diferentes responden preguntas diferentes.

## 4.3 Capa C — Scientific Memory

Conserva conocimiento normalizado ya utilizado o de alto valor.

Incluye:

- entidades resueltas;
- aliases;
- mappings;
- claims;
- relaciones;
- evidencia;
- fingerprints;
- embeddings opcionales;
- historial de resoluciones;
- conocimiento local;
- resultados de SimForge;
- snapshots.

## 4.4 Capa D — Evidence & Provenance

Cada afirmación debe conservar su procedencia.

Nunca:

```text
EGFR.optimal_temperature = 310.15
```

Sí:

```yaml
claim:
  subject: sf:protein:...
  predicate: optimal_temperature
  value: 310.15
  unit: K

provenance:
  source: ...
  source_release: ...
  source_record: ...
  retrieved_at: ...
  evidence_type: experimental
  publication: ...
```

## 4.5 Capa E — Knowledge Resolver

Convierte claims potencialmente múltiples y contradictorios en una respuesta epistemológicamente explícita.

Ejemplo:

```yaml
resolution:
  property: physiological_temperature
  status: conflicting

candidates:
  - value: 287.15
    evidence: [...]
  - value: 298.15
    evidence: [...]

recommended_action:
  mode: user_confirmation_required
```

---

# 5. Modelo de identidad universal

## 5.1 Tipos de identidad

SSKI debe distinguir:

### Identidad conceptual

“EGFR humano”.

### Identidad secuencial

Una secuencia aminoacídica exacta.

### Identidad estructural

Una realización tridimensional específica.

### Identidad química

Un grafo molecular específico.

### Identidad de especie química

Mismo compuesto padre, pero estado de protonación, carga, isotopía o tautomería definido.

### Identidad de simulación

Entidad después de preparación y parametrización.

## 5.2 Proteínas

Jerarquía recomendada:

```text
ProteinConcept
    |
    +-- ProteinIsoform
           |
           +-- ProteinSequence
                  |
                  +-- ProteinVariant
                         |
                         +-- ProteinState
```

No tratar UniProt ↔ PDB ↔ topology como equivalencias directas.

```text
ProteinSequence
     |
     | mapped_to
     v
PDB Polymer Entity
     |
     | instantiated_as
     v
PDB Chain Instance
     |
     | prepared_as
     v
PreparedStructure
     |
     | parameterized_as
     v
TopologyMolecule
```

## 5.3 Moléculas pequeñas

```text
ChemicalConcept
    |
    +-- ParentCompound
            |
            +-- ChemicalSpecies
                    |
                    +-- Stereoisomer
                            |
                            +-- Conformer
```

### Identificadores y fingerprints recomendados

- InChI;
- Standard InChIKey;
- canonical SMILES;
- isomeric SMILES;
- formal charge;
- stereochemistry;
- isotopic specification;
- molecular graph hash;
- tautomer parent;
- connectivity layer.

## 5.4 Lípidos

Los lípidos requieren modelado explícito de:

- identidad estructural;
- clase lipídica;
- cadena(s) acilo/alquilo;
- sn-position cuando se conozca;
- estado de protonación;
- composición agregada cuando la especie exacta no esté determinada;
- nomenclatura equivalente entre LIPID MAPS, SwissLipids, ChEBI y otras fuentes.

## 5.5 Membranas

Una membrana no debe ser sólo `membrane=true`.

Entidad propuesta:

```yaml
MembraneContext:
  organism:
  organelle:
  tissue:
  leaflet_asymmetry:
  composition:
    outer_leaflet: []
    inner_leaflet: []
  temperature:
  ionic_environment:
  source_claims: []
```

Debe existir también:

- `MembraneComposition`;
- `MembraneComponentFraction`;
- `BilayerModel`;
- `MembraneProteinOrientation`;
- `TransmembraneSegment`;
- `TopologyAnnotation`;
- `MembraneSimulationReference`.

---

# 6. Modelo de conocimiento y evidencia

## 6.1 Claim como unidad fundamental

La base no guarda “verdades”; guarda afirmaciones trazables.

```text
Claim = subject + predicate + object/value + context + provenance
```

Ejemplo:

```yaml
claim_id: sfclaim:...
subject: sf:protein:...
predicate: binds
object: sf:chemical:...
context:
  organism: Homo sapiens
  assay_type: biochemical
provenance:
  source: ChEMBL
  record: ...
```

## 6.2 Clases de afirmación

```text
experimental_measurement
curated_annotation
reported_observation
computational_prediction
computational_simulation
inference
user_annotation
hypothesis
recommendation
```

## 6.3 Estado epistemológico

Cada claim puede estar en:

```text
reported
supported
conflicting
superseded
retracted
uncertain
locally_validated
```

## 6.4 Contradicciones

Nunca sobrescribir automáticamente un valor incompatible.

```text
Claim A ---- CONTRADICTS ---- Claim B
```

El sistema puede resolver operacionalmente una decisión sin borrar la contradicción original.

## 6.5 Evidencia negativa

Debe existir soporte explícito para negaciones.

Ejemplo:

```text
NOT_LOCATED_IN
DOES_NOT_BIND
NO_ACTIVITY_DETECTED
```

Una ausencia de dato no equivale a evidencia negativa.

---

# 7. Fuentes científicas: catálogo inicial

El Source Registry debe tratar cada fuente como un proveedor con capacidades, versión, licencia, método de acceso, cadencia de actualización y fiabilidad contextual.

## 7.1 Proteínas, genes y secuencias

### Tier 1

- **UniProtKB / UniParc / UniRef** — identidad proteica, secuencias, anotación, cross-references.
- **NCBI Gene / RefSeq / Protein** — genes, transcriptos, proteínas, assemblies y taxonomía.
- **Ensembl** — genes, transcriptos, ortología y variantes genómicas.

### Integración

Uso primario:

- identidad;
- secuencia;
- isoformas;
- taxonomía;
- nomenclatura;
- cross-reference.

## 7.2 Estructuras macromoleculares

### Tier 1

- **wwPDB / RCSB PDB / PDBe** — estructuras experimentales y metadatos.
- **AlphaFold Protein Structure Database** — modelos predichos.
- **ModelArchive** — modelos computacionales depositados.

### Datos estructurales derivados

- assemblies;
- chain mappings;
- sequence mapping;
- ligands;
- resolution/method;
- missing residues;
- modifications;
- experimental conditions;
- predicted confidence.

## 7.3 Dominios y función

- **InterPro**;
- **Pfam** vía InterPro;
- **CATH-Gene3D**;
- **PROSITE**;
- **SMART**;
- **PANTHER**;
- **Gene Ontology**.

## 7.4 Moléculas, ligandos y metabolitos

### Tier 1

- **PubChem**;
- **ChEBI**;
- **ChEMBL**;
- **UniChem** como capa de cross-mapping químico;
- **PDB Chemical Component Dictionary (CCD)**.

### Tier 2

- DrugCentral;
- HMDB donde los términos de uso permitan la integración deseada;
- ZINC para casos compatibles con sus condiciones de acceso;
- otras fuentes especializadas.

## 7.5 Bioactividad y farmacología

- **ChEMBL**;
- **BindingDB**;
- **PubChem BioAssay**;
- **IUPHAR/BPS Guide to Pharmacology**;
- DrugCentral.

Representar:

- Ki;
- Kd;
- IC50;
- EC50;
- kon/koff cuando existan;
- assay type;
- target construct;
- species;
- conditions;
- units;
- qualifiers;
- publication.

Nunca comparar valores sin conservar contexto experimental.

## 7.6 Enzimas y reacciones

- **Rhea**;
- EC nomenclature;
- BRENDA cuando licencia y condiciones permitan el uso previsto;
- SABIO-RK;
- Reactome;
- MetaCyc/KEGG sólo conforme a términos aplicables.

## 7.7 Interacciones y complejos

- IntAct;
- Complex Portal;
- BioGRID;
- STRING con distinción clara entre evidencia, asociación e inferencia;
- Reactome.

## 7.8 Pathways

- Reactome;
- WikiPathways;
- GO-CAM;
- Rhea para reacciones bioquímicas.

## 7.9 Lípidos y membranas

### Lípidos

- **LIPID MAPS**;
- **SwissLipids**;
- ChEBI;
- Rhea.

### Proteínas de membrana

- **MemProtMD**;
- OPM;
- PDBTM/UniTmp cuando sea utilizable;
- TCDB para clasificación de transportadores;
- GPCRdb para GPCRs cuando aplique;
- IUPHAR para farmacología de receptores y canales.

### Conocimiento específico de simulación de membrana

Conservar:

- orientación;
- profundidad de inserción;
- segmentos TM;
- composición lipídica empleada;
- grosor;
- curvatura si se reporta;
- oligomerización;
- sistema experimental/biológico;
- condiciones de temperatura y salinidad;
- fuente del modelo de membrana.

## 7.10 Variantes y enfermedad

- **ClinVar**;
- dbSNP;
- gnomAD;
- Ensembl Variation;
- CIViC y otras fuentes especializadas cuando sean relevantes y licenciables.

## 7.11 Taxonomía

- NCBI Taxonomy como identidad taxonómica primaria;
- mappings complementarios de UniProt y Ensembl.

## 7.12 Literatura

- PubMed / NCBI;
- Europe PMC;
- Crossref para DOI y metadatos.

SSKI no necesita almacenar los artículos completos para utilizar sus metadatos y relaciones bibliográficas.

---

# 8. Source Registry

Cada fuente debe registrarse mediante un descriptor declarativo.

```yaml
source_id: rcsb_pdb
name: RCSB Protein Data Bank
entity_domains:
  - structure
  - protein
  - ligand

access:
  rest: true
  graphql: true
  bulk: true
  incremental: true

update_strategy:
  release_discovery: holdings_endpoint
  preferred_interval: 6h

provenance:
  version_method: source_release_or_retrieval_timestamp
  immutable_record_hash: true

license:
  redistribution: allowed
  attribution_required: true

priority:
  experimental_structure: 1
```

Campos obligatorios:

```text
source_id
name
domains
capabilities
API endpoints
bulk endpoints
release/version discovery
rate limits
license metadata
update cadence
health status
parser version
normalizer version
last successful sync
last source release observed
```

---

# 9. Contrato obligatorio de Source Adapter

Cada fuente implementa la misma interfaz conceptual.

```python
class KnowledgeSourceAdapter:
    def capabilities(self): ...
    def healthcheck(self): ...
    def current_release(self): ...
    def search(self, query): ...
    def fetch(self, external_id): ...
    def fetch_changed_since(self, cursor): ...
    def normalize(self, source_record): ...
    def provenance(self, source_record): ...
```

Opcionales:

```python
bulk_manifest()
stream_changes()
resolve_external_id()
fetch_artifact()
license_policy()
```

El núcleo de SimForge no debe contener lógica específica de UniProt, PDB o ChEMBL.

---

# 10. Estrategia de acceso: híbrida, no monolítica

## 10.1 Modo 1 — Universal index synchronization

Sincronizar únicamente información compacta y de alto valor para resolución:

- identificadores;
- aliases;
- tipo de entidad;
- fingerprints;
- checksums;
- cross-references;
- estado active/obsolete;
- release de origen.

Objetivo:

> Poder reconocer una entidad sin tener que descargar su payload completo.

## 10.2 Modo 2 — On-demand enrichment

Cuando un workflow necesite información no almacenada:

```text
query
  -> local memory
  -> source registry
  -> best source(s)
  -> fetch
  -> normalize
  -> evidence store
  -> cache
```

## 10.3 Modo 3 — Scheduled enrichment

Para entidades activas de proyectos:

```text
hot entities
  -> periodic refresh
  -> differential comparison
  -> new claims
  -> invalidate affected resolutions
```

## 10.4 Modo 4 — Bulk release ingestion cuando sea eficiente

Bulk sólo para conjuntos donde:

- el proveedor recomienda bulk;
- el dataset es razonablemente manejable;
- existen releases versionados;
- la información compacta es de uso transversal.

El archivo descargado puede ser **staging efímero**:

```text
download
 -> verify checksum
 -> parse
 -> normalize
 -> commit entities/claims
 -> record source manifest/hash
 -> delete staging payload
```

No se requiere conservar el dump original indefinidamente.

## 10.5 Modo 5 — Reproducibility artifact pinning

Si un artefacto externo fue realmente utilizado para producir resultados científicos, puede marcarse como `pinned`.

Ejemplos:

- mmCIF usado para preparar una simulación;
- FASTA usado para modelar una proteína;
- SDF usado para parametrizar un ligando;
- release específico de GO usado en un análisis.

Estos artefactos sí deben conservarse o archivarse de forma recuperable.

---

# 11. Jerarquía de almacenamiento

## 11.1 L0 — Universal Entity Index

Persistente y compacto.

Contiene:

- entity IDs;
- external IDs;
- canonical labels;
- type;
- hashes/fingerprints;
- source presence;
- obsolescence mappings.

## 11.2 L1 — Scientific Memory

Persistente.

Contiene:

- claims;
- relationships;
- evidence;
- normalized properties;
- mappings;
- resolution history;
- project/local knowledge.

## 11.3 L2 — Artifact Cache

Evictable salvo contenido `pinned`.

Contiene:

- PDB/mmCIF;
- FASTA;
- SDF/MOL;
- PAE matrices;
- ontology files;
- selected source payloads.

## 11.4 L3 — Reproducibility Archive

Inmutable.

Sólo aquello que formó parte directa de resultados científicos reproducibles.

---

# 12. Stack de persistencia recomendado

## 12.1 PostgreSQL — Source of truth operacional

Usar para:

- entities;
- identifiers;
- aliases;
- relationships;
- claims;
- evidence;
- source registry;
- synchronization cursors;
- snapshots;
- decision traces;
- local knowledge.

Extensiones útiles:

- `pgvector` para embeddings;
- índices GIN/GiST;
- JSONB para campos fuente-específicos controlados.

## 12.2 Parquet — datasets analíticos y snapshots grandes

Usar para:

- bioactivities masivas;
- mappings masivos;
- releases compactados;
- export analítico;
- historical snapshots.

Consulta con DuckDB/Polars cuando convenga.

## 12.3 Object storage — artifacts

Backend local:

- filesystem content-addressable inicialmente;
- MinIO cuando se requiera una interfaz S3;
- S3-compatible remoto en instalaciones distribuidas.

## 12.4 Graph layer

Fase inicial:

- relaciones en PostgreSQL.

No introducir Neo4j/Memgraph/Kuzu como source of truth al principio.

Cuando traversal y escala lo justifiquen:

```text
PostgreSQL source of truth
        |
        v
Graph projection
        |
        +-- Kuzu / Memgraph / Neo4j / equivalent
```

La proyección se puede reconstruir.

## 12.5 Search index

Opcional según escala:

- PostgreSQL full-text para MVP;
- OpenSearch/Elasticsearch sólo si búsquedas textuales y facetas lo requieren.

---

# 13. Content-addressable storage

Cada artefacto debe identificarse también por hash.

```text
sha256(payload)
```

Ejemplo:

```yaml
artifact_id: sfartifact:...
sha256: 64f...
media_type: chemical/x-pdbx-mmcif
source: RCSB
external_id: 9XYZ
source_release: ...
cache_state: pinned
```

Beneficios:

- deduplicación;
- integridad;
- reproducibilidad;
- detección de cambios silenciosos;
- posibilidad de recuperar la misma versión.

---

# 14. Identity Resolution Engine

## 14.1 Orden de resolución

### Nivel 1 — Identificadores exactos

Ejemplos:

```text
UniProtKB:P00533
PDB:1M17
CHEMBL:CHEMBL25
PubChemCID:2244
ChEBI:15365
```

### Nivel 2 — Fingerprints deterministas

Proteínas:

```text
SHA256(canonical_sequence)
```

Moléculas:

```text
Standard InChIKey
canonical graph hash
connectivity hash
```

### Nivel 3 — Mappings autoritativos

Ejemplos:

- UniProt ID mapping;
- SIFTS/PDB sequence mappings;
- UniChem;
- source-provided cross references.

### Nivel 4 — Similaridad computacional

Proteínas:

- sequence identity;
- alignment coverage;
- domain architecture;
- structural similarity.

Moléculas:

- substructure;
- Tanimoto fingerprints;
- stereochemistry-aware graph matching.

### Nivel 5 — Semantic candidate generation

- names;
- synonyms;
- embeddings;
- literature context.

Este nivel genera candidatos, nunca equivalencia automática.

## 14.2 Resultado de resolución

```yaml
resolution:
  input: "GLP-1R"
  entity: sf:protein:...
  confidence: 0.999
  method: authoritative_alias_plus_taxonomy
  external_ids:
    UniProtKB: P43220
```

## 14.3 Estados

```text
resolved
resolved_with_context
ambiguous
conflicting
unresolved
new_entity_candidate
```

---

# 15. Machine learning y AI dentro de SSKI

AI/ML sí forma parte de la arquitectura, pero **no es la memoria factual**.

## 15.1 Usos permitidos/recomendados

- entity candidate generation;
- synonym discovery;
- duplicate detection;
- structure/sequence similarity ranking;
- literature triage;
- relationship hypothesis generation;
- anomaly detection;
- source conflict explanation;
- natural-language query interpretation;
- source routing;
- ranking de evidencia según contexto;
- clasificación de documentos.

## 15.2 Representaciones aprendidas

Proteínas:

- sequence embeddings;
- structure embeddings;
- functional embeddings.

Moléculas:

- molecular embeddings;
- fingerprints;
- activity embeddings.

Literatura:

- document embeddings;
- claim embeddings.

## 15.3 Regla epistemológica

```text
ML similarity -> candidate / hypothesis
```

Nunca:

```text
ML similarity -> scientific fact
```

## 15.4 Aprendizaje por uso

SSKI puede aprender:

- qué fuentes responden mejor a cada pregunta;
- qué aliases aparecen frecuentemente;
- qué entidades se consultan juntas;
- qué mappings ya fueron confirmados;
- qué datos conviene precachear;
- qué resoluciones fallan con frecuencia.

Esto optimiza acceso, no reescribe hechos científicos.

---

# 16. Modelo de memoria científica

## 16.1 Global/public knowledge

Conocimiento importado o referenciado desde fuentes externas.

Namespace:

```text
public/*
```

## 16.2 Laboratory/local knowledge

Conocimiento propio del usuario/laboratorio.

Ejemplos:

- constructos;
- anotaciones manuales;
- resultados experimentales;
- interpretaciones validadas.

Namespace:

```text
local/*
```

## 16.3 Project knowledge

Sólo válido dentro de un proyecto.

```text
project/<project_id>/*
```

## 16.4 Run knowledge

Resultados específicos de una ejecución.

```text
run/<run_id>/*
```

## 16.5 Precedencia

Una anotación local puede orientar una decisión local, pero nunca sobrescribir silenciosamente el conocimiento público.

---

# 17. Actualización continua: mantener el sistema a la vanguardia

## 17.1 Knowledge Update Scheduler

Servicio permanente:

```text
Source Registry
      |
      v
Release Watcher
      |
      v
Change Detector
      |
      v
Incremental Fetch
      |
      v
Normalize + Validate
      |
      v
Commit new claims/entities
      |
      v
Invalidate affected resolutions
      |
      v
Snapshot
```

## 17.2 Cadencias

No imponer una sola frecuencia global.

Cada fuente define:

```text
release_based
weekly
daily
hourly_possible
on_demand
manual
```

El sistema consulta el release/version endpoint cuando exista y evita descargar si no hubo cambio.

## 17.3 Freshness classes

```text
HOT
  entidades de proyectos activos

WARM
  entidades consultadas recientemente

COLD
  índice global sin payload detallado
```

Política recomendada:

- HOT: refresh agresivo;
- WARM: refresh moderado;
- COLD: actualizar mappings durante releases/index sync.

## 17.4 Staleness metadata

Cada resolución debe incluir:

```yaml
freshness:
  source_release: ...
  retrieved_at: ...
  checked_at: ...
  max_age_policy: ...
  stale: false
```

---

# 18. Snapshotting científico

## 18.1 Knowledge Snapshot

Cada simulación importante debe poder referenciar un snapshot lógico.

```yaml
knowledge_snapshot:
  id: sfks:2026-09-23:000184
  created_at: ...
  ontology_version: ...
  resolver_version: ...
  source_versions:
    uniprot: ...
    pdb: ...
    chembl: ...
```

No significa copiar todos los datasets.

El snapshot registra exactamente qué versiones y claims influyeron en la decisión.

## 18.2 Decision snapshot

Más compacto aún:

```text
snapshot only the closure of knowledge actually used
```

Es decir, preservar:

- claims consultados;
- evidencia;
- source records IDs;
- source versions;
- hashes de artefactos;
- policy versions.

---

# 19. Integración con SimForge Decision Engine

El Knowledge Resolver entrega un objeto `DecisionContext`.

```yaml
DecisionContext:
  subject: sf:protein:...
  question: choose_structure

  evidence:
    experimental_structures: [...]
    predicted_structures: [...]

  constraints:
    sequence_coverage_min: 0.9

  conflicts: []
  uncertainty: ...
  source_snapshot: ...
```

El Decision Engine puede entonces aplicar una policy explícita.

```text
Knowledge -> context
Policy -> decision
Decision -> action
```

No:

```text
Knowledge -> action
```

## 19.1 DecisionTrace obligatorio

```yaml
decision_id: sfdec:...
question: structure_selection
candidates: [...]
selected: ...
policy: structure_selection:v3
knowledge_snapshot: ...
evidence_ids: [...]
uncertainty: ...
human_override: false
```

---

# 20. Protocolo de consulta

## 20.1 Consulta típica

Usuario:

```text
simforge knowledge resolve P00533
```

Proceso:

```text
1. Parse input
2. Detect namespace/identifier
3. Query Universal Entity Registry
4. If exact hit -> return entity
5. Check requested information scope
6. Check freshness
7. Query Scientific Memory
8. If insufficient/stale -> Federation Engine
9. Select authoritative source adapters
10. Fetch
11. Normalize
12. Resolve identity
13. Generate claims
14. Store provenance
15. Resolve conflicts
16. Cache useful result
17. Return KnowledgeResponse
```

## 20.2 KnowledgeResponse

```yaml
entity:
  id: sf:protein:...
  label: ...

resolution:
  status: resolved
  confidence: ...

knowledge:
  ...

freshness:
  ...

sources:
  ...

conflicts:
  ...
```

---

# 21. Protocolo de ingestión de una nueva fuente

Una nueva fuente NO entra directamente a producción.

## Fase A — Source qualification

Evaluar:

- autoridad científica;
- tipo de evidencia;
- mantenimiento activo;
- versión/release;
- API/bulk availability;
- estabilidad;
- licencia;
- redistribución;
- identificadores;
- cobertura;
- actualización;
- calidad de provenance.

## Fase B — Schema mapping

Mapear:

```text
source field -> SimForge ontology concept
```

## Fase C — Adapter implementation

Implementar contrato común.

## Fase D — Golden records

Crear casos de prueba conocidos.

## Fase E — Conflict tests

Comprobar:

- duplicados;
- obsolete IDs;
- missing values;
- contradictory mappings;
- unit normalization.

## Fase F — Shadow ingestion

Cargar sin afectar resoluciones de producción.

## Fase G — Promotion

Activar la fuente dentro de Source Registry.

---

# 22. Licenciamiento y gobernanza

Cada fuente necesita una `SourceLicensePolicy`.

```yaml
source: example
license:
  type: ...
  attribution_required: true
  redistribution_allowed: false
  cache_allowed: true
  local_use_allowed: true
```

El sistema debe poder bloquear:

- redistribución de payloads restringidos;
- publicación de datasets derivados incompatibles;
- almacenamiento permanente cuando no esté permitido.

El conocimiento normalizado también debe conservar la fuente para atribución.

---

# 23. Seguridad e integridad científica

## 23.1 Checksums

Todos los artefactos descargados que influyan en ciencia reproducible deben tener checksum.

## 23.2 TLS

Acceso a fuentes mediante HTTPS cuando esté disponible.

## 23.3 Parser sandboxing

No ejecutar contenido descargado.

## 23.4 Schema validation

Cada adapter valida payloads antes de normalizar.

## 23.5 Quarantine

Cambios inesperados de schema:

```text
source payload
  -> validation failure
  -> quarantine
  -> alert
```

Nunca interpretar silenciosamente campos desconocidos.

---

# 24. Observabilidad

Métricas mínimas:

```text
source_health
source_latency
source_error_rate
source_release_seen
last_successful_sync
entities_resolved
ambiguous_resolutions
conflicts_detected
cache_hit_rate
freshness_age
claims_ingested
claims_superseded
artifact_cache_size
```

Logs estructurados deben contener:

```text
request_id
source_id
adapter_version
release
entity_id
operation
latency
result
```

---

# 25. Resiliencia

## 25.1 Circuit breaker por fuente

Si una API falla repetidamente:

```text
OPEN
  -> use cache / alternate source
  -> retry after cooldown
```

## 25.2 Rate limiting

Cumplir límites del proveedor.

## 25.3 Retry

Sólo para errores transitorios.

No reintentar automáticamente:

- 400 schema errors;
- ambiguous identifiers;
- licensing failures.

## 25.4 Offline mode

Debe permitir:

- consultas a Scientific Memory;
- ejecución reproducible con artefactos pinned;
- uso de snapshots previos.

Debe indicar claramente cuando la frescura no pudo verificarse.

---

# 26. Membranas como dominio de primera clase

SimForge tiene un foco fuerte en sistemas de membrana; por tanto SSKI debe permitir preguntas como:

```text
¿Esta proteína es transmembranal?
¿Cuántas hélices TM tiene?
¿Cuál es su orientación?
¿En qué membrana biológica reside?
¿Qué lípidos predominan en ese contexto?
¿Qué estructuras experimentales existen?
¿Hay modelos previamente insertados en bicapa?
¿Qué oligomerización se reporta?
¿Qué composición es biológicamente defendible?
```

## 26.1 Membrane Evidence Resolver

Debe combinar, sin confundir:

- anotación de localización;
- predicción TM;
- estructura PDB;
- orientación OPM/PDBTM;
- simulaciones MemProtMD;
- literatura;
- composición lipídica de LIPID MAPS/SwissLipids;
- conocimiento local.

## 26.2 Resultado ejemplo

```yaml
membrane_context:
  status: supported

  topology:
    transmembrane_segments: 5
    evidence: [...]

  orientation:
    source: ...
    confidence: ...

  lipid_environment:
    biological_context_known: partial
    candidate_models:
      - composition: ...
        evidence: ...
```

No generar automáticamente una composición de membrana como si fuera un hecho si sólo existe evidencia indirecta.

---

# 27. Datos de simulación como conocimiento

SimForge debe ingerir sus propios resultados.

Ejemplo:

```text
public knowledge:
ligand X binds target Y
```

SimForge puede añadir:

```yaml
claim:
  subject: ligand_X
  predicate: stable_interaction_in_simulation
  object: target_Y

context:
  run_id: ...
  force_field: ...
  water_model: ...
  temperature: ...
  simulation_time: ...

provenance:
  type: simforge_md
```

Esto jamás se mezcla con una medición experimental; ambos permanecen comparables pero epistemológicamente distintos.

---

# 28. API interna propuesta

```text
GET  /knowledge/entity/{id}
POST /knowledge/resolve
POST /knowledge/query
GET  /knowledge/evidence/{claim_id}
GET  /knowledge/sources/{entity_id}
GET  /knowledge/conflicts/{entity_id}
POST /knowledge/enrich/{entity_id}
GET  /knowledge/freshness/{entity_id}
POST /knowledge/snapshot
GET  /knowledge/snapshot/{id}
```

## 28.1 Query declarativa

```yaml
subject:
  type: protein
  id: UniProtKB:P00533

request:
  - structures
  - ligands
  - domains
  - variants
  - membrane_context

freshness:
  require_current: true
```

---

# 29. CLI propuesta

```bash
simforge knowledge resolve P00533
simforge knowledge resolve 1M17
simforge knowledge resolve "eugenol"

simforge knowledge show sf:protein:...
simforge knowledge sources sf:protein:...
simforge knowledge conflicts sf:protein:...
simforge knowledge refresh sf:protein:...

simforge knowledge protein P00533 --structures --ligands --domains
simforge knowledge ligand CHEMBL25 --targets --bioactivity
simforge knowledge membrane P43220

simforge knowledge snapshot create
simforge knowledge snapshot inspect <id>

simforge knowledge status
simforge knowledge source status
```

---

# 30. Modelo de tablas mínimo

## `entity`

```text
id
entity_type
canonical_label
status
created_at
updated_at
```

## `external_identifier`

```text
entity_id
namespace
accession
version
valid_from
valid_to
source_id
```

## `entity_fingerprint`

```text
entity_id
fingerprint_type
fingerprint_value
algorithm_version
```

## `relationship`

```text
subject_id
predicate
object_id
context_json
claim_id
```

## `claim`

```text
id
subject_id
predicate
object_entity_id
value_json
context_json
claim_class
status
```

## `evidence`

```text
id
claim_id
source_id
source_record_id
source_release
publication_id
evidence_type
retrieved_at
payload_hash
```

## `source`

```text
id
name
config_json
license_json
health
last_release
last_checked
```

## `artifact`

```text
id
sha256
media_type
storage_uri
source_id
external_id
cache_policy
```

## `knowledge_snapshot`

```text
id
created_at
resolver_version
ontology_version
manifest_json
```

## `decision_trace`

```text
id
snapshot_id
question
candidates_json
selected_json
policy_version
evidence_ids
human_override
```

---

# 31. Ontología v0 requerida antes de programar

Antes del esquema SQL definitivo deben definirse formalmente:

```text
Entity
BiologicalEntity
ProteinConcept
ProteinSequence
ProteinIsoform
ProteinVariant
MacromolecularComplex
ChemicalEntity
ChemicalSpecies
Lipid
Ion
Cofactor
StructuralEntity
ExperimentalStructure
PredictedStructure
Assembly
Chain
Domain
BindingSite
MembraneContext
Reaction
EnzymeActivity
Interaction
Pathway
Claim
Evidence
Source
Snapshot
SimulationEntity
```

Y relaciones básicas:

```text
HAS_SEQUENCE
HAS_ISOFORM
HAS_VARIANT
HAS_STRUCTURE
HAS_ASSEMBLY
HAS_CHAIN
HAS_DOMAIN
HAS_BINDING_SITE
BINDS
INHIBITS
ACTIVATES
CATALYZES
PARTICIPATES_IN
INTERACTS_WITH
REQUIRES_COFACTOR
LOCATED_IN
MEMBER_OF
REPRESENTS
DERIVED_FROM
MAPS_TO
PARAMETERIZED_AS
SUPPORTED_BY
CONTRADICTS
SUPERSEDES
```

---

# 32. Directorio propuesto dentro de SimForge

```text
simforge/
  knowledge/
    __init__.py

    ontology/
      entities.py
      predicates.py
      contexts.py
      versions.py

    registry/
      entities.py
      identifiers.py
      fingerprints.py

    evidence/
      claims.py
      evidence.py
      provenance.py
      conflicts.py

    federation/
      engine.py
      routing.py
      source_registry.py
      freshness.py
      health.py

    sources/
      uniprot/
      rcsb/
      alphafold/
      modelarchive/
      interpro/
      gene_ontology/
      pubchem/
      chebi/
      unichem/
      chembl/
      bindingdb/
      rhea/
      reactome/
      lipidmaps/
      swisslipids/
      memprotmd/
      clinvar/
      ncbi/

    resolution/
      entity_resolver.py
      protein_resolver.py
      chemical_resolver.py
      structure_resolver.py
      membrane_resolver.py

    memory/
      public.py
      local.py
      project.py
      run.py

    storage/
      postgres.py
      parquet.py
      artifacts.py
      vectors.py

    snapshots/
      models.py
      builder.py
      restore.py

    api/
      models.py
      service.py

    cli/
      knowledge.py

  core/
    decision_engine.py
```

Dependencia estricta:

```text
sources -> normalization -> knowledge core
knowledge core -> KnowledgeResolver
KnowledgeResolver -> DecisionContext
DecisionContext -> DecisionEngine
```

Nunca al revés.

---

# 33. Despliegue local-first recomendado

## Perfil individual

```text
SimForge process
PostgreSQL
local content-addressable artifact directory
Parquet store
optional pgvector
```

No requiere Kubernetes.

## Perfil laboratorio

```text
SimForge API
PostgreSQL
MinIO/S3
worker pool
scheduler
optional graph/search projection
```

## Perfil institucional/HPC

```text
Knowledge API
PostgreSQL HA
S3/object storage
queue/workers
HPC-aware artifact cache
replicated read nodes
```

La semántica debe ser idéntica en los tres perfiles.

---

# 34. Estrategia de implementación por fases

## Fase 0 — Especificación, sin código productivo

Entregables:

1. `SCIENTIFIC_KNOWLEDGE_MODEL_V0.md`
2. `IDENTITY_AND_MAPPING_V0.md`
3. `EVIDENCE_AND_PROVENANCE_V0.md`
4. `SOURCE_ADAPTER_CONTRACT_V0.md`
5. `KNOWLEDGE_SNAPSHOT_V0.md`

## Fase 1 — Core de identidad

Implementar:

- Entity;
- ExternalIdentifier;
- Fingerprint;
- Source;
- Provenance;
- resolver exacto;
- PostgreSQL.

Fuentes iniciales:

```text
UniProt
RCSB PDB
PubChem
ChEBI
UniChem
```

## Fase 2 — Evidencia y estructura

Agregar:

```text
ChEMBL
InterPro
GO
AlphaFold DB
ModelArchive
Rhea
```

Implementar:

- Claim;
- Evidence;
- conflicts;
- snapshots.

## Fase 3 — Membranas

Agregar:

```text
LIPID MAPS
SwissLipids
MemProtMD
OPM/PDBTM-compatible sources
TCDB
```

Crear `MembraneResolver`.

## Fase 4 — Bioactividad y pathways

Agregar:

```text
BindingDB
Reactome
IntAct
Complex Portal
BioGRID
IUPHAR
```

## Fase 5 — Variantes y literatura

Agregar:

```text
ClinVar
dbSNP
gnomAD
PubMed
Europe PMC
```

## Fase 6 — Scientific Memory + SimForge results

Ingerir:

- docking;
- MD;
- MM/PBSA;
- umbrella/PMF;
- FEP/TI;
- QM;
- QM/MM;
- FEL;
- experimental data supplied by user.

## Fase 7 — Learned retrieval and graph intelligence

Agregar:

- embeddings;
- semantic candidate generation;
- automatic source routing;
- graph traversal optimization;
- anomaly detection;
- literature claim extraction asistida.

---

# 35. Primer vertical slice recomendado

No comenzar intentando integrar 30 fuentes.

Primer sistema end-to-end:

```text
Input
  P00533

      ↓

UniProt identity
      ↓
PDB structures
      ↓
PDB CCD ligands
      ↓
ChEMBL/PubChem chemical resolution
      ↓
InterPro domains
      ↓
Rhea/Reactome function
      ↓
KnowledgeResponse
      ↓
Snapshot
```

Debe demostrar:

- identidad;
- mapping;
- evidencia;
- source versions;
- actualización;
- cache;
- snapshot;
- DecisionContext.

Sólo después ampliar dominios.

---

# 36. Pruebas científicas obligatorias

## 36.1 Protein identity regression

- misma secuencia con IDs distintos;
- isoformas;
- accession obsoleto;
- mutante PDB;
- constructo truncado;
- cadenas repetidas.

## 36.2 Chemical identity regression

- sal vs parent compound;
- protonation state;
- tautómero;
- estereoisómero;
- molécula sin estereo definido;
- PDB CCD vs PubChem vs ChEMBL.

## 36.3 Structure mapping

- UniProt ↔ PDB residue numbering;
- missing residues;
- engineered mutations;
- insertion codes;
- multiple chains;
- biological assemblies.

## 36.4 Evidence conflict

Dos fuentes reportan valores incompatibles.

Esperado:

```text
status = conflicting
```

No overwrite.

## 36.5 Update regression

Fuente publica nueva release.

Esperado:

- nuevo release detectado;
- delta procesado;
- resoluciones afectadas invalidadas;
- snapshot anterior permanece reproducible.

---

# 37. Source priority no equivale a truth priority

No utilizar una jerarquía global tipo:

```text
PDB > AlphaFold > everything else
```

La prioridad depende de la pregunta.

Ejemplo:

```text
Question: experimental coordinates
PDB = appropriate authoritative source
AlphaFold = not experimental evidence
```

Pero:

```text
Question: full-length structural hypothesis
AlphaFold may have broader sequence coverage
```

El Knowledge Resolver evalúa **fitness-for-question**.

---

# 38. Política de “latest knowledge”

Una respuesta sólo puede llamarse “actual” si:

1. se conoce la versión/release de las fuentes relevantes;
2. la fuente fue revisada dentro de su freshness policy;
3. no existe una sincronización pendiente conocida;
4. se reporta cuándo fue verificada.

Ejemplo:

```yaml
knowledge_status:
  freshness: current
  checked_at: 2026-09-23T...
  sources:
    - source: RCSB
      release: ...
```

Si no pudo verificarse:

```text
freshness = unknown
```

Nunca asumir actualidad por cache silenciosamente.

---

# 39. Política de almacenamiento mínimo

Persistir siempre:

```text
entity identity
external IDs
fingerprints
claims used
provenance
source release
source record IDs
artifact hashes
snapshot manifests
decision traces
```

Persistir condicionalmente:

```text
normalized detailed payloads
frequently used records
active project entities
```

Evictable:

```text
large downloaded artifacts
unused API payloads
bulk staging files
```

Pinned:

```text
artifacts that directly contributed to reproducible scientific output
```

---

# 40. Protocolo de actualización por tipo de fuente

## API con release/version endpoint

```text
poll version
if unchanged: stop
if changed: incremental/full metadata sync
```

## API sin releases

```text
ETag / Last-Modified / record timestamps
```

## Bulk release

```text
manifest -> checksum -> staging -> delta -> normalized commit -> discard staging
```

## Continuously changing source

```text
hot entities -> refresh schedule
cold entities -> lazy refresh
```

---

# 41. Compatibilidad FAIR

Objetivo:

- Findable;
- Accessible;
- Interoperable;
- Reusable.

Prácticas:

- persistent internal IDs;
- machine-readable provenance;
- versioned ontologies;
- source attribution;
- standard identifiers;
- JSON/JSON-LD export;
- explicit units;
- content hashes;
- immutable snapshots.

Se puede añadir posteriormente export PROV-O/JSON-LD sin hacer de RDF el storage primario.

---

# 42. Principios de diseño finales

### P1 — Identity before data

Resolver qué entidad es antes de acumular atributos.

### P2 — Evidence before truth

Toda afirmación conoce su procedencia.

### P3 — Context before comparison

No comparar medidas científicas sin condiciones compatibles.

### P4 — Federation before replication

Acceder al ecosistema antes de copiarlo.

### P5 — Cache what is useful

Persistir lo reutilizable; recuperar el resto.

### P6 — Pin what affects reproducibility

Si influyó en un resultado científico, conservar hash, versión y preferentemente artefacto.

### P7 — Determinism before AI

IDs, hashes, mappings y estructuras químicas se resuelven determinísticamente siempre que sea posible.

### P8 — AI proposes; evidence decides

AI genera candidatos e interpretaciones, no hechos.

### P9 — Contradictions are data

No esconderlas.

### P10 — External knowledge never becomes automatic truth

El Decision Engine recibe evidencia y contexto, no valores desnudos.

### P11 — Local knowledge remains namespaced

Resultados propios no se presentan como conocimiento universal.

### P12 — Everything important is versioned

Ontology, adapters, parsers, source releases, resolver policies y decision policies.

---

# 43. Fuentes oficiales verificadas durante este diseño

Las siguientes fuentes fueron verificadas como ejemplos de acceso programático, releases o descarga al diseñar este protocolo el 2026-09-23.

## Estructuras

- RCSB PDB Web APIs Overview:  
  https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- RCSB PDB File Download Services:  
  https://www.rcsb.org/docs/programmatic-access/file-download-services
- AlphaFold Protein Structure Database FAQ / bulk downloads:  
  https://alphafold.ebi.ac.uk/faq
- ModelArchive:  
  https://modelarchive.org/

## Proteínas y dominios

- UniProt REST-linked records/help:  
  https://www.uniprot.org/help/linking_to_uniprot
- InterPro API documentation:  
  https://interpro-documentation.readthedocs.io/en/latest/api.html

## Ontologías

- Gene Ontology downloads:  
  https://geneontology.org/docs/download-ontology/
- Gene Ontology annotations/releases:  
  https://geneontology.org/docs/download-go-annotations/

## Química

- PubChem PUG REST:  
  https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- PubChem downloads:  
  https://pubchem.ncbi.nlm.nih.gov/docs/downloads
- ChEBI REST APIs:  
  https://www.ebi.ac.uk/chebi/tools
- ChEBI downloads:  
  https://www.ebi.ac.uk/chebi/downloads
- UniChem:  
  https://www.ebi.ac.uk/unichem/
- ChEMBL Web Services documentation:  
  https://github.com/chembl/GLaDOS-docs/blob/master/web-services/chembl-data-web-services.md

## Bioactividad

- BindingDB downloads:  
  https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp?all_download=yes

## Reacciones y pathways

- Reactome downloads:  
  https://reactome.org/download-data
- Reactome Content Service:  
  https://reactome.org/dev/content-service

## Lípidos y membranas

- LIPID MAPS REST:  
  https://www.lipidmaps.org/resources/rest
- LIPID MAPS downloads:  
  https://www.lipidmaps.org/databases/lmsd/download
- SwissLipids programmatic downloads:  
  https://www.swisslipids.org/api/index.php/downloadData
- MemProtMD:  
  https://memprotmd.bioch.ox.ac.uk/home/
- MemProtMD API:  
  https://memprotmd.bioch.ox.ac.uk/api/

## Variantes

- ClinVar programmatic access:  
  https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/
- ClinVar downloads/release cycle:  
  https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/

---

# 44. Decisión arquitectónica recomendada

Adoptar formalmente la siguiente arquitectura como objetivo:

```text
                     SimForge Scientific Knowledge Infrastructure

                  +----------------------------------------------+
                  |             Universal Entity Registry        |
                  +----------------------+-----------------------+
                                         |
                      +------------------+------------------+
                      |                                     |
              Federation Engine                    Scientific Memory
                      |                                     |
             Source Adapter Mesh                 Claims / Evidence
                      |                                     |
        External Scientific Databases             Local Results
                      |                                     |
                      +------------------+------------------+
                                         |
                                 Knowledge Resolver
                                         |
                                  Decision Context
                                         |
                                   Decision Engine
                                         |
                                      SimForge
```

La infraestructura no intenta contener físicamente todo el conocimiento molecular existente. Mantiene una **red universal de identidad y evidencia** suficientemente compacta para reconocer entidades globalmente y utiliza federación + caché + snapshots para tener acceso efectivo, actualizado y reproducible al resto.

---

# 45. Próximo documento

El siguiente paso recomendado no es programar adaptadores todavía.

Crear:

```text
SCIENTIFIC_KNOWLEDGE_MODEL_V0.md
```

Debe fijar formalmente:

1. entidades nucleares;
2. invariantes de identidad;
3. relaciones permitidas;
4. semántica de `Claim`;
5. semántica de `Evidence`;
6. contextos;
7. contradicciones;
8. namespaces;
9. versionado;
10. mappings UniProt ↔ PDB ↔ chemical entities ↔ topology.

Después:

```text
IDENTITY_AND_MAPPING_V0.md
EVIDENCE_AND_PROVENANCE_V0.md
SOURCE_ADAPTER_CONTRACT_V0.md
```

Sólo entonces conviene congelar el esquema SQL y comenzar el primer vertical slice.

---

## Conclusión

SSKI debe comportarse menos como una “base de datos enorme” y más como un **sistema operativo del conocimiento molecular**.

Su unidad primaria no será el archivo descargado. Será la combinación:

```text
ENTITY
  +
IDENTITY
  +
CLAIM
  +
EVIDENCE
  +
PROVENANCE
  +
VERSION
```

Los archivos, APIs, modelos y bases externas son proveedores reemplazables de esa información.

Con este diseño, SimForge puede crecer desde unas pocas entidades locales hasta un índice de escala global sin necesidad de replicar físicamente todo el ecosistema científico y sin sacrificar actualización, trazabilidad ni reproducibilidad.
