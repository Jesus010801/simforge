# ADR-004 — Content-Addressed Scientific Artifacts

**Status:** Accepted

## Decision
Physical artifacts used by a KnowledgeBundle must have cryptographic content hashes. SHA256 is the initial canonical digest.

## Consequences
Artifact identity is verifiable independently of filenames or mutable upstream URLs.
