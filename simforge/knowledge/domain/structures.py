"""Experimentally classified source records with an explicit extension path."""
from typing import Annotated, Literal
import re
from pydantic import Field, field_validator, model_validator
from .entities import Entity
from .identity import ExternalIdentifier, _ClaimRef, _EvidenceRef, _ImportProvenance, _consistent_types, _Ref
from .artifacts import ArtifactRef
from .serialization import _DomainModel, _Items, _Text, _Decimal, _Object, _nfc

_Categories = Literal['X_RAY_DIFFRACTION', 'ELECTRON_MICROSCOPY', 'SOLUTION_NMR', 'SOLID_STATE_NMR', 'NEUTRON_DIFFRACTION', 'ELECTRON_DIFFRACTION', 'FIBER_DIFFRACTION', 'OTHER_EXPERIMENTAL']
_ALIASES = {
    'x ray diffraction': 'X_RAY_DIFFRACTION', 'xray diffraction': 'X_RAY_DIFFRACTION',
    'electron microscopy': 'ELECTRON_MICROSCOPY', 'cryo em': 'ELECTRON_MICROSCOPY',
    'cryo electron microscopy': 'ELECTRON_MICROSCOPY', 'solution nmr': 'SOLUTION_NMR',
    'solid state nmr': 'SOLID_STATE_NMR', 'neutron diffraction': 'NEUTRON_DIFFRACTION',
    'electron diffraction': 'ELECTRON_DIFFRACTION', 'fiber diffraction': 'FIBER_DIFFRACTION',
    'fibre diffraction': 'FIBER_DIFFRACTION', 'other experimental': 'OTHER_EXPERIMENTAL',
}


def _label(value):
    return ' '.join(re.sub(r'[\W_]+', ' ', _nfc(value).casefold()).split())


def _contradiction(label):
    """Classify explicit contradictions only; absence is not experimental proof.

    Input consists of case-folded, punctuation-separated lexical tokens. These
    bounded phrase classes intentionally do not infer general natural language.
    """
    classes = (
        ('negated experimental origin',
         r'\b(?:non ?experimental|(?:not|no)(?: an?| any)? experimental)\b'),
        ('computational or simulated generation',
         r'\b(?:computational|simulated|simulation|in silico|'
         r'computer (?:generated|derived|modeled|modelled)|model derived)\b'),
        ('prediction',
         r'\b(?:predicted|predictions?|alphafold(?:v?[0-9]+)?|alpha fold(?: [0-9]+)?|machine learning)\b'),
        ('theoretical or comparative modelling',
         r'\b(?:theoretical|homology|comparative model(?:ing|ling)?)\b'),
        ('molecular simulation',
         r'\b(?:molecular dynamics|md (?:derived|simulation|model|structure))\b'),
    )
    for reason, pattern in classes:
        if re.search(pattern, label):
            return reason
    return None


class _ExperimentalMethod(_DomainModel):
    category: _Categories
    raw_method: _Text
    declared_origin: Literal['EXPERIMENTAL', 'COMPUTATIONAL', 'PREDICTED', 'THEORETICAL', 'UNKNOWN']
    classification_basis: Literal['SOURCE_DECLARATION', 'CURATED_ASSERTION', 'USER_ASSERTION']
    classification_claim: _ClaimRef
    provenance: _ImportProvenance

    @field_validator('category', mode='before')
    @classmethod
    def _category(cls, value):
        return _ALIASES.get(_label(value), value) if isinstance(value, str) else value

    @model_validator(mode='after')
    def _experimental(self):
        if self.declared_origin != 'EXPERIMENTAL':
            raise ValueError('explicit experimental origin is required')
        label = _label(self.raw_method)
        contradiction = _contradiction(label)
        if contradiction is not None:
            raise ValueError('method contradicts experimental origin: ' + contradiction)
        # Ambiguity is separate from contradiction and remains fail-closed.
        ambiguous = (
            r'\b(unknown|unspecified|ambiguous|unclassified|unclear|undetermined|unrecognized|'
            r'possibly|maybe|uncertain|not (?:known|determined|specified))\b'
        )
        if ('?' in self.raw_method or not label
                or label in {'na', 'n a', 'other', 'other experimental', 'model', 'method'}
                or re.search(ambiguous, label)):
            raise ValueError('unknown or ambiguous method cannot be experimental')
        known = _ALIASES.get(label)
        if known is not None and self.category != known:
            raise ValueError('experimental method category disagrees with source label')
        if known is None:
            if self.category != 'OTHER_EXPERIMENTAL':
                raise ValueError('unrecognized experimental method requires OTHER_EXPERIMENTAL')
            if self.classification_basis != 'SOURCE_DECLARATION':
                raise ValueError('experimental extension requires a source declaration')
        return self


class ExperimentalStructure(Entity):
    entity_type: Literal['experimental_structure'] = 'experimental_structure'
    origin: Literal['experimental']
    source_identifier: ExternalIdentifier
    method: _ExperimentalMethod
    artifact: ArtifactRef
    evidence_refs: Annotated[_Items[_EvidenceRef], Field(min_length=1)]
    resolution_angstrom: Annotated[_Decimal, Field(gt=0)] | None = None
    quality_metadata: _Object = Field(default_factory=dict)
    _unordered_fields = ('evidence_refs',)

    @model_validator(mode='after')
    def _source(self):
        source = self.source_identifier
        target = _Ref(entity_id=self.sf_id, entity_revision=self.entity_revision, entity_type=self.entity_type)
        if source.entity != target or source.mapping_status != 'CONFIRMED_MAPPING':
            raise ValueError('confirmed source mapping must target this structure revision')
        if source.provenance.source_record.version_status != 'PINNED':
            raise ValueError('experimental source version must be pinned')
        _consistent_types((target, self.method.classification_claim, *self.evidence_refs))
        return self
