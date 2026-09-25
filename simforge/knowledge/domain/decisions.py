"""Recorded scientific resolution, without evaluation or execution."""
from typing import Annotated, Literal
from pydantic import Field, model_validator
from .identity import EntityId, _Ref, _ClaimRef, _EvidenceRef, _consistent_types
from .artifacts import ArtifactRef, _artifact_key, _artifact_consistency
from .serialization import _DomainModel, _Text, _Token, _Symbol, _PositiveInt, _SHA256, _Items, _Decimal, _Structured


class PolicyRef(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.policy_id.value, 'policy'),)
    policy_id: EntityId
    policy_version: _Token
    definition_sha256: _SHA256


class _LiteralResult(_DomainModel):
    kind: Literal['literal']
    value: _Structured


class _ResolutionResult(_DomainModel):
    selected_entities: _Items[_Ref] = ()
    selected_artifacts: _Items[ArtifactRef] = ()
    value: _LiteralResult | None = None
    _unordered_fields = ('selected_entities', 'selected_artifacts')

    @model_validator(mode='after')
    def _present(self):
        if not (self.selected_entities or self.selected_artifacts or self.value is not None):
            raise ValueError('a resolution result cannot be empty')
        return self


class _BoundEntity(_DomainModel):
    name: _Token
    kind: Literal['entity']
    value: _Ref


class _BoundArtifact(_DomainModel):
    name: _Token
    kind: Literal['artifact']
    value: ArtifactRef


class _ScientificParameters(_DomainModel):
    """Closed Phase 1 scientific conditions, not a generic configuration bag.

    Units are encoded in field names. pH is not artificially limited to 0..14.
    New scientific parameters require a reviewed schema change. There are no
    string, mapping, list, executable, transport, or backend-valued fields.
    """
    ph: _Decimal | None = None
    temperature_kelvin: Annotated[_Decimal, Field(gt=0)] | None = None
    ionic_strength_molar: Annotated[_Decimal, Field(ge=0)] | None = None
    tolerance: Annotated[_Decimal, Field(gt=0)] | None = None
    preserve_heavy_atom_connectivity: bool | None = None
    preserve_stereochemistry: bool | None = None


class _BoundLiteral(_DomainModel):
    name: _Token
    kind: Literal['literal']
    value: _ScientificParameters

    @model_validator(mode='after')
    def _concrete(self):
        if all(getattr(self.value, name) is None for name in type(self.value).model_fields):
            raise ValueError('literal binding must specify a scientific condition')
        return self


class _Output(_DomainModel):
    name: _Token
    scientific_type: _Symbol


class _DeferredOperation(_DomainModel):
    operation_key: _Token
    target_stage: Literal['preparation']
    action: Literal['science:protonation', 'science:reconstruction', 'science:tautomer_resolution']
    method_id: _Symbol
    method_version: _Token
    method_sha256: _SHA256
    required_inputs: Annotated[_Items[_BoundEntity | _BoundArtifact | _BoundLiteral], Field(min_length=1)]
    parameters: _ScientificParameters = Field(default_factory=_ScientificParameters)
    expected_outputs: Annotated[_Items[_Output], Field(min_length=1)]
    success_criteria: Annotated[_Items[_Token], Field(min_length=1)]
    on_failure: Literal['ABORT']
    _unordered_fields = ('required_inputs', 'expected_outputs', 'success_criteria')

    @model_validator(mode='after')
    def _names(self):
        for items in (self.required_inputs, self.expected_outputs):
            names = [x.name for x in items]
            if len(set(names)) != len(names):
                raise ValueError('duplicate operation binding name')
        # A pinned scientific method identifier is data, never an execution URI.
        if not self.method_id.startswith('science:'):
            raise ValueError('method must belong to the scientific identifier namespace')
        # Preserve the existing rejection of placeholder/execution method labels.
        # This identifier check is separate from the closed parameter schema.
        parts = self.method_id.lower().replace('-', '_').split(':')
        if any(p in {'unknown', 'manual_review', 'ask_later', 'http', 'https', 'network',
                     'federation', 'gromacs', 'openmm', 'amber', 'shell', 'exec'} for p in parts):
            raise ValueError('method must identify scientific resolution, not execution or lookup')
        return self


def _ref_key(ref):
    return ref.entity_id.value, ref.entity_revision, ref.entity_type


def _require_members(members, inventory, key):
    available = {key(v) for v in inventory}
    if any(key(v) not in available for v in members):
        raise ValueError('reference missing from declared input inventory')


class Decision(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.decision_id.value, 'decision'),)
    decision_id: EntityId
    decision_revision: _PositiveInt
    question: _Text
    status: Literal['RESOLVED', 'RESOLVED_WITH_WARNING', 'DEFERRED', 'ABORT']
    policy: PolicyRef
    mandatory: bool = True
    entity_inputs: _Items[_Ref] = ()
    artifact_inputs: _Items[ArtifactRef] = ()
    claim_inputs: _Items[_ClaimRef] = ()
    evidence_inputs: _Items[_EvidenceRef] = ()
    parameters: _ScientificParameters = Field(default_factory=_ScientificParameters)
    result: _ResolutionResult | None = None
    deferred_operation: _DeferredOperation | None = None
    warnings: _Items[_Text] = ()
    rationale: _Text
    _unordered_fields = ('entity_inputs', 'artifact_inputs', 'claim_inputs', 'evidence_inputs')

    @model_validator(mode='after')
    def _state(self):
        resolved = self.status in ('RESOLVED', 'RESOLVED_WITH_WARNING')
        if (self.result is not None) != resolved:
            raise ValueError('decision state/result mismatch')
        if (self.deferred_operation is not None) != (self.status == 'DEFERRED'):
            raise ValueError('decision state/deferred operation mismatch')
        if self.status == 'RESOLVED' and self.warnings:
            raise ValueError('resolved warnings require RESOLVED_WITH_WARNING')
        if self.status == 'RESOLVED_WITH_WARNING' and not self.warnings:
            raise ValueError('warning state requires warnings')
        _consistent_types((*self.entity_inputs, *self.claim_inputs, *self.evidence_inputs))
        _artifact_consistency(self.artifact_inputs)
        if self.result is not None:
            _require_members(self.result.selected_entities, self.entity_inputs, _ref_key)
            _require_members(self.result.selected_artifacts, self.artifact_inputs, _artifact_key)
        if self.deferred_operation is not None:
            for binding in self.deferred_operation.required_inputs:
                if isinstance(binding, _BoundEntity):
                    _require_members((binding.value,), self.entity_inputs, _ref_key)
                elif isinstance(binding, _BoundArtifact):
                    _require_members((binding.value,), self.artifact_inputs, _artifact_key)
        return self


class DecisionContext(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.context_id.value, 'context'),)
    context_id: EntityId
    context_revision: _PositiveInt
    schema_version: Literal['0.1.0'] = '0.1.0'
    request_sha256: _SHA256
    entity_refs: _Items[_Ref] = ()
    artifact_refs: _Items[ArtifactRef] = ()
    claim_refs: _Items[_ClaimRef] = ()
    evidence_refs: _Items[_EvidenceRef] = ()
    decisions: _Items[Decision] = ()
    context_warnings: _Items[_Text] = ()
    _unordered_fields = ('entity_refs', 'artifact_refs', 'claim_refs', 'evidence_refs', 'decisions')

    @model_validator(mode='after')
    def _closure(self):
        _consistent_types((*self.entity_refs, *self.claim_refs, *self.evidence_refs))
        _artifact_consistency(self.artifact_refs)
        keys = [(d.decision_id.value, d.decision_revision) for d in self.decisions]
        if len(set(keys)) != len(keys):
            raise ValueError('duplicate decision identity/revision')
        operations = [d.deferred_operation.operation_key for d in self.decisions if d.deferred_operation is not None]
        if len(set(operations)) != len(operations):
            raise ValueError('duplicate deferred operation key')
        for decision in self.decisions:
            for members, inventory, key in (
                (decision.entity_inputs, self.entity_refs, _ref_key),
                (decision.artifact_inputs, self.artifact_refs, _artifact_key),
                (decision.claim_inputs, self.claim_refs, _ref_key),
                (decision.evidence_inputs, self.evidence_refs, _ref_key),
            ):
                _require_members(members, inventory, key)
        return self

    @property
    def resolved_selections(self):
        return tuple(d.result for d in self.decisions if d.result is not None)

    @property
    def warnings(self):
        return self.context_warnings + tuple(w for d in self.decisions for w in d.warnings)

    @property
    def deferred_operations(self):
        return tuple(d.deferred_operation for d in self.decisions if d.deferred_operation is not None)

    @property
    def has_mandatory_abort(self):
        return any(d.mandatory and d.status == 'ABORT' for d in self.decisions)
