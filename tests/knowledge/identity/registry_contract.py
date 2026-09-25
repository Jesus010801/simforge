"""Reusable behavioral assertions. No backend imports or private state access."""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace, FrozenInstanceError
from datetime import datetime, timezone, timedelta
from threading import Barrier
from uuid import UUID
import pytest
from simforge.knowledge.domain import Entity, EntityId, ExternalIdentifier
from simforge.knowledge.identity import (
    RegistryChange, RegistryWrite, ExternalIdentifierKey, ExternalIdentifierMatch, RegistryError,
    EntityNotFoundError, RegistryRevisionNotFoundError, EntityRevisionNotFoundError,
    RevisionConflictError, OperationConflictError, InvalidRegistryOperationError,
)

NOW = datetime(2026, 1, 2, tzinfo=timezone.utc)


def eid(n=1):
    return EntityId(value=f'sf:00000000-0000-4000-8000-{n:012x}')


def descriptor(n=1, revision=1, kind='protein', lifecycle='ACTIVE'):
    return Entity(sf_id=eid(n), entity_type=kind, entity_revision=revision, lifecycle=lifecycle)


def change(n=1, **kwargs):
    return RegistryChange(**dict(operation_id=UUID(f'10000000-0000-4000-8000-{n:012x}'),
                                 actor='curator', reason='source declaration', recorded_at=NOW, **kwargs))


def alias(n=1, revision=1, kind='protein', release='R1', **kwargs):
    data = dict(entity=dict(entity_id=eid(n), entity_revision=revision, entity_type=kind),
                namespace='uniprot', accession='P00533', mapping_status='CANDIDATE_MAPPING',
                lifecycle='ACTIVE', provenance=dict(source_record=dict(source_id='fixture',
                    record_id='P00533', source_release=release, version_status='PINNED'),
                    raw_payload_sha256='a'*64, importer_id='fixture', importer_version='1',
                    imported_at=NOW, extraction_method='fixture'))
    data.update(kwargs)
    return ExternalIdentifier(**data)


def command(n=1, op=1, expected=0, entity=None, aliases=(), **kwargs):
    return RegistryWrite(entity=entity if entity is not None else descriptor(n), aliases=aliases,
                         expected_registry_revision=expected, change=change(op), **kwargs)


def invalid(call, code, path=None):
    with pytest.raises(InvalidRegistryOperationError) as caught:
        call()
    assert caught.value.code == code
    if path is not None:
        assert caught.value.field_path == path
    return caught.value


def race(*calls):
    barrier = Barrier(len(calls))
    def run(call):
        barrier.wait(timeout=10)
        try:
            return call()
        except RegistryError as error:
            return error
    with ThreadPoolExecutor(max_workers=len(calls)) as pool:
        return tuple(pool.map(run, calls))


class RegistryBehavior:
    @pytest.mark.parametrize('kind', ['protein', 'protein_sequence', 'chemical_identity',
        'chemical_species', 'experimental_structure', 'claim', 'evidence'])
    @pytest.mark.parametrize('lifecycle', ['ACTIVE', 'DEPRECATED', 'RETRACTED'])
    def test_registration(self, registry, kind, lifecycle):
        entity = descriptor(kind=kind, lifecycle=lifecycle)
        result = registry.write(command(entity=entity))
        assert result.entity == entity and result.registry_revision == 1
        assert registry.get(eid()) == result
        assert registry.get(eid(), registry_revision=1) == result
        assert registry.get_entity_revision(eid(), 1) == entity
        assert registry.history(eid()) == (result,)

    def test_unknown_and_missing_reads(self, registry):
        for call in (lambda: registry.get(eid()), lambda: registry.history(eid()),
                     lambda: registry.get_entity_revision(eid(), 1),
                     lambda: registry.get(eid(), registry_revision=1)):
            with pytest.raises(EntityNotFoundError) as error:
                call()
            assert error.value.entity_id == eid()
        registry.write(command())
        for call, cls, field in ((lambda: registry.get(eid(), registry_revision=2), RegistryRevisionNotFoundError, 'registry_revision'),
                                 (lambda: registry.get_entity_revision(eid(), 2), EntityRevisionNotFoundError, 'entity_revision')):
            with pytest.raises(cls) as error:
                call()
            assert error.value.entity_id == eid() and getattr(error.value, field) == 2
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'absent')) == ()

    @pytest.mark.parametrize('value', [True, False, 1.0, '1', -1, 0, None])
    def test_read_revision_ingress(self, registry, value):
        invalid(lambda: registry.get_entity_revision(eid(), value), 'INVALID_ARGUMENT', ('entity_revision',))
        if value is not None:
            invalid(lambda: registry.get(eid(), registry_revision=value), 'INVALID_ARGUMENT', ('registry_revision',))

    @pytest.mark.parametrize('value', [None, 'sf:00000000-0000-4000-8000-000000000001', {}, 1])
    def test_read_identifier_ingress(self, registry, value):
        for call in (lambda: registry.get(value), lambda: registry.history(value),
                     lambda: registry.get_entity_revision(value, 1)):
            invalid(call, 'INVALID_ARGUMENT', ('entity_id',))

    def test_cas_and_history(self, registry):
        records = [registry.write(command())]
        for i in range(2, 5):
            records.append(registry.write(command(op=i, expected=i-1)))
        winner = registry.write(command(op=5, expected=4))
        with pytest.raises(RevisionConflictError) as caught:
            registry.write(command(op=6, expected=4))
        assert (caught.value.entity_id, caught.value.expected, caught.value.actual) == (eid(), 4, 5)
        assert registry.history(eid()) == (*records, winner)
        with pytest.raises(RevisionConflictError) as caught:
            registry.write(command(op=7))
        assert (caught.value.expected, caught.value.actual) == (0, 5)
        with pytest.raises(RevisionConflictError) as caught:
            registry.write(command(n=2, op=8, expected=1))
        assert caught.value.actual == 0

    @pytest.mark.parametrize('revision,kind,lifecycle,code', [
        (3, 'protein', 'ACTIVE', 'ENTITY_REVISION_SEQUENCE'),
        (1, 'chemical_identity', 'ACTIVE', 'ENTITY_TYPE_CHANGE'),
        (1, 'protein', 'RETRACTED', 'ENTITY_REVISION_REWRITE')])
    def test_invalid_transitions(self, registry, revision, kind, lifecycle, code):
        first = registry.write(command())
        error = invalid(lambda: registry.write(command(op=2, expected=1,
            entity=descriptor(revision=revision, kind=kind, lifecycle=lifecycle))), code)
        assert error.entity_id == eid() and registry.history(eid()) == (first,)

    def test_sequence_and_failed_id_reuse(self, registry):
        invalid(lambda: registry.write(command(entity=descriptor(revision=2))), 'ENTITY_REVISION_SEQUENCE')
        registry.write(command())
        registry.write(command(op=2, expected=1, entity=descriptor(revision=2)))
        invalid(lambda: registry.write(command(op=3, expected=2)), 'ENTITY_REVISION_SEQUENCE')
        assert registry.get_entity_revision(eid(), 1) == descriptor()
        assert registry.get_entity_revision(eid(), 2) == descriptor(revision=2)

    @pytest.mark.parametrize('before', ['ACTIVE', 'DEPRECATED', 'RETRACTED'])
    @pytest.mark.parametrize('after', ['ACTIVE', 'DEPRECATED', 'RETRACTED'])
    def test_all_lifecycle_transitions(self, registry, before, after):
        first = registry.write(command(entity=descriptor(lifecycle=before)))
        second = registry.write(command(op=2, expected=1, entity=descriptor(revision=2, lifecycle=after)))
        assert registry.get(eid()) == second and registry.get(eid(), registry_revision=1) == first
        assert registry.get_entity_revision(eid(), 1).lifecycle == before

    @pytest.mark.parametrize('delta', [timedelta(0), timedelta(days=-20)])
    def test_time_never_orders_history(self, registry, delta):
        first = registry.write(command())
        second = registry.write(replace(command(op=2, expected=1), change=replace(change(2), recorded_at=NOW+delta)))
        assert registry.history(eid()) == (first, second)
        assert second.registry_revision == 2
        assert second.change.recorded_at == NOW+delta

    def test_shared_state_retry_after_advancement_and_reopen(self, registry_handles):
        a, b, reopen = registry_handles
        original = command(aliases=(alias(),))
        first = a.write(original)
        head = b.write(command(op=2, expected=1))
        for handle in (b, reopen()):
            assert handle.write(original) == first
            assert handle.history(eid()) == (first, head)
            assert handle.get(eid()) == head

    @pytest.mark.parametrize('field', ['actor', 'reason', 'recorded_at', 'expected', 'entity', 'aliases'])
    def test_operation_conflict(self, registry, field):
        original = command()
        first = registry.write(original)
        if field in ('actor', 'reason', 'recorded_at'):
            value = NOW+timedelta(days=1) if field == 'recorded_at' else 'changed'
            other = replace(original, change=replace(original.change, **{field: value}))
        else:
            other = replace(original, **{'expected_registry_revision': 1} if field == 'expected' else
                            {'entity': descriptor(2)} if field == 'entity' else {'aliases': (alias(),)})
        with pytest.raises(OperationConflictError) as caught:
            registry.write(other)
        assert caught.value.operation_id == original.change.operation_id
        assert registry.history(eid()) == (first,)

    def test_normalized_retry(self, registry):
        aliases = (alias(release='R2'), alias())
        original = command(aliases=aliases)
        first = registry.write(original)
        retry = replace(original, aliases=tuple(reversed(aliases)), change=replace(original.change,
            recorded_at=NOW.astimezone(timezone(timedelta(hours=-6)))))
        assert registry.write(retry) == first
        assert registry.history(eid()) == (first,)

    def test_independent_logical_registries(self, registry_factory):
        a, b = registry_factory(), registry_factory()
        a.write(command())
        other = b.write(command(entity=descriptor(2)))  # same operation ID, different state
        assert other.registry_revision == 1
        with pytest.raises(EntityNotFoundError):
            a.get(eid(2))
        with pytest.raises(EntityNotFoundError):
            b.get(eid())

    def test_defensive_ownership(self, registry):
        aliases = [alias()]
        cmd = command(aliases=aliases)
        aliases.clear()
        first = registry.write(cmd)
        assert len(first.aliases) == 1
        with pytest.raises((FrozenInstanceError, AttributeError)):
            first.registry_revision = 22
        with pytest.raises(ValueError):
            first.aliases[0].lifecycle = 'RETRACTED'
        detached = first.aliases[0].to_plain()
        detached['provenance']['source_record']['source_release'] = 'changed'
        assert registry.get(eid()).aliases[0].provenance.source_record.source_release == 'R1'
        history = list(registry.history(eid()))
        history.clear()
        matches = list(registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533')))
        matches.clear()
        assert registry.history(eid()) == (first,)
        assert len(registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))) == 1


    def test_retry_after_unrelated_entity_write(self, registry_handles):
        a, b, reopen = registry_handles
        original = command(aliases=(alias(),))
        first = a.write(original)
        unrelated = b.write(command(n=2, op=2, aliases=(alias(2),)))
        assert reopen().write(original) == first
        assert a.history(eid()) == (first,)
        assert a.get(eid(2)) == unrelated
        assert len(a.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))) == 2

    def test_failed_cas_does_not_reserve_operation(self, registry):
        first = registry.write(command())
        with pytest.raises(RevisionConflictError):
            registry.write(command(op=2, expected=0))
        result = registry.write(command(op=2, expected=1))
        assert registry.history(eid()) == (first, result)

    def test_alias_failure_replaces_nothing(self, registry):
        first = registry.write(command(aliases=(alias(),)))
        invalid(lambda: registry.write(command(op=2, expected=1,
            entity=descriptor(revision=2), aliases=(alias(release='R2'), alias(2)))),
            'ALIAS_TARGET_MISMATCH')
        assert registry.history(eid()) == (first,)
        assert registry.get(eid()) == first
        with pytest.raises(EntityRevisionNotFoundError):
            registry.get_entity_revision(eid(), 2)
        matches = registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))
        assert len(matches) == 1 and matches[0].alias == alias()
        assert matches[0].registry_revision == 1
        assert registry.write(command(op=2, expected=1)).registry_revision == 2


class AliasBehavior:
    @pytest.mark.parametrize('bad,code', [(alias(2), 'ALIAS_TARGET_MISMATCH'),
        (alias(kind='chemical_identity'), 'ALIAS_TARGET_MISMATCH'),
        (alias(revision=2), 'ALIAS_TARGET_REVISION_MISSING')])
    def test_atomic_invalid_last_alias(self, registry, bad, code):
        invalid(lambda: registry.write(command(aliases=(alias(), bad))), code)
        with pytest.raises(EntityNotFoundError):
            registry.get(eid())
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533')) == ()
        assert registry.write(command()).registry_revision == 1

    def test_duplicate_alias_and_rollback(self, registry):
        first = registry.write(command(aliases=(alias(),)))
        invalid(lambda: registry.write(command(op=2, expected=1, aliases=(alias(), alias()))), 'DUPLICATE_ALIAS_DECLARATION')
        assert registry.get(eid()) == first
        assert registry.write(command(op=2, expected=1)).registry_revision == 2

    def test_many_to_many_and_order(self, registry):
        for n in (3, 1, 2):
            declarations = (alias(n, release='R2'), alias(n), alias(n, mapping_status='CONFIRMED_MAPPING',
                confirmation_rule=dict(rule_id='rule', version='1', definition_sha256='b'*64)))
            result = registry.write(command(n=n, op=n, aliases=declarations))
            assert result.aliases == tuple(sorted(declarations, key=lambda a: a.to_canonical_bytes()))
        matches = registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))
        assert len(matches) == 9
        assert [(x.entity.sf_id.value, x.alias.to_canonical_bytes()) for x in matches] == sorted(
            (x.entity.sf_id.value, x.alias.to_canonical_bytes()) for x in matches)
        assert {x.entity.sf_id for x in matches} == {eid(1), eid(2), eid(3)}

    def test_administrative_withdrawal(self, registry):
        original = alias()
        first = registry.write(command(aliases=(original,)))
        second = registry.write(command(op=2, expected=1))
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533')) == ()
        assert first.aliases[0] == original and first.aliases[0].lifecycle == 'ACTIVE'
        assert registry.get(eid(), registry_revision=1).aliases == (original,)
        assert second.entity == first.entity and second.aliases == ()
        assert second.change == change(2)
        retracted = alias(lifecycle='RETRACTED', release='R2')
        third = registry.write(command(op=3, expected=2, aliases=(retracted,)))
        assert third.aliases[0].provenance.source_record.source_release == 'R2'
        assert first.aliases[0].provenance.source_record.source_release == 'R1'

    @pytest.mark.parametrize('assignment_first', [False, True])
    def test_reassignment_intermediate_states(self, registry, assignment_first):
        first = registry.write(command(aliases=(alias(),)))
        withdraw = command(op=2, expected=1)
        assign = command(n=2, op=3, aliases=(alias(2),))
        key = ExternalIdentifierKey('uniprot', 'P00533')
        registry.write(assign if assignment_first else withdraw)
        assert len(registry.lookup_external(key)) == (2 if assignment_first else 0)
        registry.write(withdraw if assignment_first else assign)
        assert [x.entity.sf_id for x in registry.lookup_external(key)] == [eid(2)]
        assert registry.history(eid())[0] == first
        assert len(registry.history(eid())) == 2 and len(registry.history(eid(2))) == 1

    def test_old_revision_and_current_descriptor(self, registry):
        registry.write(command(aliases=(alias(),)))
        record = registry.write(command(op=2, expected=1, entity=descriptor(revision=2, lifecycle='RETRACTED'), aliases=(alias(), alias(revision=2))))
        matches = registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))
        assert all(x.entity == record.entity and x.registry_revision == 2 for x in matches)
        assert {x.alias.entity.entity_revision for x in matches} == {1, 2}

    @pytest.mark.parametrize('lifecycle', ['ACTIVE', 'DEPRECATED', 'RETRACTED'])
    @pytest.mark.parametrize('offset', [-1000, 1000])
    def test_no_clock_or_lifecycle_filter(self, registry, lifecycle, offset):
        declaration = alias(lifecycle=lifecycle, valid_from=NOW+timedelta(days=offset), valid_until=NOW+timedelta(days=offset+1))
        registry.write(command(aliases=(declaration,)))
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))[0].alias == declaration

    def test_source_semantics_and_evidence_preserved(self, registry):
        first = registry.write(command(aliases=(alias(),)))
        newer = alias(release='R2', mapping_status='CONFIRMED_MAPPING',
            evidence_refs=[dict(entity_id=eid(50), entity_revision=8, entity_type='evidence')],
            accession_version='100', source_record_revision='200')
        second = registry.write(command(op=2, expected=1, aliases=(alias(), newer)))
        assert second.entity.entity_revision == 1 and second.registry_revision == 2
        assert len(registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))) == 2
        assert registry.get(eid(), registry_revision=1) == first
        assert newer in second.aliases


class ConcurrencyBehavior:
    @pytest.mark.parametrize('mode', ['identical', 'creation', 'revision', 'operation', 'cross_entity_operation', 'shared_key'])
    def test_concurrent_writes(self, registry, mode, registry_overlap):
        a = command(aliases=(alias(),))
        b = a
        expected_error = None
        if mode == 'creation':
            b = replace(a, change=change(2))
            expected_error = RevisionConflictError
        elif mode == 'revision':
            registry.write(command(op=9))
            a = replace(a, expected_registry_revision=1)
            b = replace(a, change=change(2))
            expected_error = RevisionConflictError
        elif mode == 'operation':
            b = replace(a, aliases=())
            expected_error = OperationConflictError
        elif mode == 'cross_entity_operation':
            b = command(n=2, aliases=(alias(2),))
            expected_error = OperationConflictError
        elif mode == 'shared_key':
            b = command(n=2, op=2, aliases=(alias(2),))
        results = registry_overlap(lambda: registry.write(a), lambda: registry.write(b))
        errors = [r for r in results if isinstance(r, RegistryError)]
        assert len(errors) == (1 if expected_error else 0)
        if expected_error:
            assert isinstance(errors[0], expected_error)
            if expected_error is RevisionConflictError:
                assert errors[0].entity_id == a.entity.sf_id
                assert errors[0].expected == a.expected_registry_revision
                assert errors[0].actual == a.expected_registry_revision + 1
            else:
                assert errors[0].operation_id == a.change.operation_id
        if mode == 'identical':
            assert results[0] == results[1]
        expected_count = 2 if mode == 'revision' else 1
        if mode != 'cross_entity_operation':
            assert len(registry.history(eid())) == expected_count
        if mode == 'shared_key':
            assert len(registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))) == 2
            assert registry.get(eid(2)) == results[1]
        # A winner must be replayable through the journal without an extra write.
        for submitted, result in zip((a, b), results):
            if not isinstance(result, RegistryError):
                assert registry.write(submitted) == result
        if mode == 'cross_entity_operation':
            winner = next(r for r in results if not isinstance(r, RegistryError))
            loser = eid(2) if winner.entity.sf_id == eid() else eid()
            with pytest.raises(EntityNotFoundError):
                registry.get(loser)

    def test_lookup_write_boundary(self, registry, registry_overlap):
        old = (alias(), alias(release='R2'))
        new = (alias(release='R3'), alias(release='R4'), alias(release='R5'))
        previous = registry.write(command(aliases=old))
        key = ExternalIdentifierKey('uniprot', 'P00533')
        def expected(record):
            return tuple(ExternalIdentifierMatch(record.entity, record.registry_revision, a) for a in record.aliases)
        for i in range(2, 22):
            before = registry.lookup_external(key)
            assert before == expected(previous)
            target = new if i % 2 == 0 else old
            proposed = command(op=i, expected=i-1, aliases=target,
                entity=descriptor(revision=i, lifecycle='DEPRECATED' if i % 2 else 'ACTIVE'))
            result, observed = registry_overlap(lambda: registry.write(proposed),
                                                 lambda: registry.lookup_external(key))
            assert not isinstance(result, RegistryError)
            after = registry.lookup_external(key)
            assert after == expected(result)
            assert observed in (before, after)
            assert registry.get(eid()) == result
            previous = result


def registry_observation(registry):
    """Complete public view for the failure-injection scenario, backend-neutral."""
    def read(call):
        try:
            return call()
        except EntityNotFoundError as error:
            return ('entity_missing', error.entity_id)
        except RegistryRevisionNotFoundError as error:
            return ('registry_revision_missing', error.entity_id, error.registry_revision)
        except EntityRevisionNotFoundError as error:
            return ('entity_revision_missing', error.entity_id, error.entity_revision)
    entities = tuple((read(lambda: registry.get(eid(n))), read(lambda: registry.history(eid(n))),
        tuple(read(lambda: registry.get(eid(n), registry_revision=r)) for r in (1, 2, 3)),
        tuple(read(lambda: registry.get_entity_revision(eid(n), r)) for r in (1, 2, 3))) for n in (1, 2))
    matches = tuple(registry.lookup_external(ExternalIdentifierKey('uniprot', accession))
                    for accession in ('P00533', 'replacement', 'absent'))
    return entities, matches


def assert_failed_write_unchanged(registry, before, originals):
    assert registry_observation(registry) == before
    for original, result in originals:
        assert registry.write(original) == result
    assert registry_observation(registry) == before
