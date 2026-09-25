"""Additional ingress and failure-boundary regressions."""
from dataclasses import replace
import pytest
from simforge.knowledge.domain import Entity, EntityId
from simforge.knowledge.identity import ExternalIdentifierKey, EntityNotFoundError, OperationConflictError
from registry_contract import command, descriptor, alias, eid, change, invalid


@pytest.mark.parametrize('value', [None, {}, 'write', 1])
def test_bad_operation_objects(registry, value):
    invalid(lambda: registry.write(value), 'INVALID_ARGUMENT', ('command',))
    invalid(lambda: registry.lookup_external(value), 'INVALID_ARGUMENT', ('key',))


@pytest.mark.parametrize('value', [True, False, '1', 1.0, -1, None])
def test_write_revalidates_bound_command(registry, value):
    forged = command()
    object.__setattr__(forged, 'expected_registry_revision', value)
    invalid(lambda: registry.write(forged), 'INVALID_ARGUMENT', ('expected_registry_revision',))
    assert registry.write(command()).registry_revision == 1


def test_nested_change_and_key_paths_at_operation_ingress(registry):
    forged = command()
    object.__setattr__(forged.change, 'actor', '')
    invalid(lambda: registry.write(forged), 'INVALID_ARGUMENT', ('change', 'actor'))
    key = ExternalIdentifierKey('uniprot', 'x')
    object.__setattr__(key, 'namespace', 'UPPER')
    invalid(lambda: registry.lookup_external(key), 'INVALID_ARGUMENT', ('key', 'namespace'))


def test_structural_validation_precedes_replay(registry):
    original = command()
    registry.write(original)
    object.__setattr__(original.entity, 'entity_revision', True)
    invalid(lambda: registry.write(original), 'INVALID_ARGUMENT', ('entity', 'entity_revision'))
    assert registry.get(eid()).entity.entity_revision == 1


def test_operation_conflict_precedes_transition_and_alias_checks(registry):
    registry.write(command())
    with pytest.raises(OperationConflictError):
        registry.write(command(entity=descriptor(revision=3), aliases=(alias(2),)))


def test_noop_new_operation_is_new_snapshot(registry):
    first = registry.write(command())
    second = registry.write(command(op=2, expected=1))
    assert first.entity == second.entity
    assert second.registry_revision == 2
    assert registry.get_entity_revision(eid(), 1) == first.entity


def test_returned_and_ingress_objects_do_not_share_storage(registry):
    cmd = command(aliases=(alias(),))
    returned = registry.write(cmd)
    object.__setattr__(cmd.entity, 'lifecycle', 'RETRACTED')
    object.__setattr__(cmd.aliases[0].provenance.source_record, 'source_release', 'corrupted input')
    object.__setattr__(returned.entity, 'lifecycle', 'DEPRECATED')
    object.__setattr__(returned.aliases[0], 'lifecycle', 'RETRACTED')
    current = registry.get(eid())
    assert current.entity.lifecycle == 'ACTIVE'
    assert current.aliases[0].lifecycle == 'ACTIVE'
    assert current.aliases[0].provenance.source_record.source_release == 'R1'
    for record in (registry.get(eid()), registry.history(eid())[0]):
        object.__setattr__(record.entity, 'lifecycle', 'RETRACTED')
    match = registry.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))[0]
    object.__setattr__(match.alias, 'lifecycle', 'RETRACTED')
    revision = registry.get_entity_revision(eid(), 1)
    object.__setattr__(revision, 'lifecycle', 'RETRACTED')
    assert registry.get(eid()).entity.lifecycle == 'ACTIVE'
    assert registry.get_entity_revision(eid(), 1).lifecycle == 'ACTIVE'
    assert registry.write(command(aliases=(alias(),))).aliases[0].lifecycle == 'ACTIVE'


def test_storage_defect_is_not_reclassified(registry, monkeypatch):
    cmd = command()
    def defect(*args, **kwargs):
        raise RuntimeError('implementation defect')
    monkeypatch.setattr(Entity, 'model_validate', defect)
    with pytest.raises(RuntimeError):
        registry.write(cmd)


def test_operation_signature_binding(registry):
    with pytest.raises(TypeError):
        registry.get(eid(), 1)
    with pytest.raises(TypeError):
        registry.write(command(), extra=True)
    with pytest.raises(TypeError):
        registry.history()


def test_invalid_read_revision_retains_entity_id(registry):
    assert invalid(lambda: registry.get(eid(), registry_revision=0), 'INVALID_ARGUMENT').entity_id == eid()
    assert invalid(lambda: registry.get_entity_revision(eid(), 0), 'INVALID_ARGUMENT').entity_id == eid()


def _failure_scenario(creation=False):
    from simforge.knowledge.storage.memory_registry import InMemoryEntityRegistry
    registry = InMemoryEntityRegistry()
    originals = []
    commands = (command(n=2, op=20, aliases=(alias(2),)),)
    if not creation:
        commands += (command(aliases=(alias(),)),)
    for cmd in commands:
        originals.append((cmd, registry.write(cmd)))
    proposed = command(op=2, expected=0 if creation else 1,
        entity=descriptor(revision=1 if creation else 2),
        aliases=(alias(accession='replacement', revision=1 if creation else 2),))
    return registry, proposed, originals


def _trace_prepublication(creation=False, inject=None):
    """Every executed backend line occurrence up to and including publication.

    Tracing includes comprehensions, copy/validation calls, every prepared map,
    journal construction and aggregate construction. Repeated line occurrences
    are distinct boundaries. No production test hook or filesystem edits.
    """
    import inspect
    import sys
    from simforge.knowledge.storage import memory_registry as backend
    registry, proposed, originals = _failure_scenario(creation)
    source, start = inspect.getsourcelines(backend.InMemoryEntityRegistry.write)
    publication = next((start+i for i, line in enumerate(source) if 'self._state = prepared' in line), None)
    # Before remediation, publication was the final separate journal mutation.
    if publication is None:
        publication = next(start+i for i, line in enumerate(source) if 'self._operations[operation] =' in line)
    from registry_contract import registry_observation
    before = registry_observation(registry)
    boundaries = []
    active = True
    def trace(frame, event, arg):
        nonlocal active
        if active and event == 'line' and frame.f_code.co_filename == backend.__file__:
            boundaries.append((frame.f_code.co_name, frame.f_lineno))
            if len(boundaries)-1 == inject:
                raise MemoryError('injected preparation/publication failure')
            if frame.f_code is backend.InMemoryEntityRegistry.write.__code__ and frame.f_lineno == publication:
                active = False
        return trace
    sys.settrace(trace)
    try:
        result = registry.write(proposed)
    except MemoryError:
        result = None
    finally:
        sys.settrace(None)
    return boundaries, registry, proposed, originals, before, result


@pytest.mark.parametrize('creation', [False, True], ids=['revision', 'registration'])
def test_every_preparation_and_publication_failure_is_atomic(creation):
    from registry_contract import assert_failed_write_unchanged
    boundaries, _, _, _, _, _ = _trace_prepublication(creation)
    assert len(boundaries) > 20
    for index, boundary in enumerate(boundaries):
        reached, registry, proposed, originals, before, result = _trace_prepublication(creation, index)
        assert reached[-1] == boundary and result is None
        assert_failed_write_unchanged(registry, before, originals)
        committed = registry.write(proposed)
        expected = 1 if creation else 2
        assert committed.registry_revision == expected
        assert registry.write(proposed) == committed
        assert [r.registry_revision for r in registry.history(eid())] == list(range(1, expected+1))
        assert registry.get(eid()) == committed
        assert registry.get(eid(), registry_revision=expected) == committed
        assert registry.get_entity_revision(eid(), proposed.entity.entity_revision) == proposed.entity
        assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'replacement'))[0].registry_revision == expected
        print(f'closed failure boundary {creation=} {index=} {boundary=}')
