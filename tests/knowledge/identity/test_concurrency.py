"""Repeated synchronized contention exercises the shared contract assertions."""
import pytest
from registry_contract import ConcurrencyBehavior


@pytest.mark.parametrize('iteration', range(10))
def test_repeated_operation_collision(registry, iteration, registry_overlap):
    ConcurrencyBehavior().test_concurrent_writes(registry, 'cross_entity_operation', registry_overlap)


@pytest.mark.parametrize('mode', ['creation', 'revision', 'operation', 'cross_entity_operation', 'shared_key'])
def test_controlled_overlap_rejects_disabled_synchronization(registry, mode):
    from conftest import _memory_overlap
    # Invoke the unchanged semantic assertions against a disposable backend with
    # its synchronization disabled. Each relevant lost-update defect must fail.
    overlap = lambda first, second: _memory_overlap(registry, first, second, disable_synchronization=True)
    with pytest.raises(AssertionError):
        ConcurrencyBehavior().test_concurrent_writes(registry, mode, overlap)


def test_boundary_assertion_rejects_permanently_stale_lookup(registry, registry_overlap):
    class StaleLookup:
        def __init__(self):
            self.initial = None
        def write(self, command):
            return registry.write(command)
        def get(self, entity_id):
            return registry.get(entity_id)
        def lookup_external(self, key):
            if self.initial is None:
                self.initial = registry.lookup_external(key)
            return self.initial
    with pytest.raises(AssertionError):
        ConcurrencyBehavior().test_lookup_write_boundary(StaleLookup(), registry_overlap)
