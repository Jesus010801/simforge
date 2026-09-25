"""Only backend provisioning lives here; assertions remain backend-independent."""
import pytest
from simforge.knowledge.storage.memory_registry import InMemoryEntityRegistry


@pytest.fixture
def registry_factory():
    return InMemoryEntityRegistry


@pytest.fixture
def registry_handles(registry_factory):
    registry = registry_factory()
    # A memory instance IS logical state. Persistent fixtures supply two handles
    # and a reopen callable connected to that same state, without changing tests.
    return registry, registry, lambda: registry


@pytest.fixture
def registry(registry_handles):
    return registry_handles[0]


# The approved filenames include test_architecture.py, also used by the frozen
# non-package domain tests. Namespace only this suite's modules so ordinary
# repository-wide pytest collection works without changing root configuration,
# renaming an approved file, or adding an unapproved __init__.py.
class _RegistryTestModule(pytest.Module):
    def _getobj(self):
        import importlib.util
        import sys
        name = '_sski_phase2_tests_' + self.path.stem
        spec = importlib.util.spec_from_file_location(name, self.path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        except BaseException:
            sys.modules.pop(name, None)
            raise
        return module


def pytest_pycollect_makemodule(module_path, parent):
    return _RegistryTestModule.from_parent(parent, path=module_path)


def _memory_overlap(registry, first, second, *, disable_synchronization=False):
    """Backend-specific scheduling; conformance assertions see only two results.

    Hold the first call immediately before publication, then prove the second
    has attempted entry. If it is blocked, release the first; if it entered
    (negative-control backend), let it finish before releasing the first.
    Future backend fixtures can use transaction barriers for the same schedule.
    No sleeps, scheduling probability, or production synchronization hooks.
    """
    import inspect
    import sys
    from threading import Event, Thread, get_ident
    from simforge.knowledge.identity import RegistryError
    source, start = inspect.getsourcelines(InMemoryEntityRegistry.write)
    publication = next(start+i for i, line in enumerate(source) if 'self._state = prepared' in line)
    paused, release, attempted, second_done = (Event() for _ in range(4))
    second_id = [None]
    second_acquired = [False]
    original_lock = registry._lock
    results = [None, None]
    defects = []

    class ObservedLock:
        def __enter__(self):
            acquired = True if disable_synchronization else original_lock.acquire(blocking=False)
            if get_ident() == second_id[0]:
                second_acquired[0] = acquired
                attempted.set()
            if not acquired:
                original_lock.acquire()
            return self
        def __exit__(self, *args):
            if not disable_synchronization:
                original_lock.release()

    def trace(frame, event, arg):
        if event == 'line' and frame.f_code is InMemoryEntityRegistry.write.__code__ and frame.f_lineno == publication:
            paused.set()
            assert release.wait(10), 'coordinator failed to release first publication'
        return trace

    def run(index, call):
        if index == 0:
            sys.settrace(trace)
        else:
            second_id[0] = get_ident()
        try:
            results[index] = call()
        except RegistryError as error:
            results[index] = error
        except BaseException as error:
            defects.append(error)
        finally:
            sys.settrace(None)
            if index == 1:
                second_done.set()

    registry._lock = ObservedLock()
    workers = [Thread(target=run, args=(0, first)), Thread(target=run, args=(1, second))]
    try:
        workers[0].start()
        assert paused.wait(10), 'first write did not reach publication'
        workers[1].start()
        assert attempted.wait(10), 'second operation did not attempt overlapping entry'
        if second_acquired[0]:
            assert second_done.wait(10), 'unlocked second call failed to complete'
        release.set()
    finally:
        release.set()
        for worker in workers:
            if worker.ident is not None:
                worker.join(10)
        registry._lock = original_lock
    assert not any(worker.is_alive() for worker in workers)
    if defects:
        raise defects[0]
    return tuple(results)


@pytest.fixture
def registry_overlap(registry):
    """Schedule first write and a competing call across its commit boundary."""
    return lambda first, second: _memory_overlap(registry, first, second)
