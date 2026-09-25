"""Static boundaries, exact API, isolated execution, cross-process determinism."""
import ast
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest
from simforge.knowledge import domain, identity, storage
from simforge.knowledge.storage.memory_registry import InMemoryEntityRegistry

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'simforge/knowledge'
MODULES = {
    'identity.records': {'RegistryChange', 'RegistryWrite', 'RegistryRecord', 'ExternalIdentifierKey', 'ExternalIdentifierMatch'},
    'identity.registry': {'EntityRegistry'},
    'identity.errors': {'RegistryError', 'EntityNotFoundError', 'RegistryRevisionNotFoundError',
        'EntityRevisionNotFoundError', 'RevisionConflictError', 'OperationConflictError', 'InvalidRegistryOperationError'},
    'storage.memory_registry': {'InMemoryEntityRegistry'},
}
METHODS = {'write', 'get', 'get_entity_revision', 'history', 'lookup_external'}
FORBIDDEN = {'core', 'runtime', 'builders', 'executors', 'pipelines', 'workflows', 'ligand',
    'analysis', 'validators', 'adapters', 'descriptors', 'utils', 'cli', 'requests', 'httpx',
    'sqlalchemy', 'sqlite3', 'duckdb', 'psycopg', 'psycopg2', 'rdkit', 'numpy', 'scipy',
    'yaml', 'typer', 'rich', 'openmm', 'gromacs'}


def test_exact_api_and_no_convenience_interfaces():
    exports = set()
    for suffix, expected in MODULES.items():
        module = importlib.import_module('simforge.knowledge.' + suffix)
        assert set(module.__all__) == expected
        assert {name for name in vars(module) if not name.startswith('_')} == expected
        exports.update(expected)
    assert set(identity.__all__) == exports - {'InMemoryEntityRegistry'}
    assert storage.__all__ == []
    for cls in (identity.EntityRegistry, InMemoryEntityRegistry):
        assert {name for name, value in vars(cls).items() if not name.startswith('_') and callable(value)} == METHODS
    assert set(domain.__all__) == {'EntityId', 'ExternalIdentifier', 'Entity', 'Protein', 'ProteinSequence',
        'ChemicalIdentity', 'ChemicalSpecies', 'ExperimentalStructure', 'ArtifactRef', 'Claim', 'Evidence',
        'PolicyRef', 'Decision', 'DecisionContext', 'KnowledgeBundleManifest'}


def test_static_dependencies_and_no_global_state():
    allowed = {'typing', 'dataclasses', 'datetime', 'uuid', 'unicodedata', 're'}
    forbidden_calls = {'open', 'exec', 'eval', '__import__', 'import_module', 'write_text', 'write_bytes',
        'read_text', 'read_bytes', 'mkdir', 'system', 'Popen', 'uuid4', 'now', 'today'}
    for package in ('identity', 'storage'):
        for file in (BASE / package).glob('*.py'):
            tree = ast.parse(file.read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    assert all(a.name in allowed | ({'threading'} if package == 'storage' else set()) for a in node.names)
                if isinstance(node, ast.ImportFrom):
                    if node.level:
                        assert package == 'identity' and node.level == 1 and node.module in {'records', 'registry', 'errors'}
                    elif node.module == 'simforge.knowledge.domain':
                        assert all(a.name in domain.__all__ for a in node.names)
                    elif node.module.startswith('simforge.'):
                        assert package == 'storage' and node.module in {
                            'simforge.knowledge.identity.records', 'simforge.knowledge.identity.registry', 'simforge.knowledge.identity.errors'}
                    else:
                        assert node.module in allowed | ({'threading'} if package == 'storage' else set())
                if isinstance(node, ast.Call):
                    name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ''
                    assert name not in forbidden_calls, (file, node.lineno, name)
            for node in tree.body:
                if isinstance(node, (ast.Assign, ast.AnnAssign)):
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    names = {t.id for t in targets if isinstance(t, ast.Name)}
                    # Only exports and private type annotations at module level.
                    assert names <= {'__all__', '_Code', '_Path', '_T'}, (file, names)
    for file in (BASE / 'domain').glob('*.py'):
        tree = ast.parse(file.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert node.level <= 1
                assert not (node.module or '').startswith('simforge.knowledge.identity')
                assert not (node.module or '').startswith('simforge.knowledge.storage')


@pytest.mark.parametrize('entry', ['simforge.knowledge.identity', 'simforge.knowledge.storage.memory_registry'])
def test_guarded_fresh_interpreter(tmp_path, entry):
    code = r'''
import sys, json, importlib
sys.path.insert(0, sys.argv[1])
forbidden = set(json.loads(sys.argv[2]))
class Block:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in forbidden:
            raise AssertionError('forbidden import: ' + fullname)
sys.meta_path.insert(0, Block())
def audit(event, args):
    if event in {'socket.connect', 'socket.bind', 'socket.getaddrinfo', 'subprocess.Popen',
                 'os.system', 'sqlite3.connect', 'os.mkdir', 'os.remove', 'os.rename'}:
        raise AssertionError('forbidden operation: ' + event)
    if event == 'open':
        # Python source/dynamic-module loading is permitted, persistence is not.
        mode = args[1]
        flags = args[2]
        if (isinstance(mode, str) and any(c in mode for c in 'wax+')) or (isinstance(flags, int) and flags & 3):
            raise AssertionError('filesystem write')
sys.addaudithook(audit)
importlib.import_module(sys.argv[3])
if sys.argv[3].endswith('identity'):
    assert 'simforge.knowledge.storage.memory_registry' not in sys.modules
from simforge.knowledge.storage.memory_registry import InMemoryEntityRegistry
from simforge.knowledge.identity import RegistryChange, RegistryWrite, ExternalIdentifierKey
from simforge.knowledge.domain import Entity, EntityId
from uuid import UUID
from datetime import datetime, timezone
import builtins, io
def no_file(*args, **kwargs):
    raise AssertionError('unexpected filesystem operation')
builtins.open = no_file
io.open = no_file
registry = InMemoryEntityRegistry()
entity = Entity(sf_id='sf:00000000-0000-4000-8000-000000000001', entity_revision=1, entity_type='protein')
change = RegistryChange(UUID('10000000-0000-4000-8000-000000000001'), 'actor', 'reason', datetime(2026,1,2,tzinfo=timezone.utc))
command = RegistryWrite(entity, (), 0, change)
record = registry.write(command)
assert registry.write(command) == record
assert registry.get(entity.sf_id) == record
assert registry.get(entity.sf_id, registry_revision=1) == record
assert registry.get_entity_revision(entity.sf_id, 1) == entity
assert registry.history(entity.sf_id) == (record,)
assert registry.lookup_external(ExternalIdentifierKey('uniprot', 'x')) == ()
assert not (set(sys.modules) & forbidden)
assert not any(name.startswith(('simforge.knowledge.federation', 'simforge.knowledge.compiler',
    'simforge.knowledge.resolver', 'simforge.knowledge.runtime', 'simforge.knowledge.policies',
    'simforge.knowledge.importers')) for name in sys.modules)
print('isolated')
'''
    result = subprocess.run([sys.executable, '-I', '-B', '-c', code, str(ROOT), json.dumps(sorted(FORBIDDEN)), entry],
                            capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == 'isolated'


def test_parent_packages_do_not_import_backend(tmp_path):
    result = subprocess.run([sys.executable, '-I', '-B', '-c',
        "import sys; sys.path.insert(0, sys.argv[1]); import simforge.knowledge; import simforge.knowledge.storage; "
        "assert 'simforge.knowledge.storage.memory_registry' not in sys.modules; "
        "assert 'simforge.knowledge.identity' not in sys.modules", str(ROOT)], capture_output=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr


def test_hash_seed_determinism(tmp_path):
    code = r'''
import sys, json
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[1] + '/tests/knowledge/identity')
from registry_contract import command, alias, eid, change
from simforge.knowledge.storage.memory_registry import InMemoryEntityRegistry
from simforge.knowledge.identity import ExternalIdentifierKey, RevisionConflictError
r = InMemoryEntityRegistry()
for number in {'1', '2', '3'}:
    n = int(number)
    r.write(command(n=n, op=n, aliases=tuple(alias(n, release=release) for release in {'R1', 'R2'})))
r.write(command(op=4, expected=1))
records = [(x.entity.sf_id.value, x.registry_revision, x.alias.to_canonical_bytes().decode())
           for x in r.lookup_external(ExternalIdentifierKey('uniprot', 'P00533'))]
history = [(x.registry_revision, x.entity.to_canonical_bytes().decode(), [a.to_canonical_bytes().decode() for a in x.aliases]) for x in r.history(eid())]
try:
    r.write(command(op=5, expected=1))
except RevisionConflictError as error:
    failure = [error.entity_id.value, error.expected, error.actual]
print(json.dumps([records, history, failure], sort_keys=True))
'''
    outputs = []
    for seed in ('1', '777', '42'):
        result = subprocess.run([sys.executable, '-B', '-c', code, str(ROOT)],
            env=dict(os.environ, PYTHONHASHSEED=seed, PYTHONDONTWRITEBYTECODE='1'), capture_output=True, cwd=tmp_path)
        assert result.returncode == 0, result.stderr
        outputs.append(result.stdout)
    assert len(set(outputs)) == 1
