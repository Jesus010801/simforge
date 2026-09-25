"""Check cold imports in children; the pytest parent loads legacy conftest."""
import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest
from simforge.knowledge import domain
from test_serialization import samples

ROOT = Path(__file__).resolve().parents[3]
PRODUCTION = ROOT / 'simforge/knowledge/domain'
FORBIDDEN = {'core', 'runtime', 'builders', 'executors', 'pipelines', 'workflows', 'ligand', 'analysis', 'validators', 'adapters', 'descriptors', 'utils', 'cli', 'requests', 'httpx', 'sqlalchemy', 'sqlite3', 'duckdb', 'psycopg', 'psycopg2', 'rdkit', 'numpy', 'scipy', 'yaml', 'typer', 'rich', 'openmm', 'gromacs'}


def test_static_import_and_operation_boundary():
    allowed = {'__future__', 'dataclasses', 'datetime', 'decimal', 'json', 're', 'types', 'unicodedata', 'typing', 'uuid', 'hashlib', 'pydantic'}
    for file in PRODUCTION.glob('*.py'):
        tree = ast.parse(file.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                assert all(a.name.split('.')[0] in allowed for a in node.names), (file, node.lineno)
            elif isinstance(node, ast.ImportFrom) and not node.level:
                assert node.module.split('.')[0] in allowed, (file, node.lineno)
            elif isinstance(node, ast.ImportFrom):
                assert node.level == 1 and (PRODUCTION / (node.module + '.py')).exists()
            if isinstance(node, ast.Call):
                name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ''
                assert name not in {'open', 'exec', 'eval', '__import__', 'import_module', 'write_text', 'write_bytes', 'read_text', 'read_bytes', 'mkdir', 'system', 'Popen', 'uuid4', 'now', 'today'}, (file, name)


def test_exact_public_api():
    assert set(domain.__all__) == {'EntityId', 'ExternalIdentifier', 'Entity', 'Protein', 'ProteinSequence', 'ChemicalIdentity', 'ChemicalSpecies', 'ExperimentalStructure', 'ArtifactRef', 'Claim', 'Evidence', 'PolicyRef', 'Decision', 'DecisionContext', 'KnowledgeBundleManifest'}


def test_fresh_child_import_and_no_io(tmp_path):
    code = r'''
import sys, json, builtins
sys.path.insert(0, sys.argv[1])
forbidden = set(json.loads(sys.argv[2]))
class Block:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in forbidden:
            raise AssertionError('forbidden import: ' + fullname)
sys.meta_path.insert(0, Block())
def audit(event, args):
    if event in {'socket.connect', 'socket.bind', 'subprocess.Popen', 'os.system', 'sqlite3.connect'}:
        raise AssertionError('forbidden operation: ' + event)
sys.addaudithook(audit)
from simforge.knowledge import domain
assert not (set(sys.modules) & forbidden)
records = json.loads(sys.stdin.read())
def no_file(*args, **kwargs):
    raise AssertionError('unexpected file access')
builtins.open = no_file
import io
io.open = no_file
for name, payload in records:
    model = getattr(domain, name)
    value = model.from_canonical_bytes(payload.encode())
    assert value.to_canonical_bytes() == payload.encode()
assert not (set(sys.modules) & forbidden)
print('isolated')
'''
    records = [(type(v).__name__, v.to_canonical_bytes().decode()) for v in samples()]
    result = subprocess.run([sys.executable, '-I', '-B', '-c', code, str(ROOT), json.dumps(sorted(FORBIDDEN))], input=json.dumps(records), text=True, capture_output=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == 'isolated'


def test_cross_process_canonical_bytes(tmp_path):
    # Normal interpreter intentionally honors PYTHONHASHSEED. Explicit source root.
    code = r'''
import sys
sys.path.insert(0, sys.argv[1])
from decimal import Decimal, localcontext
from datetime import datetime, timezone, timedelta
from simforge.knowledge.domain import Claim
with localcontext() as ctx:
    ctx.prec = int(sys.argv[2])
    values = {key: Decimal('2.50') for key in set(['é', 'a', 'z'])}
    if sys.argv[2] == '2': values = dict(reversed(list(values.items())))
    p = dict(source_record=dict(source_id='fixture',record_id='record',version_status='UNVERSIONED'),raw_payload_sha256='a'*64,importer_id='fixture',importer_version='1',imported_at=datetime(2026,1,2,tzinfo=timezone.utc).astimezone(timezone(timedelta(hours=int(sys.argv[3])))),extraction_method='fixture')
    obj = Claim(sf_id='sf:00000000-0000-4000-8000-000000000002',entity_revision=1,subject=dict(entity_id='sf:00000000-0000-4000-8000-000000000001',entity_revision=1,entity_type='protein'),predicate='science:test',object=dict(kind='literal',value=values),assertion_provenance=[p])
    sys.stdout.buffer.write(obj.to_canonical_bytes())
'''
    outputs = []
    for seed, precision, offset in [('1','2','-6'), ('777','28','0'), ('42','60','3')]:
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONDONTWRITEBYTECODE='1')
        result = subprocess.run([sys.executable, '-B', '-c', code, str(ROOT), precision, offset], env=env, capture_output=True, cwd=tmp_path)
        assert result.returncode == 0, result.stderr
        outputs.append(result.stdout)
    assert len(set(outputs)) == 1


def test_legacy_collision_types_remain_separate():
    from runtime.artifacts import ArtifactRef
    from core.md_knowledge.evidence import Evidence
    from ligand.campaign import ChemicalIdentity
    assert domain.ArtifactRef is not ArtifactRef
    assert domain.Evidence is not Evidence
    assert domain.ChemicalIdentity is not ChemicalIdentity
