from datetime import timedelta, timezone
from decimal import Decimal, localcontext
import json
import pytest
from simforge.knowledge.domain import (
    EntityId, ExternalIdentifier, Entity, Protein, ProteinSequence, ChemicalIdentity,
    ChemicalSpecies, ExperimentalStructure, ArtifactRef, Claim, Evidence,
    PolicyRef, Decision, DecisionContext, KnowledgeBundleManifest,
)
from test_identity import uid, base, external, provenance, artifact, policy, decision, context, manifest, NOW
from test_chemistry import chemical, species
from test_structures import structure
from test_claims_evidence import claim, evidence


def samples():
    return [EntityId(value=uid()), ExternalIdentifier(**external()), Entity(**base(), entity_type='protein'),
            Protein(**base(), provenance=[provenance()]), ProteinSequence(**base(), sequence='a c d', provenance=[provenance()]),
            ChemicalIdentity(**chemical()), ChemicalSpecies(**species()), ExperimentalStructure(**structure()),
            ArtifactRef(**artifact()), Claim(**claim()), Evidence(**evidence()), PolicyRef(**policy()),
            Decision(**decision()), DecisionContext(**context()), KnowledgeBundleManifest(**manifest())]


@pytest.mark.parametrize('index', range(15))
def test_all_public_contracts_roundtrip_and_immutable(index):
    value = samples()[index]
    encoded = value.to_canonical_bytes()
    assert type(value).from_canonical_bytes(encoded).to_canonical_bytes() == encoded
    field = next(iter(type(value).model_fields))
    with pytest.raises(ValueError):
        setattr(value, field, None)
    with pytest.raises(TypeError):
        type(value).model_construct()
    with pytest.raises(TypeError):
        type(value).construct()
    with pytest.raises(TypeError):
        value.model_copy(update={})
    with pytest.raises(TypeError):
        value.copy()
    assert value.model_copy() == value
    plain = value.to_plain()
    if isinstance(plain, dict):
        plain.clear()
    assert value.to_canonical_bytes() == encoded


def value_claim(value):
    return Claim(**(claim() | dict(object=dict(kind='literal', value=value))))


def test_recursive_defensive_copy_and_nested_freeze():
    source = {'measurements': [Decimal('2.50'), {'name': 'e\u0301'}]}
    value = value_claim(source)
    before = value.to_canonical_bytes()
    source['measurements'][1]['name'] = 'changed'
    source['measurements'].append(4)
    assert value.to_canonical_bytes() == before
    with pytest.raises((AttributeError, TypeError)):
        value.object.value.data[0][1].data += ()
    with pytest.raises(ValueError):
        value.assertion_provenance[0].source_record.record_id = 'changed'


@pytest.mark.parametrize('bad', [1.0, float('nan'), float('inf'), Decimal('NaN'), Decimal('Infinity'), Decimal('-Infinity'), {1: 'x'}, {1, 2}, iter([1]), object(), lambda: None])
def test_unsupported_structured_values(bad):
    with pytest.raises(ValueError):
        value_claim(bad)


def test_cycle_rejected():
    cyclic = []
    cyclic.append(cyclic)
    with pytest.raises(ValueError):
        value_claim(cyclic)


def test_nfc_equivalence_and_key_collision():
    assert value_claim({'e\u0301': 'e\u0301'}).to_canonical_bytes() == value_claim({'é': 'é'}).to_canonical_bytes()
    with pytest.raises(ValueError):
        value_claim({'e\u0301': 1, 'é': 2})


@pytest.mark.parametrize('value', ['\ud800', {'\udfff': 1}, ['\ud800']])
def test_surrogates_rejected(value):
    with pytest.raises(ValueError):
        value_claim(value)


@pytest.mark.parametrize('source,expected', [('2.50', '2.5'), ('2.5', '2.5'), ('-0.00', '0'), ('0E+100', '0'), ('1E+3', '1000'), ('0.0000100', '0.00001'), ('123456789.012345678900', '123456789.0123456789')])
def test_decimal_golden_vectors_across_contexts(source, expected):
    outputs = []
    for precision in (2, 6, 28, 60):
        with localcontext() as ctx:
            ctx.prec = precision
            value = value_claim(Decimal(source))
            wire = json.loads(value.to_canonical_bytes())
            assert wire['payload']['object']['value'] == {'kind': 'decimal', 'value': expected}
            outputs.append(value.to_canonical_bytes())
    assert len(set(outputs)) == 1


def test_scalar_types_remain_distinct_and_equivalent_decimals_match():
    assert len({value_claim(v).to_canonical_bytes() for v in (1, Decimal('1'), True, '1')}) == 4
    assert value_claim(Decimal('2.50')).to_canonical_bytes() == value_claim(Decimal('2.5')).to_canonical_bytes()


def test_canonical_root_golden_and_unknown_envelopes():
    value = EntityId(value=uid())
    expected = b'{"model":"EntityId","payload":"sf:00000000-0000-4000-8000-000000000001","serialization_version":"sski-json-v1"}'
    assert value.to_canonical_bytes() == expected
    for bad in (expected.replace(b'sski-json-v1', b'sski-json-v2'), expected.replace(b'EntityId', b'Unknown'), b'\xef\xbb\xbf'+expected,
                b'{"model":"EntityId","model":"EntityId","payload":null,"serialization_version":"sski-json-v1"}'):
        with pytest.raises(ValueError):
            EntityId.from_canonical_bytes(bad)


def test_duplicate_normalized_json_keys_and_nonfinite_json():
    data = value_claim({'x': 1}).to_canonical_bytes()
    for bad in (data.replace(b'"x":', '"é":{"kind":"integer","value":2},"e\\u0301":'.encode()),
                data.replace(b'"value":1', b'"value":NaN')):
        with pytest.raises(ValueError):
            Claim.from_canonical_bytes(bad)


def test_ordering_dates_and_defaults():
    assert value_claim({'b': 2, 'a': 1}).to_canonical_bytes() == value_claim({'a': 1, 'b': 2}).to_canonical_bytes()
    assert value_claim([1,2]).to_canonical_bytes() != value_claim([2,1]).to_canonical_bytes()
    a = provenance()
    b = provenance() | dict(imported_at=NOW.astimezone(timezone(timedelta(hours=-6))))
    x = Protein(**base(), provenance=[a])
    y = Protein(**base(), provenance=[b])
    assert x.to_canonical_bytes() == y.to_canonical_bytes()
    assert b'2026-01-02T03:04:05.000000Z' in x.to_canonical_bytes()
    assert b'"label":null' in x.to_canonical_bytes()
    second = provenance() | dict(importer_version='2')
    assert Protein(**base(), provenance=[a,second]).to_canonical_bytes() == Protein(**base(), provenance=[second,a]).to_canonical_bytes()


def test_unknown_fields_and_malformed_tags_on_decode():
    data = json.loads(value_claim(1).to_canonical_bytes())
    data['payload']['object']['value']['kind'] = 'decimal'
    with pytest.raises(ValueError):
        Claim.from_canonical_bytes(json.dumps(data).encode())
    data = json.loads(value_claim(1).to_canonical_bytes())
    data['payload']['unknown'] = True
    with pytest.raises(ValueError):
        Claim.from_canonical_bytes(json.dumps(data).encode())


@pytest.mark.parametrize('options', [dict(strict=False), dict(from_attributes=True), dict(extra='allow'), dict(extra='ignore')])
def test_public_validation_cannot_relax_contract(options):
    data = base() | dict(entity_type='protein', entity_revision=True)
    with pytest.raises(ValueError):
        Entity.model_validate(data, **options)
    with pytest.raises(ValueError):
        Entity.model_validate_json(json.dumps(data), **options)
    with pytest.raises(ValueError):
        Entity.model_validate_strings(data, **options)


def test_all_nested_digest_fields_reject_missing_and_invalid():
    def digest_paths(value, path=()):
        if isinstance(value, dict):
            for key, item in value.items():
                if key == 'sha256' or key.endswith('_sha256'):
                    yield path + (key,)
                else:
                    yield from digest_paths(item, path + (key,))
        elif isinstance(value, list):
            for i, item in enumerate(value):
                yield from digest_paths(item, path + (i,))
    objects = samples() + [Decision(**decision('DEFERRED'))]
    for obj in objects:
        original = json.loads(obj.to_canonical_bytes())
        for path in digest_paths(original):
            for remove in (False, True):
                wire = json.loads(obj.to_canonical_bytes())
                target = wire
                for part in path[:-1]:
                    target = target[part]
                if remove:
                    del target[path[-1]]
                else:
                    target[path[-1]] = ''
                with pytest.raises(ValueError):
                    type(obj).from_canonical_bytes(json.dumps(wire).encode())


@pytest.mark.parametrize('index', range(15))
def test_legacy_io_and_pickle_entrypoints_are_disabled(index, monkeypatch):
    import builtins
    def forbidden(*args, **kwargs):
        raise AssertionError('attempted filesystem access')
    cls = type(samples()[index])
    monkeypatch.setattr(builtins, 'open', forbidden)
    with pytest.raises(TypeError):
        cls.parse_file('/must/not/be/accessed')
    with pytest.raises(TypeError):
        cls.parse_raw(b'not a pickle', allow_pickle=True)


def test_nested_digest_value_is_frozen():
    value = ArtifactRef(**artifact())
    with pytest.raises(ValueError):
        value.sha256.value = 'b'*64


@pytest.mark.parametrize('kind,data', [('object', (('\ud800', 1),)), ('object', (('é', 1), ('e\u0301', 2))), ('integer', True), ('array', [1, 2])])
def test_private_value_instances_are_revalidated(kind, data):
    from simforge.knowledge.domain.serialization import _Value
    with pytest.raises(ValueError):
        value_claim(_Value(kind, data))


def test_private_value_revalidation_normalizes_keys():
    from simforge.knowledge.domain.serialization import _Value
    assert value_claim(_Value('object', (('e\u0301', _Value('integer', 1)),))).to_canonical_bytes() == value_claim({'é': 1}).to_canonical_bytes()


class _CallerDecimal(Decimal):
    def __format__(self, spec):
        return self.rendering

    def as_tuple(self):
        raise AssertionError('subclass numeric conversion must not run')

    def is_finite(self):
        raise AssertionError('subclass finiteness must not run')


@pytest.mark.parametrize('precision', [2, 6, 28, 60])
def test_decimal_subclass_is_owned_exactly_without_context_rounding(precision):
    with localcontext() as ctx:
        ctx.prec = precision
        source = _CallerDecimal('123456789.012345678900')
        source.rendering = '7.5'
        obj = value_claim(source)
        expected = value_claim(Decimal('123456789.012345678900')).to_canonical_bytes()
        assert type(obj.object.value.data) is Decimal
        assert obj.to_canonical_bytes() == expected
        source.rendering = '9.5'
        assert obj.to_canonical_bytes() == expected
        assert Claim.from_canonical_bytes(expected) == obj


def test_datetime_subclass_cannot_change_owned_canonical_instant():
    from datetime import datetime
    class CallerDate(datetime):
        def isoformat(self, *args, **kwargs):
            return self.rendering

        def astimezone(self, *args, **kwargs):
            raise AssertionError('subclass timezone conversion must not run')

        @property
        def year(self):
            return 1999

    zone = timezone(timedelta(hours=-6))
    source = CallerDate(2026, 1, 2, 3, 4, 5, 123456, tzinfo=zone)
    source.rendering = 'not the stored instant'
    p = provenance() | dict(imported_at=source)
    obj = Protein(**base(), provenance=[p])
    expected = Protein(**base(), provenance=[provenance() | dict(
        imported_at=datetime(2026, 1, 2, 9, 4, 5, 123456, tzinfo=timezone.utc))])
    assert type(obj.provenance[0].imported_at) is datetime
    assert obj.to_canonical_bytes() == expected.to_canonical_bytes()
    source.rendering = 'changed again'
    assert obj.to_canonical_bytes() == expected.to_canonical_bytes()


def test_mutable_timezone_is_snapshotted_and_fold_preserved():
    from datetime import datetime, tzinfo
    class CallerZone(tzinfo):
        hours = 2
        def utcoffset(self, dt):
            return timedelta(hours=self.hours - dt.fold)
        def dst(self, dt):
            return timedelta(0)
    zone = CallerZone()
    source = datetime(2026, 1, 2, 3, tzinfo=zone, fold=1)
    obj = Protein(**base(), provenance=[provenance() | dict(imported_at=source)])
    before = obj.to_canonical_bytes()
    assert b'2026-01-02T02:00:00.000000Z' in before
    zone.hours = 9
    assert obj.to_canonical_bytes() == before


def test_string_subclass_behavior_is_not_retained():
    class CallerString(str):
        def __str__(self):
            return self.rendering
        def __iter__(self):
            raise AssertionError('subclass iteration must not run')
    source = CallerString('source value')
    source.rendering = 'different value'
    obj = value_claim(source)
    assert type(obj.object.value.data) is str
    expected = value_claim('source value').to_canonical_bytes()
    source.rendering = 'changed'
    assert obj.to_canonical_bytes() == expected
