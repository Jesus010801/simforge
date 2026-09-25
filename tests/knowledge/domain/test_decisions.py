import itertools
import pytest
from simforge.knowledge.domain import Decision, DecisionContext
from test_identity import decision, context, operation, ref, uid


@pytest.mark.parametrize('status', ['RESOLVED', 'RESOLVED_WITH_WARNING', 'DEFERRED', 'ABORT'])
def test_four_states_and_round_trip(status):
    value = Decision(**decision(status))
    assert Decision.from_canonical_bytes(value.to_canonical_bytes()) == value


@pytest.mark.parametrize('status,result,deferred,warnings', list(itertools.product(['RESOLVED', 'RESOLVED_WITH_WARNING', 'DEFERRED', 'ABORT'], [False, True], [False, True], [False, True])))
def test_state_matrix(status, result, deferred, warnings):
    data = decision('ABORT') | dict(status=status, result=dict(selected_entities=[ref()]) if result else None, deferred_operation=operation() if deferred else None, warnings=['warning'] if warnings else [])
    valid = result == (status in ('RESOLVED','RESOLVED_WITH_WARNING')) and deferred == (status == 'DEFERRED') and not (status == 'RESOLVED' and warnings) and not (status == 'RESOLVED_WITH_WARNING' and not warnings)
    if valid:
        assert Decision(**data).status == status
    else:
        with pytest.raises(ValueError):
            Decision(**data)


@pytest.mark.parametrize('field', ['operation_key', 'target_stage', 'action', 'method_id', 'method_version', 'method_sha256', 'required_inputs', 'expected_outputs', 'success_criteria', 'on_failure'])
def test_deferred_requires_complete_contract(field):
    data = decision('DEFERRED')
    del data['deferred_operation'][field]
    with pytest.raises(ValueError):
        Decision(**data)


@pytest.mark.parametrize('symbol', ['/bin/sh', 'sh -c x', 'https://example.org', 'science:unknown', 'science:manual_review', 'science:ask_later', 'gromacs:mdrun', 'network:lookup'])
def test_deferred_does_not_describe_execution_or_unknown(symbol):
    data = decision('DEFERRED')
    data['deferred_operation']['action'] = symbol
    with pytest.raises(ValueError):
        Decision(**data)


def test_pinned_binding_and_unique_names():
    data = decision('DEFERRED')
    data['deferred_operation']['required_inputs'][0]['value'] = ref(99)
    with pytest.raises(ValueError):
        Decision(**data)
    data = decision('DEFERRED')
    data['deferred_operation']['expected_outputs'].append(dict(name='species', scientific_type='science:other'))
    with pytest.raises(ValueError):
        Decision(**data)
    data = decision('DEFERRED')
    data['deferred_operation']['on_failure'] = 'RETRY'
    with pytest.raises(ValueError):
        Decision(**data)


def test_selections_context_closure_and_projections():
    data = context()
    value = DecisionContext(**data)
    assert value.resolved_selections and not value.has_mandatory_abort
    data['entity_refs'] = []
    with pytest.raises(ValueError):
        DecisionContext(**data)
    data = context()
    data['decisions'] = [decision('DEFERRED')]
    assert len(DecisionContext(**data).deferred_operations) == 1
    another = decision('DEFERRED') | dict(decision_id=uid(41))
    data['decisions'].append(another)
    with pytest.raises(ValueError):
        DecisionContext(**data)


def test_duplicate_decision_and_reference_type_conflicts():
    data = context()
    data['decisions'].append(decision('ABORT'))
    with pytest.raises(ValueError):
        DecisionContext(**data)
    data = context()
    data['claim_refs'] = [ref(1, 'claim')]
    with pytest.raises(ValueError):
        DecisionContext(**data)


@pytest.mark.parametrize('status', ['UNKNOWN', 'MANUAL_REVIEW', 'ready', '', 1])
def test_no_extra_decision_states(status):
    with pytest.raises(ValueError):
        Decision(**(decision('ABORT') | dict(status=status)))


@pytest.mark.parametrize('key', ['command', 'shell_command', 'executable_path', 'engine', 'engine_config', 'callback', 'api_url'])
def test_deferred_rejects_nested_execution_transport_parameters(key):
    data = decision('DEFERRED')
    data['deferred_operation']['parameters'] = {'nested': [{key: 'not scientific data'}]}
    with pytest.raises(ValueError):
        Decision(**data)
    data = decision('DEFERRED')
    data['parameters'] = {key: 'not scientific data'}
    with pytest.raises(ValueError):
        Decision(**data)


def test_deferred_literal_binding_rejects_execution_configuration():
    data = decision('DEFERRED')
    data['deferred_operation']['required_inputs'] = [dict(name='bad', kind='literal', value={'command': 'curl remote'})]
    with pytest.raises(ValueError):
        Decision(**data)


@pytest.mark.parametrize('field', ['required_inputs', 'expected_outputs', 'success_criteria'])
def test_deferred_requirements_cannot_be_empty(field):
    data = decision('DEFERRED')
    data['deferred_operation'][field] = []
    with pytest.raises(ValueError):
        Decision(**data)


def test_deferred_cannot_have_unbound_or_callable_inputs():
    data = decision('DEFERRED')
    data['deferred_operation']['required_inputs'] = [dict(name='missing', kind='lazy', value='lookup later')]
    with pytest.raises(ValueError):
        Decision(**data)
    data['deferred_operation']['required_inputs'] = [dict(name='function', kind='literal', value=lambda: 1)]
    with pytest.raises(ValueError):
        Decision(**data)


def test_resolved_literal_null_is_explicit_and_empty_result_invalid():
    data = decision()
    data['result'] = dict(value=dict(kind='literal', value=None))
    obj = Decision(**data)
    assert Decision.from_canonical_bytes(obj.to_canonical_bytes()) == obj
    data['result'] = {}
    with pytest.raises(ValueError):
        Decision(**data)


# Deliberately independent payload shapes, not a copy of production field names.
_EXECUTION_PAYLOADS = [
    {'transport': {'url': 'https://example.invalid/lookup', 'timeout': 10}},
    {'argv': ['/bin/sh', '-c', 'echo forbidden']},
    {'path': '/usr/bin/python'},
    {'environment': {'PATH': '/bin', 'OMP_NUM_THREADS': '8'}},
    {'gromacs': {'integrator': 'md', 'nsteps': 1000}},
    {'platform': 'CUDA', 'integrator': {'type': 'LangevinIntegrator', 'steps': 100}},
    {'resources': {'partition': 'gpu', 'nodes': 2, 'walltime': '01:00:00'}},
    {'conditions': [{'measurement': {'argv': ['sh', '-c', 'x']}}]},
    {'temperature_kelvin': {'value': 'https://example.invalid'}},
    {'ph': '/bin/sh'},
    {'tolerance': ['curl', 'remote']},
    {'preserve_stereochemistry': {'nested': {'cwd': '/tmp', 'stdout': 'out'}}},
    'https://example.invalid',
    ['/bin/sh', '-c', 'echo forbidden'],
]


@pytest.mark.parametrize('payload', _EXECUTION_PAYLOADS)
@pytest.mark.parametrize('location', ['decision', 'operation', 'binding'])
def test_positive_scientific_schema_rejects_execution_shapes(payload, location):
    data = decision('DEFERRED')
    if location == 'decision':
        data['parameters'] = payload
    elif location == 'operation':
        data['deferred_operation']['parameters'] = payload
    else:
        data['deferred_operation']['required_inputs'] = [dict(name='conditions', kind='literal', value=payload)]
    with pytest.raises(ValueError):
        Decision(**data)


def test_scientific_conditions_are_exact_pinned_and_roundtrip():
    from decimal import Decimal
    conditions = dict(ph=Decimal('7.40'), temperature_kelvin=Decimal('310.15'),
                      ionic_strength_molar=Decimal('0.15'), tolerance=Decimal('0.0001'),
                      preserve_heavy_atom_connectivity=True, preserve_stereochemistry=False)
    data = decision('DEFERRED')
    data['parameters'] = conditions
    data['deferred_operation']['parameters'] = conditions
    data['deferred_operation']['required_inputs'].append(dict(name='conditions', kind='literal', value=conditions))
    obj = Decision(**data)
    before = obj.to_canonical_bytes()
    conditions['ph'] = Decimal('12')
    assert obj.to_canonical_bytes() == before
    assert obj.deferred_operation.parameters.ph == Decimal('7.4')
    assert Decision.from_canonical_bytes(before) == obj
    assert b'"temperature_kelvin":"310.15"' in before


@pytest.mark.parametrize('conditions', [dict(ph=7.4), dict(ph='7.4'), dict(temperature_kelvin=0),
    dict(ionic_strength_molar=-1), dict(tolerance=0), dict(preserve_stereochemistry=1),
    dict(preserve_heavy_atom_connectivity='true')])
def test_scientific_condition_types_are_strict(conditions):
    data = decision('DEFERRED')
    data['deferred_operation']['parameters'] = conditions
    with pytest.raises(ValueError):
        Decision(**data)


@pytest.mark.parametrize('conditions', [{}, {'ph': None}])
def test_literal_scientific_binding_must_be_concrete(conditions):
    data = decision('DEFERRED')
    data['deferred_operation']['required_inputs'] = [dict(name='conditions', kind='literal', value=conditions)]
    with pytest.raises(ValueError):
        Decision(**data)


def test_scientific_condition_bounds_and_nonfinite_decimals():
    from decimal import Decimal
    for conditions in (dict(temperature_kelvin=Decimal('0')), dict(ionic_strength_molar=Decimal('-1')),
                       dict(tolerance=Decimal('0')), dict(ph=Decimal('NaN')), dict(ph=Decimal('Infinity'))):
        data = decision('DEFERRED')
        data['parameters'] = conditions
        with pytest.raises(ValueError):
            Decision(**data)
    data = decision('DEFERRED')
    data['parameters'] = dict(ph=Decimal('-1'), ionic_strength_molar=Decimal('0'))
    assert Decision(**data).parameters.ph == Decimal('-1')


@pytest.mark.parametrize('action', ['science:protonation', 'science:reconstruction', 'science:tautomer_resolution'])
def test_positive_scientific_action_vocabulary(action):
    data = decision('DEFERRED')
    data['deferred_operation']['action'] = action
    obj = Decision(**data)
    assert Decision.from_canonical_bytes(obj.to_canonical_bytes()) == obj


@pytest.mark.parametrize('method', ['science:unknown', 'science:manual_review', 'science:network',
                                   'science:gromacs', 'gromacs:mdrun', 'https://example.invalid'])
def test_pinned_method_still_rejects_placeholders_and_execution(method):
    data = decision('DEFERRED')
    data['deferred_operation']['method_id'] = method
    with pytest.raises(ValueError):
        Decision(**data)
