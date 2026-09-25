"""Immutable values and the versioned, in-memory sski-json-v1 wire contract.

Generic Pydantic JSON output is deliberately not the hashing contract.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
import json
import re
import types
import unicodedata
from typing import Annotated, ClassVar, Literal, TypeVar, Union, get_args, get_origin

from pydantic import BaseModel, BeforeValidator, ConfigDict, Field, PlainValidator, TypeAdapter, model_validator


def _nfc(value: str) -> str:
    if not isinstance(value, str):
        raise ValueError('expected a string')
    # Read the built-in string value, not subclass iteration/formatting behavior.
    value = str.__str__(value)
    if any(0xD800 <= ord(c) <= 0xDFFF for c in value):
        raise ValueError('Unicode surrogates are forbidden')
    return unicodedata.normalize('NFC', value)


def _text(value: str) -> str:
    value = _nfc(value)
    if not value.strip():
        raise ValueError('text must be nonempty')
    return value


def _token(value: str) -> str:
    value = _text(value)
    if value != value.strip() or any(unicodedata.category(c) == 'Cc' for c in value):
        raise ValueError('invalid token whitespace/control character')
    return value


def _symbol(value: str) -> str:
    value = _token(value)
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*:[A-Za-z][A-Za-z0-9_.:-]*', value):
        raise ValueError('expected a symbolic namespaced identifier, not a path or command')
    return value


def _utc(value: datetime) -> datetime:
    if not isinstance(value, datetime):
        raise ValueError('timezone-aware datetime required')
    # Base descriptors bypass subclass properties; retain no subclass methods.
    fields = {name: getattr(datetime, name).__get__(value) for name in
              ('year', 'month', 'day', 'hour', 'minute', 'second', 'microsecond', 'tzinfo', 'fold')}
    owned = datetime(**fields)
    offset = owned.utcoffset()
    if offset is None:
        raise ValueError('timezone-aware datetime required')
    # Snapshot even a mutable caller-owned tzinfo into a fixed offset, then UTC.
    return owned.replace(tzinfo=timezone(offset)).astimezone(timezone.utc)


def _finite(value: Decimal) -> Decimal:
    if not isinstance(value, Decimal):
        raise ValueError('finite Decimal required; binary floats are forbidden')
    # Tuple construction is exact and context-independent; call the base method
    # explicitly so a subclass cannot replace its numeric value or formatting.
    owned = Decimal(Decimal.as_tuple(value))
    if not owned.is_finite():
        raise ValueError('finite Decimal required; binary floats are forbidden')
    return owned


def _decimal_text(value: Decimal) -> str:
    value = _finite(value)
    if value.is_zero():
        return '0'
    result = format(value, 'f')
    if '.' in result:
        result = result.rstrip('0').rstrip('.')
    return result


_Text = Annotated[str, BeforeValidator(_text)]
_Token = Annotated[str, BeforeValidator(_token)]
_Symbol = Annotated[str, BeforeValidator(_symbol)]
_PositiveInt = Annotated[int, Field(gt=0)]
_UTC = Annotated[datetime, BeforeValidator(_utc)]
_Decimal = Annotated[Decimal, BeforeValidator(_finite)]


@dataclass(frozen=True, slots=True)
class _Value:
    """Private tree; all supported ingestion routes rebuild and validate it."""
    kind: str
    data: None | bool | int | Decimal | str | tuple[_Value, ...] | tuple[tuple[str, _Value], ...]


def _freeze(value, active=None) -> _Value:
    if isinstance(value, _Value):
        # Private value instances are still untrusted at a public ingestion boundary.
        active = set() if active is None else active
        if id(value) in active:
            raise ValueError('cyclic structured value')
        active.add(id(value))
        try:
            if value.kind == 'object':
                if type(value.data) is not tuple:
                    raise ValueError('immutable object entries required')
                pairs = {}
                for pair in value.data:
                    if type(pair) is not tuple or len(pair) != 2:
                        raise ValueError('invalid object entry')
                    key = _nfc(pair[0])
                    if key in pairs:
                        raise ValueError('NFC key collision')
                    pairs[key] = _freeze(pair[1], active)
                return _Value('object', tuple(sorted(pairs.items())))
            if value.kind == 'array':
                if type(value.data) is not tuple:
                    raise ValueError('immutable array entries required')
                return _Value('array', tuple(_freeze(v, active) for v in value.data))
            expected = {'null': type(None), 'boolean': bool, 'integer': int, 'decimal': Decimal, 'string': str}
            if value.kind not in expected or type(value.data) is not expected[value.kind]:
                raise ValueError('invalid structured scalar')
            return _freeze(value.data, active)
        finally:
            active.remove(id(value))
    if value is None:
        return _Value('null', None)
    if type(value) is bool:
        return _Value('boolean', value)
    if type(value) is int:
        return _Value('integer', value)
    if isinstance(value, Decimal):
        return _Value('decimal', _finite(value))
    if isinstance(value, str):
        return _Value('string', _nfc(value))
    if type(value) not in (dict, list, tuple):
        raise ValueError('unsupported structured value')
    active = set() if active is None else active
    if id(value) in active:
        raise ValueError('cyclic structured value')
    active.add(id(value))
    try:
        if isinstance(value, dict):
            pairs = {}
            for key, item in value.items():
                key = _nfc(key)
                if key in pairs:
                    raise ValueError('NFC key collision')
                pairs[key] = _freeze(item, active)
            return _Value('object', tuple(sorted(pairs.items())))
        return _Value('array', tuple(_freeze(item, active) for item in value))
    finally:
        active.remove(id(value))


def _freeze_object(value) -> _Value:
    result = _freeze(value)
    if result.kind != 'object':
        raise ValueError('object value required')
    return result


_Structured = Annotated[_Value, PlainValidator(_freeze)]
_Object = Annotated[_Value, PlainValidator(_freeze_object)]


def _value_wire(value: _Value):
    data = value.data
    if value.kind == 'decimal':
        data = _decimal_text(data)
    elif value.kind == 'array':
        data = [_value_wire(item) for item in data]
    elif value.kind == 'object':
        data = {key: _value_wire(item) for key, item in data}
    return {'kind': value.kind, 'value': data}


def _value_from_wire(value) -> _Value:
    if not isinstance(value, dict) or set(value) != {'kind', 'value'}:
        raise ValueError('malformed structured value tag')
    kind, data = value['kind'], value['value']
    if kind == 'decimal':
        if not isinstance(data, str) or not re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?', data):
            raise ValueError('invalid decimal representation')
        return _freeze(Decimal(data))
    if kind == 'array' and isinstance(data, list):
        return _Value('array', tuple(_value_from_wire(v) for v in data))
    if kind == 'object' and isinstance(data, dict):
        return _Value('object', tuple(sorted((k, _value_from_wire(v)) for k, v in data.items())))
    expected = {'null': type(None), 'boolean': bool, 'integer': int, 'string': str}
    if kind not in expected or type(data) is not expected[kind]:
        raise ValueError('invalid structured scalar')
    return _freeze(data)


def _wire(value):
    if isinstance(value, _Value):
        return _value_wire(value)
    if isinstance(value, _DomainModel):
        if value._scalar:
            return value.value
        return {name: _wire(getattr(value, name)) for name in type(value).model_fields}
    if isinstance(value, Decimal):
        return _decimal_text(value)
    if isinstance(value, datetime):
        return _utc(value).isoformat(timespec='microseconds').replace('+00:00', 'Z')
    if isinstance(value, tuple):
        return [_wire(v) for v in value]
    if value is None or type(value) in (str, bool, int):
        return value
    raise ValueError('unsupported canonical value')


def _bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode('utf-8')


def _sequence(value):
    if type(value) not in (list, tuple):
        raise ValueError('a list or tuple is required')
    return tuple(value)


def _unordered(value):
    encoded = [_bytes(_wire(v)) for v in value]
    if len(set(encoded)) != len(encoded):
        raise ValueError('duplicate unordered member')
    return tuple(v for _, v in sorted(zip(encoded, value), key=lambda pair: pair[0]))


_T = TypeVar('_T')
_Items = Annotated[tuple[_T, ...], BeforeValidator(_sequence)]


def _pairs(items):
    result = {}
    for key, value in items:
        key = _nfc(key)
        if key in result:
            raise ValueError('duplicate or NFC-colliding JSON key')
        result[key] = value
    return result


def _bad_constant(value):
    raise ValueError('nonfinite JSON number')


def _decode(annotation, value):
    origin = get_origin(annotation)
    if origin is Annotated:
        return _decode(get_args(annotation)[0], value)
    if annotation is _Value:
        return _value_from_wire(value)
    if origin in (Union, types.UnionType):
        for member in get_args(annotation):
            try:
                decoded = _decode(member, value)
                return TypeAdapter(member).validate_python(decoded, strict=True)
            except (ValueError, TypeError):
                continue
        raise ValueError('invalid union value')
    if origin is tuple:
        if not isinstance(value, list):
            raise ValueError('array required')
        return tuple(_decode(get_args(annotation)[0], item) for item in value)
    if isinstance(annotation, type) and issubclass(annotation, _DomainModel):
        if annotation._scalar:
            return annotation.model_validate(value)
        if not isinstance(value, dict):
            raise ValueError('model object required')
        fields = annotation.model_fields
        if set(value) - fields.keys():
            raise ValueError('unknown model fields')
        return annotation.model_validate({k: _decode(fields[k].annotation, v) for k, v in value.items()})
    if annotation is datetime:
        if not isinstance(value, str):
            raise ValueError('datetime string required')
        return _utc(datetime.fromisoformat(value.replace('Z', '+00:00')))
    if annotation is Decimal:
        if not isinstance(value, str) or not re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?', value):
            raise ValueError('decimal string required')
        return _finite(Decimal(value))
    if isinstance(value, str):
        return _nfc(value)
    return value


def _check_identity_roles(root):
    """Check explicitly declared identity roles, never infer identity from strings."""
    roles = {}
    def visit(value):
        if isinstance(value, _DomainModel):
            for key, role in value._identity_roles():
                if key in roles and roles[key] != role:
                    raise ValueError('same internal ID has conflicting entity types or contract roles')
                roles[key] = role
            for name in type(value).model_fields:
                visit(getattr(value, name))
        elif isinstance(value, tuple):
            for item in value:
                visit(item)
    visit(root)


class _DomainModel(BaseModel):
    model_config = ConfigDict(frozen=True, strict=True, extra='forbid', validate_default=True, revalidate_instances='always')
    _scalar: ClassVar[bool] = False
    _unordered_fields: ClassVar[tuple[str, ...]] = ()

    @model_validator(mode='after')
    def _sort_collections(self):
        _check_identity_roles(self)
        for name in self._unordered_fields:
            object.__setattr__(self, name, _unordered(getattr(self, name)))
        return self

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ()

    @staticmethod
    def _validation_options(options):
        # Pydantic's per-call options otherwise override a model's strict config.
        if options.get('strict') is False or options.get('from_attributes') is True:
            raise ValueError('domain validation cannot relax strictness or read arbitrary attributes')
        if options.get('extra') not in (None, 'forbid'):
            raise ValueError('domain validation must forbid extra fields')
        return dict(options, strict=True)

    @classmethod
    def model_validate(cls, obj, **kwargs):
        return super().model_validate(obj, **cls._validation_options(kwargs))

    @classmethod
    def model_validate_json(cls, json_data, **kwargs):
        return super().model_validate_json(json_data, **cls._validation_options(kwargs))

    @classmethod
    def model_validate_strings(cls, obj, **kwargs):
        return super().model_validate_strings(obj, **cls._validation_options(kwargs))

    @classmethod
    def model_construct(cls, *args, **kwargs):
        raise TypeError('unvalidated construction is not supported')

    @classmethod
    def construct(cls, *args, **kwargs):
        raise TypeError('unvalidated construction is not supported')

    @classmethod
    def parse_file(cls, *args, **kwargs):
        raise TypeError('filesystem parsing is outside the domain contract')

    @classmethod
    def parse_raw(cls, *args, **kwargs):
        raise TypeError('use from_canonical_bytes; legacy raw/pickle parsing is unsupported')

    def model_copy(self, *, update=None, deep=False):
        if update is not None:
            raise TypeError('reconstruct through validation instead of updating a copy')
        return type(self).model_validate(self)

    def copy(self, *args, **kwargs):
        raise TypeError('legacy copy is not supported')

    def to_plain(self):
        """Return a detached canonical payload (not a mutable view)."""
        return _wire(self)

    def to_canonical_bytes(self) -> bytes:
        return _bytes({'serialization_version': 'sski-json-v1', 'model': type(self).__name__, 'payload': _wire(self)})

    @classmethod
    def from_canonical_bytes(cls, data: bytes):
        if not isinstance(data, bytes) or data.startswith(b'\xef\xbb\xbf'):
            raise ValueError('UTF-8 bytes without BOM required')
        envelope = json.loads(data.decode('utf-8'), object_pairs_hook=_pairs, parse_constant=_bad_constant)
        if not isinstance(envelope, dict) or set(envelope) != {'serialization_version', 'model', 'payload'}:
            raise ValueError('invalid canonical envelope')
        if envelope['serialization_version'] != 'sski-json-v1' or envelope['model'] != cls.__name__:
            raise ValueError('unknown or mismatched serialization/model version')
        return _decode(cls, envelope['payload'])


class _SHA256(_DomainModel):
    _scalar = True
    value: str

    @model_validator(mode='before')
    @classmethod
    def _validate_digest(cls, value):
        if isinstance(value, str):
            value = {'value': value}
        if isinstance(value, dict):
            digest = value.get('value')
            if not isinstance(digest, str) or not re.fullmatch('[0-9a-fA-F]{64}', digest):
                raise ValueError('SHA256 requires exactly 64 hexadecimal characters')
            return {**value, 'value': digest.lower()}
        return value
