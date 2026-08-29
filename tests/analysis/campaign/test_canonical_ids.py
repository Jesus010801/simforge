"""Canonical scientific identity: determinism, filesystem-safety, order-independence."""
from __future__ import annotations

import random
import re

import pytest

from analysis.campaign.canonical_ids import (
    canonical_system_id, deduplicate_ids, format_replicate, sanitize_token,
    stable_disambiguator,
)

_SAFE_RE = re.compile(r"^[a-z0-9_]+$")


def test_canonical_id_is_deterministic():
    kwargs = dict(receptor="GLP-1R", partner="Semaglutide", water_model="TIP3P",
                  replicate_id="rep2")
    first = canonical_system_id(**kwargs)
    for _ in range(50):
        assert canonical_system_id(**kwargs) == first
    assert first == "glp1r__semaglutide__tip3p__rep02"


@pytest.mark.parametrize("value", [
    "GLP-1R", "SPC/E", "Exendin-4 (mut)", "  Trailing  ", "μ-opioid", "",
    "UPPER_CASE", "weird!!!chars###", "a.b.c",
])
def test_sanitize_token_is_filesystem_safe(value):
    tok = sanitize_token(value)
    assert _SAFE_RE.match(tok), f"{value!r} -> {tok!r} is not [a-z0-9_]"


def test_sanitize_token_examples():
    assert sanitize_token("GLP-1R") == "glp1r"
    assert sanitize_token("Exendin-4") == "exendin4"
    assert sanitize_token("TIP3P") == "tip3p"
    assert sanitize_token("water model") == "water_model"
    assert sanitize_token("SPC/E") == "spc_e"


def test_unknown_semantics_never_invented():
    cid = canonical_system_id(receptor=None, partner=None, water_model=None,
                              replicate_id=None)
    assert cid == "unknown__unknown__unknown__rep_unknown"
    assert sanitize_token(None) == "unknown"
    assert sanitize_token("") == "unknown"


@pytest.mark.parametrize("raw,expected", [
    ("rep1", "rep01"),
    ("r3", "rep03"),
    ("run10", "rep10"),
    ("3", "rep03"),
    (3, "rep03"),
    ("replica07", "rep07"),
    ("final", "rep_final"),
    (None, "rep_unknown"),
])
def test_format_replicate(raw, expected):
    assert format_replicate(raw) == expected


def test_format_replicate_non_numeric_is_token_not_invented_number():
    tok = format_replicate("final")
    assert not any(ch.isdigit() for ch in tok)


def test_extra_dimensions_sorted_key_order_core_fields_first():
    cid = canonical_system_id(
        receptor="R", partner="P", water_model="W", replicate_id=1,
        extra_dimensions={"temperature": "310K", "ph": "7", "salt": "150mM"},
    )
    assert cid.startswith("r__p__w__")
    assert cid.endswith("__rep01")
    # extra keys appear in sorted order: ph, salt, temperature
    assert cid.index("ph_7") < cid.index("salt_150mm") < cid.index("temperature_310k")


def test_extra_dimensions_unknown_values_dropped():
    cid = canonical_system_id(receptor="R", partner="P", water_model="W",
                              replicate_id=1, extra_dimensions={"x": "", "y": "keep"})
    assert "y_keep" in cid
    assert "x_" not in cid


def _pairs():
    return [
        ("glp1r__sema__tip3p__rep01", ["/a/traj.xtc", "digestA"]),
        ("glp1r__sema__tip3p__rep01", ["/b/traj.xtc", "digestB"]),
        ("glp1r__exendin__tip3p__rep01", ["/c/traj.xtc", "digestC"]),
    ]


def test_deduplicate_ids_distinct_suffixes_for_colliding_ids():
    mapping = deduplicate_ids(_pairs())
    assert mapping[2] == "glp1r__exendin__tip3p__rep01"          # unique -> untouched
    assert mapping[0] != mapping[1]                               # collision -> split
    assert mapping[0].startswith("glp1r__sema__tip3p__rep01__")
    assert mapping[1].startswith("glp1r__sema__tip3p__rep01__")


def test_deduplicate_ids_order_independent():
    pairs = _pairs()
    base = deduplicate_ids(pairs)
    base_id_to_tokens = {base[i]: tuple(pairs[i][1]) for i in range(len(pairs))}

    for seed in range(20):
        shuffled = pairs[:]
        random.Random(seed).shuffle(shuffled)
        m = deduplicate_ids(shuffled)
        got = {m[i]: tuple(shuffled[i][1]) for i in range(len(shuffled))}
        assert got == base_id_to_tokens, f"seed {seed} changed the id<->tokens map"


def test_deduplicate_ids_identical_tokens_still_deterministic():
    pairs = [
        ("x__y__z__rep01", ["same"]),
        ("x__y__z__rep01", ["same"]),
        ("x__y__z__rep01", ["same"]),
    ]
    m1 = deduplicate_ids(pairs)
    m2 = deduplicate_ids(list(reversed(pairs)))
    assert sorted(m1.values()) == sorted(m2.values())
    assert any(v.endswith("__dup01") for v in m1.values())
    assert len(set(m1.values())) == 3


def test_stable_disambiguator_order_free():
    assert stable_disambiguator(["a", "b", "c"]) == stable_disambiguator(["c", "a", "b"])
    assert stable_disambiguator(["a", "b"]) != stable_disambiguator(["a", "c"])
    assert len(stable_disambiguator(["a"], length=6)) == 6
