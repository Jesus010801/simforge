# core/test_structural_annotation.py
"""
Tests for the Structural Biology Annotation Layer (Phase 1).

Coverage:
    - ResidueRange parsing and validation
    - TransmembraneSegment model
    - MembraneTopologyAnnotation overlap detection
    - OrientationAnnotation geometric constraints
    - StructuralAnnotation completeness helpers
    - migrate_from_legacy_orientation backward compat
    - Full round-trip: new-format YAML → SystemState → structural_annotation
    - Full round-trip: legacy YAML → SystemState → structural_annotation (migration)
"""

from __future__ import annotations

import pytest
from pathlib import Path

from core.structural_annotation import (
    parse_residue_range,
    residues_in_range,
    TransmembraneSegment,
    BiologicalDomain,
    OrientationEvidence,
    MembraneTopologyAnnotation,
    OrientationAnnotation,
    StructuralAnnotation,
    migrate_from_legacy_orientation,
)


# ─────────────────────────────────────────────────────────────────────────────
# parse_residue_range
# ─────────────────────────────────────────────────────────────────────────────

class TestParseResidueRange:
    def test_single_range(self):
        assert parse_residue_range("1-50") == [(1, 50)]

    def test_single_residue(self):
        assert parse_residue_range("42") == [(42, 42)]

    def test_multiple_ranges(self):
        assert parse_residue_range("1-20,45-60") == [(1, 20), (45, 60)]

    def test_individual_residues(self):
        assert parse_residue_range("5,10,15") == [(5, 5), (10, 10), (15, 15)]

    def test_spaces_allowed(self):
        assert parse_residue_range("1-20, 45-60") == [(1, 20), (45, 60)]

    def test_invalid_format(self):
        with pytest.raises(ValueError, match="Formato de rango inválido"):
            parse_residue_range("abc")

    def test_invalid_format_dash_only(self):
        with pytest.raises(ValueError):
            parse_residue_range("-50")

    def test_zero_residue_rejected(self):
        with pytest.raises(ValueError, match="≥ 1"):
            parse_residue_range("0-10")

    def test_inverted_range_rejected(self):
        with pytest.raises(ValueError, match="inicio > fin"):
            parse_residue_range("50-10")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            parse_residue_range("")


# ─────────────────────────────────────────────────────────────────────────────
# residues_in_range
# ─────────────────────────────────────────────────────────────────────────────

class TestResiduesInRange:
    def test_simple_range(self):
        assert residues_in_range("1-5") == {1, 2, 3, 4, 5}

    def test_single_residue(self):
        assert residues_in_range("7") == {7}

    def test_multi_range(self):
        r = residues_in_range("1-3,10-12")
        assert r == {1, 2, 3, 10, 11, 12}


# ─────────────────────────────────────────────────────────────────────────────
# TransmembraneSegment
# ─────────────────────────────────────────────────────────────────────────────

class TestTransmembraneSegment:
    def test_plain_range(self):
        seg = TransmembraneSegment(residues="51-75")
        assert seg.residue_set() == set(range(51, 76))

    def test_with_metadata(self):
        seg = TransmembraneSegment(residues="51-75", label="TM1", helix_type="alpha")
        assert seg.label == "TM1"
        assert seg.helix_type == "alpha"

    def test_invalid_residues_rejected(self):
        with pytest.raises(ValueError):
            TransmembraneSegment(residues="abc")


# ─────────────────────────────────────────────────────────────────────────────
# BiologicalDomain
# ─────────────────────────────────────────────────────────────────────────────

class TestBiologicalDomain:
    def test_valid_domain(self):
        d = BiologicalDomain(residues="1-50", role="extracellular", label="N-loop")
        assert d.residue_set() == set(range(1, 51))

    def test_invalid_role_rejected(self):
        with pytest.raises(ValueError):
            BiologicalDomain(residues="1-50", role="cytoplasmic")


# ─────────────────────────────────────────────────────────────────────────────
# MembraneTopologyAnnotation
# ─────────────────────────────────────────────────────────────────────────────

class TestMembraneTopologyAnnotation:
    def test_no_overlap(self):
        mt = MembraneTopologyAnnotation(
            extracellular_regions=["1-50"],
            intracellular_regions=["200-250"],
            transmembrane_segments=["51-75", "90-112"],
        )
        assert mt.overlap_warnings() == []

    def test_ec_ic_overlap_detected(self):
        mt = MembraneTopologyAnnotation(
            extracellular_regions=["1-100"],
            intracellular_regions=["80-200"],
        )
        warns = mt.overlap_warnings()
        assert any("EC↔IC" in w for w in warns)

    def test_ec_tm_overlap_detected(self):
        mt = MembraneTopologyAnnotation(
            extracellular_regions=["1-60"],
            transmembrane_segments=["50-80"],
        )
        warns = mt.overlap_warnings()
        assert any("EC↔TM" in w for w in warns)

    def test_tm_segment_as_dict(self):
        mt = MembraneTopologyAnnotation(
            transmembrane_segments=[{"residues": "51-75", "label": "TM1"}],
        )
        assert isinstance(mt.transmembrane_segments[0], TransmembraneSegment)
        assert mt.transmembrane_segments[0].label == "TM1"

    def test_tm_segment_as_string(self):
        mt = MembraneTopologyAnnotation(
            transmembrane_segments=["51-75"],
        )
        assert mt.transmembrane_segments[0] == "51-75"

    def test_is_sufficient_for_orient_true(self):
        mt = MembraneTopologyAnnotation(
            extracellular_regions=["1-50"],
            intracellular_regions=["200-250"],
        )
        assert mt.is_sufficient_for_orient()

    def test_is_sufficient_for_orient_false_when_only_ec(self):
        mt = MembraneTopologyAnnotation(extracellular_regions=["1-50"])
        assert not mt.is_sufficient_for_orient()

    def test_tm_segment_count(self):
        mt = MembraneTopologyAnnotation(
            transmembrane_segments=["51-75", "90-112", "130-155"],
        )
        assert mt.tm_segment_count() == 3


# ─────────────────────────────────────────────────────────────────────────────
# OrientationAnnotation
# ─────────────────────────────────────────────────────────────────────────────

class TestOrientationAnnotation:
    def test_default_standard(self):
        o = OrientationAnnotation()
        assert o.is_standard
        assert o.extracellular_side == "+z"
        assert o.intracellular_side == "-z"

    def test_inverted(self):
        o = OrientationAnnotation(extracellular_side="-z", intracellular_side="+z")
        assert not o.is_standard

    def test_same_side_rejected(self):
        with pytest.raises(ValueError, match="opuestos"):
            OrientationAnnotation(extracellular_side="+z", intracellular_side="+z")

    def test_confidence_out_of_range(self):
        with pytest.raises(ValueError):
            OrientationAnnotation(confidence=1.5)

    def test_confidence_zero_allowed(self):
        o = OrientationAnnotation(confidence=0.0)
        assert o.confidence == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# StructuralAnnotation
# ─────────────────────────────────────────────────────────────────────────────

class TestStructuralAnnotation:
    def _make_complete(self) -> StructuralAnnotation:
        return StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
                transmembrane_segments=["31-69"],
            ),
            orientation=OrientationAnnotation(
                extracellular_side="+z",
                source="user_annotation",
                confidence=0.9,
            ),
        )

    def test_is_complete_for_orient(self):
        ann = self._make_complete()
        assert ann.is_complete_for_orient()

    def test_is_partial_without_orientation(self):
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
            ),
        )
        assert ann.is_partial_for_orient()
        assert not ann.is_complete_for_orient()

    def test_incomplete_without_topology(self):
        ann = StructuralAnnotation()
        assert not ann.is_complete_for_orient()
        assert not ann.is_partial_for_orient()

    def test_effective_confidence_complete(self):
        ann = self._make_complete()
        assert ann.effective_confidence() == 0.9

    def test_effective_confidence_partial(self):
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
            ),
        )
        assert 0.0 < ann.effective_confidence() < 1.0

    def test_effective_confidence_empty(self):
        assert StructuralAnnotation().effective_confidence() == 0.0

    def test_validation_warnings_missing_orientation(self):
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
            ),
        )
        warns = ann.validation_warnings()
        assert any("orientation" in w for w in warns)

    def test_no_warnings_when_complete(self):
        ann = self._make_complete()
        assert ann.validation_warnings() == []


# ─────────────────────────────────────────────────────────────────────────────
# migrate_from_legacy_orientation
# ─────────────────────────────────────────────────────────────────────────────

class TestMigrateFromLegacyOrientation:
    def test_full_migration(self):
        ann = migrate_from_legacy_orientation(
            extracellular_residues="1-30",
            intracellular_residues="70-100",
            tm_segments="31-69",
        )
        mt = ann.membrane_topology
        assert mt is not None
        assert mt.extracellular_regions == ["1-30"]
        assert mt.intracellular_regions == ["70-100"]
        assert mt.transmembrane_segments == ["31-69"]

    def test_partial_migration_no_tm(self):
        ann = migrate_from_legacy_orientation(
            extracellular_residues="1-30",
            intracellular_residues="70-100",
            tm_segments=None,
        )
        assert ann.membrane_topology.transmembrane_segments == []

    def test_empty_migration(self):
        ann = migrate_from_legacy_orientation(None, None, None)
        mt = ann.membrane_topology
        assert mt.extracellular_regions == []
        assert mt.intracellular_regions == []

    def test_evidence_records_migration_note(self):
        ann = migrate_from_legacy_orientation("1-30", "70-100", None)
        assert len(ann.evidence) == 1
        assert "legacy" in ann.evidence[0].notes.lower()


# ─────────────────────────────────────────────────────────────────────────────
# Round-trip: new YAML format → SystemState
# ─────────────────────────────────────────────────────────────────────────────

class TestYAMLRoundTrip:
    @pytest.fixture
    def new_format_yaml(self, tmp_path: Path) -> Path:
        content = """
project:
  name: test_membrane

components:
  - id: protein_1
    role: protein
    file: protein.pdb

structural_annotation:
  membrane_topology:
    extracellular_regions:
      - "1-50"
    intracellular_regions:
      - "200-250"
    transmembrane_segments:
      - residues: "51-75"
        label: "TM1"
  orientation:
    extracellular_side: "+z"
    intracellular_side: "-z"
    source: "user_annotation"
    confidence: 0.9

environment:
  membrane:
    enabled: true
    type: DPPC
  solvent:
    water_model: spce
  ions:
    concentration: 0.154
  temperature_K: 310.0
  duration_ns: 50.0

forcefields:
  protein: charmm36

simulation_objectives:
  - membrane_protein_dynamics
"""
        p = tmp_path / "membrane_new.yaml"
        p.write_text(content)
        return p

    @pytest.fixture
    def legacy_format_yaml(self, tmp_path: Path) -> Path:
        content = """
project:
  name: test_membrane_legacy

components:
  - id: protein_1
    role: protein
    file: protein.pdb

environment:
  membrane:
    enabled: true
    type: DPPC
    orientation:
      extracellular_residues: "1-30"
      intracellular_residues: "70-100"
      tm_segments: "31-69"
  solvent:
    water_model: spce
  ions:
    concentration: 0.154
  temperature_K: 300.0
  duration_ns: 10.0

forcefields:
  protein: opls-aa

simulation_objectives:
  - membrane_protein_dynamics
"""
        p = tmp_path / "membrane_legacy.yaml"
        p.write_text(content)
        return p

    def test_new_format_parses(self, new_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(new_format_yaml)
        assert state.structural_annotation is not None

    def test_new_format_topology(self, new_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(new_format_yaml)
        mt = state.structural_annotation.membrane_topology
        assert mt.extracellular_regions == ["1-50"]
        assert mt.intracellular_regions == ["200-250"]

    def test_new_format_tm_segment_as_object(self, new_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(new_format_yaml)
        segs = state.structural_annotation.membrane_topology.transmembrane_segments
        assert isinstance(segs[0], TransmembraneSegment)
        assert segs[0].label == "TM1"

    def test_new_format_orientation(self, new_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(new_format_yaml)
        o = state.structural_annotation.orientation
        assert o.extracellular_side == "+z"
        assert o.confidence == 0.9

    def test_new_format_is_complete(self, new_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(new_format_yaml)
        assert state.structural_annotation.is_complete_for_orient()

    def test_legacy_format_migrates(self, legacy_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(legacy_format_yaml)
        assert state.structural_annotation is not None

    def test_legacy_format_topology_preserved(self, legacy_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(legacy_format_yaml)
        mt = state.structural_annotation.membrane_topology
        assert mt.extracellular_regions == ["1-30"]
        assert mt.intracellular_regions == ["70-100"]
        assert mt.transmembrane_segments == ["31-69"]

    def test_legacy_format_no_explicit_orientation(self, legacy_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(legacy_format_yaml)
        # Legacy migration does not inject an OrientationAnnotation —
        # the user didn't specify axes, so we don't assume them.
        assert state.structural_annotation.orientation is None

    def test_legacy_is_partial_not_complete(self, legacy_format_yaml: Path):
        from core.parser import parse_yaml
        state = parse_yaml(legacy_format_yaml)
        ann = state.structural_annotation
        assert ann.is_partial_for_orient()
        assert not ann.is_complete_for_orient()
