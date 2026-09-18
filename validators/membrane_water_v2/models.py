"""Versioned typed atlas models. Geometry objects never store water coordinates."""
from dataclasses import dataclass, field
import numpy as np
from .cell import PeriodicCell, Lattice

SCHEMA_VERSION = "membrane-water-atlas/2.0"

@dataclass
class SystemGeometry:
    cell: PeriodicCell
    coordinates: np.ndarray
    radii: np.ndarray
    species: np.ndarray
    atom_names: tuple[str, ...]
    molecule_ids: np.ndarray
    radius_known: np.ndarray
    warnings: list[str] = field(default_factory=list)
    schema_version: str = SCHEMA_VERSION

@dataclass(frozen=True)
class WaterMolecule:
    uid: int
    atom_indices: tuple[int, ...]
    oxygen_index: int
    residue_number: int
    residue_name: str
    valid: bool = True

@dataclass(frozen=True)
class ExclusionModel:
    probe_radius_nm: float = .14
    hard_contact_scale: float = .70
    oxygen_radius_nm: float = .152
    boundary_tolerance_nm: float = .015
    radius_uncertainty_nm: float = .02
    schema_version: str = SCHEMA_VERSION
    def __post_init__(self):
        if self.probe_radius_nm < 0 or not 0 < self.hard_contact_scale <= 1:
            raise ValueError("Invalid exclusion model")

@dataclass
class MembraneAtlas:
    normal_axis: int
    tangent_axes: tuple[int, int]
    center_fraction: float
    lower: np.ndarray
    upper: np.ndarray
    supported: np.ndarray
    normal: np.ndarray
    warnings: list[str] = field(default_factory=list)
    schema_version: str = SCHEMA_VERSION

@dataclass
class FreeSpaceField:
    lattice: Lattice
    clearance: np.ndarray
    nearest_species: np.ndarray
    hard_excluded: np.ndarray
    nominal_free: np.ndarray
    certified_free: np.ndarray
    possible_free: np.ndarray
    uncertain: np.ndarray | None = None
    schema_version: str = SCHEMA_VERSION

@dataclass
class VoidComponent:
    id: int
    voxel_count: int
    region_ids: list[int] = field(default_factory=list)
    upper_access: bool = False
    lower_access: bool = False

@dataclass
class Region:
    id: int
    component_id: int
    voxel_count: int
    classification: str = "ambiguous"
    attributes: dict = field(default_factory=dict)
    evidence: dict = field(default_factory=dict)
    uncertainty: list[str] = field(default_factory=list)
    voxel_indices: np.ndarray | None = None
    schema_version: str = SCHEMA_VERSION

@dataclass
class Portal:
    id: int
    regions: tuple[int, int]
    face_count: int
    area_nm2: float
    max_clearance_nm: float
    certified: bool
    image_offsets: list[tuple[int, int, int]]

@dataclass
class RegionGraph:
    nodes: list[int]
    portals: list[Portal]
    schema_version: str = SCHEMA_VERSION

@dataclass
class RegionHierarchy:
    components: list[VoidComponent]
    aggregates: list[dict]
    scale_levels_nm: list[float]
    scale_components: list[dict]
    schema_version: str = SCHEMA_VERSION

@dataclass
class PathwayNetwork:
    networks: list[dict]
    schema_version: str = SCHEMA_VERSION

@dataclass
class Atlas:
    geometry: SystemGeometry
    exclusion: ExclusionModel
    field: FreeSpaceField
    membrane: MembraneAtlas
    labels: np.ndarray
    regions: list[Region]
    graph: RegionGraph
    hierarchy: RegionHierarchy
    pathways: PathwayNetwork
    warnings: list[str]
    timings: dict = field(default_factory=dict)
    schema_version: str = SCHEMA_VERSION
    components: list[VoidComponent] = field(default_factory=list)

@dataclass
class WaterMapping:
    molecule_id: int
    region_id: int | None
    classification: str
    hard_clash: bool
    membrane_interior: bool
    uncertainty: list[str]
    schema_version: str = SCHEMA_VERSION
    molecule_uid: int | None = None
    atom_indices: tuple[int, ...] = ()

@dataclass
class CleanupDecision:
    molecule_id: int
    action: str
    reason: str
    region_id: int | None
    policy: str
    schema_version: str = SCHEMA_VERSION

