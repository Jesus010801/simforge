"""Phase 3.5 — the semantic index must provide the ``System`` group that the
center / fit preprocessing steps select (``stdin = "<group>\\nSystem"`` with
``-n semantic_index.ndx``).

Before the fix, ``build_semantic_index()`` never emitted ``System``; GROMACS
then re-prompted for the second group and died with ``Cannot read from input``,
so every center/fit view built through ``simforge study`` failed.

Golden atom-list hashes below were captured from the *pre-fix* builder, so the
tests prove the existing groups are unchanged (membership and order).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, ComponentType, SemanticIndex, ViewBuildStatus,
)
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import ObservableSpec
from analysis.campaign.structure.component_detector import detect_components
from analysis.campaign.structure.index_builder import build_semantic_index
from analysis.campaign.structure.index_groups import parse_index_groups
from analysis.campaign.trajectory import requirements as reqs
from analysis.campaign.trajectory.preprocessor import (
    _identity_evidence, _plan_steps, load_view_time_index, view_identity,
)
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj, write_gro, write_pdb,
)


def _ordered(ids) -> str:
    """Order- and multiplicity-preserving hash of an atom list."""
    return hashlib.sha256(",".join(map(str, ids)).encode()).hexdigest()[:24]


def _groups(ndx) -> dict[str, tuple[int, str]]:
    return {g.name: (len(g.atom_ids), _ordered(g.atom_ids)) for g in parse_index_groups(ndx)}


def _natoms(structure: Path) -> int:
    if structure.suffix == ".gro":
        return int(structure.read_text().splitlines()[1])
    return sum(1 for l in structure.read_text().splitlines() if l[:6].strip() in ("ATOM", "HETATM"))


# Captured from the pre-fix builder (analysis/campaign/structure/index_builder.py @ Phase 3).
_GOLDEN_GRO = {
    "Water": (2000, "d51772f7b6a5b564130061de"), "Ions": (10, "e7d4d87cd548803e54199629"),
    "Receptor": (50, "dd63df3042e7e7d49daa23ff"),
    "Receptor_Backbone": (50, "dd63df3042e7e7d49daa23ff"),
    "Receptor_CA": (50, "dd63df3042e7e7d49daa23ff"),
    "Complex": (50, "dd63df3042e7e7d49daa23ff"),
    "Complex_Backbone": (50, "dd63df3042e7e7d49daa23ff"),
    "Complex_CA": (50, "dd63df3042e7e7d49daa23ff"),
}
_GOLDEN_PDB = {
    "Receptor": (120, "1c4732eaf0b754ef54c34467"),
    "Receptor_Backbone": (120, "1c4732eaf0b754ef54c34467"),
    "Receptor_CA": (120, "1c4732eaf0b754ef54c34467"),
    "Peptide": (25, "c0460c6e5ab94f9bfe36946b"),
    "Peptide_Backbone": (25, "c0460c6e5ab94f9bfe36946b"),
    "Peptide_CA": (25, "c0460c6e5ab94f9bfe36946b"),
    "Complex": (145, "4742ac2d3eddc7bc7511bd81"),
    "Complex_Backbone": (145, "4742ac2d3eddc7bc7511bd81"),
    "Complex_CA": (145, "4742ac2d3eddc7bc7511bd81"),
}
_GOLDEN_REAL = {
    "Water": (254700, "e6d4a40ff4f4ebff887baeb9"), "Ions": (610, "54d35bd580244f71e153ba4b"),
    "Membrane": (21450, "1d6003d0f841519469f8d0e1"),
    "Receptor": (6571, "b006e34cdcd447a56cd6b14e"),
    "Receptor_Backbone": (1200, "7f8b60565a58b870cd2f2daf"),
    "Receptor_CA": (400, "76633b601ab94ecaac8553e0"),
    "Complex": (6571, "b006e34cdcd447a56cd6b14e"),
    "Complex_Backbone": (1200, "7f8b60565a58b870cd2f2daf"),
    "Complex_CA": (400, "76633b601ab94ecaac8553e0"),
}


def _build(structure: Path, out: Path) -> SemanticIndex:
    det = detect_components(structure_path=structure)
    return build_semantic_index(structure_path=structure, components=det.components,
                                out_ndx=out)


@pytest.fixture
def gro(tmp_path):
    return write_gro(tmp_path / "solv.gro")


@pytest.fixture
def pdb(tmp_path):
    return write_pdb(tmp_path / "rp.pdb", {"A": 120, "B": 25})


# ═══════════════════════════════════════════════════════════════════════════════
# System group
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
@pytest.mark.parametrize("which", ["gro", "pdb"])
def test_system_emitted_with_every_atom_once_in_order(request, tmp_path, which):
    structure = request.getfixturevalue(which)
    idx = _build(structure, tmp_path / "semantic.ndx")
    groups = {g.name: g for g in parse_index_groups(idx.path)}
    assert "System" in groups and "System" in idx.group_names()
    n = _natoms(structure)
    assert groups["System"].atom_ids == tuple(range(1, n + 1))     # all, once, topology order
    meta = next(g for g in idx.groups if g.name == "System")
    assert meta.n_atoms == n


@requires_gmx
def test_system_is_deterministic(gro, tmp_path):
    a = _build(gro, tmp_path / "a.ndx")
    b = _build(gro, tmp_path / "b.ndx")
    assert Path(a.path).read_text() == Path(b.path).read_text()


@requires_gmx
@pytest.mark.parametrize("which,golden", [("gro", _GOLDEN_GRO), ("pdb", _GOLDEN_PDB)])
def test_existing_groups_unchanged(request, tmp_path, which, golden):
    idx = _build(request.getfixturevalue(which), tmp_path / "semantic.ndx")
    groups = _groups(idx.path)
    assert {k: v for k, v in groups.items() if k != "System"} == golden
    # System is appended last: positions of the pre-existing groups are unchanged
    assert list(groups)[:-1] == list(golden) and list(groups)[-1] == "System"


def test_system_is_not_a_molecular_component():
    assert "system" not in {v for k, v in vars(ComponentType).items() if not k.startswith("_")}


# ═══════════════════════════════════════════════════════════════════════════════
# View identity implications
# ═══════════════════════════════════════════════════════════════════════════════

def _evidence(tmp_path, ndx_text: str, req):
    from analysis.campaign.fingerprint import fingerprint_file
    ndx = tmp_path / "semantic.ndx"
    ndx.write_text(ndx_text)
    traj = tmp_path / "md.xtc"
    traj.write_bytes(b"traj")
    tpr = tmp_path / "md.tpr"
    tpr.write_bytes(b"tpr")
    names = [l.strip("[] \n") for l in ndx_text.splitlines() if l.startswith("[")]
    from analysis.campaign.models import SemanticIndexGroup
    sidx = SemanticIndex(path=str(ndx), groups=[SemanticIndexGroup(n, 0) for n in names])
    steps, err = _plan_steps(req, str(traj), str(tpr), None, sidx, tmp_path / "out")
    assert err is None
    fps = {str(tpr): fingerprint_file(tpr)}
    ev, groups = _identity_evidence(req, req.view_kind(), steps, [fingerprint_file(traj)],
                                    fps, str(ndx), "v")
    return view_identity(ev), groups


_BASE = "[ Receptor ]\n1 2 3\n[ Receptor_Backbone ]\n1 2\n"


def test_adding_system_does_not_change_views_that_do_not_use_it(tmp_path):
    for req in (reqs.whole_only("w"), reqs.nojump("n")):
        a, _ = _evidence(tmp_path, _BASE, req)
        b, _ = _evidence(tmp_path, _BASE + "[ System ]\n1 2 3 4 5\n", req)
        assert a == b


def test_center_fit_identity_includes_system_atom_set(tmp_path):
    from analysis.campaign.results import atom_set_hash
    for req in (reqs.fit_to("Receptor_Backbone", "f"), reqs.centered_on("Receptor", "c")):
        k5, groups = _evidence(tmp_path, _BASE + "[ System ]\n1 2 3 4 5\n", req)
        # the whole step's System comes from the -s default groups (no -n);
        # the center/fit step's System is read from the index -> atom evidence
        system = [g for g in groups if g["name"] == "System" and "atoms_sha256" in g]
        assert system and system[0]["atoms_sha256"] == atom_set_hash([1, 2, 3, 4, 5])
        k6, _ = _evidence(tmp_path, _BASE + "[ System ]\n1 2 3 4 5 6\n", req)
        assert k5 != k6


# ═══════════════════════════════════════════════════════════════════════════════
# The previously failing study path
# ═══════════════════════════════════════════════════════════════════════════════

class _ReceptorFitProbe(ObservableSpec):
    """Requests exactly the view rmsd-peptide-receptor-frame needs (receptor-backbone
    fit) on a receptor-only system, so the study path reaches preprocessing."""
    id = "probe-receptor-fit"
    required_components = (ComponentType.RECEPTOR,)

    def trajectory_requirements(self, params):
        return reqs.fit_to("Receptor_Backbone", "probe")

    def execute(self, ctx):
        v = ctx.trajectory_view
        ok = v.safe and v.path is not None
        return AnalysisResult(self.id, ctx.system.system_id,
                              AnalysisStatus.SUCCESS if ok else AnalysisStatus.FAILED,
                              message="; ".join(w.message for w in v.warnings))


@pytest.fixture
def real_study(tmp_path):
    d = tmp_path / "campaign" / "semaglutide" / "rep1"
    d.mkdir(parents=True)
    (d / "production_run.xtc").symlink_to(REAL_XTC)
    (d / "system_topol.tpr").symlink_to(REAL_TPR)
    (d / "reference.gro").symlink_to(REAL_GRO)
    return tmp_path / "campaign"


@requires_gmx
@requires_real_traj
def test_study_fit_view_builds_on_real_system(real_study):
    from analysis.campaign.orchestration.study_analyzer import run_analyze

    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in (REAL_XTC, REAL_TPR, REAL_GRO)}
    registry.register(_ReceptorFitProbe())
    from analysis.campaign.models import IntentSource
    from analysis.campaign.orchestration.study_analyzer import clear_view_cache
    from analysis.campaign.trajectory.policy import PolicyIntent

    # Phase 5: the fit reference (.gro while a .tpr exists) is ambiguous, so the
    # policy requires explicit intent — without it nothing is built.
    blocked = run_analyze(real_study, ["probe-receptor-fit"], inspect_trajectories=True)
    assert blocked.results[0].status == AnalysisStatus.FAILED
    assert "fit.reference_ambiguous" in blocked.results[0].message
    clear_view_cache()

    res = run_analyze(real_study, ["probe-receptor-fit"], inspect_trajectories=True,
                      policy_intent=PolicyIntent(IntentSource.EXPLICIT_API))
    (r,) = res.results
    assert r.status == AnalysisStatus.SUCCESS, r.message

    prov = json.loads(Path(r.provenance_path).read_text())
    view = prov["trajectory_view"]
    assert view["build_status"] == ViewBuildStatus.BUILT
    assert [o["operation"] for o in view["operations"]] == ["make_whole", "fit"]
    assert all(o["returncode"] == 0 for o in view["operations"])
    assert view["operations"][1]["stdin"] == "Receptor_Backbone\nSystem"

    ndx = Path(prov["semantic_index"]["path"])
    groups = _groups(ndx)
    assert {k: v for k, v in groups.items() if k != "System"} == _GOLDEN_REAL
    assert groups["System"][0] == 283331
    assert next(g for g in parse_index_groups(ndx) if g.name == "System").atom_ids == \
        tuple(range(1, 283332))

    manifest = json.loads(Path(view["manifest_path"]).read_text())
    assert manifest["status"] == "complete" and manifest["view_identity"] == view["cache_key"]
    sys_ev = [g for g in manifest["groups"] if g["name"] == "System" and "n_atoms" in g]
    assert sys_ev and sys_ev[0]["n_atoms"] == 283331
    assert view["time_index_ref"]["n_frames"] == 51
    assert view["time_index_ref"]["same_timeline_as_source"] is True

    from analysis.campaign.models import TrajectoryRequirements, TrajectoryView
    tv = TrajectoryView(kind=view["kind"], requirements=TrajectoryRequirements(),
                        path=view["path"], time_index_ref=view["time_index_ref"])
    assert load_view_time_index(tv).n_frames == 51

    after = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in (REAL_XTC, REAL_TPR, REAL_GRO)}
    assert after == before
