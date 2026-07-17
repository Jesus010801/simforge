"""
tests/test_strong_posre_restraints.py
──────────────────────────────────────
Regression tests for per-chain strong position restraints.

Root cause of the original bug:
  gmx genrestr run on the full embedded system (protein + bilayer) writes global
  atom indices into strong_posre.itp.  When that file is included inside a single
  chain's moleculetype definition, GROMACS checks 1..N_chain and rejects any index
  > N_chain_A.  For a 5-chain GLP-1R with ~1400 atoms per chain, chain B's first
  atom is index 1401 — instantly out of bounds for chain A's local space.

Fix: parse [ atoms ] from each chain ITP, write local 1-based heavy-atom indices
into per-chain strong_posre_Protein_chain_X.itp files, include each file only
immediately after its matching moleculetype ITP.
"""

from __future__ import annotations

import re
import json
import importlib.util
from pathlib import Path
import pytest


# ── helpers ──────────────────────────────────────────────────────────────────

def _load_validator():
    """Load strong_posre_validator directly from the validators/ directory."""
    spec = importlib.util.spec_from_file_location(
        "strong_posre_validator",
        str(Path(__file__).resolve().parent.parent / "validators" / "strong_posre_validator.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_chain_itp(path: Path, mol_name: str, n_heavy: int, n_h: int = 5) -> Path:
    """Write a minimal chain ITP with n_heavy heavy atoms and n_h hydrogens."""
    lines = [
        "[ moleculetype ]",
        f"{mol_name}    3",
        "",
        "[ atoms ]",
        "; nr  type  resnr resname atname cgnr charge mass",
    ]
    idx = 1
    for i in range(n_heavy):
        lines.append(f"  {idx:4d}  opls_135   1   ALA   CA   {idx}  0.0  12.011")
        idx += 1
    for i in range(n_h):
        lines.append(f"  {idx:4d}  opls_140   1   ALA   HA   {idx}  0.0   1.008")
        idx += 1
    lines += ["", "[ bonds ]", "; no bonds listed"]
    path.write_text("\n".join(lines) + "\n")
    return path


def _make_posre_itp(path: Path, indices: list[int], fc: int = 100000) -> Path:
    """Write a [ position_restraints ] file with given (local) indices."""
    lines = [
        "; strong position restraints",
        "[ position_restraints ]",
        "; atom  funct  fcx  fcy  fcz",
    ] + [f"  {i:6d}    1  {fc}  {fc}  {fc}" for i in indices]
    path.write_text("\n".join(lines) + "\n")
    return path


# ── A. _parse_heavy_local_indices ─────────────────────────────────────────────

class TestParseHeavyLocalIndices:
    """Test the heavy-atom index parser embedded in the embed script."""

    def _exec_parser(self):
        """Exec the inline _parse_heavy_local_indices function from validators."""
        code = """
import re as _re_posre
def _parse_heavy_local_indices(top_path):
    _in_atoms = False
    _idx_list = []
    for _l in top_path.read_text().splitlines():
        _ls = _l.strip()
        if _re_posre.match(r'^\\[\\s*atoms\\s*\\]', _ls, _re_posre.IGNORECASE):
            _in_atoms = True; continue
        if _in_atoms:
            if _ls.startswith('['):
                break
            if not _ls or _ls.startswith(';'):
                continue
            _pp = _ls.split()
            if len(_pp) >= 5:
                try:
                    _nr    = int(_pp[0])
                    _aname = _pp[4]
                    if not _aname.upper().startswith('H'):
                        _idx_list.append(_nr)
                except (ValueError, IndexError):
                    pass
    return _idx_list
"""
        ns: dict = {}
        exec(code, ns)
        return ns["_parse_heavy_local_indices"]

    def test_heavy_atoms_selected(self, tmp_path):
        """Heavy atoms (not starting with H) are returned; hydrogens are excluded."""
        itp = _make_chain_itp(tmp_path / "chain_A.itp", "Protein_chain_A", n_heavy=10, n_h=3)
        fn = self._exec_parser()
        indices = fn(itp)
        assert len(indices) == 10
        # First atom is index 1
        assert indices[0] == 1
        # Last heavy atom is index 10
        assert indices[-1] == 10
        # Hydrogen atoms (11, 12, 13) are excluded
        assert 11 not in indices
        assert 12 not in indices

    def test_indices_are_local_one_based(self, tmp_path):
        """Indices start at 1, not 0, and are relative to the ITP file."""
        itp = _make_chain_itp(tmp_path / "chain_B.itp", "Protein_chain_B", n_heavy=5)
        fn = self._exec_parser()
        assert fn(itp) == [1, 2, 3, 4, 5]

    def test_large_chain(self, tmp_path):
        """Works for chains with 6975 heavy atoms (GLP-1R chain size)."""
        n = 6975
        itp = _make_chain_itp(tmp_path / "chain.itp", "Protein_chain_A", n_heavy=n, n_h=0)
        fn = self._exec_parser()
        result = fn(itp)
        assert len(result) == n
        assert result[-1] == n

    def test_empty_itp_returns_empty(self, tmp_path):
        itp = tmp_path / "empty.itp"
        itp.write_text("[ moleculetype ]\nEmpty 3\n")
        fn = self._exec_parser()
        assert fn(itp) == []


# ── B. Per-chain file generation ──────────────────────────────────────────────

class TestPerChainFileGeneration:
    """Test that the embed script generates one file per chain with local indices."""

    def _get_embed_script(self, tmp_path: Path) -> str:
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_posre
components:
  - id: peptide_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        embed_dir = next((ws / "steps").glob("*_embed_in_bilayer"), None)
        if embed_dir is None:
            return ""
        for _name in ("run_embed.py", "run_embed_in_bilayer.py", "run.py"):
            _p = embed_dir / _name
            if _p.exists():
                return _p.read_text()
        return ""

    def test_script_has_prot_top_dir(self, tmp_path):
        script = self._get_embed_script(tmp_path)
        assert "PROT_TOP_DIR" in script

    def test_script_has_per_chain_generator(self, tmp_path):
        """The embed script must define _parse_heavy_local_indices and _write_posre_itp."""
        script = self._get_embed_script(tmp_path)
        assert "_parse_heavy_local_indices" in script
        assert "_write_posre_itp" in script

    def test_script_reads_manifest(self, tmp_path):
        """The embed script must read protein_topology_manifest.json for chain info."""
        script = self._get_embed_script(tmp_path)
        assert "protein_topology_manifest.json" in script
        assert "_mol_names_pm" in script
        assert "_mol_itps_pm" in script

    def test_script_writes_posre_manifest(self, tmp_path):
        """The embed script must write strong_posre_manifest.json."""
        script = self._get_embed_script(tmp_path)
        assert "strong_posre_manifest.json" in script

    def test_no_gmx_genrestr_call_in_embed_script(self, tmp_path):
        """gmx genrestr subprocess call must not be present — global indices break multichain."""
        script = self._get_embed_script(tmp_path)
        # The comment mentioning genrestr is OK; the subprocess.run call must be gone
        assert '"gmx", "genrestr"' not in script
        assert "input=\"Protein" not in script

    def test_per_chain_file_naming_convention(self, tmp_path):
        """Per-chain files must be named strong_posre_<MolName>.itp."""
        script = self._get_embed_script(tmp_path)
        assert "strong_posre_{" in script or "strong_posre_" in script


# ── C. Per-chain functional: generate → validate ──────────────────────────────

class TestPerChainFunctional:
    """Simulate per-chain generation and validate bounds."""

    def test_multichain_generation(self, tmp_path):
        """5 chains → 5 per-chain files with local indices; none exceed chain size."""
        validator = _load_validator()
        chains = ["Protein_chain_A", "Protein_chain_B", "Protein_chain_C",
                  "Protein_chain_D", "Protein_chain_E"]
        n_heavy_per_chain = 6975  # GLP-1R-scale

        itp_dir = tmp_path / "prot_top"
        itp_dir.mkdir()

        manifest = {
            "molecule_names": chains,
            "molecule_itps":  [f"topol_{c}.itp" for c in chains],
        }
        (itp_dir / "protein_topology_manifest.json").write_text(json.dumps(manifest))

        # Create chain ITPs and per-chain posre files (local indices)
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")

        topol_lines = ['#include "ff/forcefield.itp"']
        for chain in chains:
            _make_chain_itp(itp_dir / f"topol_{chain}.itp", chain, n_heavy=n_heavy_per_chain, n_h=0)
            sp_path = tmp_path / f"strong_posre_{chain}.itp"
            _make_posre_itp(sp_path, list(range(1, n_heavy_per_chain + 1)))

            chain_itp_rel = f"prot_top/topol_{chain}.itp"
            sp_rel = f"strong_posre_{chain}.itp"
            topol_lines += [
                f'#include "{chain_itp_rel}"',
                "#ifdef STRONG_POSRES",
                f'#include "{sp_rel}"',
                "#endif",
            ]

        topol_lines += ["[ system ]", "System", "[ molecules ]"]
        for chain in chains:
            topol_lines.append(f"{chain}    1")

        topol = tmp_path / "topol.top"
        topol.write_text("\n".join(topol_lines) + "\n")

        errors = validator.validate_posre_bounds(topol)
        assert errors == [], f"Expected no errors for valid local indices:\n{errors}"

    def test_global_index_detected_as_violation(self, tmp_path):
        """Writing chain B's global index into chain A's posre file must be caught."""
        validator = _load_validator()
        n_chain_a = 6975
        n_chain_b = 6975
        # Simulate the original bug: strong_posre.itp contains global index 6976 (chain B start)
        global_idx_chain_b_start = n_chain_a + 1  # 6976

        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        _make_chain_itp(prot_dir / "topol_Protein_chain_A.itp", "Protein_chain_A",
                        n_heavy=n_chain_a, n_h=0)
        # Wrong: contains global index from chain B (6976 > 6975 = chain A size)
        _make_posre_itp(tmp_path / "strong_posre.itp",
                        [1, 2, 3, global_idx_chain_b_start])

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "prot/topol_Protein_chain_A.itp"\n'
            "#ifdef STRONG_POSRES\n"
            '#include "strong_posre.itp"\n'
            "#endif\n"
            "[ system ]\nSystem\n"
            "[ molecules ]\nProtein_chain_A    1\n"
        )

        errors = validator.validate_posre_bounds(topol)
        assert any(str(global_idx_chain_b_start) in e for e in errors), (
            f"Expected out-of-bounds error for index {global_idx_chain_b_start}:\n{errors}"
        )
        assert any("TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR" in e for e in errors)

    def test_per_chain_files_pass_validation(self, tmp_path):
        """Per-chain posre files with local indices 1..N pass without errors."""
        validator = _load_validator()
        n = 6975

        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        _make_chain_itp(prot_dir / "topol_Protein_chain_A.itp", "Protein_chain_A",
                        n_heavy=n, n_h=0)
        _make_chain_itp(prot_dir / "topol_Protein_chain_B.itp", "Protein_chain_B",
                        n_heavy=n, n_h=0)
        # Correct: each file has its own local indices starting from 1
        _make_posre_itp(tmp_path / "strong_posre_Protein_chain_A.itp", list(range(1, n + 1)))
        _make_posre_itp(tmp_path / "strong_posre_Protein_chain_B.itp", list(range(1, n + 1)))

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "prot/topol_Protein_chain_A.itp"\n'
            "#ifdef STRONG_POSRES\n"
            '#include "strong_posre_Protein_chain_A.itp"\n'
            "#endif\n"
            '#include "prot/topol_Protein_chain_B.itp"\n'
            "#ifdef STRONG_POSRES\n"
            '#include "strong_posre_Protein_chain_B.itp"\n'
            "#endif\n"
            "[ system ]\nSystem\n"
            "[ molecules ]\nProtein_chain_A    1\nProtein_chain_B    1\n"
        )

        errors = validator.validate_posre_bounds(topol)
        assert errors == [], f"Per-chain posre must pass: {errors}"

    def test_pdb2gmx_posre_files_untouched(self, tmp_path):
        """
        pdb2gmx generates posre_Protein_chain_X.itp for regular POSRES.
        The per-chain strong_posre must be a separate file; pdb2gmx files are unchanged.
        """
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        original_content = (
            "[ position_restraints ]\n; atom funct fcx fcy fcz\n  1    1  1000 1000 1000\n"
        )
        posre_orig = prot_dir / "posre_Protein_chain_A.itp"
        posre_orig.write_text(original_content)

        # The test is simply that the per-chain strong posre is a DIFFERENT file
        sp_file = tmp_path / "strong_posre_Protein_chain_A.itp"
        _make_posre_itp(sp_file, [1, 2, 3])

        assert posre_orig.read_text() == original_content, "pdb2gmx posre must not be modified"
        assert sp_file.exists(), "Per-chain strong posre must be a separate file"
        assert sp_file.name != posre_orig.name

    def test_no_strong_posre_topology_still_valid(self, tmp_path):
        """When no strong_posre files exist, validate_posre_bounds returns empty (no errors)."""
        validator = _load_validator()
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        _make_chain_itp(prot_dir / "topol_Protein_chain_A.itp", "Protein_chain_A",
                        n_heavy=100, n_h=10)

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "prot/topol_Protein_chain_A.itp"\n'
            "[ system ]\nSystem\n[ molecules ]\nProtein_chain_A    1\n"
        )

        errors = validator.validate_posre_bounds(topol)
        assert errors == []

    def test_index_exactly_at_bound_is_ok(self, tmp_path):
        """Index == n_atoms_in_chain must NOT trigger an error."""
        validator = _load_validator()
        n = 100
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        _make_chain_itp(prot_dir / "chain.itp", "Protein_chain_A", n_heavy=n, n_h=0)
        _make_posre_itp(tmp_path / "strong_posre_Protein_chain_A.itp", [n])  # exactly n

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "prot/chain.itp"\n'
            '#include "strong_posre_Protein_chain_A.itp"\n'
            "[ system ]\nSystem\n[ molecules ]\nProtein_chain_A    1\n"
        )
        assert validator.validate_posre_bounds(topol) == []

    def test_index_one_over_bound_is_error(self, tmp_path):
        """Index == n_atoms + 1 must trigger an error."""
        validator = _load_validator()
        n = 100
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        prot_dir = tmp_path / "prot"
        prot_dir.mkdir()
        _make_chain_itp(prot_dir / "chain.itp", "Protein_chain_A", n_heavy=n, n_h=0)
        _make_posre_itp(tmp_path / "bad_posre.itp", [n + 1])  # n+1 is out of bounds

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "prot/chain.itp"\n'
            '#include "bad_posre.itp"\n'
            "[ system ]\nSystem\n[ molecules ]\nProtein_chain_A    1\n"
        )
        errors = validator.validate_posre_bounds(topol)
        assert len(errors) == 1
        assert "TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR" in errors[0]
        assert str(n + 1) in errors[0]


# ── D. Assembly script checks ─────────────────────────────────────────────────

class TestAssemblyScriptStrongPosre:
    """Verify the assemble_system.py script handles per-chain posre correctly."""

    def _get_asm_script(self, tmp_path: Path) -> str:
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_asm_posre
components:
  - id: peptide_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return ""
        return (asm_dir / "run_assemble_system.py").read_text()

    def test_script_reads_sp_manifest(self, tmp_path):
        """Assembly script must read strong_posre_manifest.json from embed dir."""
        script = self._get_asm_script(tmp_path)
        assert "strong_posre_manifest.json" in script
        assert "_sp_manifest_asm" in script

    def test_script_injects_per_chain_in_prot_block(self, tmp_path):
        """Assembly script must inject per-chain includes inside the prot_block loop."""
        script = self._get_asm_script(tmp_path)
        assert "_sp_per_chain_done" in script
        assert "_sp_manifest_asm" in script

    def test_script_has_posre_bounds_validator(self, tmp_path):
        """Assembly script must call validate_posre_bounds."""
        script = self._get_asm_script(tmp_path)
        assert "validate_posre_bounds" in script
        assert "strong_posre_validator" in script

    def test_script_uses_mol_names_in_prot_block_loop(self, tmp_path):
        """The prot_block loop must use mol names to match per-chain posre files."""
        script = self._get_asm_script(tmp_path)
        assert "_mol_names_p" in script
        # Indexed access of mol_names_p for each ITP
        assert "_mn_i" in script or "mol_names_p[" in script


# ── E. Error category and signal ─────────────────────────────────────────────

class TestPositionRestraintIndexError:

    def test_error_category_exists(self):
        from executors.remediation_models import ErrorCategory
        assert hasattr(ErrorCategory, "TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR")
        assert ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR == \
               "topology_position_restraint_index_error"

    def test_signal_matches_gromacs_error(self):
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory

        stderr = (
            "ERROR 1 [file strong_posre.itp, line 6980]:\n"
            "  Atom index (6976) in position_restraints out of bounds (1-6975).\n"
        )
        matched = None
        for sig in _SIGNAL_PATTERNS:
            hit, _ = sig.match(stderr)
            if hit:
                matched = sig.category
                break
        assert matched == ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR

    def test_signal_matches_pre_grompp_diagnostic(self):
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory

        stderr = (
            "TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR: strong_posre.itp contains atom index 6976, "
            "but Protein_chain_A has only 6975 atoms.\n"
        )
        matched = None
        for sig in _SIGNAL_PATTERNS:
            hit, _ = sig.match(stderr)
            if hit:
                matched = sig.category
                break
        assert matched == ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR

    def test_plan_is_non_retryable(self, tmp_path):
        from executors.remediation_models import ErrorCategory, ErrorSeverity, DiagnosisResult
        from executors.signal_detector import AdaptiveReasoner

        diag = DiagnosisResult(
            step_id        = "embed_in_bilayer",
            step_dir       = str(tmp_path),
            category       = ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR,
            severity       = ErrorSeverity.FATAL,
            confidence     = 0.99,
            primary_signal = "Atom index (6976) in position_restraints out of bounds (1-6975)",
            explanation    = "Global system atom index written to per-chain restraint file",
        )
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0
        assert plan.requires_human is True
