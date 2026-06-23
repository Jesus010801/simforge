# builders/workspace_builder.py

from __future__ import annotations

from pathlib import Path
import json
import shutil
import uuid
from datetime import datetime

from core.compiler_models import (
    CompilationResult,
)
from builders.step_builders.minimization_builder import (
    MinimizationBuilder,
)
from builders.builder_registry import (
    STEP_BUILDERS,
)
from core.workspace_fingerprint import (
    compute_build_signature,
    compute_builder_signature,
    compute_template_hash,
    SIMFORGE_VERSION,
)


def _find_simforge_root() -> Path:
    p = Path(__file__).parent
    for _ in range(8):
        if (p / "adapters").is_dir() and (p / "core").is_dir():
            return p
        p = p.parent
    raise RuntimeError("SimForge project root not found")


_PROT_MEMB_FILES = _find_simforge_root() / "docs" / "Prot-Memb_FILES"

# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _stamp_step(step_dir: Path, generated_at: str) -> None:
    """Add provenance fields (template_hash, generated_at) to step metadata.json."""
    meta_path = step_dir / "metadata.json"
    if not meta_path.exists():
        return
    try:
        meta = json.loads(meta_path.read_text())
        scripts = list(step_dir.glob("*.sh")) + list(step_dir.glob("*.py"))
        meta["template_hash"] = compute_template_hash(scripts)
        meta["generated_at"]  = generated_at
        meta_path.write_text(json.dumps(meta, indent=4))
    except Exception:
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# Workspace Builder
# ═══════════════════════════════════════════════════════════════════════════════

class WorkspaceBuilder:
    """
    Construye un workspace físico reproducible
    desde un CompilationResult.
    """

    def build(
        self,
        result: CompilationResult,
        output_dir: str = "simforge_runs",
        workspace_name: str | None = None,
        yaml_source: str = "",
        workspace_path: "str | Path | None" = None,
    ) -> Path:
        """
        Args:
            result:          CompilationResult from the compiler.
            output_dir:      Root directory for all workspaces.
            workspace_name:  Override the workspace directory name.
                             Defaults to inferred_system_type.
            yaml_source:     Absolute path to the original YAML config.
            workspace_path:  Explicit run directory (overrides output_dir +
                             workspace_name).  Used by the CLI to pass a
                             pre-created timestamped run directory so that
                             each compile lands in its own immutable location.
        """
        if workspace_path is not None:
            root = Path(workspace_path)
        else:
            system_name = workspace_name or (
                result.state.inferred_system_type
                or "simforge_system"
            )
            root = Path(output_dir) / system_name

        # ────────────────────────────────────────────────────────────────────
        # Core directories
        # ────────────────────────────────────────────────────────────────────

        workflow_dir = root / "workflow"

        reports_dir = root / "reports"

        metadata_dir = root / "metadata"

        steps_dir = root / "steps"

        workflow_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        reports_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        metadata_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        steps_dir.mkdir(parents=True, exist_ok=True)

        # ────────────────────────────────────────────────────────────────────
        # inputs/ — stage all external source files into the workspace
        # ────────────────────────────────────────────────────────────────────

        inputs_dir = root / "inputs"
        inputs_dir.mkdir(exist_ok=True)
        self._stage_inputs(result, inputs_dir)

        # ────────────────────────────────────────────────────────────────────
        # membrane_assets/ — stage bilayer GRO, inflategro, and FF directory
        # into the workspace so scripts need no external path at runtime.
        # ────────────────────────────────────────────────────────────────────

        membrane_assets_dir: Path | None = None
        _needs_membrane = any(
            step.step_id in ("embed_in_bilayer", "membrane_embedding")
            for step in result.execution_order
        )
        if _needs_membrane:
            membrane_assets_dir = root / "membrane_assets"
            membrane_assets_dir.mkdir(exist_ok=True)
            self._stage_membrane_assets(result, membrane_assets_dir, metadata_dir)

        # ────────────────────────────────────────────────────────────────────
        # Step folders — pass 1: create dirs + build step_dir_map
        # ────────────────────────────────────────────────────────────────────

        step_dir_map: dict[str, Path] = {
            "__workspace_root__": root.resolve(),
        }
        if membrane_assets_dir is not None:
            step_dir_map["__membrane_assets__"] = membrane_assets_dir.resolve()

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            step_dir.mkdir(exist_ok=True)
            step_dir_map[step.step_id] = step_dir.resolve()

        # ────────────────────────────────────────────────────────────────────
        # Step materialization — pass 2: builders receive full map
        # ────────────────────────────────────────────────────────────────────

        generated_at = datetime.now().isoformat(timespec="seconds")

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            builder = STEP_BUILDERS.get(step.stage.value)
            if builder is not None:
                builder.build(step, step_dir, step_dir_map)
                _stamp_step(step_dir, generated_at)

        # ────────────────────────────────────────────────────────────────────
        # Mermaid workflow
        # ────────────────────────────────────────────────────────────────────

        mermaid_path = (
            workflow_dir
            / "workflow.mmd"
        )

        mermaid_path.write_text(
            result.mermaid_graph
        )

        # ────────────────────────────────────────────────────────────────────
        # User workflow
        # ────────────────────────────────────────────────────────────────────

        workflow_txt = (
            workflow_dir
            / "workflow.txt"
        )

        lines = []

        for i, step in enumerate(
            result.user_view,
            start=1,
        ):

            lines.append(
                f"{i:02d}. {step}"
            )

        workflow_txt.write_text(
            "\n".join(lines)
        )

        # ────────────────────────────────────────────────────────────────────
        # Summary metadata
        # ────────────────────────────────────────────────────────────────────

        summary_json = (
            metadata_dir
            / "summary.json"
        )

        summary_data = {
            "system_type": (
                result.state.inferred_system_type
            ),

            "workflow_steps": (
                len(result.plan.steps)
            ),

            "blocking_issues": (
                len(result.plan.blocking_issues)
            ),

            "special_protocols": (
                result.plan.special_protocols
            ),
        }

        summary_json.write_text(
            json.dumps(
                summary_data,
                indent=4,
            )
        )

        # ────────────────────────────────────────────────────────────────────
        # Execution manifest — fuente de verdad para el executor
        # ────────────────────────────────────────────────────────────────────

        manifest_entries = []
        for i, step in enumerate(result.execution_order, start=1):
            manifest_entries.append({
                "step_id":         step.step_id,
                "dir_name":        f"{i:02d}_{step.step_id}",
                "stage":           step.stage.value,
                "step_type":       step.step_type.value,
                "automation_level": step.effective_automation_level().value,
                "blocking":        step.blocking,
                "depends_on":      step.depends_on,
            })

        compiled_at       = datetime.now().isoformat(timespec="seconds")
        builder_sig       = compute_builder_signature()
        build_sig         = compute_build_signature(yaml_source) if yaml_source else builder_sig

        manifest_data = {
            "compile_id":         str(uuid.uuid4()),
            "compiled_at":        compiled_at,
            "simforge_version":   SIMFORGE_VERSION,
            "yaml_source":        str(Path(yaml_source).resolve()) if yaml_source else "",
            "builder_signature":  builder_sig,
            "build_signature":    build_sig,
            "system_type":        result.state.inferred_system_type,
            "n_steps":            len(manifest_entries),
            "steps":              manifest_entries,
        }

        manifest_path = metadata_dir / "execution_manifest.json"
        manifest_path.write_text(json.dumps(manifest_data, indent=4))

        # ────────────────────────────────────────────────────────────────────
        # Compile report (professional scientific summary)
        # ────────────────────────────────────────────────────────────────────
        try:
            from core.report_generator import generate_compile_report
            generate_compile_report(result, root)
        except Exception:
            pass   # report failure must never block workspace creation

        return root

    # ────────────────────────────────────────────────────────────────────────
    # Input staging
    # ────────────────────────────────────────────────────────────────────────

    def _stage_membrane_assets(
        self,
        result: "CompilationResult",
        assets_dir: Path,
        metadata_dir: Path,
    ) -> None:
        """
        Copy dppc*.gro, inflategro-Jorge.pl, and oplsaa_membrane.ff/ from
        docs/Prot-Memb_FILES/ into workspace/membrane_assets/.

        Raises FileNotFoundError at compile time if any required asset is
        missing — workspaces must be self-contained before they are shipped.
        Writes metadata/membrane_assets.json listing every staged item.
        """
        # Detect which bilayer file(s) are referenced by step params.
        bilayer_files: set[str] = set()
        for step in result.execution_order:
            bf = step.params.get("bilayer_file")
            if bf:
                bilayer_files.add(bf)
        if not bilayer_files:
            bilayer_files = {"dppc512_whole.gro"}

        staged: list[dict] = []

        # ── Bilayer GRO files ─────────────────────────────────────────────
        for bilayer_file in sorted(bilayer_files):
            src = _PROT_MEMB_FILES / bilayer_file
            if not src.exists():
                raise FileNotFoundError(
                    f"Membrane asset not found: {src}\n"
                    f"Expected under {_PROT_MEMB_FILES}"
                )
            dst = assets_dir / bilayer_file
            shutil.copy2(src, dst)
            staged.append({
                "asset":       bilayer_file,
                "type":        "file",
                "source":      str(src),
                "destination": str(dst),
                "size_bytes":  src.stat().st_size,
            })

        # ── inflategro-Jorge.pl ───────────────────────────────────────────
        inflategro_src = _PROT_MEMB_FILES / "inflategro-Jorge.pl"
        if not inflategro_src.exists():
            raise FileNotFoundError(f"Membrane asset not found: {inflategro_src}")
        inflategro_dst = assets_dir / "inflategro-Jorge.pl"
        shutil.copy2(inflategro_src, inflategro_dst)
        staged.append({
            "asset":       "inflategro-Jorge.pl",
            "type":        "file",
            "source":      str(inflategro_src),
            "destination": str(inflategro_dst),
            "size_bytes":  inflategro_src.stat().st_size,
        })

        # ── oplsaa_membrane.ff/ ───────────────────────────────────────────
        ff_src = _PROT_MEMB_FILES / "oplsaa_membrane.ff"
        if not ff_src.is_dir():
            raise FileNotFoundError(f"Membrane asset not found: {ff_src}")
        ff_dst = assets_dir / "oplsaa_membrane.ff"
        if ff_dst.exists():
            shutil.rmtree(ff_dst)
        shutil.copytree(ff_src, ff_dst)
        n_ff_files = sum(1 for _ in ff_dst.rglob("*") if _.is_file())
        staged.append({
            "asset":       "oplsaa_membrane.ff",
            "type":        "directory",
            "source":      str(ff_src),
            "destination": str(ff_dst),
            "n_files":     n_ff_files,
        })

        # ── Write metadata ────────────────────────────────────────────────
        meta = {
            "staged_at":  datetime.now().isoformat(timespec="seconds"),
            "source_dir": str(_PROT_MEMB_FILES),
            "assets_dir": str(assets_dir),
            "assets":     staged,
        }
        (metadata_dir / "membrane_assets.json").write_text(
            json.dumps(meta, indent=4)
        )

    def _stage_inputs(
        self,
        result: "CompilationResult",
        inputs_dir: Path,
    ) -> None:
        """
        Copy all external source files declared in component definitions
        into workspace/inputs/.

        Raises FileNotFoundError at compile time (not runtime) if any
        declared source file is missing — workspaces must be self-contained.

        Files are stored as inputs/{component_id}{original_extension}.
        """
        for comp in result.state.components:
            if not comp.file:
                continue

            src = Path(comp.file)

            if not src.exists():
                raise FileNotFoundError(
                    f"Source file for component '{comp.id}' not found: {src}\n"
                    f"  Check the 'file:' field in your YAML config.\n"
                    f"  Path was resolved relative to the YAML directory."
                )

            ext = src.suffix or ".pdb"
            dst = inputs_dir / f"{comp.id}{ext}"

            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)