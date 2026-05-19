# core/parser.py

from pathlib import Path
import yaml

from core.models import SystemState
from core.inference import run_inference


def parse_yaml(path: str | Path) -> SystemState:
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Archivo no encontrado: {path}")

    if path.suffix not in (".yaml", ".yml"):
        raise ValueError(f"El archivo debe ser YAML: {path}")

    with open(path, "r") as f:
        raw = yaml.safe_load(f)

    # Resolver paths de componentes relativos al YAML
    base_dir = path.parent
    if "components" in raw:
        for component in raw["components"]:
            if "file" in component:
                component_path = Path(component["file"])
                if not component_path.is_absolute():
                    component["file"] = str(base_dir / component_path)

    try:
        state = SystemState(**raw)
    except Exception as e:
        raise ValueError(f"Error validando YAML:\n{e}")

    state = run_inference(state)
    return state