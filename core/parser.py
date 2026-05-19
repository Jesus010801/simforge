# core/parser.py

from pathlib import Path
import yaml

from core.models import SystemState
from core.inference import run_inference


def parse_yaml(path: str | Path) -> SystemState:
    """
    Lee un YAML y retorna un SystemState completamente inferido.

    Orden garantizado:
      1. Cargar YAML
      2. Validar estructura con Pydantic
      3. Ejecutar pipeline de inferencia biológica
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Archivo no encontrado: {path}")

    if path.suffix not in (".yaml", ".yml"):
        raise ValueError(f"El archivo debe ser YAML: {path}")

    with open(path, "r") as f:
        raw = yaml.safe_load(f)

    try:
        state = SystemState(**raw)
    except Exception as e:
        raise ValueError(f"Error validando YAML:\n{e}")

    state = run_inference(state)

    return state
