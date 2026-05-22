# builders/step_builders/_utils.py

from __future__ import annotations

import os
from pathlib import Path


def rel(from_dir: Path, to_dir: Path) -> str:
    """Relative path from from_dir to to_dir for use in shell scripts."""
    return os.path.relpath(to_dir, from_dir)
