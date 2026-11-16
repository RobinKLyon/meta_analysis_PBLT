"""Shared YAML configuration and effective-run metadata helpers."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
from typing import Any, Optional

import yaml

def load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML mapping and reject empty or malformed configurations."""

    path = Path(path)
    #check if path is absolute or relative to configs:
    if not path.is_absolute():
        path = Path(__file__).parent.parent / "configs" / path
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration must contain a mapping: {path}")
    return config


def resolve_path(value: Optional[str], base_dir: Path) -> Optional[Path]:
    """Resolve an absolute path or a path relative to the YAML file."""
    if value in (None, ""):
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def write_effective_config(
    output_root: Path,
    config: dict[str, Any],
    *,
    script: str,
    inputs: Optional[dict[str, Any]] = None,
    outputs: Optional[list[str]] = None,
) -> Path:
    """Write the effective configuration and append the current pipeline step."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    path = output_root / "config.json"
    if path.is_file():
        with path.open("r", encoding="utf-8") as handle:
            effective = json.load(handle)
    else:
        effective = {}
    effective["configuration"] = config
    effective.setdefault("pipeline_steps", []).append(
        {
            "script": script,
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "inputs": inputs or {},
            "outputs": outputs or [],
        }
    )
    path.write_text(
        json.dumps(effective, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return path
