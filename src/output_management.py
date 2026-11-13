"""Shared output-folder and configuration helpers for PBLT analyses."""

from pathlib import Path
import json


def prepare_output_root(output_root, config=None):
    """Create an output root and reject incompatible existing configurations."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    config_path = output_root / "config.json"
    if config_path.exists() and config is not None:
        existing = json.loads(config_path.read_text(encoding="utf-8"))
        existing_baseline = existing.get("parameters", {}).get("baseline_type")
        requested_baseline = config.get("parameters", {}).get("baseline_type")
        if (
            existing_baseline is not None
            and requested_baseline is not None
            and existing_baseline != requested_baseline
        ):
            raise ValueError(
                f"Output folder {output_root} already uses baseline "
                f"'{existing_baseline}', not '{requested_baseline}'. "
                "Choose a matching output folder or a new folder."
            )
    return output_root


def require_existing_output_root(output_root):
    """Validate an output root used for loading previously saved data."""
    output_root = Path(output_root)
    if not output_root.is_dir():
        raise FileNotFoundError(f"Output directory not found: {output_root}")
    config_path = output_root / "config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    return output_root
