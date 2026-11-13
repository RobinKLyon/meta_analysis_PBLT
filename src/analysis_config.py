"""Shared configuration for PBLT statistical analyses."""

from datetime import datetime
import json
from pathlib import Path


BANDS = {
    "Theta": (4, 8),
    "Alpha": (8, 12),
    "Beta": (12, 30),
    "Gamma": (30, 50),
    "HFA": (50, 150),
}
BASELINE_LABELS = {
    "div": "baseline_div",
    "zscore": "baseline_zscore",
    "log_div": "baseline_log_div",
    "raw": "baseline_raw",
    "diff": "baseline_diff",
}
OUTPUT_TYPE_CHOICES = tuple(BASELINE_LABELS.values())
BAND_DISPLAY_NAMES = {"HFA": "BGA"}


def format_band_label(band):
    """Return the publication label and frequency range for a band."""
    display_name = BAND_DISPLAY_NAMES.get(band, band)
    if band in BANDS:
        low, high = BANDS[band]
        return f"{display_name} ({low:g}-{high:g} Hz)"
    return display_name


def build_output_folder_name(baseline_type, date_label):
    """Build an output folder name independent of the outcome contrast."""
    try:
        baseline_label = BASELINE_LABELS[baseline_type]
    except KeyError as exc:
        raise ValueError(f"Unsupported baseline type: {baseline_type}") from exc
    return f"{baseline_label}_{date_label}"
OUTCOMES = {
    0: "Gain failure (0 EUR)",
    1: "Loss failure (-1 EUR)",
    2: "Gain success (+1 EUR)",
    3: "Loss success (0 EUR)",
}
ALPHA_DEFAULT = 0.05
CALCULATION_CHOICES = {
    "minus1_vs_plus1": "-1 vs +1",
    "plus1_vs_zero_and_zero_vs_minus1": "+1 vs 0 and -1 vs 0",
}
BASELINE_MAX = -2.0


def build_config(
    tmin,
    tmax,
    alpha,
    calculation,
    channels_verbose=False,
    baseline_type="div",
    take_outside=False,
    arguments=None,
):
    """Build a JSON-serializable record of the analysis parameters."""
    return {
        "algorithm": "PBLT TFR statistical analysis",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "parameters": {
            "tmin": float(tmin),
            "tmax": float(tmax),
            "alpha": float(alpha),
            "calculation": calculation,
            "calculation_label": CALCULATION_CHOICES[calculation],
            "calculation_choices": CALCULATION_CHOICES,
            "channels_verbose": bool(channels_verbose),
            "baseline_type": baseline_type,
            "take_outside": bool(take_outside),
            "baseline_max": BASELINE_MAX,
        },
        "arguments": arguments or {},
        "bands": {
            name: {"frequency_min_hz": limits[0], "frequency_max_hz": limits[1]}
            for name, limits in BANDS.items()
        },
        "outcomes": {str(code): label for code, label in OUTCOMES.items()},
    }


def save_config(config, output_root):
    """Save the effective analysis configuration under output_root."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    config_path = output_root / "config.json"
    config_path.write_text(
        json.dumps(config, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return config_path


def load_config(config_path):
    """Load a saved configuration file."""
    return json.loads(Path(config_path).read_text(encoding="utf-8"))


def find_latest_output_folder(output_root, output_type):
    """Return the most recently modified output folder for an analysis type."""
    output_root = Path(output_root)
    candidates = [
        path
        for path in output_root.glob(f"{output_type}_*")
        if path.is_dir() and (path / "config.json").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No output folder found for type '{output_type}' in {output_root}"
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)
