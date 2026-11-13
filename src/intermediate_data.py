"""Persistence helpers for integrated PBLT trial data."""

from pathlib import Path
import pandas as pd

from analysis_config import BANDS


OUTCOME_LABELS = {
    0: "gain_failure_0_eur",
    1: "loss_failure_minus1_eur",
    2: "gain_success_plus1_eur",
    3: "loss_success_0_eur",
}
REQUIRED_COLUMNS = {
    "patient",
    "position",
    "trial_index",
    "contact_index",
    "contact",
    "outcome",
    *BANDS,
}


def save_intermediate_data(integrated_data, output_root):
    """Save integrated trial-level data in one CSV per outcome and position."""
    if not integrated_data:
        return
    data = pd.concat(integrated_data, ignore_index=True)
    intermediate_root = Path(output_root) / "data" / "intermediate_data"
    intermediate_root.mkdir(parents=True, exist_ok=True)
    metadata_columns = [
        "patient", "position", "trial_index", "contact_index", "contact", "outcome"
    ]
    data_columns = metadata_columns + list(BANDS)
    missing_columns = set(data_columns) - set(data.columns)
    if missing_columns:
        raise ValueError(
            "Missing columns in integrated intermediate data: "
            + ", ".join(sorted(missing_columns))
        )
    positions = ["ventral", "dorsal"]
    if "outside" in set(data["position"].astype(str).str.lower()):
        positions.append("outside")
    for position in positions:
        position_data = data[data["position"] == position]
        position_dir = intermediate_root / position
        position_dir.mkdir(parents=True, exist_ok=True)
        for outcome, outcome_name in OUTCOME_LABELS.items():
            output_path = position_dir / f"outcome_{outcome}_{outcome_name}.csv"
            position_data[position_data["outcome"] == outcome][data_columns].to_csv(
                output_path, index=False
            )


def load_intermediate_data(output_root):
    """Load integrated trial-level data saved under an existing output root."""
    intermediate_root = Path(output_root) / "data" / "intermediate_data"
    if not intermediate_root.is_dir():
        intermediate_root = Path(output_root) / "intermediate_data"
    if not intermediate_root.is_dir():
        raise FileNotFoundError(
            f"Intermediate data directory not found in {output_root}"
        )
    data_frames = []
    positions = ["ventral", "dorsal"]
    if (intermediate_root / "outside").is_dir():
        positions.append("outside")
    for position in positions:
        position_dir = intermediate_root / position
        if not position_dir.is_dir():
            raise FileNotFoundError(f"Intermediate data directory not found: {position_dir}")
        files = sorted(position_dir.glob("outcome_*.csv"))
        if not files:
            raise FileNotFoundError(f"No intermediate outcome files found in {position_dir}")
        data_frames.extend(pd.read_csv(path) for path in files)
    data = pd.concat(data_frames, ignore_index=True)
    missing_columns = REQUIRED_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(
            "Missing columns in intermediate data: "
            + ", ".join(sorted(missing_columns))
        )
    contact_columns = ["patient", "position", "contact_index", "contact"]
    return [
        contact_data.reset_index(drop=True)
        for _, contact_data in data.groupby(contact_columns, sort=False)
    ]
