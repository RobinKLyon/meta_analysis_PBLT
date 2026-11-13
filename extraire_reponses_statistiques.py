#!/usr/bin/env python3
"""Extract response counts and ratios from the PBLT statistics CSV.

The ratio denominator is the number of distinct sites evaluated in each
position x band combination. Non-significant responses are retained in the
output to make the denominator explicit.

Exemple :
    python extraire_reponses_statistiques.py
    python extraire_reponses_statistiques.py --input path/to/results.csv
"""

from pathlib import Path
import argparse
import sys

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from analysis_config import (
    ALPHA_DEFAULT,
    BANDS,
    OUTPUT_TYPE_CHOICES,
    find_latest_output_folder,
    load_config,
)


REQUIRED_COLUMNS = {
    "patient",
    "position",
    "contact_index",
    "contact",
    "band",
    "response",
}
POSITION_ORDER = ["dorsal", "ventral"]
RESPONSE_ORDER = ["positive", "negative", "non_significant"]


def load_statistics(input_path: Path) -> pd.DataFrame:
    """Load the CSV and validate the columns required for extraction."""
    data = pd.read_csv(input_path)
    missing = REQUIRED_COLUMNS - set(data.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"Missing columns in {input_path}: {missing_text}")

    data = data.copy()
    data["position"] = data["position"].str.lower()
    unknown_positions = sorted(set(data["position"].dropna()) - set(POSITION_ORDER))
    if unknown_positions:
        raise ValueError(f"Unknown positions in {input_path}: {unknown_positions}")

    unknown_responses = sorted(set(data["response"].dropna()) - set(RESPONSE_ORDER))
    if unknown_responses:
        raise ValueError(f"Unknown responses in {input_path}: {unknown_responses}")
    return data


def extract_response_summary(data: pd.DataFrame, band_order=None) -> pd.DataFrame:
    """Aggregate responses by band and position at site level."""
    site_columns = ["patient", "position", "contact_index", "contact", "band"]
    sites = data.drop_duplicates(site_columns).copy()

    counts = (
        sites.groupby(["band", "position", "response"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=RESPONSE_ORDER, fill_value=0)
        .reset_index()
    )
    totals = (
        sites.groupby(["band", "position"], observed=False)
        .size()
        .rename("n_sites")
        .reset_index()
    )
    summary = totals.merge(counts, on=["band", "position"], how="left")
    summary = summary.rename(
        columns={
            "positive": "n_positive",
            "negative": "n_negative",
            "non_significant": "n_non_significant",
        }
    )
    for column in ["n_positive", "n_negative", "n_non_significant"]:
        summary[column] = summary[column].fillna(0).astype(int)

    summary["positive_ratio"] = summary["n_positive"] / summary["n_sites"]
    summary["negative_ratio"] = summary["n_negative"] / summary["n_sites"]
    summary["non_significant_ratio"] = (
        summary["n_non_significant"] / summary["n_sites"]
    )
    summary["positive_percent"] = 100 * summary["positive_ratio"]
    summary["negative_percent"] = 100 * summary["negative_ratio"]
    summary["non_significant_percent"] = 100 * summary["non_significant_ratio"]

    band_order = band_order or list(BANDS)
    band_order_map = {band: index for index, band in enumerate(band_order)}
    summary["_band_order"] = summary["band"].map(band_order_map).fillna(len(band_order_map))
    summary["_position_order"] = summary["position"].map(
        {position: index for index, position in enumerate(POSITION_ORDER)}
    )
    summary = (
        summary.sort_values(["_band_order", "_position_order"])
        .drop(columns=["_band_order", "_position_order"])
        .reset_index(drop=True)
    )
    return summary


def add_fisher_region_tests(summary: pd.DataFrame, alpha=0.05) -> pd.DataFrame:
    """Test whether significant-site status is associated with region."""
    summary = summary.copy()
    p_values = []
    metadata = []
    for band in summary["band"].drop_duplicates():
        band_data = summary[summary["band"] == band].set_index("position")
        table = []
        for position in POSITION_ORDER:
            row = band_data.loc[position] if position in band_data.index else None
            significant = (
                int(row["n_positive"] + row["n_negative"])
                if row is not None else 0
            )
            total = int(row["n_sites"]) if row is not None else 0
            table.append([significant, max(total - significant, 0)])
        _, p_value = fisher_exact(table, alternative="two-sided")
        p_values.append(float(p_value))
        metadata.append(band)

    tests = {
        metadata[index]: p_values[index]
        for index in range(len(metadata))
    }
    summary["fisher_p_value"] = summary["band"].map(lambda band: tests[band])
    summary["fisher_significant"] = summary["fisher_p_value"] < alpha
    return summary


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    default_outputs_root = script_dir / "outputs"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Source statistics CSV. By default, use the newest matching output folder.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Aggregated CSV to create. By default, save it in the selected output folder.",
    )
    parser.add_argument(
        "--output-type",
        choices=sorted(OUTPUT_TYPE_CHOICES),
        default="baseline_div",
        help="Baseline output type used to find the newest output folder.",
    )
    parser.add_argument(
        "--outputs-root",
        type=Path,
        default=default_outputs_root,
        help="Root directory containing typed analysis output folders.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Analysis config.json. Defaults to config.json next to the input CSV when present.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help="FDR significance threshold. Defaults to the saved config or 0.05.",
    )
    args = parser.parse_args()

    output_folder = (
        args.input.parent
        if args.input is not None
        else find_latest_output_folder(args.outputs_root, args.output_type)
    )
    input_path = args.input or output_folder / "resultats_statistiques_tous_patients.csv"
    output_path = args.output or output_folder / "resume_reponses_par_bande_region.csv"
    data = load_statistics(input_path)
    config_path = args.config or output_folder / "config.json"
    config = load_config(config_path) if config_path.exists() else None
    band_order = list(config["bands"]) if config else list(BANDS)
    configured_alpha = (
        config.get("parameters", {}).get("alpha", ALPHA_DEFAULT)
        if config else ALPHA_DEFAULT
    )
    alpha = configured_alpha if args.alpha is None else args.alpha
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be between 0 and 1, got {alpha}")
    summary = add_fisher_region_tests(
        extract_response_summary(data, band_order=band_order), alpha=alpha
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    if config and output_path.parent != config_path.parent:
        (output_path.parent / "config.json").write_text(
            config_path.read_text(encoding="utf-8"), encoding="utf-8"
        )
    print(f"Selected output folder: {output_folder.resolve()}")
    print(f"Summary created: {output_path.resolve()}")
    print(f"Extracted rows: {len(summary)}")


if __name__ == "__main__":
    main()
