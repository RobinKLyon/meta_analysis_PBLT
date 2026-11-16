"""Generate mixed-model figures from existing mixed-model result files."""

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1]))
from analyse_regression_mixte_intermediate import (
    CONTRASTS,
    plot_regional_delta_forest,
    plot_regional_deltas,
    plot_site_slopes,
)
from src.analysis_config import find_latest_output_folder


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-output", type=Path, required=False, default= find_latest_output_folder(str(Path(__file__).parents[1]) + "/outputs", "baseline_raw"))
    parser.add_argument(
        "--contrasts",
        nargs="+",
        choices=sorted(CONTRASTS),
        default=sorted(CONTRASTS),
    )
    args = parser.parse_args()

    result_root = args.input_output / "data" / "mixed_model_regression"
    if not result_root.is_dir():
        raise FileNotFoundError(f"Missing mixed-model results: {result_root}")
    figure_root = (
        args.input_output
        / "fig"
        / "summary"
        / "mixed_model_regression"
        / "regional_model"
    )
    for contrast_name in args.contrasts:
        data_dir = result_root / contrast_name
        model_results_path = data_dir / "model_results_by_band.csv"
        site_results_path = data_dir / "site_slopes_by_band.csv"
        if not model_results_path.is_file() or not site_results_path.is_file():
            raise FileNotFoundError(
                f"Missing LMM result tables for {contrast_name}: {data_dir}"
            )
        model_results = pd.read_csv(model_results_path)
        site_results = pd.read_csv(site_results_path)
        raw_dir = figure_root / "raw"
        scaled_dir = figure_root / "scaled"
        site_dir = figure_root / "site_slopes"
        plot_regional_delta_forest(
            model_results,
            contrast_name,
            raw_dir / f"{contrast_name}_regional_delta_forest.png",
        )
        plot_regional_deltas(
            model_results,
            contrast_name,
            raw_dir / f"{contrast_name}_regional_deltas.png",
        )
        plot_site_slopes(
            site_results,
            model_results,
            contrast_name,
            site_dir / f"{contrast_name}_site_slopes.png",
        )
        plot_regional_delta_forest(
            model_results,
            contrast_name,
            scaled_dir / f"{contrast_name}_regional_delta_forest.png",
            scaled=True,
        )
        plot_regional_deltas(
            model_results,
            contrast_name,
            scaled_dir / f"{contrast_name}_regional_deltas.png",
            scaled=True,
        )
    print(f"Mixed-model figures generated in: {figure_root.resolve()}")


if __name__ == "__main__":
    main()
