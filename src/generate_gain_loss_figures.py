"""Generate gain/loss category figures from saved CSV results."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analysis_config import BANDS
from plot_data import plot_bar_plot_compare_gain_loss_responses_multiple_bands


CATEGORIES = [
    "gain_only",
    "both_prefer_gain",
    "both_equal",
    "both_prefer_loss",
    "loss_only",
    "neither",
]
COMPARISONS = {
    "plus1_vs_minus1": "overall_response",
    "plus1_vs_zero": "gain_response",
    "minus1_vs_zero": "loss_response",
}


def clean_generated_outputs(output: Path) -> None:
    data_dir = output / "data" / "gain_loss_responses"
    for path in data_dir.glob("electrode_responses_*.csv"):
        path.unlink()
    for path in data_dir.glob("outside_significant_*.csv"):
        path.unlink()
    filtered_dir = data_dir / "filtered_electrode_lists"
    if filtered_dir.is_dir():
        shutil.rmtree(filtered_dir)

    summary_dir = (
        output / "fig" / "summary" / "plus1_vs_zero_and_zero_vs_minus1"
    )
    if not summary_dir.exists():
        return
    for path in summary_dir.rglob("*.png"):
        if (
            "gain_loss_site_responses" in path.name
            or "merged_gain_loss_site_responses" in path.name
        ):
            path.unlink()
    for path in summary_dir.glob("bands_*"):
        if path.is_dir():
            shutil.rmtree(path)


def ensure_plot_source(output: Path, bands: list[str]) -> Path:
    """Create the statistical CSV when it is absent but intermediates exist."""
    data_dir = output / "data" / "gain_loss_responses"
    source = data_dir / "pairwise_comparison_results_by_band.csv"
    if source.is_file():
        return source

    config_path = output / "config.json"
    intermediate_dir = output / "data" / "intermediate_data"
    if not config_path.is_file() or not intermediate_dir.is_dir():
        raise FileNotFoundError(
            f"Missing plot data CSV: {source}. "
            "Run analyse_pblt_tous_patients.py first to create intermediate data."
        )
    config = json.loads(config_path.read_text(encoding="utf-8"))
    parameters = config.get("parameters", {})
    analysis_script = Path(__file__).with_name("analyse_pblt_tous_patients.py")
    command = [
        sys.executable,
        str(analysis_script),
        "--load-intermediate-data",
        "--output",
        str(output),
        "--calculation",
        "plus1_vs_zero_and_zero_vs_minus1",
        "--baseline_type",
        str(parameters.get("baseline_type", "div")),
        "--tmin",
        str(parameters.get("tmin", 0.0)),
        "--tmax",
        str(parameters.get("tmax", 0.8)),
        "--bands",
        *bands,
        "--no-plots",
    ]
    if parameters.get("take_outside", False):
        command.append("--take_outside")
    subprocess.run(command, check=True)
    if not source.is_file():
        raise FileNotFoundError(
            f"Statistical generation completed but did not create {source}"
        )
    return source


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Analysis output folder containing data/gain_loss_responses.",
    )
    parser.add_argument(
        "--bands",
        nargs="+",
        choices=list(BANDS),
        default=list(BANDS),
        help="Frequency bands to plot.",
    )
    parser.add_argument(
        "--plot-outside",
        choices=("all", "exclude", "only"),
        default="all",
        help="Plot all regions, exclude outside, or plot outside only.",
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        choices=("dorsal", "ventral", "outside"),
        default=("dorsal", "ventral", "outside"),
        help="Regions to include in CSV exports and figures.",
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        choices=CATEGORIES,
        default=CATEGORIES,
        help="Merged response categories to export and display.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove previously generated electrode CSVs and category figures first.",
    )
    args = parser.parse_args()

    data_dir = args.output / "data" / "gain_loss_responses"
    if args.clean:
        clean_generated_outputs(args.output)
    source = ensure_plot_source(args.output, args.bands)

    data = pd.read_csv(source)
    required = {
        "band",
        "position",
        "positive_merged_response_category",
        "negative_merged_response_category",
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {source}: {sorted(missing)}")

    data = data[
        data["band"].isin(args.bands) & data["position"].isin(args.regions)
    ].copy()
    if args.plot_outside == "exclude":
        data = data[data["position"] != "outside"].copy()
    elif args.plot_outside == "only":
        data = data[data["position"] == "outside"].copy()
    if data.empty:
        raise ValueError("No rows remain after applying band and region filters.")
    effective_regions = [
        region for region in args.regions if region in set(data["position"])
    ]

    export_dir = data_dir / "filtered_electrode_lists"
    export_dir.mkdir(parents=True, exist_ok=True)
    for comparison, response_column in COMPARISONS.items():
        comparison_data = data.copy()
        for direction in ("positive", "negative"):
            category_column = f"{direction}_merged_response_category"
            direction_data = comparison_data[
                comparison_data[response_column].eq(direction)
                & comparison_data[category_column].isin(args.categories)
            ].copy()
            direction_data.insert(5, "comparison", comparison)
            direction_data.to_csv(
                export_dir
                / f"{comparison}_{direction}_bands_{'_'.join(args.bands)}_"
                f"regions_{'_'.join(effective_regions)}.csv",
                index=False,
            )

    figure_dir = (
        args.output
        / "fig"
        / "summary"
        / "plus1_vs_zero_and_zero_vs_minus1"
        / f"bands_{'_'.join(args.bands)}"
        / f"regions_{'_'.join(effective_regions)}"
        / f"categories_{'_'.join(args.categories)}"
    )
    for sign in ("positive", "negative"):
        plot_data = data.rename(
            columns={f"{sign}_merged_response_category": "response_category"}
        )
        plot_data = plot_data[plot_data["response_category"].isin(args.categories)]
        plot_bar_plot_compare_gain_loss_responses_multiple_bands(
            plot_data,
            figure_dir / sign,
            BANDS=BANDS,
            metric="percent",
            categories=CATEGORIES,
            bands_for_final_bar_plot=args.bands,
            output_name=f"merged_gain_loss_site_responses_{sign}_by_region",
            response_sign=sign,
            denominator_results=data,
        )

    print(f"Plot data source: {source.resolve()}")
    print(f"Filtered electrode CSVs: {export_dir.resolve()}")
    print(f"Figures created in: {figure_dir.resolve()}")


if __name__ == "__main__":
    main()
