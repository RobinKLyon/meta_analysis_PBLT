"""Run the gain/loss response plotting script."""

from src.plot_reponses_gain_loss import load_summary, plot_summary
from src.analysis_config import (
    BANDS,
    OUTPUT_TYPE_CHOICES,
    find_latest_output_folder,
    load_config,
)
import argparse
from pathlib import Path

def main() -> None:
    script_dir = Path(__file__).resolve().parent
    default_outputs_root = script_dir / "outputs"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Summary CSV. By default, use the newest matching output folder.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where figures are saved. By default, save them in the selected output folder.",
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
        help="Analysis config.json. Defaults to config.json next to the summary when present.",
    )
    args = parser.parse_args()

    output_folder = (
        args.summary.parent
        if args.summary is not None
        else find_latest_output_folder(args.outputs_root, args.output_type)
    )
    summary_path = args.summary or output_folder / "resume_reponses_par_bande_region.csv"
    output_dir = args.output_dir or output_folder / "figures_reponses"
    summary = load_summary(summary_path)
    config_path = args.config or output_folder / "config.json"
    config = load_config(config_path) if config_path.exists() else None
    bands = list(config["bands"]) if config else list(BANDS)
    plot_summary(
        summary,
        "percent",
        output_dir / "reponses_gain_loss_ratios.png",
        bands=bands,
    )
    plot_summary(
        summary,
        "count",
        output_dir / "reponses_gain_loss_comptes.png",
        bands=bands,
    )
    if config and output_dir != config_path.parent:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "config.json").write_text(
            config_path.read_text(encoding="utf-8"), encoding="utf-8"
        )
    print(f"Selected output folder: {output_folder.resolve()}")
    print(f"Figures created in: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
