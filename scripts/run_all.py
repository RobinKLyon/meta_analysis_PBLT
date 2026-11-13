"""Run the PBLT analysis stages in dependency order."""

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from src.analysis_config import find_latest_output_folder
from src.configuration import load_yaml


def artifact_is_complete(step: str, output: Path, artifacts: dict[str, Path]) -> bool:
    """Return whether a stage has produced the files required by later stages."""
    if step == "lmm":
        contrasts = (
            "gain_plus1_vs_zero",
            "loss_zero_vs_minus1",
            "minus1_vs_plus1",
        )
        return all(
            (artifacts[step] / contrast / filename).is_file()
            for contrast in contrasts
            for filename in ("model_results_by_band.csv", "site_slopes_by_band.csv")
        )
    return artifacts.get(step, Path()).is_dir()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/run_all.yaml"))
    parser.add_argument("--input-output", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_yaml(args.config)
    root = Path(__file__).parents[1]
    scripts = Path(__file__).parent
    output = args.input_output.resolve() if args.input_output else None
    steps = config.get("steps", {})

    if output is None and steps.get("intermediate", True):
        command = [
            sys.executable,
            str(scripts / "01_generate_intermediate.py"),
            "--config",
            str(root / "configs" / "default.yaml"),
        ]
        if args.overwrite:
            command.append("--overwrite")
        subprocess.run(command, cwd=root, check=True)
        output = find_latest_output_folder(root / "outputs", "baseline_raw")
    if output is None:
        output = find_latest_output_folder(root / "outputs", "baseline_raw")

    stage_commands = [
        ("comparisons", "02_run_comparisons.py", "comparisons.yaml"),
        ("response_plots", "03_generate_response_plots.py", "default.yaml"),
        ("lmm", "04_run_lmm.py", "lmm.yaml"),
        ("lmm_figures", "05_generate_lmm_figures.py", None),
        ("temporal_response_categories", "06_temporal_response_categories.py", None),
        ("tfr_response_categories", "07_plot_tfr_response_categories.py", None),
        ("insula_csv", "08_generate_insula_plot_csv.py", None),
    ]
    existing_artifacts = {
        "comparisons": output / "data" / "gain_loss_responses",
        "response_plots": output / "fig" / "summary",
        "lmm": output / "data" / "mixed_model_regression",
        "lmm_figures": output / "fig" / "summary" / "mixed_model_regression",
        "temporal_response_categories": output / "data" / "temporal_response_categories",
        "tfr_response_categories": output / "data" / "tfr_power_by_response_category",
        "insula_csv": output / "data" / "gain_loss_responses" / "mni_site_csv",
    }
    for step, filename, config_name in stage_commands:
        if not steps.get(step, True):
            continue
        if (
            not args.overwrite
            and config.get("reuse_existing", True)
            and artifact_is_complete(step, output, existing_artifacts)
        ):
            print(f"Skipping existing step '{step}': {existing_artifacts[step]}")
            continue
        command = [sys.executable, str(scripts / filename)]
        if config_name:
            command.extend(["--config", str(root / "configs" / config_name)])
        if filename not in {"05_generate_lmm_figures.py", "06_temporal_response_categories.py",
                            "07_plot_tfr_response_categories.py"}:
            command.extend(["--input-output", str(output)])
        else:
            command.extend(["--input-output", str(output)])
        if filename == "08_generate_insula_plot_csv.py":
            command.extend(["--band", str(config.get("insula_band", "HFA"))])
        if args.overwrite and filename in {"02_run_comparisons.py", "03_generate_response_plots.py",
                                            "04_run_lmm.py"}:
            command.append("--overwrite")
        subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    main()
