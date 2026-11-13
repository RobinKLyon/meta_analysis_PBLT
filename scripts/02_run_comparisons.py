"""Run independent site-level event comparisons on intermediate data."""

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from src.analysis_config import find_latest_output_folder
from src.configuration import load_yaml, resolve_path, write_effective_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/comparisons.yaml"))
    parser.add_argument("--input-output", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_yaml(args.config)
    input_output = args.input_output
    if input_output is None:
        root = resolve_path("./outputs", args.config.parent)
        input_output = find_latest_output_folder(
            root, config.get("input", {}).get("latest_output_type", "baseline_raw")
        )
    input_output = input_output.resolve()
    alpha = str(config.get("statistics", {}).get("alpha", 0.05))
    script = Path(__file__).parents[1] / "analyse_pblt_tous_patients.py"
    for calculation in ("minus1_vs_plus1", "plus1_vs_zero_and_zero_vs_minus1"):
        command = [
            sys.executable, str(script), "--output", str(input_output),
            "--load-intermediate-data", "--calculation", calculation,
            "--alpha", alpha, "--no-plots",
        ]
        subprocess.run(command, check=True)
    write_effective_config(
        input_output,
        config,
        script=Path(__file__).name,
        inputs={"input_output": str(input_output)},
        outputs=["data/gain_loss_responses", "resultats_statistiques_tous_patients.csv"],
    )


if __name__ == "__main__":
    main()
