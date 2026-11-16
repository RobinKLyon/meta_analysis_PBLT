"""Run the configured per-band mixed-model regression analysis."""

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from src.configuration import load_yaml, write_effective_config
from src.analysis_config import find_latest_output_folder

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("lmm.yaml"))
    parser.add_argument("--input-output", type=Path, required=False, default= find_latest_output_folder(str(Path(__file__).parents[1]) + "/outputs", "baseline_raw"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_yaml(args.config)
    model = config.get("model", {})

    command = [
        sys.executable,
        str(Path(__file__).parents[1] / "src" / "analyse_regression_mixte_intermediate.py"),
        "--output", str(args.input_output),
        "--model", str(model.get("preset", "regional_patient")),
        "--p-adjust", str(model.get("p_adjust", "fdr_bh")),
    ]
    subprocess.run(command, check=True)
    write_effective_config(
        args.input_output,
        config,
        script=Path(__file__).name,
        inputs={"input_output": str(args.input_output.resolve())},
        outputs=["data/mixed_model_regression"],
    )


if __name__ == "__main__":
    main()
