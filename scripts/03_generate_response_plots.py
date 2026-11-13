"""Generate response barplots from saved comparison CSV files."""

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from src.configuration import load_yaml, write_effective_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/default.yaml"))
    parser.add_argument("--input-output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_yaml(args.config)
    command = [
        sys.executable,
        str(Path(__file__).parents[1] / "generate_gain_loss_figures.py"),
        "--output", str(args.input_output),
    ]
    if config.get("intermediate", {}).get("include_outside", False):
        command.extend(["--plot-outside", "all"])
    else:
        command.extend(["--plot-outside", "exclude"])
    subprocess.run(command, check=True)
    write_effective_config(
        args.input_output,
        config,
        script=Path(__file__).name,
        inputs={"input_output": str(args.input_output.resolve())},
        outputs=["fig/summary", "data/gain_loss_responses"],
    )


if __name__ == "__main__":
    main()
