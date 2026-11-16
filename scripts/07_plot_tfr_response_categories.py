"""Generate TFR/t-value plots grouped by response category."""

import argparse
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from src.analysis_config import find_latest_output_folder
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-output", type=Path, required=False, default= find_latest_output_folder(str(Path(__file__).parents[1]) + "/outputs", "baseline_raw"))
    parser.add_argument("--sign", choices=("positive", "negative", "both"), default="both")
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--regions", nargs="+", default=None)
    args = parser.parse_args()
    command = [
            sys.executable,
            str(Path(__file__).parents[1] / "src" / "plot_tfr_power_by_response_category.py"),
            "--output", str(args.input_output), "--sign", args.sign,
        ]
    if args.bands:
        command.extend(["--bands", *args.bands])
    if args.regions:
        command.extend(["--regions", *args.regions])
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
