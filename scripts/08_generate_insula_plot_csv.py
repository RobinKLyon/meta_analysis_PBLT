"""Generate MNI/Insula Plotter CSV files for a selected frequency band."""

import argparse
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).parents[1]))
from src.analysis_config import find_latest_output_folder
print(sys.path)
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-output", type=Path, required=False, default= find_latest_output_folder(str(Path(__file__).parents[1]) + "/outputs", "baseline_raw"))
    parser.add_argument("--band", required=True)
    parser.add_argument("--sign", choices=("positive", "negative"), default="positive")
    args = parser.parse_args()
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parents[1] / "src" / "plot_sites_mni_3d.py"),
            "--output", str(args.input_output),
            "--bands", args.band,
            "--sign", args.sign,
        ],
        check=True,
    )
    path = (
        args.input_output / "data" / "gain_loss_responses" / "mni_site_csv"
        / f"mni_sites_{args.sign}_{args.band}_with_t_value.csv"
    )
    print(f"Insula Plotter CSV: {path.resolve()}")


if __name__ == "__main__":
    main()
