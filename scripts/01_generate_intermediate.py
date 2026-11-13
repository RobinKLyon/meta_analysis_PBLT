"""Generate reusable intermediate data from epoch-level TFR files."""

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1]))
from src.analysis_config import find_latest_output_folder
from src.configuration import load_yaml, resolve_path, write_effective_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/default.yaml"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_yaml(args.config)
    dataset = config.get("dataset", {})
    intermediate = config.get("intermediate", {})
    output_root = resolve_path(intermediate.get("output_root"), args.config.parent)
    if output_root is None:
        raise ValueError("intermediate.output_root is required.")
    if output_root.exists() and any(output_root.iterdir()) and not (
        args.overwrite or config.get("output", {}).get("overwrite", False)
    ):
        raise FileExistsError(
            f"Output directory is not empty: {output_root}. Use --overwrite."
        )
    command = [
        sys.executable,
        str(Path(__file__).parents[1] / "analyse_pblt_tous_patients.py"),
        "--intermediate-only",
        "--baseline_type",
        str(intermediate.get("baseline_type", "raw")),
        "--tmin",
        str(intermediate.get("time_min", 0.0)),
        "--tmax",
        str(intermediate.get("time_max", 0.8)),
    ]
    if output_root.name.lower() != "outputs":
        command.extend(["--output", str(output_root)])
    if dataset.get("tfr_root"):
        command.extend(["--data-root", str(resolve_path(dataset["tfr_root"], args.config.parent))])
    if dataset.get("roi_file"):
        command.extend(["--roi-file", str(resolve_path(dataset["roi_file"], args.config.parent))])
    patients = intermediate.get("patients")
    if patients:
        command.extend(["--patients", *map(str, patients)])
    if intermediate.get("include_outside", False):
        command.append("--take_outside")
    subprocess.run(command, check=True)
    if output_root.name.lower() == "outputs":
        output_root = find_latest_output_folder(output_root, f"baseline_{intermediate.get('baseline_type', 'raw')}")
    write_effective_config(
        output_root,
        config,
        script=Path(__file__).name,
        inputs={"config": str(args.config.resolve())},
        outputs=["data/intermediate_data", "config.json"],
    )


if __name__ == "__main__":
    main()
