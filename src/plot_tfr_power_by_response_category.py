"""Plot frequency-band power over time for gain/loss response categories."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from analysis_config import BANDS
from behavorial_organize_data import get_pos_events, prepare_outcomes
from electrodes_lebelling_masks import (
    bipolar_to_monopolar,
    classify_contacts_insula,
    get_chanlist_from_info,
    monopolar_channels_from_bipolar,
)
from load_tfr_data import get_hdr_info, load_patient, load_tfr_data_wya
from signal_processing import integrate_frequency
from paths_data import path_to_PBLT, path_to_julich_output


CATEGORIES = [
    "loss_only",
    "gain_only",
    "both",
    "both_prefer_loss",
    "both_prefer_gain",
    "neither",
]
CATEGORY_COLORS = {
    "loss_only": "#A80000",
    "gain_only": "#16803C",
    "both": "#A67C52",
    "both_prefer_loss": "#F4A6A6",
    "both_prefer_gain": "#8CCB9B",
    "neither": "#BDBDBD",
}


def load_config(output: Path) -> dict:
    path = output / "config.json"
    if not path.is_file():
        raise FileNotFoundError(f"Missing config.json: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def build_category_lookup(data: pd.DataFrame, sign: str) -> dict:
    column = f"{sign}_merged_response_category"
    required = {"patient", "position", "contact_index", "band", column}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing category columns: {sorted(missing)}")
    return {
        (
            str(row.patient),
            str(row.position),
            int(row.contact_index),
            str(row.band),
        ): str(getattr(row, column))
        for row in data.itertuples(index=False)
    }


def event_t_values(epochs: np.ndarray, outcomes: np.ndarray, event_code: int) -> np.ndarray:
    """Test event power against outcome-zero power at each post-event time."""
    event_epochs = epochs[outcomes == event_code]
    zero_epochs = epochs[outcomes == 0]
    if len(event_epochs) < 2 or len(zero_epochs) < 2:
        return np.full(epochs.shape[1], np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        return ttest_ind(
            event_epochs,
            zero_epochs,
            axis=0,
            equal_var=False,
            nan_policy="omit",
        ).statistic


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bands", nargs="+", choices=list(BANDS), default=list(BANDS))
    parser.add_argument(
        "--regions",
        nargs="+",
        choices=("dorsal", "ventral", "outside"),
        default=("dorsal", "ventral"),
    )
    parser.add_argument("--sign", choices=("positive", "negative", "both"), default="both")
    parser.add_argument("--categories", nargs="+", choices=CATEGORIES, default=CATEGORIES)
    parser.add_argument("--tmin", type=float, default=-1.0)
    parser.add_argument("--tmax", type=float, default=2.0)
    parser.add_argument("--integration-tmin", type=float, default=None)
    parser.add_argument("--integration-tmax", type=float, default=None)
    args = parser.parse_args()

    config = load_config(args.output)
    parameters = config.get("parameters", {})
    tmin = float(args.tmin)
    tmax = float(args.tmax)
    integration_tmin = float(
        parameters.get("tmin", 0.0)
        if args.integration_tmin is None
        else args.integration_tmin
    )
    integration_tmax = float(
        parameters.get("tmax", 0.8)
        if args.integration_tmax is None
        else args.integration_tmax
    )
    root = Path(path_to_PBLT)
    roi_path = Path(path_to_julich_output) / "PBLT_ROIs.csv"
    source = args.output / "data" / "gain_loss_responses" / "pairwise_comparison_results_by_band.csv"
    if not source.is_file():
        raise FileNotFoundError(f"Missing category CSV: {source}")
    category_data = pd.read_csv(source)
    signs = ("positive", "negative") if args.sign == "both" else (args.sign,)
    category_lookups = {
        sign: build_category_lookup(category_data, sign) for sign in signs
    }
    output_dir = args.output / "fig" / "summary" / "tfr_power_by_response_category"
    csv_dir = args.output / "data" / "tfr_power_by_response_category"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)

    traces = []
    patients = sorted(category_data["patient"].astype(str).unique())
    roi_data = pd.read_csv(roi_path, sep=";")
    for patient in patients:
        header, tfr_path = load_patient(patient, root)
        info = get_hdr_info(header)
        data = load_tfr_data_wya(
            tfr_path,
            info["N_trial"],
            info["N_chan"],
            info["N_freqs"],
            info["N_time"],
        )
        times = np.asarray(info["times"])
        time_mask = (times >= tmin) & (times <= tmax)
        if time_mask.sum() < 2:
            raise ValueError(f"Too few time points for {patient} in {tmin}-{tmax}s")
        band_signals = integrate_frequency(data, info["freqs"])
        events = get_pos_events(patient, root)
        _, _, outcomes = prepare_outcomes(events, data.shape[0])
        channels = get_chanlist_from_info(info)
        mono = monopolar_channels_from_bipolar(channels)
        first, _ = bipolar_to_monopolar(channels, mono)
        classification = classify_contacts_insula(
            roi_data, first, patient.removesuffix("_PBLT")
        )
        for contact_index, row in enumerate(classification.itertuples(index=False)):
            position = str(row.classification).lower()
            if position not in args.regions:
                continue
            for band in args.bands:
                band_rows = category_data[
                    (category_data["patient"].astype(str) == patient)
                    & (category_data["position"] == position)
                    & (category_data["contact_index"] == contact_index)
                    & (category_data["band"] == band)
                ]
                if band_rows.empty:
                    continue
                signal = band_signals[band][:, time_mask, contact_index]
                traces.append(
                    {
                        "patient": patient,
                        "position": position,
                        "contact_index": contact_index,
                        "band": band,
                        "times": times[time_mask],
                        "epochs": signal,
                        "outcomes": np.asarray(outcomes),
                        **{
                            sign: category_lookups[sign].get(
                                (patient, position, contact_index, band), "neither"
                            )
                            for sign in signs
                        },
                    }
                )

    outcome_codes = {"+1": 2, "-1": 1}
    for sign in signs:
        for band in args.bands:
            for category in args.categories:
                figure_dir = output_dir / sign / band / category
                data_dir = csv_dir / sign / band / category
                figure_dir.mkdir(parents=True, exist_ok=True)
                data_dir.mkdir(parents=True, exist_ok=True)
                for region in args.regions:
                    figure, axis = plt.subplots(figsize=(10, 5))
                    csv_rows = []
                    plotted = False
                    for event, outcome_code in outcome_codes.items():
                        site_t_values = []
                        site_ids = set()
                        for item in traces:
                            item_category = item[sign]
                            category_matches = (
                                item_category == category
                                or (
                                    category == "both"
                                    and item_category.startswith("both_")
                                )
                            )
                            if (
                                item["position"] != region
                                or item["band"] != band
                                or not category_matches
                            ):
                                continue
                            outcome_mask = item["outcomes"] == outcome_code
                            if np.any(outcome_mask):
                                site_t = event_t_values(
                                    item["epochs"], item["outcomes"], outcome_code
                                )
                                if np.any(np.isfinite(site_t)):
                                    site_t_values.append(site_t)
                                    site_ids.add(
                                        (item["patient"], item["contact_index"])
                                    )
                        if not site_t_values:
                            continue
                        plotted = True
                        site_values = np.vstack(site_t_values)
                        values = np.nanmean(site_values, axis=0)
                        counts = np.sum(np.isfinite(site_values), axis=0)
                        sem = np.zeros_like(values)
                        valid_sem = counts > 1
                        if np.any(valid_sem):
                            with np.errstate(invalid="ignore", divide="ignore"):
                                sem[valid_sem] = (
                                    np.nanstd(
                                        site_values[:, valid_sem], axis=0, ddof=1
                                    )
                                    / np.sqrt(counts[valid_sem])
                                )
                        time_values = next(
                            item["times"]
                            for item in traces
                            if item["position"] == region
                            and item["band"] == band
                        )
                        color = "#16803C" if event == "+1" else "#A80000"
                        axis.plot(
                            time_values,
                            values,
                            color=color,
                            label=f"Event {event} (n sites={len(site_ids)})",
                        )
                        axis.fill_between(
                            time_values,
                            values - sem,
                            values + sem,
                            color=color,
                            alpha=0.2,
                            linewidth=0,
                        )
                        for index, (time, value) in enumerate(
                            zip(time_values, values)
                        ):
                            csv_rows.append(
                                {
                                    "sign": sign,
                                    "region": region,
                                    "band": band,
                                    "category": category,
                                    "event": event,
                                    "time": time,
                                    "mean_t_value": value,
                                    "sem_t_value": sem[index],
                                    "lower_sem": values[index] - sem[index],
                                    "upper_sem": values[index] + sem[index],
                                    "n_sites": len(site_ids),
                                }
                            )
                    if not plotted:
                        plt.close(figure)
                        continue
                    axis.axvline(0, color="black", linewidth=0.8, alpha=0.5)
                    axis.axvspan(
                        integration_tmin,
                        integration_tmax,
                        color="grey",
                        alpha=0.2,
                        linewidth=0,
                        label="Integration window",
                    )
                    axis.set_title(
                        f"{region.title()} — {band} — {category} ({sign})"
                    )
                    axis.set_xlabel("Time (s)")
                    axis.set_ylabel("t-value (event vs outcome 0)")
                    axis.grid(alpha=0.25)
                    axis.legend(loc="best")
                    figure.tight_layout()
                    figure.savefig(
                        figure_dir / f"{region}_power_over_time.png",
                        dpi=250,
                        bbox_inches="tight",
                    )
                    plt.close(figure)
                    pd.DataFrame(csv_rows).to_csv(
                        data_dir / f"{region}_power_over_time.csv",
                        index=False,
                    )

    print(f"Category source: {source.resolve()}")
    print(f"Figures: {output_dir.resolve()}")
    print(f"Trace CSVs: {csv_dir.resolve()}")


if __name__ == "__main__":
    main()
