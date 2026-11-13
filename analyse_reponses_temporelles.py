"""Count response categories while varying the end of the integration window."""

import argparse
import json
from pathlib import Path
import sys
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

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
from paths_data import path_to_PBLT, path_to_julich_output
from signal_processing import integrate_frequency


CATEGORIES = [
    "gain_only",
    "both_prefer_gain",
    "both_equal",
    "both_prefer_loss",
    "loss_only",
]


def load_config(output: Path) -> dict:
    config_path = output / "config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing config.json: {config_path}")
    return json.loads(config_path.read_text(encoding="utf-8"))


def build_analysis_signature(args: argparse.Namespace) -> dict:
    return {
        "bands": list(args.bands),
        "regions": list(args.regions),
        "window_start": args.window_start,
        "window_end": args.window_end,
        "step": args.step,
        "interval_size": args.interval_size,
        "n_points": args.n_points,
        "mode": args.mode,
        "alpha": args.alpha,
    }


def load_cached_results(
    data_dir: Path, signature: dict
) -> Optional[pd.DataFrame]:
    result_path = data_dir / "temporal_response_by_site_and_window.csv"
    metadata_path = data_dir / "analysis_parameters.json"
    if not result_path.is_file() or not metadata_path.is_file():
        return None
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if metadata != signature:
        return None
    return pd.read_csv(result_path)


def save_analysis_metadata(data_dir: Path, signature: dict) -> None:
    (data_dir / "analysis_parameters.json").write_text(
        json.dumps(signature, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def smooth_values(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(values, kernel, mode="same")


def compare(values: np.ndarray, first: int, second: int) -> tuple[float, float]:
    x = values[first]
    y = values[second]
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    test = mannwhitneyu(x, y, alternative="two-sided")
    return float(np.median(x) - np.median(y)), float(test.pvalue)


def classify_window_rows(rows: list[dict], alpha: float) -> pd.DataFrame:
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    for comparison in ("gain", "loss", "overall"):
        p_column = f"p_raw_{comparison}"
        corrected_column = f"p_fdr_{comparison}"
        result[corrected_column] = np.nan
        for (_, _, _, band), indices in result.groupby(
            ["patient", "position", "contact_index", "band"]
        ).groups.items():
            valid = result.loc[indices, p_column].notna()
            valid_indices = result.loc[indices].index[valid]
            if len(valid_indices):
                result.loc[valid_indices, corrected_column] = multipletests(
                    result.loc[valid_indices, p_column], method="fdr_bh"
                )[1]
        result[f"significant_{comparison}"] = (
            result[corrected_column].notna()
            & (result[corrected_column] < alpha)
        )
    result["overall_response"] = np.select(
        [
        result["significant_overall"] & (result["overall_difference"] > 0),
        result["significant_overall"] & (result["overall_difference"] < 0),
        ],
        ["positive", "negative"],
        default="non_significant",
    )
    result["gain_positive"] = (
        result["significant_gain"] & (result["gain_difference"] > 0)
    )
    result["gain_negative"] = (
        result["significant_gain"] & (result["gain_difference"] < 0)
    )
    result["loss_positive"] = (
        result["significant_loss"] & (result["loss_difference"] > 0)
    )
    result["loss_negative"] = (
        result["significant_loss"] & (result["loss_difference"] < 0)
    )
    for sign in ("positive", "negative"):
        both = result[f"gain_{sign}"] & result[f"loss_{sign}"]
        result[f"{sign}_category"] = np.select(
            [
                result[f"gain_{sign}"] & ~result[f"loss_{sign}"],
                both & (result["overall_response"] == "positive"),
                both & (result["overall_response"] == "non_significant"),
                both & (result["overall_response"] == "negative"),
                result[f"loss_{sign}"] & ~result[f"gain_{sign}"],
            ],
            [
                "gain_only",
                "both_prefer_gain",
                "both_equal",
                "both_prefer_loss",
                "loss_only",
            ],
            default="neither",
        )
    return result


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
    parser.add_argument("--window-start", type=float, default=-1.0)
    parser.add_argument("--window-end", type=float, default=3.0)
    parser.add_argument(
        "--step",
        type=float,
        default=0.25,
        help="Distance between consecutive interval starts or cumulative endpoints.",
    )
    parser.add_argument(
        "--interval-size",
        "--interval_size",
        dest="interval_size",
        type=float,
        default=0.25,
        help="Duration of each consecutive integration interval in seconds.",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=None,
        help="Number of temporal windows. Overrides the generated sampling grid.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=1,
        help="Moving-average convolution width in plotted points; 1 disables smoothing.",
    )
    parser.add_argument(
        "--mode",
        choices=("interval", "cumulative", "both"),
        default="both",
        help="Analyze consecutive intervals, cumulative 0-to-tmax windows, or both.",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    if args.interval_size <= 0 or args.step <= 0:
        raise ValueError("--interval-size and --step must be positive.")
    if args.n_points is not None and args.n_points < 1:
        raise ValueError("--n-points must be at least 1.")
    if args.smooth_window < 1:
        raise ValueError("--smooth-window must be at least 1.")

    load_config(args.output)
    signature = build_analysis_signature(args)
    data_dir = args.output / "data" / "temporal_response_categories"
    figure_dir = args.output / "fig" / "summary" / "temporal_response_categories"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    results = load_cached_results(data_dir, signature)
    if results is not None:
        print(f"Reusing cached temporal response data: {data_dir.resolve()}")

    roi_data = pd.read_csv(Path(path_to_julich_output) / "PBLT_ROIs.csv", sep=";")
    patients = sorted(
        pd.read_csv(
            args.output
            / "data"
            / "gain_loss_responses"
            / "pairwise_comparison_results_by_band.csv"
        )["patient"].astype(str).unique()
    )
    if args.n_points is not None:
        interval_starts = np.linspace(
            args.window_start,
            args.window_end - args.interval_size,
            args.n_points,
        )
        cumulative_ends = np.linspace(0.0, args.window_end, args.n_points)
        cumulative_ends = cumulative_ends[cumulative_ends > 0]
    else:
        interval_starts = np.arange(
            args.window_start, args.window_end - args.interval_size / 2, args.step
        )
        cumulative_ends = np.arange(
            args.step, args.window_end + args.step / 2, args.step
        )
    window_specs = []
    if args.mode in ("interval", "both"):
        window_specs.extend(
            ("interval", start, start + args.interval_size)
            for start in interval_starts
        )
        print(f"Analyzing {len(interval_starts)} consecutive intervals of size {args.interval_size}s")
    if args.mode in ("cumulative", "both"):
        window_specs.extend(
            ("cumulative_0_to_tmax", 0.0, end)
            for end in cumulative_ends
        )
    rows = []
    root = Path(path_to_PBLT)
    if results is None:
      for patient in tqdm(patients):
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
        events = get_pos_events(patient, root)
        _, _, outcomes = prepare_outcomes(events, data.shape[0])
        channels = get_chanlist_from_info(info)
        mono = monopolar_channels_from_bipolar(channels)
        first, _ = bipolar_to_monopolar(channels, mono)
        classification = classify_contacts_insula(
            roi_data, first, patient.removesuffix("_PBLT")
        )
        band_signals = integrate_frequency(data, info["freqs"])
        for contact_index, contact_row in enumerate(
            classification.itertuples(index=False)
        ):
            position = str(contact_row.classification).lower()
            if position not in args.regions:
                continue
            for band in args.bands:
                signal = band_signals[band][:, :, contact_index]
                for window_mode, window_start, window_end in window_specs:
                    mask = (times >= window_start) & (times <= window_end)
                    if mask.sum() < 2:
                        continue
                    
                    try:
                        integrated = np.trapezoid(signal[:, mask], times[mask], axis=1)
                    except Exception as e:
                        integrated = np.trapz(signal[:, mask], times[mask], axis=1)
                        continue
                    values = {
                        code: integrated[outcomes == code]
                        for code in (0, 1, 2, 3)
                    }
                    gain_difference, gain_p = compare(values, 2, 0)
                    loss_difference, loss_p = compare(values, 1, 3)
                    overall_difference, overall_p = compare(values, 2, 1)
                    rows.append(
                        {
                            "patient": patient,
                            "position": position,
                            "contact_index": contact_index,
                            "band": band,
                            "window_start": window_start,
                            "window_end": window_end,
                            "window_mode": window_mode,
                            "interval_size": window_end - window_start,
                            "n_timepoints": int(mask.sum()),
                            "gain_difference": gain_difference,
                            "loss_difference": loss_difference,
                            "overall_difference": overall_difference,
                            "p_raw_gain": gain_p,
                            "p_raw_loss": loss_p,
                            "p_raw_overall": overall_p,
                        }
                    )

        results = classify_window_rows(rows, args.alpha)
        if results.empty:
            raise ValueError("No temporal response results were generated.")
        results.to_csv(
            data_dir / "temporal_response_by_site_and_window.csv", index=False
        )
        save_analysis_metadata(data_dir, signature)
    for window_mode in results["window_mode"].unique():
        mode_data_dir = data_dir / window_mode
        mode_figure_dir = figure_dir / window_mode
        mode_data_dir.mkdir(parents=True, exist_ok=True)
        mode_figure_dir.mkdir(parents=True, exist_ok=True)
        mode_results = results[results["window_mode"] == window_mode]
        mode_results.to_csv(
            mode_data_dir / "temporal_response_by_site_and_window.csv",
            index=False,
        )
        for band in args.bands:
            band_data_dir = mode_data_dir / band
            band_data_dir.mkdir(parents=True, exist_ok=True)
            mode_results[mode_results["band"] == band].to_csv(
                band_data_dir / "temporal_response_by_site_and_window.csv",
                index=False,
            )

    signs = ("positive", "negative") if args.sign == "both" else (args.sign,)
    count_rows = []
    for sign in signs:
        category_column = f"{sign}_category"
        selected = results[results[category_column].isin(CATEGORIES)]
        counts = (
            selected.groupby(
                [
                    "window_mode",
                    "window_start",
                    "window_end",
                    "band",
                    "position",
                    category_column,
                ]
            )
            .agg(n_sites=("contact_index", "nunique"))
            .reset_index()
            .rename(columns={category_column: "category"})
        )
        grid = (
            results[
                ["window_mode", "window_start", "window_end"]
            ].drop_duplicates()
            .merge(
                pd.DataFrame(
                    {
                        "band": args.bands,
                        "key": 1,
                    }
                ),
                how="cross",
            )
            .merge(
                pd.DataFrame(
                    {
                        "position": args.regions,
                        "key": 1,
                    }
                ),
                on="key",
            )
            .merge(
                pd.DataFrame({"category": CATEGORIES, "key": 1}),
                on="key",
            )
            .drop(columns="key")
        )
        counts = grid.merge(
            counts,
            on=[
                "window_mode",
                "window_start",
                "window_end",
                "band",
                "position",
                "category",
            ],
            how="left",
        )
        counts["n_sites"] = counts["n_sites"].fillna(0).astype(int)
        counts["sign"] = sign
        count_rows.append(counts)
        for band in args.bands:
            band_counts = counts[counts["band"] == band]
            for mode in band_counts["window_mode"].unique():
                mode_band_counts = band_counts[band_counts["window_mode"] == mode]
                mode_band_counts = mode_band_counts.copy()
                mode_band_counts["n_sites_smoothed"] = np.nan
                for (position, category), indices in mode_band_counts.groupby(
                    ["position", "category"]
                ).groups.items():
                    ordered = mode_band_counts.loc[indices].sort_values(
                        ["window_start", "window_end"]
                    )
                    mode_band_counts.loc[
                        ordered.index, "n_sites_smoothed"
                    ] = smooth_values(
                        ordered["n_sites"].to_numpy(dtype=float),
                        args.smooth_window,
                    )
                band_data_dir = data_dir / mode / band
                band_data_dir.mkdir(parents=True, exist_ok=True)
                mode_band_counts.to_csv(
                    band_data_dir / f"{sign}_response_category_counts.csv",
                    index=False,
                )
        for region in args.regions:
            for band in args.bands:
                for mode in counts["window_mode"].unique():
                    region_data = counts[
                        (counts["position"] == region)
                        & (counts["band"] == band)
                        & (counts["window_mode"] == mode)
                    ]
                    band_figure_dir = figure_dir / mode / band
                    band_figure_dir.mkdir(parents=True, exist_ok=True)
                    fig, axis = plt.subplots(figsize=(10, 6))
                    for category in CATEGORIES:
                        line = region_data[region_data["category"] == category]
                        x_values = (
                            line["window_end"]
                            if mode == "cumulative_0_to_tmax"
                            else line["window_start"]
                        )
                        y_values = line["n_sites"].to_numpy(dtype=float)
                        smoothed = smooth_values(y_values, args.smooth_window)
                        line = line.copy()
                        line["n_sites_smoothed"] = smoothed
                        axis.plot(
                            x_values,
                            line["n_sites_smoothed"],
                            marker="o",
                            label=category,
                        )
                    axis.set(
                        title=f"{region.title()} — {band} — "
                        f"{sign} response categories by {mode}",
                        xlabel="Integration end t (s)"
                        if mode == "cumulative_0_to_tmax"
                        else "Integration interval start (s)",
                        ylabel="Number of responding electrodes",
                    )
                    axis.grid(alpha=0.25)
                    axis.legend(loc="best")
                    fig.tight_layout()
                    fig.savefig(
                        band_figure_dir / f"{sign}_{region}_counts.png",
                        dpi=250,
                        bbox_inches="tight",
                    )
                    plt.close(fig)
    pd.concat(count_rows, ignore_index=True).to_csv(
        data_dir / "response_category_counts_by_tmax.csv", index=False
    )
    print(f"Temporal response data: {data_dir.resolve()}")
    print(f"Temporal response figures: {figure_dir.resolve()}")


if __name__ == "__main__":
    main()
