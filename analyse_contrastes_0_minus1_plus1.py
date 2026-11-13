#!/usr/bin/env python3
"""Analyze gain and loss contrasts separately.

For every site and frequency band, this script computes:

    gain: (+1 gain success) - (0 gain failure)
    loss: (0 loss success) - (-1 loss failure)

The two contrast distributions are compared with a two-sided Mann-Whitney
test for each site and frequency band.
"""

from pathlib import Path
import argparse
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from analysis_config import BANDS, OUTPUT_TYPE_CHOICES, find_latest_output_folder
from intermediate_data import load_intermediate_data


GAIN_ZERO_OUTCOME = 0
LOSS_MINUS_ONE_OUTCOME = 1
GAIN_PLUS_ONE_OUTCOME = 2
LOSS_ZERO_OUTCOME = 3


def compute_site_contrasts(site_data, alpha):
    """Compute gain/loss means and a comparison for one site and band."""
    rows = []
    for band in BANDS:
        gain_zero = site_data.loc[site_data["outcome"] == GAIN_ZERO_OUTCOME, band]
        gain_plus_one = site_data.loc[
            site_data["outcome"] == GAIN_PLUS_ONE_OUTCOME, band
        ]
        loss_minus_one = site_data.loc[
            site_data["outcome"] == LOSS_MINUS_ONE_OUTCOME, band
        ]
        loss_zero = site_data.loc[site_data["outcome"] == LOSS_ZERO_OUTCOME, band]

        gain_zero_values = gain_zero.dropna().to_numpy(dtype=float)
        gain_plus_one_values = gain_plus_one.dropna().to_numpy(dtype=float)
        loss_minus_one_values = loss_minus_one.dropna().to_numpy(dtype=float)
        loss_zero_values = loss_zero.dropna().to_numpy(dtype=float)

        gain_zero_mean = (
            float(np.mean(gain_zero_values)) if len(gain_zero_values) else np.nan
        )
        gain_plus_one_mean = (
            float(np.mean(gain_plus_one_values))
            if len(gain_plus_one_values)
            else np.nan
        )
        loss_minus_one_mean = (
            float(np.mean(loss_minus_one_values))
            if len(loss_minus_one_values)
            else np.nan
        )
        loss_zero_mean = (
            float(np.mean(loss_zero_values)) if len(loss_zero_values) else np.nan
        )

        gain_contrast = (
            gain_plus_one_mean - gain_zero_mean
            if np.isfinite(gain_plus_one_mean) and np.isfinite(gain_zero_mean)
            else np.nan
        )
        loss_contrast = (
            loss_zero_mean - loss_minus_one_mean
            if np.isfinite(loss_zero_mean) and np.isfinite(loss_minus_one_mean)
            else np.nan
        )

        # Build trial-level contrast distributions without mixing gain and loss
        # outcome categories.
        gain_contrast_values = (
            gain_plus_one_values[:, None] - gain_zero_values[None, :]
            if len(gain_plus_one_values) and len(gain_zero_values)
            else np.array([])
        ).ravel()
        loss_contrast_values = (
            loss_zero_values[:, None] - loss_minus_one_values[None, :]
            if len(loss_zero_values) and len(loss_minus_one_values)
            else np.array([])
        ).ravel()
        if len(gain_contrast_values) and len(loss_contrast_values):
            p_value = float(
                mannwhitneyu(
                    gain_contrast_values,
                    loss_contrast_values,
                    alternative="two-sided",
                ).pvalue
            )
        else:
            p_value = np.nan

        rows.append(
            {
                "patient": site_data["patient"].iloc[0],
                "position": site_data["position"].iloc[0],
                "contact_index": int(site_data["contact_index"].iloc[0]),
                "contact": str(site_data["contact"].iloc[0]),
                "band": band,
                "n_gain_zero": len(gain_zero_values),
                "n_gain_plus1": len(gain_plus_one_values),
                "n_loss_minus1": len(loss_minus_one_values),
                "n_loss_zero": len(loss_zero_values),
                "mean_gain_zero": gain_zero_mean,
                "mean_gain_plus1": gain_plus_one_mean,
                "mean_loss_minus1": loss_minus_one_mean,
                "mean_loss_zero": loss_zero_mean,
                "gain_contrast_plus1_minus_zero": gain_contrast,
                "loss_contrast_zero_minus_minus1": loss_contrast,
                "contrast_difference_gain_minus_loss": (
                    gain_contrast - loss_contrast
                    if np.isfinite(gain_contrast) and np.isfinite(loss_contrast)
                    else np.nan
                ),
                "p_value": p_value,
                "significant": bool(np.isfinite(p_value) and p_value < alpha),
            }
        )
    return rows


def analyze_intermediate_data(integrated_data, alpha):
    """Compute gain/loss contrasts for every site and frequency band."""
    results = []
    for site_data in integrated_data:
        results.extend(compute_site_contrasts(site_data, alpha))
    return pd.DataFrame(results)


def summarize_contrasts(results):
    """Summarize mean contrasts and standard errors across all sites."""
    summary = (
        results.groupby("band", sort=False)
        .agg(
            n_sites=("gain_contrast_plus1_minus_zero", "count"),
            mean_gain_contrast=("gain_contrast_plus1_minus_zero", "mean"),
            sem_gain_contrast=("gain_contrast_plus1_minus_zero", "sem"),
            mean_loss_contrast=("loss_contrast_zero_minus_minus1", "mean"),
            sem_loss_contrast=("loss_contrast_zero_minus_minus1", "sem"),
            n_significant=("significant", "sum"),
        )
        .reindex(list(BANDS))
        .reset_index()
    )
    summary["significant_percent"] = (
        100 * summary["n_significant"] / summary["n_sites"]
    )
    return summary


def save_contrast_results(results, output_root, alpha, source_config):
    """Save site results, band summary, plot, report, and configuration."""
    output_root = Path(output_root)
    output_dir = output_root / "data" / "contrast_gain_vs_loss"
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "contrast_results.csv", index=False)

    summary = summarize_contrasts(results)
    summary.to_csv(output_dir / "contrast_summary_by_band.csv", index=False)

    x = np.arange(len(summary))
    width = 0.36
    figure, axis = plt.subplots(figsize=(11, 6))
    axis.bar(
        x - width / 2,
        summary["mean_gain_contrast"],
        width,
        yerr=summary["sem_gain_contrast"],
        capsize=4,
        label="Gain: +1 - 0",
        color="#4C78A8",
    )
    axis.bar(
        x + width / 2,
        summary["mean_loss_contrast"],
        width,
        yerr=summary["sem_loss_contrast"],
        capsize=4,
        label="Loss: 0 - (-1)",
        color="#F58518",
    )
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xticks(x)
    axis.set_xticklabels(summary["band"])
    axis.set_xlabel("Frequency band")
    axis.set_ylabel("Mean integrated-power contrast across sites")
    axis.set_title("Gain and loss contrasts across all sites")
    axis.legend()
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_dir / "contrast_summary_by_band.png", dpi=300, bbox_inches="tight")
    plt.close(figure)

    report = [
        "REPORT — gain and loss contrasts",
        "=" * 48,
        "Gain contrast: gain +1 success - gain 0 failure.",
        "Loss contrast: loss 0 success - loss -1 failure.",
        "Test: two-sided Mann-Whitney test between trial-level gain and loss contrasts.",
        f"Alpha: {alpha}.",
    ]
    (output_dir / "contrast_report.txt").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )

    config = {
        "algorithm": "PBLT gain versus loss contrast analysis",
        "source_output": str(output_root.resolve()),
        "source_config": source_config,
        "parameters": {
            "alpha": float(alpha),
            "gain_zero_outcome": GAIN_ZERO_OUTCOME,
            "gain_plus_one_outcome": GAIN_PLUS_ONE_OUTCOME,
            "loss_minus_one_outcome": LOSS_MINUS_ONE_OUTCOME,
            "loss_zero_outcome": LOSS_ZERO_OUTCOME,
            "test": "two-sided Mann-Whitney on separate gain/loss pairwise contrasts",
        },
        "bands": source_config.get("bands", {}),
    }
    (output_dir / "config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return output_dir


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--output-type",
        choices=sorted(OUTPUT_TYPE_CHOICES),
        default="baseline_div",
        help="Baseline output type used to find the newest output folder.",
    )
    parser.add_argument("--outputs-root", type=Path, default=script_dir / "outputs")
    parser.add_argument("--alpha", type=float, default=None)
    args = parser.parse_args()

    output_root = args.output or find_latest_output_folder(
        args.outputs_root, args.output_type
    )
    source_config_path = output_root / "config.json"
    if not source_config_path.is_file():
        raise FileNotFoundError(f"Source configuration not found: {source_config_path}")
    source_config = json.loads(source_config_path.read_text(encoding="utf-8"))
    alpha = (
        float(args.alpha)
        if args.alpha is not None
        else float(source_config["parameters"]["alpha"])
    )

    integrated_data = load_intermediate_data(output_root)
    results = analyze_intermediate_data(integrated_data, alpha)
    result_dir = save_contrast_results(results, output_root, alpha, source_config)
    n_sites = results[["patient", "position", "contact_index"]].drop_duplicates().shape[0]
    print(f"Contrast results saved in: {result_dir.resolve()}")
    print(f"Analyzed sites: {n_sites}")


if __name__ == "__main__":
    main()
