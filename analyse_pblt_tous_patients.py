#!/usr/bin/env python3
"""
Analysis of TFR responses by frequency band and ventral/dorsal insular contacts
for all PBLT patients.

The script reuses functions from the project modules. The integrate_frequency
and integrate_time functions are taken from notebook cell 27 because they are
not present in an importable Python module.
paths_data, load_tfr_data, behavorial_organize_data, electrodes_lebelling_masks,
signal_processing.

Example:
    python analyse_pblt_tous_patients.py
    python analyse_pblt_tous_patients.py --patients GRE_2014_AUDe
    python analyse_pblt_tous_patients.py --tmin 0 --tmax 0.8 --alpha 0.05

Results are written to <--output>/ventral and <--output>/dorsal, with a
detailed CSV per contact and a global report.
"""

from pathlib import Path
from tqdm import tqdm
import argparse
import json
import sys
import traceback

from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
from matplotlib.ticker import FuncFormatter

matplotlib.use("Agg")  # Disable graphical windows for batch execution.
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from paths_data import path_to_PBLT, path_to_julich_output, PBLT_ROI_csv_name
from behavorial_organize_data import patient_names_from_root
from signal_processing import (
    compare_outcomes,
    compare_gain_loss_outcomes,
)
from plot_data import plot_outcome_comparison ,plot_bar_plot_compare_gain_loss_responses , plot_bar_plot_compare_gain_loss_responses_multiple_bands
from analysis_config import (
    ALPHA_DEFAULT,
    BANDS,
    BASELINE_MAX,
    CALCULATION_CHOICES,
    build_config,
    build_output_folder_name,
    save_config,
)
from pblt_pipeline import run_patient
from intermediate_data import save_intermediate_data, load_intermediate_data
from output_management import prepare_output_root, require_existing_output_root

 

def run_statistics_and_plots(
    integrated_data, output_root, tmin, tmax, alpha, calculation
):
    """Run statistics and contact plots after intermediate data are saved."""
    patient_results = []
    for values_df in integrated_data:
        patient = values_df["patient"].iloc[0]
        position = values_df["position"].iloc[0]
        contact_idx = int(values_df["contact_index"].iloc[0])
        contact_name = values_df["contact"].iloc[0]
        folder = position
        out_dir_fig = output_root / "fig" / "per_unit" / calculation / folder
        out_dir_data = output_root / "data" / folder
        if calculation == "minus1_vs_plus1":
            result_df = compare_outcomes(
                values_df,
                patient,
                contact_name,
                contact_idx,
                position,
                out_dir_data,
                tmin,
                tmax,
                alpha,
            )
            plot_outcome_comparison(
                values_df,
                result_df,
                patient,
                contact_name,
                contact_idx,
                position.title(),
                out_dir_fig,
                tmin,
                tmax,
                alpha,
                bands=list(BANDS.keys()),
            )
        elif calculation == "plus1_vs_zero_and_zero_vs_minus1":
            result_df = compare_gain_loss_outcomes(
                values_df, patient, contact_name, contact_idx, position, alpha
            )
        else:
            raise ValueError(f"Unsupported calculation: {calculation}")
        patient_results.append(result_df)

    return pd.concat(patient_results, ignore_index=True) if patient_results else pd.DataFrame()



def write_report(results, output_root, alpha, calculation_label, bands_for_final_bar_plot=None):
    """Save tables and a text report by position and frequency band."""
    output_root.mkdir(parents=True, exist_ok=True)
    if results.empty:
        (output_root / "compte_rendu.txt").write_text(
            "No statistical result was produced. Check the data, contact "
            "classification, and loading errors.\n",
            encoding="utf-8"
        )
        return

    results.to_csv(output_root / "resultats_statistiques_tous_patients.csv", index=False)
    for folder, position in [
        ("ventral", "Ventral"),
        ("dorsal", "Dorsal"),
        ("outside", "Outside"),
    ]:
        sub = results[results["position"] == folder]
        data_dir = output_root / "data" / folder
        data_dir.mkdir(parents=True, exist_ok=True)
        sub.to_csv(data_dir / "resultats_statistiques.csv", index=False)

    summary = (
        results.groupby(["band", "position"], as_index=False)
        .agg(
            n_sites=("contact_index", "size"),
            n_positive=("response", lambda values: (values == "positive").sum()),
            n_negative=("response", lambda values: (values == "negative").sum()),
        )
    )
    summary["positive_percent"] = 100 * summary["n_positive"] / summary["n_sites"]
    summary["negative_percent"] = 100 * summary["n_negative"] / summary["n_sites"]
    fisher_rows = []
    for band in BANDS:
        band_results = results[results["band"] == band]
        fisher_row = {"band": band}
        for sign in ("positive", "negative"):
            contingency = []
            for position in ("dorsal", "ventral"):
                position_results = band_results[
                    band_results["position"] == position
                ]
                n_significant = int(
                    (position_results["response"] == sign).sum()
                )
                n_not_significant = int(
                    len(position_results) - n_significant
                )
                contingency.append([n_significant, n_not_significant])
            odds_ratio, p_value = fisher_exact(
                contingency, alternative="two-sided"
            )
            fisher_row[f"fisher_{sign}_odds_ratio"] = odds_ratio
            fisher_row[f"fisher_{sign}_p_value"] = p_value
            fisher_row[f"fisher_{sign}_significant"] = p_value < alpha
        fisher_rows.append(fisher_row)
    summary = summary.merge(pd.DataFrame(fisher_rows), on="band", how="left")
    summary_dir = output_root / "data" / "summary" / "minus1_vs_plus1"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_dir / "minus1_vs_plus1_by_band_and_region.csv", index=False)
    _write_response_barplot(
        summary,
        output_root / "fig" / "summary" / "minus1_vs_plus1"
        / "minus1_vs_plus1_by_band_and_region.png",
        title="-1 vs +1 responses by frequency band and insular region",
        positive_subtitle="Positive responses (+1 > -1)",
        negative_subtitle="Negative responses (-1 > +1)",
    )

    lines = [
        "REPORT — statistical responses of insular contacts",
        "=" * 66,
        f"Comparison: {calculation_label}.",
        f"Test: two-sided Mann-Whitney; FDR correction across each contact's bands; alpha={alpha}.",
        "Positive response: median difference (outcome 2 - outcome 1) > 0 and FDR p < alpha.",
        "Negative response: median difference (outcome 2 - outcome 1) < 0 and FDR p < alpha.",
        "All other contact/band combinations are counted as non-significant.",
        "",
    ]
    for folder, position in [
        ("ventral", "Ventral"),
        ("dorsal", "Dorsal"),
        ("outside", "Outside"),
    ]:
        sub = results[results["position"] == folder]
        lines += [f"\n{position.upper()}", "-" * len(position)]
        if sub.empty:
            lines.append("No contact analyzed.")
            continue
        for band in BANDS:
            band_df = sub[sub["band"] == band]
            if band_df.empty:
                continue
            positive = int((band_df["response"] == "positive").sum())
            negative = int((band_df["response"] == "negative").sum())
            nonsig = int((band_df["response"] == "non_significant").sum())
            n_contacts = int(band_df["contact"].nunique())
            n_patients = int(band_df["patient"].nunique())
            lines.append(
                f"{band}: {positive} positive, {negative} negative, "
                f"{nonsig} non-significant (evaluated contacts={len(band_df)}, "
                f"distinct contacts={n_contacts}, patients={n_patients})."
            )
    lines += [
        "",
        "Note: counts are contact x band tests, not necessarily unique contacts "
        "responding in at least one band. P-values are corrected separately "
        "within each contact, as in the notebook.",
    ]
    (output_root / "compte_rendu.txt").write_text("\n".join(lines), encoding="utf-8")


def _write_response_barplot(
    summary,
    output_path,
    title,
    positive_subtitle,
    negative_subtitle,
):
    """Save a positive/negative response barplot for each band and region."""
    bands = list(BANDS)
    colors = {
        "dorsal": "#4C78A8",
        "ventral": "#F58518",
        "outside": "#777777",
    }
    light_colors = {
        "dorsal": "#A9C4E0",
        "ventral": "#F8C28D",
        "outside": "#C7C7C7",
    }
    figure, axis = plt.subplots(figsize=(12, 7))
    x_positions = np.arange(len(bands))
    width = 0.34
    max_value = 0.0

    present_positions = [
        position for position in ("dorsal", "ventral", "outside")
        if position in set(summary["position"])
    ]
    offsets = np.linspace(
        -width * (len(present_positions) - 1) / 2,
        width * (len(present_positions) - 1) / 2,
        len(present_positions),
    )
    for position, offset in zip(present_positions, offsets):
        subset = summary[summary["position"] == position].set_index("band")
        positive_values = [
            float(subset.loc[band, "positive_percent"])
            if band in subset.index else 0.0
            for band in bands
        ]
        negative_values = [
            float(subset.loc[band, "negative_percent"])
            if band in subset.index else 0.0
            for band in bands
        ]
        max_value = max(max_value, *(positive_values or [0]), *(negative_values or [0]))
        bars_x = x_positions + offset
        if position == "outside":
            positive_colors = [light_colors[position]] * len(bands)
            negative_colors = [light_colors[position]] * len(bands)
        else:
            positive_colors = [
                colors[position]
                if bool(subset.loc[band, "fisher_positive_significant"])
                else light_colors[position]
                for band in bands
            ]
            negative_colors = [
                colors[position]
                if bool(subset.loc[band, "fisher_negative_significant"])
                else light_colors[position]
                for band in bands
            ]
        axis.bar(bars_x, positive_values, width=width, color=positive_colors,
                 label=position.title())
        axis.bar(bars_x, [-value for value in negative_values], width=width,
                 color=negative_colors)

        for x_value, band, positive, negative in zip(
            bars_x, bands, positive_values, negative_values
        ):
            total = int(subset.loc[band, "n_sites"]) if band in subset.index else 0
            positive_count = (
                int(subset.loc[band, "n_positive"]) if band in subset.index else 0
            )
            negative_count = (
                int(subset.loc[band, "n_negative"]) if band in subset.index else 0
            )
            axis.annotate(
                f"{positive_count}/{total}\n{positive:.1f}%",
                (x_value, positive), ha="center", va="bottom", fontsize=8,
                xytext=(0, 3), textcoords="offset points",
            )
            axis.annotate(
                f"{negative_count}/{total}\n{negative:.1f}%",
                (x_value, -negative), ha="center", va="top", fontsize=8,
                xytext=(0, -3), textcoords="offset points",
            )

    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xticks(x_positions)
    axis.set_xticklabels(bands)
    axis.set_xlabel("Frequency band")
    axis.set_ylabel("Responding sites (%)")
    axis.set_title(title, pad=28)
    axis.text(0.5, 1.01, positive_subtitle, transform=axis.transAxes,
              ha="center", va="bottom", fontsize=10)
    axis.text(0.5, 0.975, negative_subtitle, transform=axis.transAxes,
              ha="center", va="top", fontsize=10)
    axis.text(
        0.5, 0.945,
        "Dark bars: Fisher test significant; light bars: not significant",
        transform=axis.transAxes, ha="center", va="top", fontsize=9,
        color="#444444",
    )
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{abs(value):g}%"))
    axis.legend(title="Region")
    axis.grid(axis="y", alpha=0.25)
    limit = max_value * 1.25 if max_value else 1
    axis.set_ylim(-limit, limit)
    figure.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)



def write_gain_loss_report(
    results,
    output_root,
    alpha,
    calculation_label,
    bands_for_final_bar_plot=None,
    overall_results=None,
    plot_region="all",
    save_figures=True,
):
    """Save gain/loss classifications and regional Fisher-tested plots."""
    output_root = Path(output_root)
    contrast_dir = output_root / "data" / "gain_loss_responses"
    figure_dir = output_root / "fig" / "summary" / "plus1_vs_zero_and_zero_vs_minus1"
    contrast_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    if plot_region not in {"all", "include", "exclude", "only"}:
        raise ValueError(f"Unknown plot_region: {plot_region}")
    selected_bands = (
        list(bands_for_final_bar_plot)
        if bands_for_final_bar_plot is not None
        else list(BANDS)
    )
    unknown_bands = set(selected_bands) - set(BANDS)
    if unknown_bands:
        raise ValueError(f"Unknown frequency bands: {sorted(unknown_bands)}")
    if selected_bands != list(BANDS) or plot_region != "all":
        band_label = "_".join(selected_bands)
        figure_dir = figure_dir / f"bands_{band_label}" / f"outside_{plot_region}"
    figure_dir.mkdir(parents=True, exist_ok=True)
    results = results[results["band"].isin(selected_bands)].copy()
    overall_results = overall_results[
        overall_results["band"].isin(selected_bands)
    ].copy() if overall_results is not None else None
    if plot_region == "exclude":
        results = results[results["position"] != "outside"].copy()
        if overall_results is not None:
            overall_results = overall_results[
                overall_results["position"] != "outside"
            ].copy()
    elif plot_region == "only":
        results = results[results["position"] == "outside"].copy()
        if overall_results is not None:
            overall_results = overall_results[
                overall_results["position"] == "outside"
            ].copy()
    results["gain_positive"] = (
        results["significant_gain"]
        & (results["median_gain_plus1_minus_zero"] > 0)
    )
    results["gain_negative"] = (
        results["significant_gain"]
        & (results["median_gain_plus1_minus_zero"] < 0)
    )
    results["loss_positive"] = (
        results["significant_loss"]
        & (results["median_loss_minus1_minus_zero"] > 0)
    )
    results["loss_negative"] = (
        results["significant_loss"]
        & (results["median_loss_minus1_minus_zero"] < 0)
    )
    if overall_results is None:
        raise ValueError(
            "Overall -1 vs +1 results are required to split dual gain/loss responses."
        )
    overall_columns = [
        "patient",
        "position",
        "contact_index",
        "band",
        "response",
    ]
    overall = overall_results[overall_columns].rename(
        columns={"response": "overall_response"}
    )
    results = results.merge(
        overall,
        on=["patient", "position", "contact_index", "band"],
        how="left",
        validate="one_to_one",
    )
    results.to_csv(contrast_dir / "gain_loss_results_by_band.csv", index=False)
    results["gain_response"] = np.select(
        [
            results["gain_positive"],
            results["gain_negative"],
        ],
        ["positive", "negative"],
        default="non_significant",
    )
    results["loss_response"] = np.select(
        [
            results["loss_positive"],
            results["loss_negative"],
        ],
        ["positive", "negative"],
        default="non_significant",
    )
    results["minus1_vs_plus1_response"] = results["overall_response"]
    results.to_csv(contrast_dir / "pairwise_comparison_results_by_band.csv", index=False)
    electrode_columns = [
        "patient",
        "position",
        "contact_index",
        "contact",
        "band",
    ]
    response_exports = [
        ("plus1_vs_minus1", "overall_response"),
        ("plus1_vs_zero", "gain_response"),
        ("minus1_vs_zero", "loss_response"),
    ]
    combined_exports = []
    for comparison_name, response_column in response_exports:
        export_columns = electrode_columns + [
            response_column,
            "p_fdr_gain" if comparison_name == "plus1_vs_zero" else
            "p_fdr_loss" if comparison_name == "minus1_vs_zero" else
            "p_fdr_overall" if "p_fdr_overall" in results else "overall_response",
        ]
        export_columns = list(dict.fromkeys(
            column for column in export_columns if column in results.columns
        ))
        comparison_results = results[export_columns].copy()
        comparison_results = comparison_results.rename(
            columns={response_column: "response"}
        )
        comparison_results.insert(5, "comparison", comparison_name)
        comparison_results.to_csv(
            contrast_dir / f"electrode_responses_{comparison_name}.csv",
            index=False,
        )
        for direction in ("positive", "negative"):
            comparison_results[
                comparison_results["response"] == direction
            ].to_csv(
                contrast_dir
                / f"electrode_responses_{comparison_name}_{direction}.csv",
                index=False,
            )
        combined_exports.append(comparison_results)
    pd.concat(combined_exports, ignore_index=True).to_csv(
        contrast_dir / "electrode_responses_all_comparisons.csv",
        index=False,
    )

    site_columns = ["patient", "position", "contact_index", "contact"]
    site_results = (
        results.groupby(site_columns, as_index=False)
        .agg(
            gain_positive=("gain_positive", "any"),
            gain_negative=("gain_negative", "any"),
            loss_positive=("loss_positive", "any"),
            loss_negative=("loss_negative", "any"),
        )
    )
    site_results["positive_response_category"] = np.select(
        [
            site_results["gain_positive"] & site_results["loss_positive"],
            site_results["gain_positive"],
            site_results["loss_positive"],
        ],
        ["both", "gain_only", "loss_only"],
        default="neither",
    )
    site_results["negative_response_category"] = np.select(
        [
            site_results["gain_negative"] & site_results["loss_negative"],
            site_results["gain_negative"],
            site_results["loss_negative"],
        ],
        ["both", "gain_only", "loss_only"],
        default="neither",
    )
    site_results.to_csv(contrast_dir / "gain_loss_site_responses.csv", index=False)
    categories = [
        "gain_only",
        "both_prefer_gain",
        "both_equal",
        "both_prefer_loss",
        "loss_only",
        "neither",
    ]
    for sign in ("positive", "negative"):
        category_column = f"{sign}_response_category"
        site_results.groupby(["position", category_column]).size().unstack(
            fill_value=0
        ).reindex(columns=categories, fill_value=0).to_csv(
            contrast_dir / f"gain_loss_site_responses_{sign}_by_region.csv"
        )

    comparison_specs = {
        "gain": {
            "difference": "median_gain_plus1_minus_zero",
            "significant": "significant_gain",
            "label": "+1 vs 0",
            "filename": "gain_plus1_vs_zero_by_band_and_region",
        },
        "loss": {
            "difference": "median_loss_minus1_minus_zero",
            "significant": "significant_loss",
            "label": "-1 vs 0",
            "filename": "loss_minus1_vs_zero_by_band_and_region",
        },
    }
    bands = selected_bands
    colors = {
        "dorsal": "#4C78A8",
        "ventral": "#F58518",
        "outside": "#777777",
    }
    light_colors = {
        "dorsal": "#A9C4E0",
        "ventral": "#F8C28D",
        "outside": "#C7C7C7",
    }
    for comparison, spec in comparison_specs.items():
        direction = np.where(
            results[spec["significant"]] & (results[spec["difference"]] > 0),
            "positive",
            np.where(
                results[spec["significant"]] & (results[spec["difference"]] < 0),
                "negative",
                "non_significant",
            ),
        )
        directed_results = results.assign(direction=direction)
        outside_significant = directed_results[
            (directed_results["position"] == "outside")
            & (directed_results["direction"] != "non_significant")
        ].copy()
        outside_significant.to_csv(
            contrast_dir
            / f"outside_significant_{comparison}_sites_by_band.csv",
            index=False,
        )
        summary = (
            directed_results.groupby(["band", "position"], as_index=False)
            .agg(
                n_sites=("contact_index", "size"),
                n_positive=("direction", lambda values: (values == "positive").sum()),
                n_negative=("direction", lambda values: (values == "negative").sum()),
            )
        )
        summary["positive_percent"] = 100 * summary["n_positive"] / summary["n_sites"]
        summary["negative_percent"] = 100 * summary["n_negative"] / summary["n_sites"]
        fisher_rows = []
        for band in bands:
            band_results = directed_results[directed_results["band"] == band]
            fisher_row = {"band": band}
            for sign in ("positive", "negative"):
                contingency = []
                for position in ("dorsal", "ventral"):
                    position_results = band_results[
                        band_results["position"] == position
                    ]
                    n_significant = int(
                        (position_results["direction"] == sign).sum()
                    )
                    n_not_significant = int(
                        len(position_results) - n_significant
                    )
                    contingency.append([n_significant, n_not_significant])
                odds_ratio, p_value = fisher_exact(
                    contingency, alternative="two-sided"
                )
                fisher_row[f"fisher_{sign}_odds_ratio"] = odds_ratio
                fisher_row[f"fisher_{sign}_p_value"] = p_value
                fisher_row[f"fisher_{sign}_significant"] = p_value < alpha
            fisher_rows.append(fisher_row)
        fisher_summary = pd.DataFrame(fisher_rows)
        summary = summary.merge(fisher_summary, on="band", how="left")
        summary.to_csv(
            contrast_dir / f"{spec['filename']}.csv",
            index=False,
        )

        figure, axis = plt.subplots(figsize=(12, 7))
        x_positions = np.arange(len(bands))
        width = 0.34
        max_value = 0.0
        present_positions = [
            position for position in ("dorsal", "ventral", "outside")
            if position in set(summary["position"])
        ]
        offsets = np.linspace(
            -width * (len(present_positions) - 1) / 2,
            width * (len(present_positions) - 1) / 2,
            len(present_positions),
        )
        for position, offset in zip(present_positions, offsets):
            subset = summary[summary["position"] == position].set_index("band")
            positive_values = [
                float(subset.loc[band, "positive_percent"])
                if band in subset.index else 0.0
                for band in bands
            ]
            negative_values = [
                float(subset.loc[band, "negative_percent"])
                if band in subset.index else 0.0
                for band in bands
            ]
            max_value = max(max_value, *(positive_values or [0]), *(negative_values or [0]))
            bars_x = x_positions + offset
            if position == "outside":
                positive_colors = [light_colors[position]] * len(bands)
                negative_colors = [light_colors[position]] * len(bands)
            else:
                positive_colors = [
                    colors[position]
                    if bool(subset.loc[band, "fisher_positive_significant"])
                    else light_colors[position]
                    for band in bands
                ]
                negative_colors = [
                    colors[position]
                    if bool(subset.loc[band, "fisher_negative_significant"])
                    else light_colors[position]
                    for band in bands
                ]
            axis.bar(
                bars_x, positive_values, width=width,
                color=positive_colors,
                label=position.title(),
            )
            axis.bar(
                bars_x, [-value for value in negative_values],
                width=width,
                color=negative_colors,
            )
            for x_value, band, positive, negative in zip(
                bars_x, bands, positive_values, negative_values
            ):
                total = int(subset.loc[band, "n_sites"]) if band in subset.index else 0
                positive_count = (
                    int(subset.loc[band, "n_positive"])
                    if band in subset.index else 0
                )
                negative_count = (
                    int(subset.loc[band, "n_negative"])
                    if band in subset.index else 0
                )
                axis.annotate(
                    f"{positive_count}/{total}\n{positive:.1f}%",
                    (x_value, positive), ha="center", va="bottom",
                    fontsize=8, xytext=(0, 3), textcoords="offset points",
                )
                axis.annotate(
                    f"{negative_count}/{total}\n{negative:.1f}%",
                    (x_value, -negative), ha="center", va="top",
                    fontsize=8, xytext=(0, -3), textcoords="offset points",
                )
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set_xticks(x_positions)
        axis.set_xticklabels(bands)
        axis.set_xlabel("Frequency band")
        axis.set_ylabel("Responding sites (%)")
        axis.set_title(f"{spec['label']} responses by frequency band and insular region", pad=28)
        axis.text(
            0.5, 1.01, f"Positive responses ({spec['label']}, first outcome > second outcome)",
            transform=axis.transAxes, ha="center", va="bottom", fontsize=10,
        )
        axis.text(
            0.5, 0.975, f"Negative responses ({spec['label']}, first outcome < second outcome)",
            transform=axis.transAxes, ha="center", va="top", fontsize=10,
        )
        axis.text(
            0.5, 0.945,
            "Dark bars: significant regional Fisher test; light bars: non-significant",
            transform=axis.transAxes, ha="center", va="top", fontsize=9,
            color="#444444",
        )
        axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{abs(value):g}%"))
        axis.legend(title="Region")
        axis.grid(axis="y", alpha=0.25)
        limit = max_value * 1.25 if max_value else 1
        axis.set_ylim(-limit, limit)
        figure.tight_layout()
        if save_figures:
            figure_path = (figure_dir / f"{spec['filename']}.png").resolve()
            figure.savefig(
                figure_path,
                dpi=300,
                bbox_inches="tight",
            )
        plt.close(figure)

    for sign in ("positive", "negative"):
        band_category_results = results.copy()
        both = (
            band_category_results[f"gain_{sign}"]
            & band_category_results[f"loss_{sign}"]
        )
        category_values = np.select(
            [
                band_category_results[f"gain_{sign}"]
                & ~band_category_results[f"loss_{sign}"],
                both & (band_category_results["overall_response"] == "positive"),
                both & (band_category_results["overall_response"] == "non_significant"),
                both & (band_category_results["overall_response"] == "negative"),
                band_category_results[f"loss_{sign}"]
                & ~band_category_results[f"gain_{sign}"],
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
        band_category_results["response_category"] = category_values
        band_category_results[f"{sign}_merged_response_category"] = category_values
        band_category_results.to_csv(
            contrast_dir / f"gain_loss_site_responses_{sign}_by_band.csv",
            index=False,
        )
        sign_figure_dir = figure_dir / sign
        if save_figures:
            plot_bar_plot_compare_gain_loss_responses_multiple_bands(
                band_category_results,
                sign_figure_dir,
                BANDS=BANDS,
                metric="percent",
                categories=categories,
                bands_for_final_bar_plot=bands_for_final_bar_plot,
                output_name=f"gain_loss_site_responses_{sign}_by_region_all_bands",
                response_sign=sign,
            )
            plot_bar_plot_compare_gain_loss_responses_multiple_bands(
                band_category_results,
                sign_figure_dir,
                BANDS=BANDS,
                metric="percent",
                categories=categories,
                bands_for_final_bar_plot=bands_for_final_bar_plot,
                output_name=f"merged_gain_loss_site_responses_{sign}_by_region_all_bands",
                response_sign=sign,
            )

    merged = results.copy()
    for sign in ("positive", "negative"):
        both = merged[f"gain_{sign}"] & merged[f"loss_{sign}"]
        merged[f"{sign}_merged_response_category"] = np.select(
            [
                merged[f"gain_{sign}"] & ~merged[f"loss_{sign}"],
                both & (merged["overall_response"] == "positive"),
                both & (merged["overall_response"] == "non_significant"),
                both & (merged["overall_response"] == "negative"),
                merged[f"loss_{sign}"] & ~merged[f"gain_{sign}"],
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
    merged.to_csv(
        contrast_dir / f"merged_comparisons_{plot_region}_by_band.csv",
        index=False,
    )
    merged.to_csv(
        contrast_dir / "pairwise_comparison_results_by_band.csv",
        index=False,
    )

    (contrast_dir / "report.txt").write_text(
        "\n".join(
            [
                "REPORT — gain/loss response categories",
                "=" * 48,
                f"Comparison: {calculation_label}.",
                "Gain response: significant +1 vs 0 gain comparison in at least one band.",
                "Loss response: significant 0 vs -1 loss comparison in at least one band.",
                "Electrode significance is determined by Benjamini-Hochberg FDR.",
                "Fisher's exact test compares significant versus non-significant "
                "electrodes between dorsal and ventral regions for each band.",
                f"Alpha: {alpha}.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_output_readme(output_root, calculation):
    """Document the output structure and the analyses generated by the run."""
    output_root = Path(output_root)
    readme = """# PBLT statistical analysis outputs

The output contains the integrated data and the selected first-stage analysis.
The regional mixed-model script can subsequently compute all three contrasts:

- `gain_plus1_vs_zero`: gain comparison between 0 and +1;
- `loss_zero_vs_minus1`: loss comparison between -1 and 0;
- `minus1_vs_plus1`: comparison between -1 and +1 outcomes.

The output-folder name is based on the baseline method, not on the selected
contrast. For example, `baseline_div_YYYYMMDD` uses division baseline
correction.

## Folders

- `config.json`: parameters and command arguments used for the analysis.
- `data/intermediate_data/`: integrated trial data, organized by insular region.
- `data/dorsal/` and `data/ventral/`: statistical result tables for `-1 vs +1`.
- `data/gain_loss_responses/`: tables for the separate gain/loss comparisons.
- `data/mixed_model_regression/`: global mixed-model results and site slopes.
- `data/summary/`: summary tables used by the frequency-band barplots.
- `fig/per_unit/<analysis>/`: one figure per analyzed contact and region.
- `fig/summary/<analysis>/`: summary barplots by frequency band and region.
- `fig/summary/mixed_model_regression/`: regional and site-slope figures.

Positive bars are displayed above zero and negative bars below zero. Negative
axis tick labels show absolute percentages for readability.
Electrode significance is determined by the FDR-adjusted electrode tests.
Fisher's exact test then compares responding versus non-responding electrode
counts between dorsal and ventral regions. Existing gain/loss bars are dark
when this regional Fisher test is significant and light otherwise.
"""
    (output_root / "README.md").write_text(readme, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patients", nargs="*", default=None,
                        help="Optional patient list (names without the _PBLT suffix).")
    parser.add_argument(
        "--calculation",
        choices=sorted(CALCULATION_CHOICES),
        default="minus1_vs_plus1",
        help="First-stage statistical comparison to run.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output directory. By default, outputs/<baseline>_<date> is used.",
    )
    parser.add_argument(
        "--data-root",
        default=None,
        help="Root directory containing the epoch-level TFR patient data.",
    )
    parser.add_argument(
        "--roi-file",
        default=None,
        help="ROI table used to classify contacts as dorsal, ventral, or outside.",
    )
    parser.add_argument(
        "--load-intermediate-data",
        action="store_true",
        help="Load integrated data from an existing output directory instead of reprocessing TFR data.",
    )
    parser.add_argument("--tmin", type=float, default=0.0, help="Integration window start (s).")
    parser.add_argument("--tmax", type=float, default=0.8, help="Integration window end (s).")
    parser.add_argument("--alpha", type=float, default=ALPHA_DEFAULT, help="Significance threshold.")
    parser.add_argument(
        "--channels_verbose",
        action="store_true",
        help="Display bipolar and monopolar channel information.",
    )
    parser.add_argument(
        "--take_outside",
        action="store_true",
        help="Include contacts classified as outside the insula.",
    )
    parser.add_argument('--baseline_type', choices=['div', 'zscore', 'log_div' ,"raw" , "diff"], default='div',
                        help="Baseline correction method: division (div), z-score (zscore), or log division (log_div).")
    parser.add_argument(
        "--bands", "--bffbp",
        dest="bffbp",
        nargs="+",
        choices=list(BANDS),
        default=list(BANDS),
        help="Frequency bands to include in the generated summary figures.",
    )
    parser.add_argument(
        "--plot-outside",
        choices=("all", "exclude", "only"),
        default="all",
        help="Show all regions, exclude outside, or show outside only in summary figures.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Compute and save statistical CSV results without generating figures.",
    )
    parser.add_argument(
        "--intermediate-only",
        action="store_true",
        help="Stop after writing intermediate data and configuration.",
    )
    
    args = parser.parse_args()

    calculation_label = CALCULATION_CHOICES[args.calculation]
    date_label = datetime.now().strftime("%m%d%H%M")
    if args.load_intermediate_data and not args.output:
        raise ValueError(
            "--output is required when --load-intermediate-data is enabled."
        )
    root = Path(args.data_root) if args.data_root else Path(path_to_PBLT)
    roi_csv = (
        Path(args.roi_file)
        if args.roi_file
        else Path(path_to_julich_output) / PBLT_ROI_csv_name
    )
    output_root = (
        Path(args.output)
        if args.output
        else Path("outputs")
        / build_output_folder_name(args.baseline_type, date_label)
    )
    config = build_config(
        tmin=args.tmin,
        tmax=args.tmax,
        alpha=args.alpha,
        calculation=args.calculation,
        channels_verbose=args.channels_verbose,
        baseline_type=args.baseline_type,
        take_outside=args.take_outside,
    )
    if args.load_intermediate_data:
        require_existing_output_root(output_root)
        existing_config_path = output_root / "config.json"
        existing_config = json.loads(
            existing_config_path.read_text(encoding="utf-8")
        )
        parameters = existing_config["parameters"]
        args.tmin = float(parameters["tmin"])
        args.tmax = float(parameters["tmax"])
        args.alpha = float(parameters["alpha"])
        args.baseline_type = parameters.get("baseline_type", args.baseline_type)
        args.take_outside = bool(parameters.get("take_outside", False))
        calculation_label = CALCULATION_CHOICES[args.calculation]
        parameters["calculation"] = args.calculation
        parameters["calculation_label"] = calculation_label
        parameters["calculation_choices"] = CALCULATION_CHOICES
        all_integrated_data = load_intermediate_data(output_root)
        failures = []
        patients = sorted(
            {frame["patient"].iloc[0] for frame in all_integrated_data}
        )
        config = existing_config
        config["arguments"] = {
            "patients_requested": args.patients,
            "patients_used": patients,
            "calculation": args.calculation,
            "output": str(output_root.resolve()),
            "load_intermediate_data": True,
            "tmin": args.tmin,
            "tmax": args.tmax,
            "alpha": args.alpha,
            "channels_verbose": args.channels_verbose,
            "baseline_type": args.baseline_type,
            "take_outside": args.take_outside,
            "bands_for_final_bar_plot": args.bffbp,
            "plot_outside": args.plot_outside,
        }
        config_path = save_config(config, output_root)
        print(f"Configuration updated: {config_path.resolve()}")
        print(
            f"Loaded intermediate integrated data for "
            f"{len(all_integrated_data)} contacts."
        )
    else:
        prepare_output_root(output_root, config)
        patients = args.patients if args.patients else patient_names_from_root(root)
        if not patients:
            raise RuntimeError(f"No patient found in {root}")
        config["arguments"] = {
            "patients_requested": args.patients,
            "patients_used": patients,
            "calculation": args.calculation,
            "output": str(output_root.resolve()),
            "load_intermediate_data": False,
            "tmin": args.tmin,
            "tmax": args.tmax,
            "alpha": args.alpha,
            "channels_verbose": args.channels_verbose,
            "baseline_type": args.baseline_type,
            "take_outside": args.take_outside,
            "bands_for_final_bar_plot": args.bffbp,
            "plot_outside": args.plot_outside,
        }
        config_path = save_config(config, output_root)
        print(f"Configuration saved to: {config_path.resolve()}")

        all_integrated_data = []
        failures = []
        for patient in tqdm(patients if len(patients) > 1 else patients, desc="Patients", unit="patient"):
            try:
                integrated_data = run_patient(
                    patient,
                    root,
                    roi_csv,
                    args.tmin,
                    args.tmax,
                    args.channels_verbose,
                    args.baseline_type,
                    args.take_outside,
                )
                if integrated_data:
                    all_integrated_data.extend(integrated_data)
                    print(
                        f"[OK] {patient}: "
                        f"{len(integrated_data)} integrated contacts prepared."
                    )
            except Exception as exc:
                failures.append((patient, repr(exc)))
                print(f"[ERROR] {patient}: {exc}")
                traceback.print_exc()

        save_intermediate_data(all_integrated_data, output_root)
        print(
            f"Intermediate integrated data saved for "
            f"{len(all_integrated_data)} contacts."
        )
        if args.intermediate_only:
            config["arguments"]["intermediate_only"] = True
            save_config(config, output_root)
            print(f"Intermediate-only analysis completed: {output_root.resolve()}")
            return
    results = run_statistics_and_plots(
        all_integrated_data,
        output_root,
        args.tmin,
        args.tmax,
        args.alpha,
        args.calculation,
    )
    if args.calculation == "minus1_vs_plus1":
        write_report(results, output_root, args.alpha, calculation_label , bands_for_final_bar_plot=args.bffbp)
    else:
        overall_results = pd.concat(
            [
                compare_outcomes(
                    values_df,
                    values_df["patient"].iloc[0],
                    values_df["contact"].iloc[0],
                    int(values_df["contact_index"].iloc[0]),
                    values_df["position"].iloc[0],
                    output_root / "data" / "mixed_model_regression",
                    args.tmin,
                    args.tmax,
                    args.alpha,
                )
                for values_df in all_integrated_data
            ],
            ignore_index=True,
        )
        (output_root / "data" / "mixed_model_regression").mkdir(
            parents=True, exist_ok=True
        )
        overall_results.to_csv(
            output_root / "data" / "mixed_model_regression" / "minus1_vs_plus1_by_band.csv",
            index=False,
        )
        write_gain_loss_report(
            results,
            output_root,
            args.alpha,
            calculation_label,
            bands_for_final_bar_plot=args.bffbp,
            overall_results=overall_results,
            plot_region=args.plot_outside,
            save_figures=not args.no_plots,
        )
    write_output_readme(output_root, args.calculation)
    if failures:
        (output_root / "patients_en_erreur.txt").write_text(
            "\n".join(f"{p}: {e}" for p, e in failures), encoding="utf-8"
        )
    print(f"\nCompleted. Results: {output_root.resolve()}")
    print(
        f"Patients processed with integrated data: "
        f"{len(patients) - len(failures)} / {len(patients)}"
    )
    if failures:
        print(f"Patients with errors: {len(failures)} (see patients_en_erreur.txt)")

if __name__ == "__main__":
    main()
