#!/usr/bin/env python3
"""Fit the regional mixed model to integrated intermediate PBLT data.

For every contrast and frequency band, the model is:

    response_value ~ outcome_binary * ventral
                     + outcome_binary:C(site_id) + C(site_id)
                     + (1 | patient)

The outcome-by-region interaction directly tests whether the outcome slope
differs between dorsal and ventral regions. A second site-slope model adds
outcome-by-site interactions to derive one fitted outcome slope per site.

Contrasts:
    gain: outcome 0 versus outcome +1
    loss: outcome -1 versus outcome 0
    overall: outcome -1 versus outcome +1
"""

from pathlib import Path
import argparse
import json
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm, ttest_1samp
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from analysis_config import BANDS, find_latest_output_folder, format_band_label
from intermediate_data import load_intermediate_data
CONTRASTS = {
    "gain_plus1_vs_zero": {
        "outcomes": (0, 2),
        "labels": ("0", "+1"),
        "description": "gain: +1 versus 0",
    },
    "loss_zero_vs_minus1": {
        "outcomes": (1, 3),
        "labels": ("-1", "0"),
        "description": "loss: 0 versus -1",
    },
    "minus1_vs_plus1": {
        "outcomes": (1, 2),
        "labels": ("-1", "+1"),
        "description": "overall: +1 versus -1",
    },
}
P_ADJUST_METHODS = {
    "fdr_bh": "Benjamini-Hochberg FDR",
    "holm": "Holm",
    "bonferroni": "Bonferroni",
    "none": "No correction",
}
MODEL_PRESETS = {
    "regional": {
        "formula": "response_value ~ outcome_binary * ventral + C(site_id)",
        "description": "Regional interaction with site fixed effects",
    },
    "regional_patient": {
        "formula": "response_value ~ outcome_binary * ventral + C(site_id)",
        "description": "Regional interaction with site fixed effects and patient random intercept",
    },
    "main_effects": {
        "formula": "response_value ~ outcome_binary + ventral + C(site_id)",
        "description": "Outcome and region main effects without interaction",
    },
    "outcome_only": {
        "formula": "response_value ~ outcome_binary + C(site_id)",
        "description": "Outcome effect with site fixed effects",
    },
}


def prepare_trials(integrated_data, contrast_name):
    """Combine site data and encode the selected task contrast as 0/1."""
    contrast = CONTRASTS[contrast_name]
    frames = []
    for site_data in integrated_data:
        subset = site_data[site_data["outcome"].isin(contrast["outcomes"])].copy()
        if subset.empty:
            continue
        outcome_to_binary = {
            contrast["outcomes"][0]: 0,
            contrast["outcomes"][1]: 1,
        }
        subset["outcome_binary"] = subset["outcome"].map(outcome_to_binary)
        subset["ventral"] = (subset["position"] == "ventral").astype(int)
        subset["site_id"] = (
            subset["patient"].astype(str)
            + "|"
            + subset["position"].astype(str)
            + "|"
            + subset["contact_index"].astype(str)
        )
        frames.append(subset)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def normal_two_sided_p_value(estimate, standard_error):
    """Compute a Wald p-value, returning NaN for unavailable uncertainty."""
    if not np.isfinite(estimate) or not np.isfinite(standard_error):
        return np.nan
    if standard_error <= 0:
        return 0.0 if estimate != 0 else np.nan
    return float(2 * norm.sf(abs(estimate / standard_error)))


def adjust_p_values(values, method):
    """Apply the selected correction while preserving unavailable values."""
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    adjusted = np.full(values.shape, np.nan)
    if valid.any():
        if method == "none":
            adjusted[valid] = values[valid]
        else:
            adjusted[valid] = multipletests(values[valid], method=method)[1]
    return adjusted


def fit_band_model(model_data, band, model_name="regional"):
    """Fit the regional mixed model and derive one slope per site."""
    selected = model_data[
        [
            "site_id",
            "patient",
            "position",
            "ventral",
            "outcome_binary",
            band,
        ]
    ].dropna()
    selected = selected.rename(columns={band: "response_value"})
    if (
        selected["site_id"].nunique() < 2
        or selected["outcome_binary"].nunique() < 2
        or selected["ventral"].nunique() < 2
    ):
        raise ValueError("At least two sites, both outcomes, and both regions are required.")

    if model_name not in MODEL_PRESETS:
        raise ValueError(f"Unknown model preset: {model_name}")
    region_formula = MODEL_PRESETS[model_name]["formula"]
    site_formula = (
        "response_value ~ outcome_binary * ventral + "
        "outcome_binary:C(site_id) + C(site_id)"
    )
    fitted = None
    fit_method = None
    fit_errors = []
    def fit_model(formula, groups):
        model = smf.mixedlm(
            formula,
            selected,
            groups=groups,
            re_formula="1",
        )
        fitted = None
        fit_method = None
        fit_errors = []
        for method in ("lbfgs", "powell", "cg"):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    candidate = model.fit(
                        reml=True, method=method, maxiter=1000, disp=False
                    )
                if getattr(candidate, "converged", False):
                    fitted = candidate
                    fit_method = method
                    break
                fit_errors.append(f"{method}: non-converged")
            except (np.linalg.LinAlgError, ValueError, RuntimeError) as exc:
                fit_errors.append(f"{method}: {exc}")
        if fitted is None:
            raise RuntimeError(
                "No optimizer converged (" + "; ".join(fit_errors) + ")"
            )
        return fitted, fit_method

    region_fitted, region_fit_method = fit_model(
        region_formula, selected["patient"]
    )
    site_fitted, site_fit_method = fit_model(site_formula, selected["patient"])
    try:
        interaction_p = float(region_fitted.pvalues.get("outcome_binary:ventral", np.nan))
        regional_params = region_fitted.params
        regional_covariance = region_fitted.cov_params()
        outcome_name = "outcome_binary"
        interaction_name = "outcome_binary:ventral"
        dorsal_estimate = float(regional_params.get(outcome_name, np.nan))
        ventral_estimate = dorsal_estimate + float(
            regional_params.get(interaction_name, np.nan)
        )
        dorsal_se = float(
            np.sqrt(regional_covariance.loc[outcome_name, outcome_name])
        )
        ventral_variance = (
            regional_covariance.loc[outcome_name, outcome_name]
            + regional_covariance.loc[interaction_name, interaction_name]
            + 2 * regional_covariance.loc[outcome_name, interaction_name]
        )
        ventral_se = float(np.sqrt(max(ventral_variance, 0)))
        dorsal_p = normal_two_sided_p_value(dorsal_estimate, dorsal_se)
        ventral_p = normal_two_sided_p_value(ventral_estimate, ventral_se)
        interaction_estimate = float(region_fitted.params.get("outcome_binary:ventral", np.nan))
        interaction_standard_error = float(
            region_fitted.bse.get("outcome_binary:ventral", np.nan)
        )
    except (AttributeError, TypeError):
        interaction_p = np.nan
        dorsal_estimate = np.nan
        ventral_estimate = np.nan
        dorsal_se = np.nan
        ventral_se = np.nan
        dorsal_p = np.nan
        ventral_p = np.nan
        interaction_estimate = np.nan
        interaction_standard_error = np.nan

    site_rows = []
    site_metadata = selected.drop_duplicates("site_id")[
        ["site_id", "patient", "position", "ventral"]
    ]
    for site in site_metadata.itertuples(index=False):
        prediction_data = pd.DataFrame(
            {
                "site_id": [site.site_id, site.site_id],
                "ventral": [site.ventral, site.ventral],
                "outcome_binary": [0, 1],
            }
        )
        predictions = np.asarray(site_fitted.predict(prediction_data), dtype=float)
        site_rows.append(
            {
                "band": band,
                "site_id": site.site_id,
                "patient": site.patient,
                "position": site.position,
                "site_slope": float(predictions[1] - predictions[0]),
            }
        )

    region_model_summary = region_fitted.summary().as_text()
    return {
        "band": band,
        "interaction_p_value": interaction_p,
        "interaction_estimate": interaction_estimate,
        "interaction_standard_error": interaction_standard_error,
        "dorsal_delta": dorsal_estimate,
        "dorsal_delta_se": dorsal_se,
        "dorsal_delta_p_value": dorsal_p,
        "ventral_delta": ventral_estimate,
        "ventral_delta_se": ventral_se,
        "ventral_delta_p_value": ventral_p,
        "interaction_ci_low": (
            interaction_estimate - 1.96 * interaction_standard_error
            if np.isfinite(interaction_estimate)
            and np.isfinite(interaction_standard_error)
            else np.nan
        ),
        "interaction_ci_high": (
            interaction_estimate + 1.96 * interaction_standard_error
            if np.isfinite(interaction_estimate)
            and np.isfinite(interaction_standard_error)
            else np.nan
        ),
        "n_trials": len(selected),
        "n_sites": selected["site_id"].nunique(),
        "response_scale_factor": (
            1.0 / float(selected["response_value"].std(ddof=1))
            if float(selected["response_value"].std(ddof=1)) > 0
            else np.nan
        ),
        "region_fit_method": region_fit_method,
        "region_summary": region_model_summary,
        "site_fit_method": site_fit_method,
        "site_summary": site_fitted.summary().as_text(),
        "site_slopes": pd.DataFrame(site_rows),
    }


def analyze_contrast(
    integrated_data,
    contrast_name,
    p_adjust_method="fdr_bh",
    model_name="regional",
):
    """Fit one model per frequency band for one gain or loss contrast."""
    model_data = prepare_trials(integrated_data, contrast_name)
    model_rows = []
    site_slope_frames = []
    errors = []
    for band in BANDS:
        try:
            model_info = fit_band_model(model_data, band, model_name=model_name)
            site_slope_frames.append(model_info.pop("site_slopes"))
            model_rows.append(model_info)
        except (ValueError, np.linalg.LinAlgError, RuntimeError) as exc:
            errors.append({"band": band, "error": str(exc)})

    model_results = pd.DataFrame(model_rows)
    if not model_results.empty:
        model_results["p_adjust_method"] = p_adjust_method
        model_results["model_preset"] = model_name
        valid = model_results["interaction_p_value"].notna()
        corrected = adjust_p_values(
            model_results["interaction_p_value"].to_numpy(), p_adjust_method
        )
        model_results["interaction_p_fdr"] = corrected
        for column in ("dorsal_delta_p_value", "ventral_delta_p_value"):
            valid_delta = model_results[column].notna()
            corrected_delta = adjust_p_values(
                model_results[column].to_numpy(), p_adjust_method
            )
            model_results[column.replace("_p_value", "_p_fdr")] = corrected_delta
        model_results = model_results.drop(
            columns=[
            ],
            errors="ignore",
        )
    site_results = (
        pd.concat(site_slope_frames, ignore_index=True)
        if site_slope_frames
        else pd.DataFrame()
    )
    return model_results, site_results, errors


def significance_stars(p_value):
    """Return significance stars, or no annotation for non-significant values."""
    if not np.isfinite(p_value) or p_value >= 0.05:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def plot_regional_delta_forest(
    model_results, contrast_name, output_path, scaled=False
):
    """Plot the regional slope difference with 95% confidence intervals."""
    available = model_results.dropna(
        subset=["interaction_estimate", "interaction_ci_low", "interaction_ci_high"]
    )
    bands = [band for band in BANDS if band in set(available["band"])]
    if not bands:
        return
    figure, axis = plt.subplots(
        1, 1, figsize=(8.5, max(4.5, 0.9 * len(bands) + 2.5))
    )
    y_positions = np.arange(len(bands))[::-1]
    values = available.set_index("band").loc[bands]
    if scaled:
        values = values.copy()
        scale_by_band = model_results.set_index("band")["response_scale_factor"]
        values["scale_factor"] = values.index.to_series().map(scale_by_band)
        for column in (
            "interaction_estimate",
            "interaction_ci_low",
            "interaction_ci_high",
        ):
            values[column] = values[column] * values["scale_factor"]
    estimates = values["interaction_estimate"].to_numpy(dtype=float)
    ci_low = values["interaction_ci_low"].to_numpy(dtype=float)
    ci_high = values["interaction_ci_high"].to_numpy(dtype=float)
    axis.errorbar(
        estimates,
        y_positions,
        xerr=[
            estimates - ci_low,
            ci_high - estimates,
        ],
        fmt="o",
        color="#333333",
        ecolor="#333333",
        elinewidth=2,
        capsize=5,
        markersize=7,
    )
    x_low = float(np.nanmin(ci_low))
    x_high = float(np.nanmax(ci_high))
    x_span = x_high - x_low or 1.0
    axis.set_xlim(x_low - 0.08 * x_span, x_high + 0.22 * x_span)
    for y_position, band in zip(y_positions, bands):
        p_value = values.loc[band, "interaction_p_fdr"]
        axis.text(
            ci_high[bands.index(band)] + 0.02 * x_span,
            y_position,
            f"p={p_value:.3g}",
            va="center",
            fontsize=10,
        )
    axis.axvline(0, color="black", linewidth=0.9)
    axis.set_xlabel(
        "Δ ventral − Δ dorsal (regional outcome slope"
        + (" scaled by band" if scaled else "")
        + ", 95% CI)"
    )
    tick_labels = bands
    if scaled:
        tick_labels = [
            f"{format_band_label(band)}\n"
            f"(scale × {values.loc[band, 'scale_factor']:.3g})"
            for band in bands
        ]
    axis.set_yticks(y_positions, tick_labels)
    axis.set_title(
        f"Regional slope difference: {CONTRASTS[contrast_name]['description']}"
    )
    axis.grid(axis="x", alpha=0.25)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_regional_deltas(
    model_results, contrast_name, output_path, scaled=False
):
    """Plot dorsal and ventral regional outcome slopes side by side."""
    required = ["dorsal_delta", "ventral_delta", "dorsal_delta_se", "ventral_delta_se"]
    available = model_results.dropna(subset=required).copy()
    bands = [band for band in BANDS if band in set(available["band"])]
    if not bands:
        return
    if scaled:
        scales = model_results.set_index("band")["response_scale_factor"]
        available["scale_factor"] = available["band"].map(scales)
        for column in ("dorsal_delta", "ventral_delta", "dorsal_delta_se", "ventral_delta_se"):
            available[column] = available[column] * available["scale_factor"]
    values = available.set_index("band").loc[bands]
    figure, axis = plt.subplots(
        1, 1, figsize=(8.5, max(4.5, 0.9 * len(bands) + 2.5))
    )
    y_centers = np.arange(len(bands))[::-1]
    colors = {"dorsal": "#4C78A8", "ventral": "#F58518"}
    for region, offset, estimate_column, se_column in (
        ("dorsal", 0.12, "dorsal_delta", "dorsal_delta_se"),
        ("ventral", -0.12, "ventral_delta", "ventral_delta_se"),
    ):
        estimates = values[estimate_column].to_numpy(dtype=float)
        errors = 1.96 * values[se_column].to_numpy(dtype=float)
        axis.errorbar(
            estimates,
            y_centers + offset,
            xerr=errors,
            fmt="o",
            color=colors[region],
            ecolor=colors[region],
            capsize=4,
            linewidth=2,
            label=region.capitalize(),
        )
        p_column = f"{region}_delta_p_fdr"
        for band_index, band in enumerate(bands):
            p_value = values.loc[band, p_column]
            estimate = values.loc[band, estimate_column]
            interval_low = estimate - 1.96 * values.loc[band, se_column]
            interval_high = estimate + 1.96 * values.loc[band, se_column]
            star = (
                significance_stars(p_value)
                if interval_low > 0 or interval_high < 0
                else ""
            )
            if star:
                axis.text(
                    values.loc[band, estimate_column]
                    + 1.96 * values.loc[band, se_column]
                    + 0.02 * (values[estimate_column].max() - values[estimate_column].min() or 1),
                    y_centers[band_index] + offset,
                    star,
                    color=colors[region],
                    va="center",
                    fontweight="bold",
                )
    interaction_stars = [
        significance_stars(values.loc[band, "interaction_p_fdr"])
        for band in bands
    ]
    x_max = float(
        (values["dorsal_delta"] + 1.96 * values["dorsal_delta_se"]).max()
    )
    x_max = max(
        x_max,
        float(
            (values["ventral_delta"] + 1.96 * values["ventral_delta_se"]).max()
        ),
    )
    x_min = float(
        (values["dorsal_delta"] - 1.96 * values["dorsal_delta_se"]).min()
    )
    x_min = min(
        x_min,
        float(
            (values["ventral_delta"] - 1.96 * values["ventral_delta_se"]).min()
        ),
    )
    x_span = x_max - x_min or 1.0
    for band_index, star in enumerate(interaction_stars):
        if star:
            axis.text(
                x_max + 0.06 * x_span,
                y_centers[band_index],
                star,
                color="black",
                va="center",
                fontweight="bold",
            )
    axis.axvline(0, color="black", linewidth=0.9)
    tick_labels = bands
    if scaled:
        tick_labels = [
            f"{format_band_label(band)}\n"
            f"(scale × {values.loc[band, 'scale_factor']:.3g})"
            for band in bands
        ]
    axis.set_yticks(y_centers, tick_labels)
    axis.set_xlabel(
        "Regional outcome slope"
        + (" scaled by band" if scaled else "")
        + " (95% CI)"
    )
    axis.set_title(
        f"Regional outcome slopes: {CONTRASTS[contrast_name]['description']}"
    )
    axis.legend(frameon=False)
    axis.grid(axis="x", alpha=0.25)
    axis.set_xlim(x_min - 0.05 * x_span, x_max + 0.18 * x_span)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_site_slopes(
    site_results, model_results, contrast_name, output_path, scaled=False
):
    """Plot site slopes in one independently scaled subfigure per band."""
    if site_results.empty:
        return
    bands = [band for band in BANDS if band in set(site_results["band"])]
    if not bands:
        return
    values = site_results.copy()
    scale_by_band = model_results.set_index("band")["response_scale_factor"]
    values["scale_factor"] = values["band"].map(scale_by_band)
    if scaled:
        values["site_slope"] = values["site_slope"] * values["scale_factor"]

    colors = {"dorsal": "#4C78A8", "ventral": "#F58518"}
    figure, axis = plt.subplots(
        1,
        len(bands),
        figsize=(max(4.0 * len(bands), 10), 5.8),
        squeeze=False,
    )
    axes = axis.ravel()
    rng = np.random.default_rng(42)
    p_values_by_region = {}
    for region in ("dorsal", "ventral"):
        raw_p_values = []
        for band in bands:
            region_values = values.loc[
                (values["band"] == band) & (values["position"] == region),
                "site_slope",
            ].dropna()
            raw_p_values.append(
                float(ttest_1samp(region_values, 0).pvalue)
                if len(region_values) > 1
                else np.nan
            )
        adjusted = np.full(len(bands), np.nan)
        valid = np.isfinite(raw_p_values)
        if valid.any():
            adjusted[valid] = multipletests(
                np.asarray(raw_p_values)[valid], method="fdr_bh"
            )[1]
        p_values_by_region[region] = adjusted

    interaction_values = model_results.set_index("band")[
        "interaction_p_fdr"
    ]
    for band_index, (axis, band) in enumerate(zip(axes, bands)):
        plot_values = []
        for position, region in enumerate(("dorsal", "ventral"), start=1):
            region_values = values.loc[
                (values["band"] == band) & (values["position"] == region),
                "site_slope",
            ].dropna()
            plot_values.append(region_values.to_numpy(dtype=float))
            axis.scatter(
                position + rng.uniform(-0.09, 0.09, len(region_values)),
                region_values,
                color=colors[region],
                alpha=0.65,
                s=22,
                zorder=2,
            )
        axis.boxplot(
            plot_values,
            positions=[1, 2],
            widths=0.5,
            patch_artist=True,
            boxprops={"facecolor": "#A6A6A6", "alpha": 0.35},
            medianprops={"color": "black", "linewidth": 1.2},
            whiskerprops={"color": "#555555"},
            capprops={"color": "#555555"},
            zorder=3,
        )
        non_empty_values = [item for item in plot_values if len(item)]
        all_values = (
            np.concatenate(non_empty_values)
            if non_empty_values
            else np.array([], dtype=float)
        )
        if len(all_values):
            value_range = float(np.ptp(all_values)) or 1.0
            y_min = float(np.min(all_values)) - 0.10 * value_range
            y_max = float(np.max(all_values)) + 0.25 * value_range
            axis.set_ylim(y_min, y_max)
            for position, region in enumerate(("dorsal", "ventral")):
                region_values = plot_values[position]
                if len(region_values):
                    star = significance_stars(
                        p_values_by_region[region][band_index]
                    )
                    if star:
                        axis.text(
                            position + 1,
                            float(np.max(region_values)) + 0.10 * value_range,
                            star,
                            color=colors[region],
                            ha="center",
                            va="bottom",
                            fontsize=15,
                            fontweight="bold",
                        )
            interaction_p = interaction_values.get(band, np.nan)
            interaction_star = significance_stars(interaction_p)
            if np.isfinite(interaction_p):
                annotation = (
                    f"{interaction_star} p={interaction_p:.3g}"
                    if interaction_star
                    else f"p={interaction_p:.3g}"
                )
                axis.text(
                    1.5,
                    float(np.max(all_values)) + 0.17 * value_range,
                    annotation,
                    color="black",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold" if interaction_star else "normal",
                )
        scale_text = (
            f"\nscale × {scale_by_band.loc[band]:.3g}" if scaled else ""
        )
        axis.set_xticks([1, 2], ["Dorsal", "Ventral"])
        axis.set_title(f"{format_band_label(band)}{scale_text}")
        axis.axhline(0, color="black", linewidth=0.9)
        axis.grid(axis="y", alpha=0.25)
        if band_index == 0:
            axis.set_ylabel(
                "Site-specific outcome slope"
                + (" scaled by band" if scaled else "")
            )
    figure.suptitle(
        f"Site-specific outcome slopes: "
        f"{CONTRASTS[contrast_name]['description']}"
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def save_results(
    output_root,
    contrast_name,
    model_results,
    site_results,
    errors,
):
    """Save model tables, summaries, plots, and model text summaries."""
    output_root = Path(output_root)
    data_dir = output_root / "data" / "mixed_model_regression" / contrast_name
    figure_root = (
        output_root
        / "fig"
        / "summary"
        / "mixed_model_regression"
        / "regional_model"
    )
    raw_dir = figure_root / "raw"
    scaled_dir = figure_root / "scaled"
    site_dir = figure_root / "site_slopes"
    data_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    scaled_dir.mkdir(parents=True, exist_ok=True)
    site_dir.mkdir(parents=True, exist_ok=True)
    model_results.drop(
        columns=["summary", "region_summary", "site_summary"],
        errors="ignore",
    ).to_csv(
        data_dir / "model_results_by_band.csv", index=False
    )
    site_results.to_csv(data_dir / "site_slopes_by_band.csv", index=False)
    model_text = []
    for row in model_results.itertuples(index=False):
        model_text.extend(
            [
                f"Band: {row.band}",
                "Region-interaction model:",
                row.region_summary,
                "Site-slope model:",
                row.site_summary,
                "",
            ]
        )
    if errors:
        model_text.extend(["Model errors:"])
        model_text.extend(f"{item['band']}: {item['error']}" for item in errors)
    (data_dir / "model_summaries.txt").write_text(
        "\n".join(model_text) + "\n", encoding="utf-8"
    )
    (data_dir / "p_value_notes.txt").write_text(
        "interaction_p_fdr: FDR-adjusted p-value across the frequency "
        "bands for this contrast.\n"
        "dorsal_delta_p_fdr and ventral_delta_p_fdr: FDR-adjusted "
        "p-values for the dorsal and ventral outcome slopes against zero.\n"
        "The regional model is response_value ~ outcome_binary * ventral "
        "+ C(site_id) + (1 | patient).\n"
        "The site-slope model adds outcome_binary:C(site_id) to this formula. "
        "Site-specific slopes include the site interaction and the regional "
        "ventral effect where applicable.\n"
        "Colored figure stars are shown only when the corresponding regional "
        "95% confidence interval excludes zero. Black stars indicate a "
        "significant FDR-adjusted regional interaction.\n",
        encoding="utf-8",
    )
    plot_regional_delta_forest(
        model_results,
        contrast_name,
        raw_dir / f"{contrast_name}_regional_delta_forest.png",
    )
    plot_regional_delta_forest(
        model_results,
        contrast_name,
        scaled_dir / f"{contrast_name}_regional_delta_forest.png",
        scaled=True,
    )
    plot_regional_deltas(
        model_results,
        contrast_name,
        raw_dir / f"{contrast_name}_regional_deltas.png",
        scaled=False,
    )
    plot_regional_deltas(
        model_results,
        contrast_name,
        scaled_dir / f"{contrast_name}_regional_deltas.png",
        scaled=True,
    )
    plot_site_slopes(
        site_results,
        model_results,
        contrast_name,
        site_dir / f"{contrast_name}_site_slopes.png",
        scaled=False,
    )


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--output-type",
        default="baseline_div",
        help="Output prefix used when selecting the latest baseline folder.",
    )
    parser.add_argument(
        "--outputs-root", type=Path, default=script_dir / "outputs"
    )
    parser.add_argument(
        "--p-adjust",
        choices=sorted(P_ADJUST_METHODS),
        default="fdr_bh",
        help="Multiple-testing correction for LMM p-values.",
    )
    parser.add_argument(
        "--model",
        choices=sorted(MODEL_PRESETS),
        default="regional_patient",
        help="LMM formula preset.",
    )
    args = parser.parse_args()
    output_root = args.output or find_latest_output_folder(
        args.outputs_root, args.output_type
    )
    integrated_data = load_intermediate_data(output_root)
    all_errors = {}
    for contrast_name in CONTRASTS:
        model_results, site_results, errors = analyze_contrast(
            integrated_data,
            contrast_name,
            p_adjust_method=args.p_adjust,
            model_name=args.model,
        )
        save_results(
            output_root,
            contrast_name,
            model_results,
            site_results,
            errors,
        )
        all_errors[contrast_name] = errors

    config_path = output_root / "config.json"
    if config_path.is_file():
        config = json.loads(config_path.read_text(encoding="utf-8"))
        config.setdefault("additional_analyses", {})["mixed_model_regression"] = {
            "region_interaction_model": (
                "response_value ~ outcome_binary * ventral + C(site_id) "
                "+ (1 | patient)"
            ),
            "site_slope_model": (
                "response_value ~ outcome_binary * ventral + "
                "outcome_binary:C(site_id) + C(site_id) "
                "+ (1 | patient)"
            ),
            "contrasts": {
                key: value["description"] for key, value in CONTRASTS.items()
            },
            "random_patient_intercept": True,
            "regional_test": "outcome_binary:ventral interaction",
            "multiple_testing_correction": (
            P_ADJUST_METHODS[args.p_adjust]
            ),
            "p_adjust_method": args.p_adjust,
            "model_preset": args.model,
            "model_description": MODEL_PRESETS[args.model]["description"],
        }
        config_path.write_text(
            json.dumps(config, indent=2, ensure_ascii=True) + "\n",
            encoding="utf-8",
        )
    print(f"Completed mixed-model regression analysis: {output_root.resolve()}")
    for contrast_name, errors in all_errors.items():
        if errors:
            print(f"{contrast_name}: {len(errors)} band(s) could not be fitted.")


if __name__ == "__main__":
    main()
