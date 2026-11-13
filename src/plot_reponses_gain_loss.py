#!/usr/bin/env python3
"""Generate plots of positive and negative responses by region.

The script reads the summary produced by extraire_reponses_statistiques.py and
saves two figures: one as a percentage of sites and one as a site count. For
each band, dorsal bars are on the left and ventral bars on the right; positive
responses are above zero and negative responses are below zero.

Exemple :
    python run_plot_reponses_gain_loss.py
    python run_plot_reponses_gain_loss.py --summary path/to/summary.csv
"""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

try:
    from .analysis_config import BANDS, format_band_label, load_config
except ImportError:
    from analysis_config import BANDS, format_band_label, load_config

POSITION_LABELS = {
    "dorsal": "Dorsal",
    "ventral": "Ventral",
    "outside": "Outside",
}
COLORS = {
    "dorsal": "#4C78A8",
    "ventral": "#F58518",
    "outside": "#777777",
}
LIGHT_COLORS = {
    "dorsal": "#A9C4E0",
    "ventral": "#F8C28D",
    "outside": "#C7C7C7",
}


def load_summary(summary_path: Path) -> pd.DataFrame:
    required = {
        "band",
        "position",
        "n_sites",
        "n_positive",
        "n_negative",
        "positive_percent",
        "negative_percent",
        "fisher_p_value",
        "fisher_significant",
    }
    summary = pd.read_csv(summary_path)
    missing = required - set(summary.columns)
    if missing:
        raise ValueError(
            f"Missing columns in summary {summary_path}: "
            + ", ".join(sorted(missing))
        )
    return summary


def plot_summary(
    summary: pd.DataFrame,
    value_kind: str,
    output_path: Path,
    bands=None,
) -> None:
    """Create one figure with a pair of positions for each band."""
    if value_kind == "percent":
        positive_column = "positive_percent"
        negative_column = "negative_percent"
        ylabel = "Responding sites (%)"
        title = "Gain/loss responses by frequency band and region"
        positive_subtitle = "Positive responses (gain > loss)"
        negative_subtitle = "Negative responses (gain < loss)"
        label_format = "{count}/{total}\n{value:.1f}%"
    elif value_kind == "count":
        positive_column = "n_positive"
        negative_column = "n_negative"
        ylabel = "Number of responding sites"
        title = "Number of gain/loss responses by frequency band and region"
        positive_subtitle = "Positive responses (gain > loss)"
        negative_subtitle = "Negative responses (gain < loss)"
        label_format = "{value:.0f}\n({count}/{total})"
    else:
        raise ValueError(f"Unknown value type: {value_kind}")

    bands = bands or list(BANDS)
    figure, axis = plt.subplots(figsize=(12, 7))
    band_positions = range(len(bands))
    width = 0.34
    positions = [
        position for position in ("dorsal", "ventral", "outside")
        if position in set(summary["position"])
    ]
    offsets = {
        position: offset
        for position, offset in zip(
            positions,
            [float(value) for value in np.linspace(
                -width * (len(positions) - 1) / 2,
                width * (len(positions) - 1) / 2,
                len(positions),
            )],
        )
    }
    bar_positions = {
        position: [index + offsets[position] for index in band_positions]
        for position in positions
    }

    max_value = 0.0
    for position in positions:
        position_label = POSITION_LABELS[position]
        subset = summary[summary["position"] == position].set_index("band")
        positive_values = [
            float(subset.loc[band, positive_column]) if band in subset.index else 0.0
            for band in bands
        ]
        negative_values = [
            float(subset.loc[band, negative_column]) if band in subset.index else 0.0
            for band in bands
        ]
        max_value = max(max_value, *(positive_values or [0]), *(negative_values or [0]))
        x_values = bar_positions[position]
        bar_colors = [
            COLORS[position]
            if position != "outside" and bool(subset.loc[band, "fisher_significant"])
            else LIGHT_COLORS[position]
            for band in bands
        ]
        axis.bar(
            x_values,
            positive_values,
            width=width,
            color=bar_colors,
            label=position_label,
        )
        axis.bar(
            x_values,
            [-value for value in negative_values],
            width=width,
            color=bar_colors,
        )

        for x_value, band, positive, negative in zip(
            x_values, bands, positive_values, negative_values
        ):
            total = int(subset.loc[band, "n_sites"]) if band in subset.index else 0
            positive_count = int(subset.loc[band, "n_positive"]) if band in subset.index else 0
            negative_count = int(subset.loc[band, "n_negative"]) if band in subset.index else 0
            axis.annotate(
                label_format.format(value=positive, count=positive_count, total=total),
                (x_value, positive),
                ha="center",
                va="bottom",
                fontsize=8,
                xytext=(0, 3),
                textcoords="offset points",
            )
            axis.annotate(
                label_format.format(value=negative, count=negative_count, total=total),
                (x_value, -negative),
                ha="center",
                va="top",
                fontsize=8,
                xytext=(0, -3),
                textcoords="offset points",
            )

    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xticks(list(band_positions))
    axis.set_xticklabels([format_band_label(band) for band in bands])
    axis.set_ylabel(ylabel)
    axis.set_xlabel("Frequency band")
    axis.set_title(title, pad=28)
    axis.text(
        0.5,
        1.01,
        positive_subtitle,
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=10,
    )
    axis.text(
        0.5,
        0.975,
        negative_subtitle,
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=10,
    )
    axis.text(
        0.5,
        0.945,
        "Dark bars: Fisher test significant; light bars: not significant",
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
    )
    axis.yaxis.set_major_formatter(
        FuncFormatter(
            lambda value, _: (
                f"{abs(value):g}%" if value_kind == "percent" else f"{abs(value):g}"
            )
        )
    )
    axis.legend(title="Region")
    axis.grid(axis="y", alpha=0.25)
    limit = max_value * 1.25 if max_value else 1
    axis.set_ylim(-limit, limit)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
