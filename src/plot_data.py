import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from analysis_config import format_band_label
def plot_bar_plot_compare_gain_loss_responses(percentages, figure_dir, metric=" \\% "):
    axis = percentages.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 6),
        color=["#4C78A8", "#F58518", "#54A24B", "#BDBDBD"],
    )

    # Ajouter la valeur dans chaque segment coloré
    for container in axis.containers:
        if metric == " \\% ":
            labels = [
                f"{bar.get_height():.1f}%" if bar.get_height() > 0 else ""
                for bar in container
            ]
        else:
            labels = [
                f"{bar.get_height():.0f}" if bar.get_height() > 0 else ""
                for bar in container
            ]
        axis.bar_label(
            container,
            labels=labels,
            label_type="center",
            fontsize=9,
            color="white",
        )

    axis.set_xlabel("Insular subregion")
    axis.set_ylabel("Sites (%)" if metric in {" \\% ", "percent"} else "Sites (count)")
    axis.set_title("Gain/loss response categories by insular subregion")
    axis.legend(
        title="Response category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    axis.set_xticklabels(
        [label.title() for label in percentages.index], rotation=0
    )

    if metric == " \\% ":
        axis.set_ylim(0, 100)

    axis.grid(axis="y", alpha=0.25)
    axis.set_axisbelow(True)

    axis.figure.tight_layout()
    axis.figure.savefig(
        figure_dir / f"gain_loss_site_responses_by_region{metric}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(axis.figure)
   
def plot_bar_plot_compare_gain_loss_responses_multiple_bands(
    results,
    figure_dir,
    BANDS,
    metric=" \\% ",
    categories=None,
    bands_for_final_bar_plot=None,
    output_name=None,
    response_sign=None,
    denominator_results=None,

):

    if categories is None:
        categories = [
            "gain_only",
            "both_prefer_gain",
            "both_equal",
            "both_prefer_loss",
            "loss_only",
            "neither",
        ]
    if bands_for_final_bar_plot is not None:
        bands = [band for band in bands_for_final_bar_plot if band in BANDS]
    else:
        bands = list(BANDS.keys())

    color_by_category = {
        "gain_only": "#16803C",
        "both_prefer_gain": "#8CCB9B",
        "both_equal": "#A67C52",
        "both_prefer_loss": "#F4A6A6",
        "loss_only": "#A80000",
        "neither": "#BDBDBD",
    }
    positions = [
        position for position in ("dorsal", "ventral", "outside")
        if position in set(results["position"])
    ]

    figure, axis = plt.subplots(figsize=(12, 7))

    x_positions = np.arange(len(bands))
    width = 0.34

    offsets = np.linspace(
        -width * (len(positions) - 1) / 2,
        width * (len(positions) - 1) / 2,
        len(positions),
    )
    for position, offset in zip(positions, offsets):
        counts_by_band = []
        totals_by_band = []

        for band in bands:
            subset = results[
                results["band"].str.lower().eq(band.lower())
                & results["position"].eq(position)
            ]
            denominator_subset = (
                denominator_results[
                    denominator_results["band"].str.lower().eq(band.lower())
                    & denominator_results["position"].eq(position)
                ]
                if denominator_results is not None
                else subset
            )
            counts_by_band.append(
                subset["response_category"]
                .value_counts()
                .reindex(categories, fill_value=0)
                .astype(int)
                .to_numpy()
            )
            totals_by_band.append(int(len(denominator_subset)))

        counts_by_band = np.asarray(counts_by_band, dtype=int)
        totals_by_band = np.asarray(totals_by_band, dtype=int)
        if metric in {" \\% ", "percent"}:
            values_by_category = np.divide(
                counts_by_band * 100,
                totals_by_band[:, None],
                out=np.zeros_like(counts_by_band, dtype=float),
                where=totals_by_band[:, None] > 0,
            )
        else:
            values_by_category = counts_by_band.astype(float)
        bars_x = x_positions + offset
        bottoms = np.zeros(len(bands))

        for i, category in enumerate(categories):
            heights = values_by_category[:, i]

            bars = axis.bar(
                bars_x,
                heights,
                width=width,
                bottom=bottoms,
                color=color_by_category[category],
                label=category if position == "dorsal" else "_nolegend_",
                edgecolor="white",
                linewidth=0.5,
            )

            labels = [
                f"{int(count)}/{int(total_sites)}" if count > 0 else ""
                for count, total_sites in zip(
                    counts_by_band[:, i], totals_by_band
                )
            ]

            axis.bar_label(
                bars,
                labels=labels,
                label_type="center",
                fontsize=8,
                color="white",
            )

            bottoms += heights

    axis.set_xticks(x_positions)
    axis.set_xticklabels(bands)
    axis.set_xlabel("Frequency band")
    axis.set_ylabel("Sites (%)" if metric in {" \\% ", "percent"} else "Sites (count)")
    sign_label = f" ({response_sign} responses)" if response_sign else ""
    axis.set_title(
        "Gain/loss response categories by frequency band and insular subregion"
        + sign_label,
        pad=15,
    )

    if metric in {" \\% ", "percent"}:
        axis.set_ylim(0, 100)

    # Legend of the six sign-specific response categories.
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(
        handles,
        labels,
        title="Response category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    axis.grid(axis="y", alpha=0.25)
    axis.set_axisbelow(True)

    # Indicate the region below each bar group.
    for x in x_positions:
        for position, offset in zip(positions, offsets):
            axis.text(
                x + offset,
                -0.08,
                {"dorsal": "D", "ventral": "V", "outside": "O"}[position],
                transform=axis.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=9,
            )

    figure.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    if metric == " \\% ":
        metric_ = "percent"
    if metric == " (count) ":
        metric_ = "count"
    filename = output_name or f"gain_loss_site_responses_by_region_all_bands{metric_}"
    figure_path = (Path(figure_dir) / f"{filename}.png").resolve()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        figure_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)

def plot_tfr(tfr, times, freqs, trial= None, channel = None, title=None, cnorm=None, cmap="viridis", horizontal_line_time=None , save_file_path=None, save_file=False,show_fig=True):
    plt.figure(figsize=(10, 5))
    plt.pcolormesh(times, freqs, tfr.T, shading="auto", cmap=cmap, norm=cnorm)
    if horizontal_line_time is not None:
        for hline in horizontal_line_time:
            plt.axvline(x=hline, color='r', linestyle='--', label=f'Temps {hline}s')
    plt.yscale("log")
    plt.xlabel("Temps (s)")
    plt.ylabel("Frequency (Hz)")
    if title:
        plt.title(title)
    else:
        plt.title(f"TFR — essai {trial}, canal {channel}")
    plt.colorbar(label="Valeur TFR")
    plt.tight_layout()
    if save_file:
        if save_file_path is None:
            save_file_path = f"TFR_trial_{trial}_channel_{channel}.png"
        plt.savefig(save_file_path, dpi=300)
    if show_fig:
        plt.show()
    plt.close()

def plot_outcome_comparison(df, result_df, patient, contact_name, contact_idx, position,
                            output_dir, tmin, tmax, alpha=0.05 , bands=None):

    # One plot per contact with all bands; no intermediate plots.
    fig, axes = plt.subplots(1, max(1, len(bands)), figsize=(4 * max(1, len(bands)), 4.8), squeeze=False)
    axes = axes.ravel()
    for ax, band in zip(axes, bands):
        x = df.loc[df["outcome"] == 1, band].dropna().to_numpy()
        y = df.loc[df["outcome"] == 2, band].dropna().to_numpy()
        if len(x) and len(y):
            ax.boxplot([x, y], positions=[1, 2], widths=0.5)
            rng = np.random.default_rng(42)
            ax.scatter(1 + rng.uniform(-0.08, 0.08, len(x)), x, alpha=0.65, s=14)
            ax.scatter(2 + rng.uniform(-0.08, 0.08, len(y)), y, alpha=0.65, s=14)
        row = result_df.loc[result_df["band"] == band]
        p_fdr = row["p_fdr"].iloc[0] if len(row) else np.nan
        marker = "*" if np.isfinite(p_fdr) and p_fdr < alpha else "ns"
        band_label = format_band_label(band)
        ax.set_title(f"{band_label}\nFDR p={p_fdr:.3g} ({marker})" if np.isfinite(p_fdr) else f"{band_label}\np unavailable")
        ax.set_xticks([1, 2], ["Outcome 1\n(-1 €)", "Outcome 2\n(+1 €)"])
        ax.set_ylabel("Integrated power (normalized TFR)")
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(
        f"{patient} — {position} — contact {contact_name}\n"
        f"Outcome 1 (-1 EUR) vs outcome 2 (+1 EUR), window {tmin:g}-{tmax:g} s"
    )
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_contact = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(contact_name))
    fig.savefig(output_dir / f"{patient}_{position}_contact{contact_idx}_{safe_contact}.png",
                dpi=160, bbox_inches="tight")
    plt.close(fig)
    return result_df
