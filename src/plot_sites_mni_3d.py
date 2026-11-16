"""Plot response-category sites in MNI space with Julich parcel backgrounds."""

import argparse
from pathlib import Path
import re
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from scipy.stats import ttest_ind
from paths_data import path_to_julich_output


CATEGORY_COLORS = {
    "gain_only": "#16803C",
    "both_prefer_gain": "#8CCB9B",
    "both_equal": "#A67C52",
    "both_prefer_loss": "#F4A6A6",
    "loss_only": "#A80000",
    "neither": "#BDBDBD",
}
EXPORT_CATEGORY_COLORS = {
    "loss_only": CATEGORY_COLORS["loss_only"],
    "gain_only": CATEGORY_COLORS["gain_only"],
    "both": CATEGORY_COLORS["both_equal"],
    "both_prefer_gain": CATEGORY_COLORS["both_prefer_gain"],
    "both_prefer_loss": CATEGORY_COLORS["both_prefer_loss"],
    "neither": CATEGORY_COLORS["neither"],
}
PARCEL_COLORS = [
    "#8DD3C7",
    "#FFFFB3",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
]


def normalize_contact(value: object) -> str:
    text = str(value).strip().upper().replace(" ", "")
    match = re.fullmatch(r"([A-Z]+['’‘]?)0*(\d+)", text)
    if not match:
        return text
    prefix = match.group(1).replace("’", "'").replace("‘", "'")
    return f"{prefix}{int(match.group(2))}"


def extract_parcel_id(value: object) -> str:
    match = re.search(r"\b(Ia|Id)(\d+)\b", str(value))
    if not match:
        raise ValueError(f"Unable to extract a Juelich parcel ID from {value!r}")
    return f"{match.group(1)}{match.group(2)}"


def load_coordinates(path: Path, origin_path: Optional[Path] = None) -> pd.DataFrame:
    data = pd.read_csv(path, sep=";")
    required = {"PAT", "CHAN", "X", "Y", "Z", "Julich_insula_region"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing coordinate columns in {path}: {sorted(missing)}")
    data = data.rename(
        columns={
            "PAT": "patient",
            "CHAN": "contact",
            "Julich_insula_region": "parcel",
        }
    )
    data["contact_key"] = data["contact"].map(normalize_contact)
    for column in ("X", "Y", "Z"):
        data[column] = pd.to_numeric(
            data[column].astype(str).str.replace(",", ".", regex=False),
            errors="coerce",
        )
    data = data.dropna(subset=["X", "Y", "Z"])
    data["parcel"] = data["parcel"].replace(
        to_replace=r"^\s*(?:nan|none)?\s*$", value=np.nan, regex=True
    )
    known = data.dropna(subset=["parcel"]).copy()
    for index, row in data[data["parcel"].isna()].iterrows():
        if known.empty:
            raise ValueError("No known Juelich parcel is available for assignment.")
        distances = np.linalg.norm(
            known[["X", "Y", "Z"]].to_numpy()
            - row[["X", "Y", "Z"]].to_numpy(dtype=float),
            axis=1,
        )
        data.loc[index, "parcel"] = known.iloc[int(np.argmin(distances))]["parcel"]
    data["parcel_id"] = data["parcel"].map(extract_parcel_id)
    if data["parcel_id"].isna().any():
        raise ValueError("Every coordinate must have a Juelich parcel ID.")
    if origin_path is not None and origin_path.is_file():
        origin = pd.read_excel(origin_path)
        origin_required = {"PAT", "CHAN", "X", "Y", "Z"}
        missing_origin = origin_required - set(origin.columns)
        if missing_origin:
            raise ValueError(
                f"Missing coordinate columns in {origin_path}: {sorted(missing_origin)}"
            )
        origin = origin.rename(columns={"PAT": "patient", "CHAN": "contact"})
        origin["contact_key"] = origin["contact"].map(normalize_contact)
        for column in ("X", "Y", "Z"):
            origin[column] = pd.to_numeric(origin[column], errors="coerce")
        known_keys = set(zip(data["patient"], data["contact_key"]))
        missing_origin = origin[
            ~origin.apply(
                lambda row: (row["patient"], row["contact_key"]) in known_keys,
                axis=1,
            )
        ].dropna(subset=["X", "Y", "Z"])
        if not missing_origin.empty:
            known = data[
                ["X", "Y", "Z", "parcel", "parcel_id", "Insula_DV",
                 "Julich_insula_region_rank"]
            ].copy()
            additions = []
            for row in missing_origin.itertuples(index=False):
                distances = np.linalg.norm(
                    known[["X", "Y", "Z"]].to_numpy()
                    - np.array([row.X, row.Y, row.Z], dtype=float),
                    axis=1,
                )
                nearest = known.iloc[int(np.argmin(distances))]
                additions.append(
                    {
                        "patient": row.patient,
                        "contact": row.contact,
                        "contact_key": row.contact_key,
                        "X": row.X,
                        "Y": row.Y,
                        "Z": row.Z,
                        "parcel": nearest["parcel"],
                        "parcel_id": nearest["parcel_id"],
                        "Insula_DV": nearest["Insula_DV"],
                        "Julich_insula_region_rank": nearest[
                            "Julich_insula_region_rank"
                        ],
                    }
                )
            data = pd.concat([data, pd.DataFrame(additions)], ignore_index=True)
    return data.drop_duplicates(["patient", "contact_key"], keep="first")


def load_export_coordinates(path: Path) -> pd.DataFrame:
    data = pd.read_excel(path)
    required = {"PAT", "CHAN", "X", "Y", "Z"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing coordinate columns in {path}: {sorted(missing)}")
    data = data.rename(
        columns={"PAT": "patient", "CHAN": "contact"}
    )[["patient", "contact", "X", "Y", "Z"]]
    data["contact_key"] = data["contact"].map(normalize_contact)
    for column in ("X", "Y", "Z"):
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=["X", "Y", "Z"])
    data = data.drop_duplicates(["patient", "contact_key"], keep="first")
    if len(data) != 83:
        raise ValueError(f"Expected 83 export sites in {path}, found {len(data)}.")
    return data


def load_plus1_minus1_t_values(
    output_root: Path, bands: list[str]
) -> pd.DataFrame:
    intermediate_root = output_root / "data" / "intermediate_data"
    files = sorted(intermediate_root.rglob("outcome_*.csv"))
    if not files:
        raise FileNotFoundError(
            f"No intermediate outcome CSVs found under {intermediate_root}"
        )
    values = pd.concat((pd.read_csv(path) for path in files), ignore_index=True)
    required = {"patient", "contact", "outcome", *bands}
    missing = required - set(values.columns)
    if missing:
        raise ValueError(
            f"Missing columns for +1/-1 t-values: {sorted(missing)}"
        )
    values["contact_key"] = values["contact"].map(normalize_contact)
    rows = []
    for (patient, contact_key, contact), site_values in values.groupby(
        ["patient", "contact_key", "contact"], sort=False
    ):
        for band in bands:
            plus_one = site_values.loc[site_values["outcome"] == 2, band].dropna()
            minus_one = site_values.loc[site_values["outcome"] == 1, band].dropna()
            t_value = (
                float(
                    ttest_ind(
                        plus_one,
                        minus_one,
                        equal_var=False,
                        nan_policy="omit",
                    ).statistic
                )
                if len(plus_one) >= 2 and len(minus_one) >= 2
                else np.nan
            )
            rows.append(
                {
                    "patient": patient,
                    "contact_key": contact_key,
                    "contact": contact,
                    "band": band,
                    "t_value": t_value,
                }
            )
    result = pd.DataFrame(rows)
    if result["t_value"].isna().any():
        missing_sites = result[result["t_value"].isna()]
        raise ValueError(
            "Unable to compute +1/-1 t-values for some sites: "
            f"{missing_sites[['patient', 'contact', 'band']].to_dict('records')}"
        )
    return result


def export_site_csvs(
    response: pd.DataFrame,
    t_values: pd.DataFrame,
    origin_path: Path,
    output_dir: Path,
    bands: list[str],
    sign: str,
) -> None:
    coordinates = load_export_coordinates(origin_path)
    category_column = f"{sign}_merged_response_category"
    response_lookup = response.merge(
        t_values,
        on=["patient", "contact_key", "band"],
        how="left",
        validate="one_to_one",
    ).set_index(
        ["patient", "contact_key", "band"]
    ).to_dict("index")
    output_dir.mkdir(parents=True, exist_ok=True)
    for band in bands:
        rows = []
        for site_number, site in enumerate(
            coordinates.sort_values(["patient", "contact_key"]).itertuples(index=False),
            start=1,
        ):
            site_result = response_lookup.get(
                (site.patient, site.contact_key, band), {}
            )
            category = site_result.get(category_column, "neither")
            category = "both" if category == "both_equal" else category
            if category not in EXPORT_CATEGORY_COLORS:
                category = "neither"
            rows.append(
                {
                    "site_id": f"S{site_number:02d}",
                    "label": f"{site.patient}_{site.contact}",
                    "x": site.X,
                    "y": site.Y,
                    "z": site.Z,
                    "color": EXPORT_CATEGORY_COLORS[category],
                    "t_value": site_result.get("t_value", np.nan),
                }
            )
        pd.DataFrame(rows).to_csv(
            output_dir / f"mni_sites_{sign}_{band}_with_t_value.csv",
            index=False,
        )


def add_2d_hull(axis, points, color, alpha=0.12):
    if len(points) < 3:
        return
    hull = ConvexHull(points[:, :2])
    polygon = points[hull.vertices]
    axis.fill(
        polygon[:, 0],
        polygon[:, 1],
        color=color,
        alpha=0.28,
        linewidth=0,
    )


def add_3d_hull(axis, points, color, alpha=0.08):
    if len(points) < 4:
        return
    hull = ConvexHull(points)
    for simplex in hull.simplices:
        triangle = points[simplex]
        axis.plot_trisurf(
            triangle[:, 0],
            triangle[:, 1],
            triangle[:, 2],
            color=color,
            alpha=0.18,
            linewidth=0,
        )


def parcel_background(data: pd.DataFrame):
    parcels = list(data.groupby("parcel", dropna=False))
    colors = {
        parcel: PARCEL_COLORS[index % len(PARCEL_COLORS)]
        for index, (parcel, _) in enumerate(parcels)
    }
    return parcels, colors


def parcel_label(parcel: object) -> str:
    match = re.search(r"\b(Id\d+|Ia\d+)", str(parcel))
    return match.group(1) if match else str(parcel)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bands", nargs="+", required=True)
    parser.add_argument(
        "--sign", choices=("positive", "negative"), default="positive"
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        choices=("dorsal", "ventral", "outside"),
        default=("dorsal", "ventral"),
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        choices=list(CATEGORY_COLORS),
        default=list(CATEGORY_COLORS),
    )
    parser.add_argument(
        "--coordinates",
        type=Path,
        default=Path(path_to_julich_output) / "PBLT_ROIs.csv",
    )
    parser.add_argument(
        "--origin-coordinates",
        type=Path,
        default=Path(path_to_julich_output) / "PBLT_ORIGIN_INSULA.xlsx",
    )
    args = parser.parse_args()

    source = (
        args.output
        / "data"
        / "gain_loss_responses"
        / f"gain_loss_site_responses_{args.sign}_by_band.csv"
    )
    if not source.is_file():
        raise FileNotFoundError(f"Missing response-category CSV: {source}")
    response = pd.read_csv(source)
    category_column = f"{args.sign}_merged_response_category"
    required = {"patient", "contact", "position", "band", category_column}
    missing = required - set(response.columns)
    if missing:
        raise ValueError(f"Missing response columns in {source}: {sorted(missing)}")
    response["contact_key"] = response["contact"].map(normalize_contact)
    t_values = load_plus1_minus1_t_values(args.output, args.bands)
    export_site_csvs(
        response,
        t_values,
        args.origin_coordinates,
        args.output / "data" / "gain_loss_responses" / "mni_site_csv",
        args.bands,
        args.sign,
    )
    coordinates = load_coordinates(args.coordinates, args.origin_coordinates)
    response = response[response["band"].isin(args.bands)].copy()
    response["plot_category"] = response[category_column].where(
        response[category_column].isin(CATEGORY_COLORS), "neither"
    )
    response = response[
        ["patient", "contact_key", "band", "position", "plot_category"]
    ].drop_duplicates(["patient", "contact_key", "band"])
    merged = response.merge(
        coordinates,
        on=["patient", "contact_key"],
        how="inner",
    )
    merged = merged[merged["position"].isin(args.regions)].copy()
    if merged.empty:
        raise ValueError("No MNI sites remain after applying region filters.")
    if len(merged) != len(response):
        missing_sites = response.merge(
            coordinates[["patient", "contact_key"]],
            on=["patient", "contact_key"],
            how="left",
            indicator=True,
        )
        missing_sites = missing_sites[missing_sites["_merge"] == "left_only"]
        raise ValueError(
            "Coordinate data do not cover every barplot site. Missing sites: "
            f"{missing_sites[['patient', 'contact_key']].to_dict('records')}"
        )
    base_dir = args.output / "fig" / "summary" / "mni_sites" / args.sign
    data_base = args.output / "data" / "gain_loss_responses"
    for include_neither in (True, False):
        output_dir = base_dir / (
            "with_neither" if include_neither else "without_neither"
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        plot_data = (
            merged
            if include_neither
            else merged[merged["plot_category"] != "neither"]
        )
        plot_data.to_csv(
            data_base
            / f"mni_site_categories_{args.sign}_"
            f"{'with_neither' if include_neither else 'without_neither'}.csv",
            index=False,
        )
        count_table = (
            plot_data.groupby(["band", "position", "plot_category"])
            .size()
            .rename("n_sites")
            .reset_index()
        )
        count_table.to_csv(
            data_base
            / f"mni_site_category_counts_{args.sign}_"
            f"{'with_neither' if include_neither else 'without_neither'}.csv",
            index=False,
        )
        for band in args.bands:
            band_data = plot_data[plot_data["band"] == band]
            if band_data.empty:
                continue
            parcels, parcel_colors = parcel_background(coordinates)

            figure = plt.figure(figsize=(20, 7))
            axis_3d = figure.add_subplot(131, projection="3d")
            axis_xz = figure.add_subplot(132)
            axis_yz = figure.add_subplot(133)
            for parcel, parcel_data in parcels:
                points = parcel_data[["X", "Y", "Z"]].to_numpy()
                color = parcel_colors[parcel]
                add_3d_hull(axis_3d, points, color)
                add_2d_hull(axis_xz, points[:, [0, 2]], color)
                add_2d_hull(axis_yz, points[:, [1, 2]], color)
                label = f"{parcel_label(parcel)} ({parcel_data['Insula_DV'].iloc[0]})"
                centroid = points.mean(axis=0)
                axis_3d.text(*centroid, label, fontsize=7, color="black")
                axis_xz.text(centroid[0], centroid[2], label, fontsize=7)
                axis_yz.text(centroid[1], centroid[2], label, fontsize=7)

            for category in args.categories:
                sites = band_data[band_data["plot_category"] == category]
                if sites.empty:
                    continue
                color = CATEGORY_COLORS[category]
                axis_3d.scatter(
                    sites["X"],
                    sites["Y"],
                    sites["Z"],
                    s=75,
                    color=color,
                    edgecolor="black",
                    linewidth=0.8,
                    label=f"{category} (n={len(sites)})",
                )
                axis_xz.scatter(
                    sites["X"],
                    sites["Z"],
                    s=75,
                    color=color,
                    edgecolor="black",
                    linewidth=0.8,
                    label=category,
                )
                axis_yz.scatter(
                    sites["Y"],
                    sites["Z"],
                    s=75,
                    color=color,
                    edgecolor="black",
                    linewidth=0.8,
                    label=category,
                )

            for site in band_data.itertuples(index=False):
                site_label = f"{site.patient}:{site.contact} {site.parcel_id}"
                axis_3d.text(
                    site.X,
                    site.Y,
                    site.Z,
                    site_label,
                    fontsize=5,
                    color="black",
                )
                axis_xz.text(site.X, site.Z, site.parcel_id, fontsize=6)
                axis_yz.text(site.Y, site.Z, site.parcel_id, fontsize=6)

            axis_3d.set(xlabel="MNI X", ylabel="MNI Y", zlabel="MNI Z")
            axis_3d.set_title(f"{band} — 3D MNI")
            axis_xz.set(xlabel="MNI X", ylabel="MNI Z", title=f"{band} — XZ projection")
            axis_yz.set(ylabel="MNI Z", xlabel="MNI Y", title=f"{band} — YZ projection")
            for axis in (axis_xz, axis_yz):
                axis.grid(alpha=0.2)
                axis.set_aspect("equal", adjustable="datalim")
            handles, labels = axis_3d.get_legend_handles_labels()
            axis_3d.legend(handles, labels, loc="best", fontsize=8)
            axis_xz.legend(loc="best", fontsize=7)
            axis_yz.legend(loc="best", fontsize=7)
            figure.suptitle(
                f"MNI response sites — {band} — {args.sign} — "
                f"{'with neither' if include_neither else 'without neither'} "
                "(parcel background: Julich insula)"
            )
            figure.tight_layout()
            figure.savefig(output_dir / f"{band}_mni_sites.png", dpi=250)
            plt.close(figure)

    print(f"Matched site data: {args.output / 'data' / 'gain_loss_responses'}")
    print(f"MNI figures: {base_dir.resolve()}")


if __name__ == "__main__":
    main()
