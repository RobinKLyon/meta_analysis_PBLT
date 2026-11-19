"""Shared MNI ROI-coordinate lookup for PBLT analysis.

Matches contacts by (PAT, CHAN), handles French decimal commas and keeps
X/Y/Z uppercase across analysis and intermediate exports.
"""

from pathlib import Path
import re
import unicodedata
import warnings

import numpy as np
import pandas as pd

COORDINATE_COLUMNS = ("X", "Y", "Z")


def normalize_patient(value):
    """Canonical patient key, e.g. GRE_2014_AUDe_PBLT -> GRE_2014_AUDE."""
    if pd.isna(value):
        return ""
    return str(value).strip().upper().removesuffix("_PBLT")


def normalize_contact(value):
    """Canonical electrode name, independent of zero padding.

    Examples:
        X'02 -> X'2 ; X'2 -> X'2
        X02  -> X2   ; X2  -> X2

    The apostrophe remains significant: X2 != X'2.
    """
    if pd.isna(value):
        return ""
    contact = (
        unicodedata.normalize("NFKC", str(value))
        .strip()
        .replace("’", "'")
        .replace("′", "'")
        .replace("＇", "'")
        .upper()
    )
    # An electrode identifier is a letter prefix, an optional prime, and
    # a numeric index.  Convert only the numeric index to an integer so
    # leading zeroes no longer prevent (PAT, CHAN) matching.
    match = re.fullmatch(r"([A-Z]+)('?)0*(\d+)", contact)
    if match:
        prefix, prime, index = match.groups()
        return f"{prefix}{prime}{int(index)}"
    # Preserve unknown label forms rather than silently conflating them.
    return contact


def parse_mni_coordinate(value):
    """Parse a CSV value such as '-29,877' into a real Python float."""
    if pd.isna(value):
        return np.nan
    value = str(value).strip().replace("\u00a0", "").replace(" ", "")
    if not value:
        return np.nan
    return pd.to_numeric(value.replace(",", "."), errors="coerce")


def _canonical_column_name(value):
    return "".join(char.casefold() for char in str(value) if char.isalnum())


def build_roi_coordinate_lookup(roi_csv):
    """Return {(NORMALIZED_PAT, NORMALIZED_CHAN): {'X': float, ...}}."""
    roi_csv = Path(roi_csv)
    rois = pd.read_csv(roi_csv, sep=";", dtype=str, encoding="utf-8-sig")
    col_names = {_canonical_column_name(col): col for col in rois.columns}

    def require_col(*candidates):
        for candidate in candidates:
            col = col_names.get(_canonical_column_name(candidate))
            if col is not None:
                return col
        raise ValueError(
            f"ROI CSV {roi_csv}: missing column from {candidates}; "
            f"columns available: {list(rois.columns)}"
        )

    patient_col = require_col("PAT")
    contact_col = require_col("CHAN")
    cols = {
        axis: require_col(axis, f"MNI_{axis}", f"{axis}_MNI")
        for axis in COORDINATE_COLUMNS
    }

    lookup = {}
    for _, row in rois.iterrows():
        key = (normalize_patient(row[patient_col]), normalize_contact(row[contact_col]))
        if not all(key):
            continue
        coordinates = {
            axis: parse_mni_coordinate(row[column])
            for axis, column in cols.items()
        }
        if key not in lookup:
            lookup[key] = coordinates
            continue
        # Repeated rows are permitted only if their existing values agree.
        for axis in COORDINATE_COLUMNS:
            old, new = lookup[key][axis], coordinates[axis]
            if pd.isna(old):
                lookup[key][axis] = new
            elif pd.notna(new) and not np.isclose(old, new, rtol=0, atol=1e-5):
                raise ValueError(
                    f"Conflicting {axis} coordinates in {roi_csv} "
                    f"for {key}: {old} versus {new}"
                )

    n_complete = sum(
        all(pd.notna(coords[axis]) for axis in COORDINATE_COLUMNS)
        for coords in lookup.values()
    )
    print(
        f"[ROI MNI] {roi_csv.resolve()} | "
        f"{n_complete}/{len(lookup)} distinct contacts with complete X/Y/Z"
    )
    return lookup


def add_roi_coordinates_to_frames(integrated_data, roi_csv):
    """Attach correct X/Y/Z to each contact DataFrame, preserving existing values."""
    lookup = build_roi_coordinate_lookup(roi_csv)
    enriched = []
    missing = []
    for original in integrated_data:
        frame = original.copy()
        if frame.empty:
            enriched.append(frame)
            continue
        patient = frame["patient"].iloc[0]
        contact = frame["contact"].iloc[0]
        key = (normalize_patient(patient), normalize_contact(contact))
        coords = lookup.get(key)
        for axis in COORDINATE_COLUMNS:
            existing = (
                frame[axis].map(parse_mni_coordinate)
                if axis in frame.columns
                else pd.Series(np.nan, index=frame.index, dtype=float)
            )
            if coords is not None and pd.notna(coords[axis]):
                frame[axis] = float(coords[axis])
            else:
                frame[axis] = existing

        if not frame[list(COORDINATE_COLUMNS)].notna().all().all():
            missing.append(key)
        enriched.append(frame)

    if missing:
        unique = sorted(set(missing))
        preview = ", ".join(f"{pat}/{chan}" for pat, chan in unique[:10])
        suffix = " ..." if len(unique) > 10 else ""
        warnings.warn(
            f"ROI coordinates missing or incomplete for {len(unique)} "
            f"contact(s): {preview}{suffix} (source: {Path(roi_csv).resolve()})",
            RuntimeWarning,
        )
    else:
        print(f"[OK] X/Y/Z coordinates complete for {len(enriched)} contact DataFrames")
    return enriched
