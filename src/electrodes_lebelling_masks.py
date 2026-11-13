import pandas as pd
import re
import numpy as np
from pathlib import Path


def classify_contacts_insula(ROIS, first_contacts, patient):
    """
    Retourne la classification anatomique des contacts :
    Dorsal, Ventral ou Outside.

    ROIS doit contenir PAT, CHAN, Insula_DV.
    first_contacts : liste des contacts monopolaire.
    patient : identifiant du patient dans ROIS['PAT'].
    """

    def normalize(name):
        # Uniformise les apostrophes, mais ne les supprime pas.
        name = (
            str(name).strip().upper()
            .replace("’", "'")
            .replace("‘", "'")
            .replace(" ", "")
        )

        # Préserve la distinction entre X02 et X'02.
        match = re.fullmatch(r"([A-Z]+['’‘]?)0*(\d+)", name)
        if match:
            prefix, number = match.groups()
            return f"{prefix}{int(number)}"

        return name

    # Sélection du patient
    rois_patient = ROIS[ROIS["PAT"] == patient].copy()
    # Correspondance contact normalisé -> classification
    lookup = {}

    for _, row in rois_patient.iterrows():
        key = normalize(row["CHAN"])
        value = row.get("Insula_DV", None)

        if pd.isna(value) or str(value).strip() == "":
            value = "Outside"

        # Si un contact est présent plusieurs fois dans le CSV ROI,
        # on conserve sa première classification.
        lookup.setdefault(key, value)

    # Une classification pour chaque contact, dans le même ordre
    results = []

    for contact in first_contacts:
        key = normalize(contact)
        classification = lookup.get(key, "Outside")
        results.append({
            "contact": contact,
            "classification": classification
        })

    return pd.DataFrame(results)


def augment_roi_with_origin_insula(rois, origin_path):
    """Add source-insula contacts missing from the Julich ROI CSV.

    Contacts listed in PBLT_ORIGIN_INSULA are authoritative insular contacts.
    A missing dorsal/ventral label inherits the label of the nearest known
    insular contact from the same patient in MNI coordinate space.
    """
    origin_path = Path(origin_path)
    if not origin_path.is_file():
        return rois
    origin = pd.read_excel(origin_path)
    required = {"PAT", "CHAN", "X", "Y", "Z"}
    if not required.issubset(origin.columns):
        raise ValueError(
            f"Missing columns in origin insula file: {sorted(required - set(origin.columns))}"
        )
    rois = rois.copy()
    if "Insula_DV" not in rois.columns:
        rois["Insula_DV"] = np.nan

    def key(frame):
        return list(zip(frame["PAT"].astype(str), frame["CHAN"].astype(str)))

    existing_keys = set(key(rois))
    missing = origin[
        ~origin.apply(
            lambda row: (str(row["PAT"]), str(row["CHAN"])) in existing_keys,
            axis=1,
        )
    ]
    if not missing.empty:
        additions = pd.DataFrame(index=missing.index, columns=rois.columns)
        for column in rois.columns:
            if column in missing.columns:
                additions[column] = missing[column].values
        rois = pd.concat([rois, additions], ignore_index=True)

    origin_keys = set(key(origin))
    for index, row in rois.iterrows():
        contact_key = (str(row["PAT"]), str(row["CHAN"]))
        if contact_key not in origin_keys:
            continue
        if str(row.get("Insula_DV", "")).strip() not in {"", "nan", "None"}:
            continue
        same_patient = rois[
            (rois["PAT"].astype(str) == str(row["PAT"]))
            & rois["Insula_DV"].isin(["Dorsal", "Ventral"])
        ].copy()
        if same_patient.empty:
            continue
        coordinates = pd.to_numeric(
            row[["X", "Y", "Z"]].astype(str).str.replace(",", ".", regex=False),
            errors="coerce",
        )
        known_coordinates = same_patient[["X", "Y", "Z"]].apply(
            lambda column: pd.to_numeric(
                column.astype(str).str.replace(",", ".", regex=False),
                errors="coerce",
            )
        )
        valid_known = known_coordinates.notna().all(axis=1)
        if not bool(coordinates.notna().all()) or not bool(valid_known.any()):
            continue
        same_patient = same_patient.loc[valid_known]
        known_coordinates = known_coordinates.loc[valid_known]
        distances = np.linalg.norm(
            known_coordinates.to_numpy() - coordinates.to_numpy(dtype=float),
            axis=1,
        )
        nearest = same_patient.iloc[int(np.argmin(distances))]
        rois.loc[index, "Insula_DV"] = nearest["Insula_DV"]
    return rois

import re


def monopolar_channels_from_bipolar(chanlist):
    """Extract the unique monopolar endpoints represented in bipolar channels."""
    contact_pattern = r"[A-Z]+['’‘]?0*\d+"
    endpoints = []
    for channel in chanlist:
        matches = re.findall(contact_pattern, str(channel).strip().upper().replace(" ", ""))
        if len(matches) != 2:
            raise ValueError(f"Invalid bipolar channel name: {channel}")
        endpoints.extend(matches)
    return np.unique(endpoints)


def get_chanlist_from_info(info):
    chanlist = info["channels"][0][0]
    return [str(c[0][0]) for c in chanlist]


def bipolar_to_monopolar(chanlist, monopolar_list):
    """
    Parameters
    ----------
    chanlist : list[str]
        Example: ['E02E01', "A'02A'01"]
    monopolar_list : iterable[str]
        Monopolar contact names, for example ROIS['CHAN'].values

    Returns
    -------
    first_contacts, second_contacts : np.ndarray
        Monopolar contacts corresponding to both endpoints.
    """
    

    def normalize(name):
        # Preserve the apostrophe: X02 and X'02 are distinct.
        name = (
            str(name).strip().upper()
            .replace("’", "'")
            .replace("‘", "'")
            .replace(" ", "")
        )

        match = re.fullmatch(r"([A-Z]+')?0*(\d+)", name)
        if not match:
            return name

        prefix = match.group(1) or ""
        number = match.group(2)
        return f"{prefix}{int(number)}"


    first_contacts, second_contacts = [], []

    # A contact can be X02, X'02, or A'02.
    contact_pattern = r"[A-Z]+['’‘]?0*\d+"
    bipolar_pattern = rf"({contact_pattern})({contact_pattern})"

    for ch in chanlist:
        
        ch_clean = str(ch).strip().upper().replace(" ", "")
        match = re.fullmatch(bipolar_pattern, ch_clean)

        if not match:
            raise ValueError(f"Unrecognized bipolar channel name: {ch}")

        c1, c2 = match.groups()
        key1, key2 = normalize(c1), normalize(c2)


        first_contacts.append(key1)
        second_contacts.append(key2)

    return np.array(first_contacts), np.array(second_contacts)
