"""Core patient-level PBLT processing."""

import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from analysis_config import BASELINE_MAX
from behavorial_organize_data import (
    get_pos_events,
    prepare_outcomes,
    eliminate_outside_and_doublons_sites,
)
from electrodes_lebelling_masks import (
    get_chanlist_from_info,
    classify_contacts_insula,
    bipolar_to_monopolar,
    monopolar_channels_from_bipolar,
    augment_roi_with_origin_insula,
)
from load_tfr_data import load_patient, get_hdr_info, load_tfr_data_wya
from signal_processing import integrate_time, integrate_frequency, get_baselined_data


def run_patient(
    patient, root, roi_csv, tmin, tmax, channels_verbose, baseline_type,
    take_outside=False
):
    """Load, baseline, integrate and classify one patient's contacts."""
    print(f"\n=== Patient: {patient} ===")
    header, tfr_path = load_patient(patient, root)
    info = get_hdr_info(header)
    n_trial, n_chan, n_freqs, n_time, times, freqs = (
        info[key]
        for key in ["N_trial", "N_chan", "N_freqs", "N_time", "times", "freqs"]
    )
    data = np.asarray(
        load_tfr_data_wya(tfr_path, n_trial, n_chan, n_freqs, n_time),
        dtype=np.float32,
    )
    if data.ndim != 4:
        raise ValueError(f"TFR is not 4D for {patient}: {data.shape}")
    if data.shape[0] != n_trial:
        raise ValueError(
            f"Declared N_trial ({n_trial}) != data.shape[0] ({data.shape[0]})"
        )

    times = np.asarray(times)
    data_baseline = get_baselined_data(
        data,
        times,
        float(times[0]),
        BASELINE_MAX,
        baseline_type=baseline_type,
    )
    events_df = get_pos_events(patient, root)
    _, _, outcomes = prepare_outcomes(events_df, data.shape[0])
    band_signals = integrate_frequency(data_baseline, freqs)

    bipolar_chanlist = get_chanlist_from_info(info)
    rois = pd.read_csv(roi_csv, delimiter=";")
    rois = augment_roi_with_origin_insula(
        rois,
        Path(roi_csv).with_name("PBLT_ORIGIN_INSULA.xlsx"),
    )
    if channels_verbose:
        print(f"Bipolar contacts: {bipolar_chanlist}")
        patient_channels = rois[rois["PAT"] == patient]["CHAN"].values
        print(f"Monopolar channels: {patient_channels}", len(patient_channels),
              "insula channels for", patient)
    monopolar_chanlist = monopolar_channels_from_bipolar(bipolar_chanlist)
    first_contacts, second_contacts = bipolar_to_monopolar(
        bipolar_chanlist, monopolar_chanlist
    )
    classification = classify_contacts_insula(
        ROIS=rois,
        first_contacts=first_contacts,
        patient=patient.removesuffix("_PBLT"),
    )
    if not take_outside:
        classification = eliminate_outside_and_doublons_sites(classification)

    integrated_data = []
    position_specs = [("Ventral", "ventral"), ("Dorsal", "dorsal")]
    if take_outside:
        position_specs.append(("Outside", "outside"))
    for position, folder in position_specs:
        group = classification[classification["classification"] == position]
        for ch_idx, ch_name in group["contact"].items():
            if not isinstance(ch_idx, (int, np.integer)) or ch_idx < 0 or ch_idx >= data.shape[-1]:
                warnings.warn(f"{patient}: invalid contact index {ch_idx}; skipped.")
                continue
            values_df = integrate_time(
                band_signals, times, outcomes, tmin=tmin, tmax=tmax,
                channel=int(ch_idx)
            )
            integrated_data.append(values_df.assign(
                trial_index=values_df.index,
                patient=patient,
                position=folder,
                contact_index=int(ch_idx),
                contact=str(ch_name),
            ))
    return integrated_data
