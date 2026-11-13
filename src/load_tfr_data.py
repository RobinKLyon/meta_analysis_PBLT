
import json
from pathlib import Path
from scipy.io import loadmat
import numpy as np
def get_sampling_rate(json_path, patient):
    """Retourne la fréquence d'échantillonnage du patient en Hz."""
    patient_no_pblt = patient.removesuffix("_PBLT")  # Supprime le suffixe "_PBLT" si présent
    with open(Path(json_path), "r", encoding="utf-8") as f:
        config = json.load(f)
    return config["patients"][patient_no_pblt]["sampling_frequency"]["sampling_frequency_hz"]
def load_patient(patient, dest_dir):
    dest_dir = Path(dest_dir)

    mat_files = list(dest_dir.rglob(f"{patient}*TFR.mat"))
    wya_files = list(dest_dir.rglob(f"{patient}*TFR.wya"))

    if not mat_files or not wya_files:
        raise FileNotFoundError(f"Fichiers TFR manquants pour {patient}")

    header = loadmat(mat_files[0])
    tfr_path = wya_files[0]

 

    return header, tfr_path

def get_hdr_info(header):
    hdr = header["hdr"]
    hdr_tfr = hdr[0, 0]
    times = hdr_tfr["timelist"].ravel()
    freqs = hdr_tfr["freqlist"].ravel()
    dimsize = hdr["dimsiz"][0,0][0] #trials , temps , fréquences, channels
    # Reconstruire le tableau : essais × temps × fréquences × canaux
    N_trial =dimsize[0]
    N_chan = dimsize[3]
    N_freqs= dimsize[2]
    N_time = dimsize[1]
    assert len(times) == N_time, "La longueur de times ne correspond pas à N_time"
    assert len(freqs) == N_freqs, "La longueur de freqs ne correspond pas à N_freqs"
    return { "N_trial": N_trial, "N_chan": N_chan, "N_freqs": N_freqs, "N_time": N_time, "dimsize": dimsize, "times": times, "freqs": freqs , "channels": header["hdr"]["chanlist"]}   
  
def load_tfr_data_wya(wya_path, N_trial, N_chan, N_freqs, N_time):
    # Charger les données binaires
    data = np.fromfile(wya_path, dtype=np.float32)
    # Reconstruire le tableau : essais × temps × fréquences × canaux
    data = data.reshape(N_trial, N_time, N_freqs, N_chan, order="F")
    return data



