import pandas as pd
from pathlib import Path
import numpy as np

from load_tfr_data import get_sampling_rate
# Codes des stimuli (Picture)
PICTURE_CODES = {
    12: ("gain", "up", 1),
    17: ("gain", "up", 1),
    22: ("gain", "down", 2),
    27: ("gain", "down", 2),
    32: ("loss", "up", 1),
    37: ("loss", "up", 1),
    42: ("loss", "down", 2),
    47: ("loss", "down", 2),
}

# Codes de feedback
GAIN_CORRECT = {13, 18, 23, 28}
GAIN_INCORRECT = {14, 19, 24, 29}

LOSS_INCORRECT = {33, 38, 43, 48}
LOSS_CORRECT = {34, 39, 44, 49}

FEEDBACK_CODES = (
    GAIN_CORRECT | GAIN_INCORRECT |
    LOSS_INCORRECT | LOSS_CORRECT
)

# Valurs finals des 4 catégories :
# 0 = échec pendant un essai de gain
# 1 = échec pendant un essai de perte
# 2 = réussite pendant un essai de gain
# 3 = réussite pendant un essai de perte
OUTCOME_LABELS = {
    0: "gain_failure_0_euro",
    1: "loss_failure_minus_1_euro",
    2: "gain_success_plus_1_euro",
    3: "loss_success_0_euro",
}

OUTCOME_REWARDS = {
    0: 0,
    1: -1,
    2: 1,
    3: 0,
}
def classify_outcome(row):
    """Retourne le code de catégorie 0, 1, 2 ou 3."""
    fb = row["feedback_code"]
    condition = row["condition"]

    if pd.isna(fb):
        return np.nan

    fb = int(fb)

    if condition == "gain":
        if fb in GAIN_INCORRECT:
            return 0  # Échec gain : 0 €
        if fb in GAIN_CORRECT:
            return 2  # Réussite gain : +1 €

    elif condition == "loss":
        if fb in LOSS_INCORRECT:
            return 1  # Échec perte : -1 €
        if fb in LOSS_CORRECT:
            return 3  # Réussite perte : 0 €

    return np.nan  # Code inattendu ou incohérent



def prepare_outcomes(events_df, n_tfr):
    """Reconstruit les essais et aligne les feedbacks chronologiques sur les TFR."""
    trials = reconstruct_trials(events_df).copy()
    trials["outcome_category"] = trials.apply(classify_outcome, axis=1)
    trials["outcome_label"] = trials["outcome_category"].map(OUTCOME_LABELS)
    trials["reward_euro"] = trials["outcome_category"].map(OUTCOME_REWARDS)
    trials["feedback_says_correct"] = trials.apply(expected_feedback_status, axis=1)
    trials["response_feedback_agree"] = (
        trials["response_matches_picture"] == trials["feedback_says_correct"]
    )

    feedback_trials = (
        trials.dropna(subset=["feedback_code"])
        .sort_values("feedback_time")
        .reset_index(drop=True)
        .copy()
    )
    if len(feedback_trials) != n_tfr:
        raise ValueError(
            f"Alignement impossible : {n_tfr} TFR, mais "
            f"{len(feedback_trials)} essais avec feedback. "
            "Vérifier les événements .pos avant d'associer les essais."
        )
    feedback_trials.insert(0, "tfr_trial_index", np.arange(n_tfr))
    outcome_array = (
        feedback_trials["outcome_category"].fillna(-1).astype(int).to_numpy()
    )
    return trials, feedback_trials, outcome_array

def patient_names_from_root(root):
    """Déduit les patients des sous-dossiers contenant un fichier PBLT .pos."""
    root = Path(root)
    found = set()
    if not root.exists():
        raise FileNotFoundError(f"Répertoire PBLT introuvable : {root}")
    for pos in root.glob("*/*.pos"):
        name = pos.stem
        if name.endswith("_PBLT"):
            name = name[:-5]
        if name:
            found.add(name)
    return sorted(found)

def expected_feedback_status(row):
    if pd.isna(row["feedback_code"]):
        return None

    fb = int(row["feedback_code"])

    if row["condition"] == "gain":
        if fb in GAIN_CORRECT:
            return True
        if fb in GAIN_INCORRECT:
            return False

    if row["condition"] == "loss":
        if fb in LOSS_CORRECT:
            return True
        if fb in LOSS_INCORRECT:
            return False

    return None
def reconstruct_trials(events_df):    
    # Ordre chronologique
    events_df = events_df.sort_values("time").reset_index(drop=True)

    picture_idx = np.flatnonzero(events_df["code"].isin(PICTURE_CODES))
    feedback_idx = np.flatnonzero(events_df["code"].isin(FEEDBACK_CODES))

    trials = []

    # Un trial par feedback : le feedback est la référence
    for f_idx in feedback_idx:
        feedback = events_df.iloc[f_idx]

        # Chercher le Picture précédent le plus proche
        previous_pictures = picture_idx[picture_idx < f_idx]

        if len(previous_pictures):
            p_idx = previous_pictures[-1]
            picture = events_df.iloc[p_idx]
            picture_time = picture["time"]
            picture_code = int(picture["code"])

            # Appuis entre le Picture et le feedback
            between = events_df.iloc[p_idx + 1:f_idx]
            button_events = between[between["code"].isin([1, 2])]

            response_codes = button_events["code"].astype(int).tolist()
            response = response_codes[-1] if response_codes else np.nan
            response_time = (
                button_events.iloc[-1]["time"]
                if len(button_events) else np.nan
            )
            reaction_time = (
                response_time - picture_time
                if not np.isnan(response_time) else np.nan
            )

            condition, direction, expected_response = PICTURE_CODES.get(
                picture_code, (None, None, None)
            )

            response_matches_picture = (
                response == expected_response
                if pd.notna(response) and expected_response is not None
                else np.nan
            ) # verify that the response matches the expected response (patient is correct or not)

        else:
            # Aucun Picture antérieur
            p_idx = None
            picture_time = np.nan
            picture_code = np.nan
            response_codes = []
            response = np.nan
            response_time = np.nan
            reaction_time = np.nan
            condition = None
            direction = None
            expected_response = np.nan
            response_matches_picture = np.nan

        trials.append({
            "picture_time": picture_time,
            "picture_code": picture_code,
            "feedback_time": feedback["time"],
            "feedback_code": int(feedback["code"]),
            "response": response,
            "reaction_time": reaction_time,
            "n_button_presses": len(response_codes),
            "response_codes": response_codes,
            "response_time": response_time,
            "condition": condition,
            "direction": direction,
            "expected_response": expected_response,
            "response_matches_picture": response_matches_picture,
        })

    trials_df = pd.DataFrame(trials)



    return trials_df
def normalize_contact_name(name):
    """Uniformise les apostrophes sans fusionner les contacts distincts."""
    return str(name).strip().replace("’", "'").replace("‘", "'")

def eliminate_outside_and_doublons_sites(classification):
    classification = classification[
        classification["classification"].isin(["Ventral", "Dorsal"])
    ].copy()

    # Éviter d'analyser plusieurs fois le même contact
    # sans fusionner X02 et X'02.
    classification["contact_key"] = classification["contact"].map(
        normalize_contact_name
    )


    classification = classification.drop_duplicates(
        subset=["contact_key", "classification"],
        keep="first"
    )



    return classification
def get_pos_events(patient, root):
    """Charge les événements comportementaux du patient."""
    pos_path = get_pos_file_path(patient, root)
    config_path = Path(root) / "config.json"
    sampling_rate = get_sampling_rate(config_path, patient)
    return load_pos_file(pos_path, sampling_rate=sampling_rate)

def load_pos_file(pos_file_path, sampling_rate=None):
    """Charge les événements .pos : colonne 0 = temps, colonne 1 = code."""
    if not pos_file_path.exists():
        raise FileNotFoundError(f"Fichier introuvable : {pos_file_path}")

    # Le fichier possède une ligne d'en-tête, comme dans la cellule existante.
    events = np.loadtxt(pos_file_path)

    if events.ndim == 1:
        events = events.reshape(1, -1)

    if events.shape[1] < 2:
        raise ValueError("Le fichier .pos doit contenir au moins 2 colonnes.")

    return pd.DataFrame({
        "time": events[:, 0]/sampling_rate if sampling_rate else events[:, 0],
        "code": events[:, 1].astype(int),
        "event_row": np.arange(len(events))
    })
def get_pos_file_path(patient, path_to_PBLT):
    """
    Retourne le chemin du fichier .pos pour un patient donné.
    """
    patient = patient.removesuffix("_PBLT")  # Supprime le suffixe "_PBLT" si présent
    pos_file_path = Path(path_to_PBLT) / patient / f"{patient}_PBLT.pos"
    return pos_file_path