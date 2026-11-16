
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from analysis_config import ALPHA_DEFAULT, BANDS

# 1. Intégrer chaque TFR dans les bandes fréquentielles
def integrate_frequency(data, freqs, bands=BANDS):
    # data : trials × times × freqs × channels
    data = np.asarray(data)
    result = {}
    for band, (fmin, fmax) in bands.items():
        mask = (freqs >= fmin) & (freqs <= fmax)
        if mask.sum() < 2:
            raise ValueError(f"Fréquences insuffisantes pour {band}")

        # Intégration fréquentielle : trials × times × channels
        try:
            result[band] = np.trapz(
                data[:, :, mask, :],
                x=np.asarray(freqs)[mask],
                axis=2,
            )/(fmax - fmin)  # Normalisation par la largeur de bande
        except Exception as e:
            result[band] = np.trapezoid(
                data[:, :, mask, :],
                x=np.asarray(freqs)[mask],
                axis=2,
            )/(fmax - fmin)  # Normalisation par la largeur de bande

    return result

# 3. Intégrer dans le temps et sélectionner/agréger les contacts
def integrate_time(band_signals, times, outcomes, tmin=0, tmax=0.8,
                   channel=0):
    """
    Intègre les signaux dans le temps.
    Retourne une valeur par essai et par bande.
    channel=0 : contact choisi.
    channel=None : moyenne des contacts.
    """
    times = np.asarray(times)
    outcomes = np.asarray(outcomes)
    mask = (times >= tmin) & (times <= tmax)

    if mask.sum() < 2:
        raise ValueError("Il faut au moins deux points temporels dans la fenêtre.")

    values = {}
    for band, signal in band_signals.items():
        # signal : trials × times × channels
        selected = signal[:, mask, :]

        # Intégrale temporelle : trials × channels
        try:
            integral = np.trapz(selected, x=times[mask], axis=1)
        except Exception as e:
            integral = np.trapezoid(selected, x=times[mask], axis=1)
    
        # Une valeur par essai
        values[band] = (
            integral[:, channel] if channel is not None
            else np.nanmean(integral, axis=1)
        )

    return pd.DataFrame(values).assign(outcome=outcomes)


def baseline_mask(times, t_min_baseline, t_max_baseline):
    return (times < t_max_baseline) & (times >= t_min_baseline)
def baseline_mean(data, times, t_min_baseline, t_max_baseline):
    baseline_times = baseline_mask(times, t_min_baseline, t_max_baseline)
    return np.mean(data[:, baseline_times, :, :], axis=1)
def baseline_std(data, times, t_min_baseline, t_max_baseline):
    baseline_times = baseline_mask(times, t_min_baseline, t_max_baseline)
    return np.std(data[:, baseline_times, :, :], axis=1)
def baseline_zscore(data, times, t_min_baseline, t_max_baseline):
    baseline_mean_values = baseline_mean(data, times, t_min_baseline, t_max_baseline)
    baseline_std_values = baseline_std(data, times, t_min_baseline, t_max_baseline)
    return (data - baseline_mean_values[:, np.newaxis, :, :]) / baseline_std_values[:, np.newaxis, :, :]
def baseline_div(data, times, t_min_baseline, t_max_baseline):
    baseline_mean_values = baseline_mean(data, times, t_min_baseline, t_max_baseline)
    return data / baseline_mean_values[:, np.newaxis, :, :]
def log_baseline_div(data, times, t_min_baseline, t_max_baseline):
    baseline_mean_values = baseline_mean(data, times, t_min_baseline, t_max_baseline)
    return np.log(data / baseline_mean_values[:, np.newaxis, :, :])



def get_baselined_data(data, times, t_baseline_min, t_baseline_max, baseline_type='div'):
    """
    Apply baseline correction to the data using the specified method.
    """
    if baseline_type == 'div':
        return baseline_div(data, times, t_baseline_min, t_baseline_max)
    elif baseline_type == 'zscore':
        return baseline_zscore(data, times, t_baseline_min, t_baseline_max)
    elif baseline_type == 'log_div':
        return log_baseline_div(data, times, t_baseline_min, t_baseline_max)
    elif baseline_type == 'raw':
        return data
    elif baseline_type == 'diff':
        baseline_mean_values = baseline_mean(data, times, t_baseline_min, t_baseline_max)
        return data - baseline_mean_values[:, np.newaxis, :, :]
    else:
        raise ValueError(f"Unknown baseline_type: {baseline_type}")

def compare_outcomes(df, patient, contact_name, contact_idx, position,
                     output_dir, tmin, tmax, alpha=ALPHA_DEFAULT):
    """
    Compare outcome 1 (-1 €) à outcome 2 (+1 €) pour chaque bande.
    Une réponse positive signifie outcome 2 > outcome 1 ; négative signifie
    outcome 2 < outcome 1. Significance is based on FDR-adjusted p-values.
    """
    bands = [b for b in BANDS if b in df.columns]
    results = []
    for band in bands:
        x = df.loc[df["outcome"] == 1, band].dropna().to_numpy()  # outcome 1
        y = df.loc[df["outcome"] == 2, band].dropna().to_numpy()  # outcome 2
        p = mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) and len(y) else np.nan
        results.append({
            "patient": patient,
            "position": position,
            "contact_index": int(contact_idx),
            "contact": str(contact_name),
            "band": band,
            "n_outcome_1": len(x),
            "n_outcome_2": len(y),
            "median_outcome_1": float(np.median(x)) if len(x) else np.nan,
            "median_outcome_2": float(np.median(y)) if len(y) else np.nan,
            "median_difference_2_minus_1": (
                float(np.median(y) - np.median(x)) if len(x) and len(y) else np.nan
            ),
            "p_raw": p,
        })

    valid_indices = [i for i, r in enumerate(results) if np.isfinite(r["p_raw"])]
    if valid_indices:
        corrected = multipletests(
            [results[i]["p_raw"] for i in valid_indices], method="fdr_bh"
        )[1]
        for i, p_corr in zip(valid_indices, corrected):
            results[i]["p_fdr"] = float(p_corr)
    for r in results:
        r.setdefault("p_fdr", np.nan)
        r["significant"] = bool(np.isfinite(r["p_fdr"]) and r["p_fdr"] < alpha)
        diff = r["median_difference_2_minus_1"]
        r["response"] = (
            "positive" if r["significant"] and diff > 0 else
            "negative" if r["significant"] and diff < 0 else
            "non_significant"
        )

    
    result_df = pd.DataFrame(results)
    return result_df


def compare_gain_loss_outcomes(
    df, patient, contact_name, contact_idx, position, alpha=ALPHA_DEFAULT
):
    """Compare gain (+1 - 0) and loss (-1 - 0) outcomes by band."""
    comparisons = {
        "gain": (2, 0, "gain_plus1_minus_zero"),
        "loss": (1, 3, "loss_minus1_minus_zero"),
    }
    results = []
    for band in [b for b in BANDS if b in df.columns]:
        row = {
            "patient": patient,
            "position": position,
            "contact_index": int(contact_idx),
            "contact": str(contact_name),
            "band": band,
        }
        for name, (higher_outcome, lower_outcome, prefix) in comparisons.items():
            higher = df.loc[df["outcome"] == higher_outcome, band].dropna().to_numpy()
            lower = df.loc[df["outcome"] == lower_outcome, band].dropna().to_numpy()
            p_value = (
                mannwhitneyu(higher, lower, alternative="two-sided").pvalue
                if len(higher) and len(lower)
                else np.nan
            )
            row.update({
                f"n_{prefix}": len(higher),
                f"n_{prefix}_reference": len(lower),
                f"median_{prefix}": (
                    float(np.median(higher) - np.median(lower))
                    if len(higher) and len(lower)
                    else np.nan
                ),
                f"p_raw_{name}": p_value,
            })
        results.append(row)

    for comparison in comparisons:
        p_column = f"p_raw_{comparison}"
        corrected_column = f"p_fdr_{comparison}"
        valid = [i for i, row in enumerate(results) if np.isfinite(row[p_column])]
        corrected = (
            multipletests([results[i][p_column] for i in valid], method="fdr_bh")[1]
            if valid
            else []
        )
        for i, p_value in zip(valid, corrected):
            results[i][corrected_column] = float(p_value)
        for row in results:
            row.setdefault(corrected_column, np.nan)
            row[f"significant_{comparison}"] = bool(
                np.isfinite(row[corrected_column]) and row[corrected_column] < alpha
            )
    return pd.DataFrame(results)