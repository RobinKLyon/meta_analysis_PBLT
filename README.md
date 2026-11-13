# PBLT Meta-Analysis

Analysis scripts for PBLT TFR responses by frequency band, insular region,
response category, and MNI site.

## Contents

- `analyse_pblt_tous_patients.py`: first-stage statistical analysis.
- `analyse_regression_mixte_intermediate.py`: configurable LMM analysis.
- `analyse_reponses_temporelles.py`: temporal response-category analysis.
- `plot_sites_mni_3d.py`: MNI site and parcel visualizations.
- `plot_tfr_power_by_response_category.py`: TFR plots by response category.
- `src/`: reusable loading, preprocessing, statistics, and plotting modules.
- `README_SCRIPTS.md`: detailed script documentation.
- `README_LMM.md`: LMM documentation.

## Data and configuration

Raw TFR data, atlas spreadsheets, and generated outputs are intentionally not
included in this repository. Configure local paths in `src/paths_data.py` and
keep generated results under `outputs/`.

## Example

Run the MNI visualization from the project directory:

```powershell
python .\plot_sites_mni_3d.py `
    --output .\outputs\baseline_raw_YYYYMMDD `
    --bands Theta Alpha Beta Gamma HFA `
    --sign positive
```

See `README_SCRIPTS.md` for the complete commands and expected input files.
