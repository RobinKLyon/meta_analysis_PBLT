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
- `configs/`: user-editable YAML configurations for the pipeline.
- `scripts/`: independent stage entry points and `run_all.py`.

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
See `scripts/README.md` for the refactored stage-by-stage workflow.

## Pipeline refactoring

The project is being organized as independent pipeline stages:

1. generate reusable intermediate data from epoch-level TFR files;
2. run configurable site-level event comparisons;
3. generate response barplots from comparison CSV files;
4. fit regional mixed models, separately for each frequency band;
5. generate mixed-model figures;
6. run the temporal response-category bonus analysis;
7. generate TFR/t-value category figures;
8. export CSV files for the Insula Plotter.

The YAML files in `configs/` are the starting point for this organization.
Each completed stage should preserve its effective parameters in the output
`config.json`. Existing scripts remain available while their responsibilities
are being separated, so previous commands remain usable during the transition.
