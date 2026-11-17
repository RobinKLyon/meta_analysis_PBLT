 PBLT Meta-Analysis

Analysis pipeline for investigating time-frequency representations (TFR) of intracranial EEG responses during reward- and loss-based decision-making tasks, with a particular focus on the human insula.

The pipeline provides statistical analyses and visualizations of PBLT responses according to:

- Frequency bands
- Insular anatomical regions
- Response categories
- MNI electrode locations

## 1. Experimental task

The experimental paradigm consists of two types of decision-making games:

**Reward game**

- Reward onset: `+1`
- No reward onset: `0`

**Loss game**

- No loss onset: `0`
- Loss onset: `-1`

These events define the different outcome conditions used to investigate neural responses.

**Current analysis scope:** The input dataset currently contains time-frequency responses associated with reward- and loss game events onsets.

## 2. Input data organization

The analysis requires two types of input data:

1. Patient-level electrophysiological data.
2. Anatomical information about the insular recording sites (`PBLT_ROIs.csv`).

### 2.1. Patient-level data

Patient data must be organized in a common directory containing one subfolder per patient.

Each patient folder contains three files sharing the same patient identifier:

```text
data/
│
├── GRE_2014_AUDe/
│   ├── GRE_2014_AUDe_PBLT.pos
│   ├── GRE_2014_AUDe_PBLT_TFR.mat
│   └── GRE_2014_AUDe_PBLT_TFR.wya
│
├── GRE_2014_DINg/
│   ├── GRE_2014_DINg_PBLT.pos
│   ├── GRE_2014_DINg_PBLT_TFR.mat
│   └── GRE_2014_DINg_PBLT_TFR.wya
│
├── GRE_2015_STAe/
│   ├── GRE_2015_STAe_PBLT.pos
│   ├── GRE_2015_STAe_PBLT_TFR.mat
│   └── GRE_2015_STAe_PBLT_TFR.wya
│
└── ...
```

The naming convention is:

`PATIENT_ID_PBLT_EXTENSION`

For example, `GRE_2014_AUDe` is the patient identifier shared by all three files.

**Required file types:**


| Extension | File example                 |
| ----------- | ------------------------------ |
| `.pos`    | `GRE_2014_AUDe_PBLT.pos`     |
| `.mat`    | `GRE_2014_AUDe_PBLT_TFR.mat` |
| `.wya`    | `GRE_2014_AUDe_PBLT_TFR.wya` |

The three files must be kept together in the corresponding patient folder.

### 2.2. Anatomical ROI information: `PBLT_ROIs.csv`

**Automatic anatomical annotation**

The GitHub repository [`<span>JULICH_ATLAS_PARTIALISATION AUTOMATIQUE</span>`](https://github.com/) provides an automated anatomical annotation pipeline.

If the following columns are missing from the CSV file, this pipeline allows users to automatically generate and add them:


| Column                                   | Description                                                                      |
| ------------------------------------------ | ---------------------------------------------------------------------------------- |
| `<span>Julich_insula_region</span>`      | Anatomical insular region assigned using the Julich atlas, including hemisphere. |
| `<span>Julich_insula_region_rank</span>` | Rank associated with the assigned Julich insular region.                         |
| `<span>Insula_DV</span>`                 | Dorsal or ventral classification of the recording site.                          |

These columns are used by the PBLT Meta-Analysis pipeline to perform region-specific analyses and distinguish dorsal and ventral insular recording sites.

**Note:** These columns are required for the anatomical analyses. If they are already present in `<span>PBLT_ROIs.csv</span>`, no additional annotation step is necessary.


`PBLT_ROIs.csv` contains the anatomical information associated with the insular recording sites.

It is a **semicolon-separated CSV file** (`;`) containing 82 recording sites and 8 columns.

Example:


| PAT           | CHAN |       X |      Y |     Z | Julich_insula_region    | Julich_insula_region_rank | Insula_DV |
| --------------- | ------ | --------: | -------: | ------: | ------------------------- | --------------------------: | ----------- |
| GRE_2014_AUDe | Q02  |  36,324 | 15,468 | 12,13 | Area Id6 (Insula) right |                         3 | Dorsal    |
| GRE_2014_DINg | Z'8  | -34,569 | 25,691 | 7,928 | Area Id7 (Insula) left  |                         2 | Dorsal    |

**Column descriptions:**


| Column                      | Description                                                                      |
| ----------------------------- | ---------------------------------------------------------------------------------- |
| `PAT`                       | Patient identifier. Must match the corresponding patient folder name.            |
| `CHAN`                      | Electrode contact / recording channel identifier.                                |
| `X`, `Y`, `Z`               | MNI coordinates of the recording site (mm).                                      |
| `Julich_insula_region`      | Anatomical insular region assigned using the Julich atlas, including hemisphere. |
| `Julich_insula_region_rank` | Rank associated with the assigned Julich insular region.                         |
| `Insula_DV`                 | Dorsal or ventral classification of the recording site (`Dorsal` / `Ventral`).   |

**Important formatting conventions:**

- The CSV uses `;` as the column separator.
- Decimal coordinates use `,` as the decimal separator.
- Patient identifiers and channel names must be consistent with those used in the electrophysiological data.
- Each row corresponds to one recording site.

When loading the file with Python/Pandas, use:

```python
import pandas as pd

rois = pd.read_csv(
    "PBLT_ROIs.csv",
    sep=";",
    decimal=","
)
```

### 2.3. Expected input directory structure

The complete input directory should follow this organization:

```text
project/
│
├── data/
│   ├── GRE_2014_AUDe/
│   │   ├── GRE_2014_AUDe_PBLT.pos
│   │   ├── GRE_2014_AUDe_PBLT_TFR.mat
│   │   └── GRE_2014_AUDe_PBLT_TFR.wya
│   │
│   ├── GRE_2014_DINg/
│   │   └── ...
│   │
│   ├── GRE_2015_STAe/
│   │   └── ...
│   │
│   └── PBLT_ROIs.csv
│
├── configs/
├── scripts/
├── src/
├── outputs/
└── README.md
```

The patient-level files and the ROI table must be accessible to the pipeline through the paths specified in the configuration files.

---

# 3. Pipeline entry points

Run these commands from the repository root with the `meta_env` environment activated. All stages preserve the existing statistical and plotting implementations; the entry points only select inputs, parameters, and outputs.

Relative paths in YAML files are resolved relative to the YAML file. Therefore, the default `../outputs` value points to the repository-level `outputs/` directory when using the configurations in `configs/`.



## 1. Intermediate data

```powershell
python .\scripts\01_generate_intermediate.py `
    --config .\configs\default.yaml
```

This loads epoch-level TFR files, integrates the configured bands and time
window, classifies contacts, and writes reusable CSV files. Use `--overwrite`
to reuse an explicitly populated output location.

## 2. Comparisons

```powershell
python .\scripts\02_run_comparisons.py `
    --config .\configs\comparisons.yaml `
    --input-output .\outputs\baseline_raw_YYYYMMDD
```

The input is optional; when omitted, the newest compatible output is selected.
The YAML defines the event mapping and correction metadata. The current
implementation produces the existing `-1/+1`, `+1/0`, and `0/-1` results.

## 3. Response plots

```powershell
python .\scripts\03_generate_response_plots.py `
    --config .\configs\default.yaml `
    --input-output .\outputs\baseline_raw_YYYYMMDD
```

This stage reads saved comparison CSVs and creates simple, merged, positive,
and negative response plots. It does not reload TFR data.

## 4. Mixed models

```powershell
python .\scripts\04_run_lmm.py `
    --config .\configs\lmm.yaml `
    --input-output .\outputs\baseline_raw_YYYYMMDD
```

Models are currently fitted independently for each frequency band. The
selected model preset, contrasts, p-value correction, and source output are
recorded in `config.json`.

## 5. Mixed-model figures

```powershell
python .\scripts\05_generate_lmm_figures.py `
    --input-output .\outputs\baseline_raw_YYYYMMDD
```

This stage reads the saved model tables and regenerates the regional raw
slopes, regional forest plots, scaled reference plots, and site-slope
figures without refitting the models.

## 6. Temporal response categories

```powershell
python .\scripts\06_temporal_response_categories.py `
    --input-output .\outputs\baseline_raw_YYYYMMDD `
    --mode interval
```

Use `interval`, `cumulative`, or `both`. Cached results are reused when the
requested parameters match.

## 7. TFR response-category plots

```powershell
python .\scripts\07_plot_tfr_response_categories.py `
    --input-output .\outputs\baseline_raw_YYYYMMDD `
    --sign both
```

This retains the existing event/site t-value calculation and creates the
frequency-band category traces.

## 8. Insula Plotter CSV

```powershell
python .\scripts\08_generate_insula_plot_csv.py `
    --input-output .\outputs\baseline_raw_YYYYMMDD `
    --band HFA `
    --sign positive
```

The output contains the 83 sites and is written under
`data\gain_loss_responses\mni_site_csv\`.

## Run all

```powershell
python .\scripts\run_all.py `
    --config .\configs\run_all.yaml
```

Use `--input-output` to start from an existing intermediate-data output.
Existing outputs are never overwritten unless `--overwrite` is supplied.

## Contents

- `scripts`: step by step scripts running, main scripts are

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
