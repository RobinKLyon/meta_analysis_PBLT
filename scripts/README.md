# Pipeline entry points

Run these commands from the repository root with the `meta_env` environment
activated. All stages preserve the existing statistical and plotting
implementations; the entry points only select inputs, parameters, and outputs.
Relative paths in YAML files are resolved relative to the YAML file. Therefore
the default `../outputs` value points to the repository-level `outputs\`
directory when using the configurations in `configs\`.

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
