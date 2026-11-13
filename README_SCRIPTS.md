# PBLT analysis scripts

Run commands from:

```powershell
cd C:\Users\Robin\Documents\THESE_CODE\PROJET_META_ANALYSE\Data_manadgement
```

## Which script does what?

### `analyse_pblt_tous_patients.py`

This is the first-stage entry point. It loads patient TFR files, applies the
selected baseline, integrates the configured frequency bands and time window,
labels contacts as dorsal or ventral, and saves reusable trial-level CSV files
under `data\intermediate_data\`.

It also runs the requested first-stage contact statistics:

```powershell
python .\analyse_pblt_tous_patients.py `
    --baseline_type div `
    --calculation minus1_vs_plus1 `
    --output .\outputs\baseline_div_YYYYMMDD
```

For separate gain and loss contact statistics:

```powershell
python .\analyse_pblt_tous_patients.py `
    --baseline_type div `
    --calculation plus1_vs_zero_and_zero_vs_minus1 `
    --output .\outputs\baseline_div_YYYYMMDD
```

The output folder is named by baseline (`baseline_div`, `baseline_raw`,
`baseline_zscore`, `baseline_log_div` or `baseline_diff`), so every contrast
can be computed from the same folder. If `--output` is omitted, a new
baseline-based folder is created automatically.

To reuse existing intermediate CSV files without loading raw TFR data:

```powershell
python .\analyse_pblt_tous_patients.py `
    --output .\outputs\baseline_div_10211803 `
    --load-intermediate-data `
    --calculation minus1_vs_plus1
```

The command above explicitly generates the overall `+1` versus `-1`
comparison. To run it in the second existing output, use:

```powershell
python .\analyse_pblt_tous_patients.py `
    --output .\outputs\baseline_raw_10211948 `
    --load-intermediate-data `
    --calculation minus1_vs_plus1
```

The resulting plot is
`fig\summary\minus1_vs_plus1\minus1_vs_plus1_by_band_and_region.png`.
Electrode significance comes from the electrode-level FDR tests. A separate
Fisher exact test is then performed for positive electrodes and for negative
electrodes, comparing each sign against non-significant electrodes between
dorsal and ventral regions. Each bar uses the Fisher result for its own sign.

An existing folder must contain `config.json`, and its saved baseline must
match the requested analysis. This prevents accidentally mixing data produced
with different preprocessing.

To include all contacts classified as outside the insula, create a separate
output folder and add `--take_outside`:

```powershell
python .\analyse_pblt_tous_patients.py `
    --baseline_type div `
    --take_outside `
    --calculation plus1_vs_zero_and_zero_vs_minus1 `
    --output .\outputs\baseline_div_all_sites_1023
```

This generates both existing gain/loss figures:

- `fig\summary\plus1_vs_zero_and_zero_vs_minus1\gain_plus1_vs_zero_by_band_and_region.png`
- `fig\summary\plus1_vs_zero_and_zero_vs_minus1\loss_zero_vs_minus1_by_band_and_region.png`

The figures include dorsal, ventral and outside sites. Fisher tests remain
regional tests between dorsal and ventral; outside bars are displayed in grey
and are not included in those two-region Fisher tests.

Gain/loss category figures are now separated by response sign and saved only
as percentages:

- `fig\summary\plus1_vs_zero_and_zero_vs_minus1\positive\`: positive
  gain/loss response categories;
- `fig\summary\plus1_vs_zero_and_zero_vs_minus1\negative\`: negative
  gain/loss response categories.

The corresponding site-by-band data are saved as
`gain_loss_site_responses_positive_by_band.csv` and
`gain_loss_site_responses_negative_by_band.csv`. The `both` response is
subdivided using the FDR-tested global `-1 vs +1` contrast into
`both_prefer_loss`, `both_equal`, and `both_prefer_gain`. The complete
category order is `gain_only`, `both_prefer_gain`, `both_equal`,
`both_prefer_loss`, `loss_only`, and `neither`.

The complete merged pairwise dataset is
`data\gain_loss_responses\pairwise_comparison_results_by_band.csv`. It contains
the FDR-tested results for `+1 vs 0`, `-1 vs 0`, and `-1 vs +1`, including the
three response labels and the merged positive/negative categories. The
corresponding merged category plots are saved as
`merged_gain_loss_site_responses_positive_by_region_all_bandspercent.png` and
`merged_gain_loss_site_responses_negative_by_region_all_bandspercent.png`.

Summary generation can be restricted without recomputing intermediate data:

```powershell
python .\analyse_pblt_tous_patients.py `
    --load-intermediate-data `
    --output .\outputs\baseline_div_all_sites_1023 `
    --calculation plus1_vs_zero_and_zero_vs_minus1 `
    --bands HFA Gamma `
    --plot-outside exclude
```

Use `--bands Theta Alpha Beta Gamma HFA` for all bands,
`--plot-outside all` to display dorsal, ventral and outside, or
`--plot-outside only` to generate outside-only figures. Filtered figures are
stored below
`fig\summary\plus1_vs_zero_and_zero_vs_minus1\bands_<bands>\outside_<mode>\`.
Each stacked segment is labelled `responding_sites/total_sites` and its
percentage, for example `12/40` and `30.0%`.

### `generate_gain_loss_figures.py`

This is the plotting-only entry point. It never loads TFR data and never
recomputes statistics: it reads
`data\gain_loss_responses\pairwise_comparison_results_by_band.csv`.
The generated category figures display counts only, with each segment labelled
`responding_sites/total_sites` (for example `12/40`), without percentages.

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_div_all_sites_1023 `
    --bands Theta Alpha Beta Gamma HFA `
    --plot-outside all
```

Use `--plot-outside exclude` or `--plot-outside only` to control outside
contacts. Figures are written below
`fig\summary\plus1_vs_zero_and_zero_vs_minus1\bands_<bands>\outside_<mode>\`.
To produce CSV results without figures, add `--no-plots` to the analysis
command. The analysis script also writes
`outside_significant_gain_sites_by_band.csv` and
`outside_significant_loss_sites_by_band.csv` in
`data\gain_loss_responses\`; these contain the significant outside sites for
each comparison and their direction.

The analysis also exports electrode lists for each comparison:

- `electrode_responses_plus1_vs_minus1_positive.csv`
- `electrode_responses_plus1_vs_minus1_negative.csv`
- `electrode_responses_plus1_vs_zero_positive.csv`
- `electrode_responses_plus1_vs_zero_negative.csv`
- `electrode_responses_minus1_vs_zero_positive.csv`
- `electrode_responses_minus1_vs_zero_negative.csv`

Each row identifies the patient, region, electrode, and frequency band. The
combined file `electrode_responses_all_comparisons.csv` contains the same
information with an additional `comparison` column.

### `analyse_regression_mixte_intermediate.py`

This is the second-stage inferential analysis. It reads the intermediate CSVs
and automatically computes all three contrasts in the selected output folder:

- `gain_plus1_vs_zero`;
- `loss_minus1_vs_zero`;
- `minus1_vs_plus1`.

```powershell
python .\analyse_regression_mixte_intermediate.py `
    --output .\outputs\baseline_div_10211803
```

To select the newest folder of a baseline type:

```powershell
python .\analyse_regression_mixte_intermediate.py `
    --output-type baseline_div
```

After the intermediate data have been generated, the same command can be
rerun to regenerate all LMM CSV results and figures without reloading raw TFR
files:

```powershell
python .\analyse_regression_mixte_intermediate.py `
    --output .\outputs\baseline_raw_11031234 `
    --p-adjust none
```

This writes the three contrasts (`gain_plus1_vs_zero`,
`loss_zero_vs_minus1`, and `minus1_vs_plus1`) under
`data\mixed_model_regression\`. Regional delta plots and forest plots are
written under
`fig\summary\mixed_model_regression\regional_model\raw\` and
`fig\summary\mixed_model_regression\regional_model\scaled\`; site-slope plots
are written under
`fig\summary\mixed_model_regression\regional_model\site_slopes\`.

To regenerate the gain/loss barplots after the statistical CSVs exist:

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral `
    --categories gain_only both_prefer_gain both_equal both_prefer_loss loss_only neither `
    --clean
```

The plotting script reads
`data\gain_loss_responses\pairwise_comparison_results_by_band.csv`; it does
not reload raw TFR files. Use `--regions dorsal ventral outside` when the
dataset was generated with `--take_outside`.

### `plot_tfr_power_by_response_category.py`

This script reloads the TFR data and integrates the raw power in each
frequency band without baseline normalization. For every site and time point,
it computes a Welch t-value comparing the event epochs (`+1` or `-1`) with
the outcome-zero epochs at the same time point. For every site category, the
event t-value traces are then averaged across sites and shown with a SEM
envelope.

```powershell
python .\plot_tfr_power_by_response_category.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral `
    --sign both `
    --categories loss_only gain_only both both_prefer_loss both_prefer_gain neither `
    --tmin -1 --tmax 2 `
    --integration-tmin 0 --integration-tmax 1.5 `
```

The summary figures are saved in
`fig\summary\tfr_power_by_response_category\<sign>\<band>\<category>`, and the
corresponding CSV files in
`data\tfr_power_by_response_category\<sign>\<band>\<category>`. Each category
folder contains one plot and CSV per region. CSV files contain the event label, `mean_t_value`, `sem_t_value`,
`lower_sem`, `upper_sem`, and `n_sites`. The integration window is shown as a
grey shaded region.

### `analyse_reponses_temporelles.py`

This analysis evaluates the effect of the integration end time. For each
site, band, and consecutive time window, it integrates raw power and applies the
same two-sided Mann-Whitney comparisons as the gain/loss barplot analysis:
`+1 vs 0`, `-1 vs 0`, and `+1 vs -1`. FDR correction is applied across the
tested intervals for each site and band. The `neither` category is excluded
from the plots. Missing categories remain explicitly plotted at zero.

```powershell
python .\analyse_reponses_temporelles.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral `
    --window-start -1 --window-end 3 `
    --step 0.25 --interval-size 0.25 `
    --smooth-window 3 --mode both
```

The site/window results are saved in
`data\temporal_response_categories\temporal_response_by_site_and_window.csv`.
Counts are saved in `response_category_counts_by_tmax.csv`, and interval
figures are written to `fig\summary\temporal_response_categories`. With
`--mode both`, consecutive intervals and cumulative windows `[0,t]`, with
`t = 0.25, 0.5, ..., 3 s`, are generated. Cumulative figures use `t` (the
window end) on the x-axis. Detailed outputs are organized by mode
and frequency band:

```text
data\temporal_response_categories\<mode>\<band>\
fig\summary\temporal_response_categories\<mode>\<band>\
```

`--interval-size` controls the duration of consecutive windows, `--step`
controls their spacing and cumulative endpoint spacing, and `--n-points` can
replace the regular grid with a fixed number of windows. `--smooth-window`
applies a moving-average convolution to the plotted counts; `1` disables
smoothing. Each band directory contains the site/window CSV and the
positive/negative category count CSVs, including `n_timepoints`,
`interval_size`, `n_sites`, and `n_sites_smoothed`.

The script writes `analysis_parameters.json` and reuses the saved site/window
CSV when the requested analysis parameters match. Changing the interval size,
sampling step, number of points, mode, bands, regions, or alpha automatically
triggers a new calculation.

### `plot_sites_mni_3d.py`

This script matches response categories with the MNI coordinates in
`PBLT_ROIs.csv` and creates one figure per band. Each figure contains a 3D
MNI scatter plot plus XZ and YZ projections. Julich insula parcel points are
shown as a translucent colored background; response categories are overlaid
with distinct high-contrast colors. The same 82 sites as the gain/loss barplot
source are displayed (57 dorsal and 25 ventral per band), including `neither`
in gray. Missing parcel labels are assigned to the nearest known Jülich parcel.
Every plotted site is exported with a `parcel_id` and annotated with its
`IdX`/`IaX` parcel identifier in the 3D, XZ, and YZ views. Parcel labels (`Id5`,
`Id6`, etc.) and their
dorsal/ventral status are shown in each projection. Categories are read from
`gain_loss_site_responses_<sign>_by_band.csv`, the same per-band source used
by the merged gain/loss response plot; `both_prefer_gain`, `both_equal`, and
`both_prefer_loss` therefore remain separate.

```powershell
python .\plot_sites_mni_3d.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --sign positive `
    --regions dorsal ventral `
    --categories gain_only both_prefer_gain both_equal both_prefer_loss loss_only neither
```

The command generates both output modes:

- `fig\summary\mni_sites\<sign>\with_neither\`: all sites, including gray
  `neither` points;
- `fig\summary\mni_sites\<sign>\without_neither\`: only responding sites.

The corresponding matched tables and category count tables are saved in
`data\gain_loss_responses\` with `_with_neither` and `_without_neither`
suffixes. The counts can be compared directly with
`gain_loss_site_responses_<sign>_by_band.csv`.

The same command also writes one six-column visualization CSV per band under
`data\gain_loss_responses\mni_site_csv\`. Each file contains exactly the 83
sites from `PBLT_ORIGIN_INSULA.xlsx` with the schema
`site_id,label,x,y,z,color,t_value`. `t_value` is Welch's t-statistic for
`+1 - (-1)` at the corresponding site and frequency band. `both_equal` is
exported as the brown `both` color, and sites absent from the statistical
category table are exported in gray as `neither`. The files use the suffix
`_with_t_value.csv`.

### `extraire_reponses_statistiques.py`

This extracts site response counts and ratios from
`resultats_statistiques_tous_patients.csv`. It is optional and is useful for
the response-ratio barplots:

```powershell
python .\extraire_reponses_statistiques.py `
    --output-type baseline_div
```

The default input and output are respectively
`resultats_statistiques_tous_patients.csv` and
`resume_reponses_par_bande_region.csv` in the selected output folder.

The summary also contains two-sided Fisher exact tests for every band and
direction. The FDR-adjusted electrode-level tests first determine which
electrodes are responding. Fisher then compares the responding versus
non-responding electrode counts between dorsal and ventral regions. Fisher
significance uses the configured alpha threshold; it is not another FDR
correction.

### `run_plot_reponses_gain_loss.py`

This creates the count and percentage figures from the extracted response
summary. Run it after `extraire_reponses_statistiques.py`:

```powershell
python .\run_plot_reponses_gain_loss.py `
    --output-type baseline_div
```

Dark bars indicate a significant Fisher region association; light bars
indicate a non-significant association. This applies separately to positive
and negative bars in both the percentage and count figures.

### `analyse_contrastes_0_minus1_plus1.py`

This is an optional descriptive gain-versus-loss comparison at site level. It
uses the same intermediate data and writes its results to
`data\contrast_gain_vs_loss\`. It is not required for the regional mixed model:

```powershell
python .\analyse_contrastes_0_minus1_plus1.py `
    --output .\outputs\baseline_div_10211803
```

## Recommended execution order

For a new analysis:

1. Run `analyse_pblt_tous_patients.py` on raw data.
2. Run `analyse_regression_mixte_intermediate.py` on the resulting folder.
3. Optionally run `extraire_reponses_statistiques.py`, then
   `run_plot_reponses_gain_loss.py`.
4. Optionally run `analyse_contrastes_0_minus1_plus1.py` for the separate
   descriptive gain/loss comparison.

The reusable implementation is in `src\`. Entry-point scripts only parse
arguments and orchestrate the workflow. In particular,
`src\intermediate_data.py` owns intermediate CSV persistence and
`src\output_management.py` validates output-folder reuse.

## Output structure

```text
outputs\<baseline_folder>\
  config.json
  data\
    intermediate_data\
    mixed_model_regression\<contrast>\
  fig\
    per_unit\
    summary\
```

The mixed-model figures are organized under
`fig\summary\mixed_model_regression\regional_model\raw`,
`scaled` and `site_slopes`.
