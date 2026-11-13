# Gain/loss electrode lists and barplots

Run these commands from:

```powershell
cd C:\Users\Robin\Documents\THESE_CODE\PROJET_META_ANALYSE\Data_manadgement
```

## 1. Generate statistical CSV results

Use the analysis script to recompute the statistical results without creating
figures:

```powershell
python .\analyse_pblt_tous_patients.py `
    --load-intermediate-data `
    --output .\outputs\baseline_div_all_sites_1023 `
    --calculation plus1_vs_zero_and_zero_vs_minus1 `
    --take_outside `
    --no-plots
```

The source file for plotting is:

```text
outputs\baseline_div_all_sites_1023\data\gain_loss_responses\pairwise_comparison_results_by_band.csv
```

## 2. Export electrode lists and generate barplots

The plotting script reads the source CSV only. It does not reload TFR data or
recompute statistics.

If `pairwise_comparison_results_by_band.csv` is missing but the output folder
contains `data\intermediate_data\` and `config.json`, the plotting script
automatically runs the statistical CSV step with the saved baseline and time
window before plotting.

Example with all bands, all regions, and all categories:

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_div_all_sites_1023 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral outside `
    --categories gain_only both_prefer_gain both_equal both_prefer_loss loss_only neither `
    --clean
```

The `--clean` option removes only the generated electrode-list CSV directory,
the legacy electrode-list CSVs (`electrode_responses_*.csv`,
`outside_significant_*.csv`) and generated category figures. It does not
delete the statistical source CSV.

## Available filters

Bands:

```text
Theta Alpha Beta Gamma HFA
```

Regions:

```text
dorsal ventral outside
```

Categories:

```text
gain_only
both_prefer_gain
both_equal
both_prefer_loss
loss_only
neither
```

For example, to export and plot only HFA dorsal and outside sites for the two
dual-response categories:

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_div_all_sites_1023 `
    --bands HFA `
    --regions dorsal outside `
    --categories both_prefer_gain both_prefer_loss `
    --clean
```

The bar height remains a percentage of all sites in each
band and region. Segment labels show counts only, in the form `n/total`.
Filtering categories changes which segments are displayed but does not change
the denominator.

## Output locations

Filtered electrode lists are saved under:

```text
outputs\baseline_div_all_sites_1023\data\gain_loss_responses\filtered_electrode_lists\
```

The files are separated by comparison and direction:

```text
plus1_vs_minus1_positive_*.csv
plus1_vs_minus1_negative_*.csv
plus1_vs_zero_positive_*.csv
plus1_vs_zero_negative_*.csv
minus1_vs_zero_positive_*.csv
minus1_vs_zero_negative_*.csv
```

Figures are saved under:

```text
outputs\baseline_div_all_sites_1023\fig\summary\plus1_vs_zero_and_zero_vs_minus1\
bands_<bands>\regions_<regions>\categories_<categories>\
```

Two figures are generated in this directory:

```text
positive\merged_gain_loss_site_responses_positive_by_region.png
negative\merged_gain_loss_site_responses_negative_by_region.png
```

## Raw non-baselined dataset

To generate a raw analysis over 0 to 1.5 seconds, then run the mixed model
and barplots:

```powershell
python .\analyse_pblt_tous_patients.py `
    --baseline_type raw `
    --calculation plus1_vs_zero_and_zero_vs_minus1 `
    --output .\outputs\baseline_raw_11031234 `
    --tmin 0 `
    --tmax 1.5 `
    --no-plots

python .\analyse_regression_mixte_intermediate.py `
    --output .\outputs\baseline_raw_11031234 `
    --p-adjust none

python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral `
    --clean
```

The plotting script automatically regenerates the pairwise CSV if it is
missing and intermediate data plus `config.json` are available.

## Regenerate plots after the analyses

Once the statistical analysis has finished, regenerate the gain/loss barplots
directly from the saved pairwise CSV:

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_raw_11031234 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral `
    --categories gain_only both_prefer_gain both_equal both_prefer_loss loss_only neither `
    --clean
```

For the division-baseline dataset, replace the output path:

```powershell
python .\generate_gain_loss_figures.py `
    --output .\outputs\baseline_div_all_sites_1023 `
    --bands Theta Alpha Beta Gamma HFA `
    --regions dorsal ventral outside `
    --categories gain_only both_prefer_gain both_equal both_prefer_loss loss_only neither `
    --clean
```

The barplots are written below:

```text
fig\summary\plus1_vs_zero_and_zero_vs_minus1\
```

### Regenerate LMM and forest plots

After intermediate trial CSVs have been generated, run the mixed-model script.
It recomputes the three contrasts and regenerates the regional delta and forest
figures:

```powershell
python .\analyse_regression_mixte_intermediate.py `
    --output .\outputs\baseline_raw_11031234
```

For the division-baseline dataset:

```powershell
python .\analyse_regression_mixte_intermediate.py `
    --output .\outputs\baseline_div_all_sites_1023 `
    --p-adjust none
```

The LMM CSV results are saved under:

```text
data\mixed_model_regression\
```

The forest plots and regional delta plots are saved under:

```text
fig\summary\mixed_model_regression\regional_model\raw\
fig\summary\mixed_model_regression\regional_model\scaled\
```

The site-slope figures are saved under:

```text
fig\summary\mixed_model_regression\regional_model\site_slopes\
```

The LMM script uses the intermediate data already saved in
`data\intermediate_data\`; it does not reload raw TFR files. If the
intermediate data are missing, run the first-stage analysis before launching
the LMM command.
