# PBLT statistical analysis outputs

The output contains the integrated data and the selected first-stage analysis.
The regional mixed-model script can subsequently compute all three contrasts:

- `gain_plus1_vs_zero`: gain comparison between 0 and +1;
- `loss_zero_vs_minus1`: loss comparison between -1 and 0;
- `minus1_vs_plus1`: comparison between -1 and +1 outcomes.

The output-folder name is based on the baseline method, not on the selected
contrast. For example, `baseline_div_YYYYMMDD` uses division baseline
correction.

## Folders

- `config.json`: parameters and command arguments used for the analysis.
- `data/intermediate_data/`: integrated trial data, organized by insular region.
- `data/dorsal/` and `data/ventral/`: statistical result tables for `-1 vs +1`.
- `data/gain_loss_responses/`: tables for the separate gain/loss comparisons.
- `data/mixed_model_regression/`: global mixed-model results and site slopes.
- `data/summary/`: summary tables used by the frequency-band barplots.
- `fig/per_unit/<analysis>/`: one figure per analyzed contact and region.
- `fig/summary/<analysis>/`: summary barplots by frequency band and region.
- `fig/summary/mixed_model_regression/`: regional and site-slope figures.

Positive bars are displayed above zero and negative bars below zero. Negative
axis tick labels show absolute percentages for readability.
Electrode significance is determined by the FDR-adjusted electrode tests.
Fisher's exact test then compares responding versus non-responding electrode
counts between dorsal and ventral regions. Existing gain/loss bars are dark
when this regional Fisher test is significant and light otherwise.
