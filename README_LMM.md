# Methodology: regional linear mixed model

## Scientific question

The analysis tests whether the response to the outcome differs between dorsal
and ventral insular sites. It is repeated independently for every frequency
band and for each configured outcome contrast:

- gain: outcome `0` versus `+1`;
- loss: outcome `-1` versus `0`;
- overall: outcome `-1` versus `+1`.

The two outcome `0` categories are kept separate because they belong to
different experimental tasks.

## Trial preparation

The intermediate files contain one row per integrated trial and site. For a
given contrast, only its two outcome categories are retained and recoded as:

```text
outcome_binary = 0 for the first category
outcome_binary = 1 for the second category
```

The regional indicator is:

```text
ventral = 0 for dorsal sites
ventral = 1 for ventral sites
```

Each site is identified by patient, region and contact index:

```text
site_id = patient | position | contact_index
```

The response variable is the integrated TFR value for one frequency band.
Because this is continuous, the analysis uses a linear mixed model rather than
a binomial model.

## Fitted model

For every frequency band, the regional inferential model is:

```text
response_value ~ outcome_binary * ventral + C(site_id) + (1 | patient)
```

The expanded interpretation is:

```text
response = intercept
         + outcome effect
         + ventral baseline effect
         + outcome-by-ventral interaction
         + site baseline effect
         + patient random intercept
         + residual error
```

The model is fitted with `statsmodels` using REML. The script tries the
optimizers `lbfgs`, `powell`, then `cg`, and keeps the first converged fit.
If no optimizer converges, that band is recorded in the model error output
instead of silently being treated as a successful analysis.

`C(site_id)` controls for baseline differences between contacts. The patient
random intercept accounts for repeated observations within patients.

For the site-level visualization, a second model is fitted with the additional
site-specific outcome interactions:

```text
response_value ~ outcome_binary * ventral
                 + outcome_binary:C(site_id) + C(site_id)
                 + (1 | patient)
```

This second model is used only to derive fitted site slopes. It is not used
for the regional interaction p-value, because its interaction coefficient is
conditional on the reference site after site-specific outcome interactions
have been added.

## Estimated effects

Because dorsal sites are the reference level (`ventral = 0`):

```text
delta dorsal = beta_outcome_binary
delta ventral = beta_outcome_binary + beta_outcome_binary:ventral
```

The regional difference is the interaction coefficient:

```text
delta ventral - delta dorsal
    = beta_outcome_binary:ventral
```

For a site `s`, its fitted outcome slope is obtained from the site-slope model:

```text
site slope(s) =
    beta_outcome_binary
    + beta_outcome_binary:C(site_id)[s]
    + beta_outcome_binary:ventral * ventral(s)
```

The site term is zero for the reference site. The site-slope boxplots show
these fitted slopes separately for dorsal and ventral contacts, with the
individual site estimates overlaid.
Each subfigure also shows the regional interaction p-value in black between
the two boxplots, with black significance stars when the FDR-adjusted test is
significant.

The standard error of the ventral delta includes the covariance between the
outcome and interaction coefficients:

```text
Var(delta ventral)
  = Var(beta_outcome)
  + Var(beta_interaction)
  + 2 Cov(beta_outcome, beta_interaction)
```

All confidence intervals use the Wald approximation:

```text
estimate +/- 1.96 * standard_error
```

## Multiple testing and annotations

For each contrast, Benjamini-Hochberg FDR correction is applied independently
across the five frequency bands to:

1. the regional interaction p-values;
2. the dorsal delta p-values;
3. the ventral delta p-values.

Only FDR-adjusted p-values are exported in the result CSV. In the regional
delta figures, colored stars are shown only when the corresponding 95%
confidence interval excludes zero. A black star between dorsal and ventral
rows indicates a significant FDR-adjusted regional interaction. No text is
shown for nonsignificant results.

## Frequency-band scaling

The scaled regional and site-slope figures divide each band-specific estimate
and uncertainty by:

```text
response_scale_factor = 1 / SD(response_value)
```

The response standard deviation is computed from the non-missing trials used
for that band and contrast. Scaling preserves the zero reference and the sign
of every estimate. It changes visualization units only; the CSV estimates,
confidence intervals and p-values remain unscaled.

## Outputs

For every contrast, the regression stage writes:

```text
data/mixed_model_regression/<contrast>/model_results_by_band.csv
data/mixed_model_regression/<contrast>/model_summaries.txt
data/mixed_model_regression/<contrast>/p_value_notes.txt
```

The figures are written to:

```text
fig/summary/mixed_model_regression/regional_model/
```

The figures are organized in `raw/`, `scaled/` and `site_slopes/`. Band labels
include their frequency ranges, with `BGA (50-150 Hz)` used instead of `HFA`.
The three regional figure types are the interaction forest plot, the raw
regional delta plot and the scaled regional delta plot. The site-slope
figures are stored separately:

```text
site_slopes/<contrast>_site_slopes.png
```

The scale factor is printed below each frequency-band label in every scaled
regional figure. Site-slope boxplots are intentionally kept on their original
scale.
