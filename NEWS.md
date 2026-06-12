# easyRasch2jmv 2.0.0

Major consistency and documentation release: every analysis was reviewed
for consistent UX, analytical options, and documentation detail. The
complete record is in
[CHANGELOG-2.0.0.md](https://github.com/pgmj/easyRasch2jmv/blob/main/CHANGELOG-2.0.0.md);
highlights:

- **New analysis: Partial Gamma Local Dependence** (`locdepgamma`),
  testing item pairs for residual association in both rest-score
  directions, with significance/magnitude filters and top-N display.
- **Q3**: new per-pair table with simulation-based "above"/"below"
  flagging (per-pair HDCI intervals) alongside the global cutoff;
  observed Q3 precision raised from 2 to 4 decimals.
- **Targeting**: automatic CML-to-MML fallback when response categories
  are sparse (ports `easyRasch2::RMtargeting()`).
- **ICC plot rebuilt as a ggplot** (ports `easyRasch2::RMitemCatProb()`),
  with the model auto-selected (RM/PCM) instead of always PCM.
- **Module-wide standardizations**: Observed/Expected column labels with
  "Expected range" super-titles; "overfit"/"underfit" misfit vocabulary
  in the infit and restscore analyses; HDCI width defaults unified at
  99%; BH p-adjustment hardcoded; seeds always applied (reproducible by
  default); raw values with jamovi number formatting throughout.
- **Robustness**: all simulations now require >= 20 successful
  iterations and a 50% success rate, with the dominant failure reason
  reported -- degenerate cutoffs from near-total simulation failure are
  no longer possible.
- **Sparse-data alerts**: every analysis warns when response categories
  have fewer than 3 observations; DIF analyses check within each group
  and point to the tileplot.
- **Documentation**: sample-size and missing-data handling notes
  everywhere (including where statistics and locations use different
  samples), table footnotes for abbreviations and classification rules,
  and explanatory messages instead of silent returns when too few items
  are selected.

# easyRasch2jmv 1.0.0

- Fixes based on Jamovi Module Audit Report
  - Added missing `ggrepel` dependency in DESCRIPTION/Imports
  - Bumped version number to 1.x
  - Table structure modifications for lrdif, reliability, and locdepq3

# easyRasch2jmv 0.5.1

## Targeting plot

- The bottom panel (item threshold locations) now optionally shows
  horizontal **confidence intervals** around each threshold, mirroring
  the behaviour of `easyRasch2::RMtargeting()`. New options:
  * **Show CIs around item threshold locations** (default: on).
  * **Confidence level (%)** (default: 95).
  Polytomous items get CI bars dodged per threshold so they don't
  overlap. The x-axis auto-expands to include the CI endpoints, and
  the caption gains a note describing the CI width. SEs come from
  `eRm`: `erm_out$se.beta` for dichotomous, `thresholds(fit)$se.thresh`
  for polytomous.

## Visual consistency with easyRasch2

- All ggplot output now uses three shared internal theme helpers
  (`er2_axis_margins()`, `er2_plot_caption()`, and `er2_caption()`)
  that mirror what `easyRasch2` itself applies to its R-package plots.
  Net effect for jamovi users:
  * Every plot has the same extra breathing room around the x and y
    axis titles.
  * Every figure caption renders left-aligned, italic, at 10 pt with a
    "Note. " prefix.
  * Long captions are wrapped at 90 characters via `strwrap()` so they
    no longer run off the right edge of the plot (this fixes a
    cut-off caption that was visible on the conditional-infit plot in
    earlier versions).
- Figure captions changed their prefix from "Note: ..." (colon) to the
  APA-conventional "Note. ". The body text of each
  caption is unchanged.

# easyRasch2jmv 0.5.0

## Q3 residual correlations

- Added a **per-pair Q3 simulation plot** to the Q3 analysis (visible
  when *Compute simulation-based cutoff* is enabled), mirroring the
  design of the conditional-infit MSQ plot: one row per item pair with
  a `ggdist` dot cloud of the simulated null distribution, the
  per-pair simulation median as a black dot, and the observed Q3 from
  the mirt fit overlaid as an orange diamond. A dashed reference line
  at 0 (Q3 = 0 = local independence) makes it easy to read pairs that
  fall above the simulated null cloud.
- New **Number of item pairs to plot** option (default 10). Pairs are
  ranked by their *deviance from the simulated null*
  (`|observed Q3 − median(simulated Q3 per pair)|`), with the
  most-deviant pair at the top of the plot. Set higher to see more
  pairs; set very high (e.g. 500) to see every pair.
- The underlying simulation function now retains and aggregates
  per-pair Q3 residuals across iterations (in addition to the
  pre-existing `mean` / `max` scalars). The global cutoff scalar
  (`p99` of `max − mean`) is unchanged.

## Partial gamma DIF

- Added a **p-value sign.** star-string column (`""` / `"."` / `"*"` /
  `"**"` / `"***"`) sourced from `iarm::partgam_DIF()`, displayed
  immediately after the adjusted p-value. Matches the convention used
  by item-restscore output in this module.
- Retitled `Adjusted p (BH)` → `Adj. p-value (BH)` for header
  consistency with item-restscore.

## Item-restscore

- The `Abs. difference` column is now a *signed* `Difference`
  (observed minus expected). Positive values indicate
  over-discrimination (often associated with local dependence);
  negative values indicate under-discrimination (often associated
  with multidimensionality or noise). Sorting by *Sort by absolute
  difference* still puts the largest-magnitude misfits at the top
  but now lets you see the direction at a glance.

## Decimal formatting (jamovi `format: zto`)

- **Sum-score-to-Logit table** — `Logit score` column now uses jamovi's
  `zto` format for consistent decimal-place display.
- **Q3 residual correlation matrix** — the entire matrix now uses
  `zto` formatting.
- **Item-restscore** — the `Location` column now uses `zto` formatting.

## Renaming alignment with easyRasch2 0.8.0

- All internal calls to `easyRasch2::RM*()` were updated for the
  upstream renaming in easyRasch2 0.8.0 (e.g. `RMpartgamDIF()` →
  `RMdifGamma()`, `RMitemrestscore()` → `RMitemRestscore()`, etc.).
  Jamovi UI labels and option names are unchanged.

# easyRasch2jmv 0.4.1

Only minor changes/fixes:

- Slight modification of the dynamic CFA cutoff to always rely on the .scaled
  metrics rather than using the .robust metrics when available. In practice this
  makes little difference, since the important thing is to use the same metric from 
  the simulations as from the observed values to evaluate model fit. 
  The change is for consistency and to make sure that one gets fit values from 
  each iteration in the simulation, even with small samples. You may  
  note different actual metrics reported compared to the earlier version, but 
  model flagging should be the same.
- Removed the color legend for PCA loadings/locations plot.

# easyRasch2jmv 0.4.0

- New function for testing unidimensionality with Confirmatory Factor Analysis
  (CFA) using `lavaan` and simulation-based cutoff values for model fit metrics.

# easyRasch2jmv 0.3.1

- New function for LRT-based DIF analysis of categorical DIF variables.
- New function for PCA (principal component analysis) of standardized Rasch model
  residuals with (experimental) simulation based cutoff value for the largest 
  PCA eigenvalue.
- Partial Gamma DIF has an added option for a tileplot of response data grouped
  by the DIF variable.

# easyRasch2jmv 0.3.0 

New functions:

- Four reliability metrics: Person Separation Index (PSI) via
  eRm::SepRel(); empirical reliability via mirt::empirical_rxx(); Cronbach's 
  alpha (including bootstrap confidence intervals); and
  Relative Measurement Uncertainty (RMU) computed from mirt plausible
  values using the Bignardi, Kievit & Bürkner (2025) split-half
  correlation method.
- Transformation table from ordinal sum score to logit score, with either WLE or
  EAP scores.
- Non-parametric bootstrap of item-restscore associations,
  recommended for use with large sample sizes (n > 800). 
  See <https://pgmj.github.io/rasch_itemfit/> for more details.
- Conditional item infit with missing data, using multiple imputation. Since conditional 
  item infit needs complete data, incomplete responses are discarded. This function
  uses `mice` to impute n datasets with complete responses, then calculates infit from
  all datasets and pools the results using Rubin's rules. Optionally also uses
  simulation to determine appropriate cutoff for interpreting infit.
  
Fixes/modifications:

- Partial gamma DIF bug fix: no longer allows numeric type variable as DIF variable.
- New data validation to better handle labelled ordinal data from SPSS, etc.
- Conditional item infit now clearly states it only uses complete response data,
  and also notes the complete response sample size in the table/figure footnote/caption.

# easyRasch2jmv 0.2.4

- Targeting (Wright map) figure added.

# easyRasch2jmv 0.2.3.2

- Updated description of module
- Fix for clearWith showLegend in item threshold probability figures.
- Added references

# easyRasch2jmv 0.2.3

- Now available from Jamovi Library
- Five functions implemented. Beyond the three listed previously:
  - Partial gamma DIF
  - Item threshold probability figures

# easyRasch2jmv 0.1.0

- Initial release with 3 functions
  - Item-restscore
  - Q3 residual correlations
  - Conditional item infit
