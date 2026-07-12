# easyRasch2jmv 3.0.0

- **CICC: corrected and extended total-score grouping.** The grouping
  option previously mislabelled as "equal-width theta intervals" is in
  fact quantile-based grouping of the **total score** (approximately
  equal numbers of respondents per group; all CICC grouping happens on
  the total-score scale, following Buchardt, Christensen & Jensen,
  2023). The option is now labelled accordingly, a new **equal-width
  total-score intervals** method was added (via the easyRasch2
  development version), and the number-of-intervals setting is greyed
  out when "each total score separately" is selected (it never applied
  there). The error-band option is now described as what it is — a
  model-based band around the expected curve, complementary to the
  error bars on observed means, not a replacement — and the note
  explains how to read it. Buchardt et al. (2023) added to the
  references. (Manual score cutpoints are available in
  easyRasch2::RMitemICCPlot() via `score_breaks`; deliberately not
  exposed in the GUI.)
- **Larger plot text for package-drawn figures**: every figure rendered
  directly by easyRasch2 (CICC, item probability curves, targeting,
  both CFA figures, the Q3 heatmap, the LR-DIF locations figure, and
  the response-distribution tileplots) now has its text enlarged to the
  module's base size 15, matching the module's own plots, while each
  figure's specific theming is preserved.
- **The module now installs easyRasch2 from GitHub** (development
  version, `Remotes: pgmj/easyRasch2`) rather than CRAN, picking up
  unrounded data-frame output, the all-NA-respondent crash fixes,
  seed-reproducible RMU, and the LD gamma SE/CI columns without
  waiting for the next CRAN release.
- **Andersen LR-test DIF now uses the easyRasch2 R package directly**
  (`easyRasch2::RMdifLR()`, which deliberately remains based on
  `eRm::LRtest()`): per-group location/threshold tables, the MaxDiff
  flagging, and the per-group locations figure are all identical to the
  R package. The response-distribution tileplot is drawn by
  `easyRasch2::RMplotTile()`. This completes the migration: **all 15
  analyses now delegate their computations to easyRasch2**, and the
  module's vendored simulation and estimation code is fully retired.
- **Partial gamma DIF now uses the easyRasch2 R package directly**
  (`easyRasch2::RMdifGamma()` / `RMdifGammaCutoff()`). The simulated
  expected ranges' generating model moves from eRm-MLE person parameters
  to the WLE theta pool, shifting the ranges slightly; observed gammas,
  SEs, and p-values are unchanged (same `iarm` statistic). The
  response-distribution tileplot is now drawn by
  `easyRasch2::RMplotTile()` (the function the module's version was
  ported from). The simulated-gamma figure keeps the module's rendering,
  including the observed gamma's 95% Wald CI segment.
- **Partial gamma local dependence now uses the easyRasch2 R package
  directly** (`easyRasch2::RMlocdepGamma()`), and gains the
  simulation-based machinery its DIF sibling already had:
  - **New: simulation-based per-pair expected ranges** (parametric
    bootstrap under the fitted model, no true LD by construction), with
    Lower/Upper columns and above-range flags on both direction tables.
  - **New: a per-pair figure** (`easyRasch2::RMlocdepGammaPlot()`):
    simulated partial-gamma distributions as dot clouds with observed
    values as orange diamonds, ranked by deviation from the simulated
    null.
  - The SE and 95% CI columns come from the easyRasch2 development
    version's data output (with an older easyRasch2 they are sourced
    from the same underlying `iarm::partgam_LD()` call, identical by
    construction).
- The module's vendored simulation code is now fully retired: with all
  simulation-based analyses delegated to easyRasch2, the only retained
  helper is the observed-CFA fallback used when the CFA cutoff
  simulation fails.
- **Residual PCA now uses the easyRasch2 R package directly**
  (`easyRasch2::RMdimResidualPCA()` / `RMdimResidualPCACutoff()`).
  Standardized residuals move from eRm's `itemfit()` to the native
  CML/WLE residuals, and the variance partition now **retains extreme
  scorers** (WLE person estimates are finite there; previously only
  non-extreme cases entered the partition). Eigenvalues, loadings, the
  partition percentages, and the simulated cutoff all shift accordingly.
- **Dynamic CFA fit-index cutoffs now use the easyRasch2 R package
  directly** (`easyRasch2::RMdimCFACutoff()` / `RMdimCFA()` /
  `RMdimCFAPlot()`); the simulated datasets are now generated from
  psychotools CML + WLE (was eRm), shifting the cutoffs slightly. Two
  additions:
  - **New: standardized-loadings table and figure** — each item's
    observed one-factor loading against its simulated expected range,
    flagged above/below, pointing to the items driving
    multidimensionality.
  - The fit-index figure is now drawn by `RMdimCFAPlot()`.
- **New analysis: Item Characteristic Curves (CICC)** — conditional item
  characteristic curves via `easyRasch2::RMitemICCPlot()`: model-expected
  item score curves with observed class-interval averages overlaid
  (Kreiner-style graphical item fit), with options for interval
  construction (equal-width theta intervals or per raw score), number of
  intervals, confidence intervals or error bands, and a minimum
  observations-per-interval filter. An optional DIF variable draws the
  observed averages separately per group, with partial-gamma DIF
  annotations per item.
- **Item probability curves now use the easyRasch2 R package directly**:
  probabilities come from `easyRasch2::RMitemCatProb()` (CML via
  psychotools instead of eRm; curve positions shift very slightly for
  polytomous data). The polytomous faceted plot is drawn by the package;
  the dichotomous joint-ICC view (all items in one panel) remains a
  module-specific presentation of the package-computed curves.
- **Targeting plot now uses the easyRasch2 R package directly**: the
  Wright map is drawn by `easyRasch2::RMtargeting()` and the threshold
  table comes from `easyRasch2::RMitemParameters()` (both using the same
  CML-with-MML-sparse-fallback estimator selection as before).
  Consequences:
  - The **person histogram now uses WLE person locations including
    respondents with extreme (minimum/maximum) scores**, which the
    previous eRm-MLE estimates extrapolated; the person distribution
    shifts accordingly. Threshold standard errors shift by a few percent
    (psychotools vs eRm covariance).
  - **New: the threshold table gains Lower/Upper Wald CI columns**
    (shown when CIs are enabled), matching the intervals drawn in the
    figure.
- **Reliability now uses the easyRasch2 R package directly**
  (`easyRasch2::RMreliability()`), replacing the module's adapted
  implementation. Three substantive changes:
  - The **"Empirical" row is replaced by "Marginal"** — the model-based
    marginal reliability (Green, 1984; CML test information integrated
    over the estimated latent distribution). The previous
    `mirt::empirical_rxx()` estimate was the EAP-twin of the PSI and
    therefore largely redundant; the marginal coefficient is the
    complementary model-based view, and a large PSI-vs-Marginal gap now
    serves as an off-target diagnostic. The theta-estimator option now
    affects only the RMU plausible values.
  - **PSI moves from `eRm::SepRel()` to the native WLE-based person
    separation index** (CML item parameters via psychotools; min/max
    scorers excluded). Values shift, noticeably for samples with many
    extreme scores.
  - **The bootstrap option now yields HDCIs for Cronbach's alpha, PSI,
    and Marginal reliability** (respondents resampled, all three
    recomputed natively per resample), not just alpha. The option kept
    its internal name (saved analyses continue to work) but was retitled
    accordingly.
- **Sum score to logit transformation now uses the easyRasch2 R package
  directly** (`easyRasch2::RMscoreSE()`), replacing the module's verbatim
  port of the same WLE solver. Point estimates and standard errors are
  unchanged (the 2.0.1 release had already aligned the module's WLE SEM
  with `RMscoreSE()`); the internal CML fit moves from eRm to
  psychotools, and the module's now-unused ported WLE helpers
  (`utils-theta.R`) were removed.
- **Item-restscore and bootstrap item-restscore now use the easyRasch2
  R package directly**, replacing the module's adapted implementations.
  - Item-restscore delegates to `easyRasch2::RMitemRestscore()`: same
    gamma statistics and BH-adjusted p-values (conditional, hence
    engine-invariant), but the **Rel. location column shifts slightly**
    (CML thresholds via psychotools with a WLE person-mean reference,
    replacing eRm-MLE). Values pass through the R package's data-frame
    output (unrounded with the easyRasch2 development version).
  - Bootstrap item-restscore runs the bootstrap via
    `easyRasch2::RMitemRestscoreBoot(output = "raw")` — numerically
    identical draws and classifications with the same seed — with the
    percentage table and violin plot built from the raw per-iteration
    data (percentages remain unrounded). The dichotomous per-iteration
    refit moves from `eRm::RM` to `psychotools::pcmodel`, which is
    substantially faster; the default number of iterations was raised
    from 200 to 250, matching the other simulation-based analyses.
  - Both analyses drop respondents with no responses at all up front
    (reported in the note), as the Q3 and infit analyses do.
- **Conditional item infit (standard and multiple-imputation) now uses
  the easyRasch2 R package directly**, replacing the module's adapted
  implementation. The per-analysis model fitting moves from eRm (MLE
  person parameters) to conditional maximum likelihood via psychotools
  with Warm's weighted likelihood (WLE) person estimates. Consequences:
  - Observed infit MSQ values are engine-invariant (they match the
    previous eRm route to ~1e-6), but the **Rel. location column shifts
    slightly** (the person-mean reference is now WLE instead of eRm-MLE)
    and the simulation-based cutoffs shift slightly (the generating theta
    pool is now WLE). Results are numerically identical to
    `easyRasch2::RMitemInfit()` / `RMitemInfitCutoff()` /
    `RMitemInfitMI()` / `RMitemInfitCutoffMI()` with the same seed and
    iterations.
  - The bootstrap is substantially faster per iteration, so the standard
    analysis' default iterations were raised from 200 to 250, matching
    the other simulation-based analyses (the small-sample exception —
    around 100 iterations can outperform more; Johansson, 2025 — still
    applies and remains documented in the notes). The MI analysis'
    default of 500 total iterations is unchanged.
  - The MI analysis keeps its jamovi-specific imputation layer — the
    imputation-method choice, retry logic, and **auxiliary variables**
    (not available in the R package, which takes a ready-made mids
    object) — and passes an items-only mids object to easyRasch2.
  - Table values pass through the R package's data-frame output. The
    module installs the easyRasch2 development version (GitHub), which
    returns unrounded values; with easyRasch2 1.0.0 they arrive rounded
    (infit and pooled SE to 3 decimals, Rel. location to 2).
  - The HDCI width option's maximum changed from 100 to 99.9 in both
    analyses (consistent with the Q3 analysis).
- **Q3 residual analysis now uses the easyRasch2 R package directly**
  (easyRasch2 1.0.0 is on CRAN), replacing the module's adapted
  implementation. Estimation moves from MML (mirt) to conditional maximum
  likelihood item parameters (psychotools) with Warm's weighted likelihood
  (WLE) person estimates — true to the Rasch tradition and finite at
  extreme scores. Results are numerically identical to
  `easyRasch2::RMlocdepQ3()` / `RMlocdepQ3Cutoff()` / `RMlocdepQ3Plot()`
  run with the same seed and iterations. Consequences:
  - Q3 values shift slightly relative to earlier module versions
    (CML/WLE vs MML Q3 typically correlate > 0.99); simulation-based
    cutoffs shift accordingly. Saved analyses will show updated numbers
    when re-run.
  - The parametric bootstrap is substantially faster per iteration
    (roughly 3x for polytomous, 14x for dichotomous data), so the
    default number of iterations was raised from 100 to 250, matching
    the other simulation-based analyses.
  - **New output: Q3 heatmap** — a lower-triangle tile plot of the
    observed Q3 matrix with a diverging fill centred on the mean
    off-diagonal Q3, and pairs above the global dynamic cut-off outlined
    in black. Shown when the cutoff simulation is enabled.
  - The HDCI width option's maximum changed from 100 to 99.9 (the
    package requires a width strictly below 100%).
  - Data with no complete response rows are no longer rejected:
    CML/WLE estimation handles incomplete response patterns directly,
    so the previous "no complete cases" stop was removed. Respondents
    with no responses at all on the selected items are excluded and
    reported in the note.

# easyRasch2jmv 2.0.1

- **"Misfit" column renamed "Flagged"** in the item-restscore, bootstrap
  item-restscore, conditional infit, and MI conditional infit tables, so
  the flag column has a consistent header across the module (the Q3,
  partial gamma DIF, CFA, residual-PCA and LR-DIF tables already used
  "Flagged"). The cell content is unchanged — still "overfit"/"underfit"
  for the restscore/infit analyses and "above"/"below" for Q3 — so the
  direction of misfit remains visible; only the header and its footnote
  wording changed.
- **Item-restscore**: removed the "p-value sign." (significance stars)
  column. It is redundant given the exact adjusted p-value is reported
  and the Flagged column already marks significant deviations; the
  footnote dropped the star-coding legend accordingly.
- **Sum score to logit transformation**: the WLE standard error is now
  the information-based `1 / sqrt(I(theta))` evaluated at each estimate
  (as in catR / TAM and `easyRasch2::RMscoreSE()`), replacing the iarm
  "expected SEM" used previously. Point estimates are unchanged; the
  standard errors differ — they are larger at the score extremes, which
  is the more accurate behaviour, and Warm's bias correction keeps the
  lowest and highest scores finite. The WLE solver and grand-mean-zero
  threshold centring are now shared internal helpers ported from
  `easyRasch2`, and the output reproduces `RMscoreSE()` exactly. Only the
  score-to-logit WLE path is affected — no other analysis reports a
  per-score WLE standard error.

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
- **Robustness**: simulations require >= 20 successful iterations before
  producing cutoffs (degenerate cutoffs from near-total failure are no
  longer possible), with graceful degradation and the dominant failure
  reason reported; a caveat is shown when fewer than 100 succeed. There
  is no success-rate gate, so legitimate small-sample runs are not
  blocked.
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
