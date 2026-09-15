# easyRasch2jmv 3.2.1

Responses to the 3.2.0 audit. No results change.

## Person Change

- Two computations that sat outside `changeCache` now use it. The pooled
  critical value behind the conditional figure's reference line is a second
  `RMpersonChange()` enumeration, and it moved onto the cache's cold path.
  The retest SD simulation (`RMretestSD()`, 500 iterations by default plus an
  optional 500-iteration bootstrap) gained its own cache element, keyed on the
  options that move it. Both used to rerun on every option change, so
  toggling a *Save to dataset* checkbox with both boxes ticked cost a full
  recomputation. With both ticked on 150 respondents and 6 items, that toggle
  goes from 6.4 s to 0.7 s.
- The item pairing table is built and labelled as soon as variables are
  dropped in, rather than after the enumeration finishes. It exists to be
  checked before the figure it validates, and it used to appear last and go
  blank again on every option change.
- The retest SD table's four rows are created up front for the same reason.

## Reliability

- The conditional precision figure is built once in the analysis and stored,
  rather than redrawn by the render function. jamovi calls that function on
  every resize and every export, so with *Bootstrap confidence band* on, each
  one recomputed the band: 2.6 s at the default 200 iterations on 50
  respondents and 5 items, 25 s at the 2000-iteration maximum. Redrawing is
  now 0.1 s. The interval width and the seed also stopped invalidating the
  cached curve while the band is off, which is the only thing they reach.

## Saved file size

- Figures held in results state are stored built rather than as `ggplot`
  objects. A `ggplot` is far larger than the figure it describes: in
  ggplot2 4.0.3 an empty one is 102 KB compressed, and the three person-fit
  maps came to 287, 287 and 752 KB. jmvcore warns above 500 KB for a single
  element, and all of it is written into the saved `.omv`. Built, the same
  figures are 12 to 16 KB. Measured on 500 respondents and 20 items, the
  person-fit cache drops from 1022 KB to 62 KB, the person-change cache from
  1181 KB to 39 KB, and the reliability curve from 374 KB to 32 KB.
- `requiresData` removed from 26 of the module's 27 figures. It makes jamovi
  read the whole dataset before rendering, on every redraw and export, and
  only the tree-based DIF figure reads it.

## Labels

- Five checkboxes added in 3.2.0 were retitled to name what they turn on
  rather than the act of turning it on, matching the rest of the module:
  *Marginal reliability as a reference line*, *Respondent distribution behind
  the curve*, *Region reaching the reliability benchmark*, *Value labels as
  category labels*, *Retest SD estimated from these data*.

# easyRasch2jmv 3.2.0

Built against easyRasch2 1.3.1. Adds a longitudinal analysis, folds the new
conditional-precision curve into Reliability, and changes two results: marginal
reliability and the bottom panel of the targeting plot. Also carries the
jamovi library-audit items that were prepared for 3.1.1, which was never
released.

## New analysis: Person Change

- **Person Change** tests, for each respondent, whether their location on the
  latent variable moved between two occasions by more than measurement error
  allows (easyRasch2::RMpersonChange()). It sits under *Person-level* beside
  Person Parameters and Person Fit.

  The two occasions go in two variable boxes and are **paired by their order in
  the boxes, not by name**. A pairing table is shown above the figure for that
  reason: a mispaired analysis runs without complaint and produces a
  plausible-looking result. Data must be in wide format, both occasions on the
  same row.

- **Critical values are enumerated exactly, and there is no setting for them.**
  The sum score is a sufficient statistic, so a respondent who moves to an
  extreme score produces a large change and a large standard error together and
  the ratio is damped. Both tails of the null are pulled in, so the critical
  value sits below 1.96, and further below on short scales. Enumerating it
  needs no simulation, no iteration count and no random seed, and the analysis
  repeats identically.

- Occasion-to-occasion fluctuation can be treated as noise rather than change
  by switching the null and supplying a per-occasion retest SD. A retest SD can
  also be **estimated from the two occasions** (easyRasch2::RMretestSD()),
  valid only when they are a stability study with no expected change. It is
  reported rather than applied, so the value in force is always the one shown
  in the options.

- Critical values can be computed **per respondent** instead of pooled. An
  extra figure then plots each respondent's critical value against their own
  location, which shows the values pulling in toward the ends of the scale
  where fewer scores remain attainable.

- The per-respondent results table is off by default, with a flagged-only
  filter. Person locations, the change, its standard error, the change index,
  the p-value and the classification can be saved to the dataset.

- **Not available in jamovi: the single-respondent case.** It needs item
  parameters from an external calibration, which there is no way to supply
  here. Use easyRasch2::RMpersonChange() in R.

## Reliability

- **Marginal reliability changes formula and its values move upward**, more so
  on short scales. It is now the latent-density-weighted mean of the
  conditional reliability rather than Green's (1984) subtractive coefficient,
  which could fall below zero and was floored there. The row is renamed
  *Marginal (curve mean)*.

- New optional **conditional precision curve**
  (easyRasch2::RMreliabilityCurve()): the standard error of measurement, test
  information, or conditional reliability across the latent scale, with the
  respondent distribution behind it. A benchmark shades the region reaching a
  given reliability and reports the share of respondents inside it. The
  accompanying summary reports the superseded Green coefficient alongside the
  new one, so results from earlier versions can be reconciled, together with
  the sample's mean and SD person location, the mean and SD of the latent
  distribution fitted to them, and the average SEM and test information.

- The grouped heading over the confidence bounds now names the width in force,
  as `95% HDCI` rather than a bare `HDCI`.

## Targeting Plot

- **The bottom panel now draws response-category bands**, showing where each
  category is the most likely answer, with the threshold estimates and their
  confidence intervals beneath. Categories that are never most likely, and
  threshold reversals, are marked in red. Estimates are unchanged. The previous
  dot-and-whisker panel is still available under *Bottom panel*.

- The band palette, its start and end points, and the spacing between item rows
  are adjustable. Where the selected variables carry value labels and agree on
  them, the categories are labelled with those instead of their scores.

## Figures

- Captions are slightly larger and the same size everywhere. They were 10 pt on
  figures the module drew and 9 pt on figures drawn by easyRasch2; both are now
  10.5 pt.

## Menu

- The analyses are grouped into subgroups instead of one flat list: Item fit,
  Dimensionality, Local dependence, DIF, Item probability curves, Targeting,
  and Person-level.

## Dynamic CFA Fit-Index Cutoffs

- The standardized loadings table is now built before the analysis runs, one
  row per selected item with the item names already filled in. Previously it
  appeared as an empty placeholder and restructured itself once the
  simulation finished.
- Item names that are not valid R names, such as `Item 1` or `3 months`, no
  longer break the observed-fit fallback used when the cutoff simulation
  fails. lavaan cannot parse such names, so the fit runs under placeholder
  names and the item names are restored afterwards. The same fix is in
  easyRasch2, where these names made the cutoff simulation itself fail.

## Simulation settings

- The iteration settings that had no upper limit now have one, matching the
  three analyses that already did. Iterations are capped at 5000 in Bootstrap
  Item-Restscore, Dynamic CFA Fit-Index Cutoffs, Conditional Item Infit,
  Conditional Item Infit (Multiple Imputation), Martin-Löf Test, Partial
  Gamma DIF, and Person Fit. In Reliability, plausible-value draws are capped
  at 10000, RMU iterations at 500, and bootstrap iterations at 5000.

  Each of these drives a sequential simulation loop, so run time is linear in
  the value and one extra digit was enough to make jamovi look like it had
  stopped responding. An analysis saved with a value above a new limit will
  reopen at the limit.

## References and labels

- Every entry in the analysis reference lists now carries a link. Six
  references are new, for Person Change and the reliability curve.
- Five variable-box and section labels are title-cased to match the rest of
  the module: `Subscale 1 Items` and `Subscale 2 Items` in Martin-Löf Test,
  `Auxiliary Variables (optional)` in Conditional Item Infit (Multiple
  Imputation), and `Save to Dataset` in Person Fit and Person Parameters.

# easyRasch2jmv 3.1.0

Conditional Item Infit and the two local dependence analyses now flag on the
multiplicity-corrected bootstrap p-value rather than on the expected range,
following Johansson (2026), <https://doi.org/10.31234/osf.io/7pqz4_v2> and
easyRasch2 1.2.0. The remaining simulation-based analyses are unchanged and
follow when their own studies are complete.

## Changed defaults, Conditional Item Infit and both local dependence analyses

- **Bootstrap p-values** default to on. They still require the simulation,
  which stays off by default, so a fresh analysis does not run a bootstrap
  until asked. New for *Partial Gamma Local Dependence*, which had no such
  option, together with its **Multiple-comparison correction** setting.
- **Number of simulation iterations** 250 to 400 in all three analyses, the
  calibrated floor.
- **HDCI width** 99 to 95. The range is a description of where a fitting
  item or pair is expected to fall, not a decision rule.

  The familywise cost of the range is worse for pairs than for items, since
  pairs grow quadratically: a 95% range implies about 37% over the nine items
  of a nine-item scale and about 84% over its 36 pairs. Turning bootstrap
  p-values off now says so in both local dependence analyses.

## Data coded from 1

- Items coded with 1 as the lowest category are now recoded automatically by
  subtracting 1, instead of being refused with "requires items scored
  starting at 0". Every analysis states that it happened.

  The test is deliberately narrow: the shift applies only when *every*
  item's lowest observed response is exactly 1, which is the one arrangement
  where 1-based coding is the obvious reading, and all items are then
  shifted identically. Ragged minima (some items starting at 1, others
  higher), a minimum above 1, and anything else still produce the same error
  as before, and empty categories are handled downstream exactly as they are
  today. A uniform shift is lossless, verified against the same data coded
  from 0.

## Conditional ICC

- The note below the plot now explains the grouping *rule* and no longer
  states how many groups were formed. The requested number is not always the
  realised one: quantile groups merge where total scores tie at a boundary,
  and equal-width intervals can end up with no respondents in them. From
  easyRasch2 1.2.0 the figure caption reports the grouping actually used, so
  the note points there rather than restating a number that could disagree
  with it.

## Table notes

- Turning bootstrap p-values off now states the familywise error rate the
  expected range implies, computed from the width and item count in use
  (`1 - width^k`, so 37% for a 95% range over nine items).
- The iteration caveat is two-tier: below 400 the correction is mildly
  liberal, between 400 and 1000 error rates are calibrated but decisions
  remain seed-dependent.
- The advice that around 100 iterations could improve detection with small
  samples is withdrawn. That advantage came from an expected range that had
  not converged and carried an inflated familywise error rate.
- **Every simulation-based analysis now reports iterations that were lost**,
  not just runs that ended up short. A simulated dataset is discarded when it
  cannot be refitted, usually because an item ended up with an unused
  response category, and until now only Conditional Item Infit said so. A run
  of 400 that delivers 391 previously looked like a run of 400.

Requires easyRasch2 1.2.0 or later for the plot's interval to match the
table's. Against 1.1.x the analysis still runs, but the dot plot draws a
fixed 99.9% whisker regardless of the width setting.

## Partial Gamma Local Dependence

- **New: bootstrap p-values**, with the same *Multiple-comparison correction*
  choices as *Yen's Q3*. With expected ranges on they are the default, and the
  asymptotic adjusted p-value and significance columns give way to them. Turn
  them off to flag against the expected range as before.
- **New column: Gamma pair (max)**, the larger of the pair's two rest-score
  directions. It is the statistic that is tested, so a pair can be flagged
  while the Partial gamma shown in one table sits inside the range.

- The note under each table said a pair is flagged when its partial gamma
  falls outside the expected range, and that the two tables together test
  both rest-score directions. Since easyRasch2 1.2.0 each pair is tested once,
  on the larger of its two directions, and carries the same result in both
  tables. The notes now say so.

# easyRasch2jmv 3.0.0

The migration release: **all analyses now delegate their
computations to the `easyRasch2` R package**, and the module's vendored
estimation and simulation code is retired. Estimation moves to conditional
maximum likelihood (CML) via psychotools, with Warm's weighted likelihood (WLE)
person estimates throughout, replacing the earlier eRm and mirt
engines. Overall consequences:

- Every analysis is numerically identical to its easyRasch2
  counterpart run with the same seed and iterations.
- Observed conditional statistics are essentially unchanged (infit
  matches the previous engine to ~1e-6; gamma statistics and their
  p-values are identical; CML/WLE Q3 typically correlates > 0.99 with
  the previous MML Q3), but **relative item locations, person
  distributions, and simulation-based expected ranges shift slightly**
  — the person reference and the generating theta pool are now WLE
  instead of eRm-MLE. Saved analyses show updated numbers when re-run.
- WLE estimates are finite at extreme scores, so extreme scorers are
  now retained where they previously were extrapolated or excluded
  (the targeting person histogram and the residual-PCA variance
  partition).
- The parametric bootstraps are substantially faster per iteration,
  and **default iterations are unified at 250** across the
  simulation-based analyses. The HDCI width maximum is 99.9%
  throughout. Option changes that do not affect a simulation or
  bootstrap (sorting, filters, p-value settings, plot options, the DIF
  tileplot) no longer rerun it, in every simulation-based analysis —
  including the MI conditional infit analysis, where a sort-order
  change previously reran the imputation and simulation pipeline.
- Table values pass through the package's unrounded data-frame output
  (jamovi's number formatting applies).
- Respondents with no responses on any selected item are excluded up
  front and reported in the note; the Q3 analysis no longer requires
  complete response rows (CML handles incomplete patterns directly).

New functionality:

- **New analysis: Tree-Based DIF** — model-based recursive partitioning
  for DIF via `easyRasch2::RMdifTree()` (psychotree Rasch/PCM trees;
  Strobl et al., 2015): the sample is split wherever item parameters
  are unstable along one or more covariates, so groups need not be
  pre-specified — several covariates at once, continuous covariates
  with data-driven cutpoints, and interactions as nested splits. Each
  split's items get an effect size (Mantel-Haenszel on the ETS Delta
  scale for dichotomous data, partial gamma for polytomous; Henninger
  et al., 2023, 2025) classified A/B/C, with optional iterative
  purification, optional pruning of all-negligible splits, an
  adjustable classification alpha, and an optional p-value adjustment
  for the partial-gamma classification. Output: the tree figure with
  per-node item profiles, and the per-split effect-size table. The
  output notes that the A/B/C boundaries are conventions rather than
  sample-calibrated values, pointing to the Partial Gamma DIF analysis
  with simulation-based cutoffs for a calibrated test.
- **New analysis: Martin-Löf Test** — likelihood-ratio test of
  unidimensionality against an a priori two-subscale partition of the
  items, generalised to polytomous models (Christensen, Bjorner,
  Kreiner & Petersen, 2002), via `easyRasch2::RMdimMartinLof()`. The
  p-value comes from Monte Carlo simulation under the unidimensional
  null (Christensen & Kreiner, 2007; the asymptotic chi-square is
  biased toward conservatism), with optional Besag-Clifford sequential
  stopping for faster runs when the null is compatible with the data.
  Alongside the test, the correlation between the two subscales' WLE
  person estimates is reported with a 95% CI. Two default figures: the Monte
  Carlo null distribution with the observed statistic marked, and the
  observed-vs-expected subscore cross-table residual heatmap
  (`RMdimMartinLofResiduals()`, with an optional minimum-expected-count
  filter for sparse cells). The output states prominently that the
  partition must be a priori: testing a split suggested by the same
  data (e.g. the first residual-PCA contrast) invalidates the p-value.
  Complete cases only (at least 30 required).
- **New analysis: Person Fit** — per-respondent person-fit statistics
  with Monte-Carlo resampled p-values, via `easyRasch2::RMpersonFit()`:
  conditional infit and outfit MSQ (conditional on the total score, so
  no biased person estimate enters and partial missingness is handled
  directly) and the standardized log-likelihood lz, each shown as a
  person-fit map (statistic vs person location, flagged respondents
  highlighted) — all three by default, individually selectable.
  Significance comes from resampling under the fitted model (the
  asymptotic nulls are unreliable; Müller, 2020; Sinharay, 2016), with
  a choice of flagging direction for the MSQ statistics (two-sided, or
  underfit-only for more power against careless responding) and an
  adjustable flagging alpha (default 0.05; a per-person screening flag,
  so about alpha of respondents are flagged by chance under fit). The
  flag, the statistics, and their p-values can be **saved as variables
  in the dataset** — e.g. to filter out aberrant respondents before
  rerunning other analyses. Extreme scorers cannot be assessed and get
  empty cells. A summary table reports assessed/extreme counts and
  flagged counts overall and per statistic. The default of 500
  resampling iterations (rather than the module-wide 250) reflects that
  per-person p-values have resolution 1/iterations.
- **New analysis: Person Parameters** — estimates each respondent's
  location on the latent variable (theta, in logits) with its standard
  error of measurement, by WLE (Warm's weighted likelihood, the
  default; finite at extreme scores) or EAP (normal prior estimated
  from the data by marginal maximum likelihood), via
  `easyRasch2::RMpersonParameters()`. The estimates can be **saved as
  new variables in the dataset** — theta, SEM, sum score, number of
  items answered, and an extreme-score flag — for use in other jamovi
  analyses or for export; respondents with partially missing responses
  are retained (theta estimated from the items they answered), and the
  saved columns stay row-aligned with the spreadsheet. The displayed
  output is a compact summary table (mean, SD, median, MAD, IQR, range
  of theta, mean SEM, extreme-score counts) and a theta histogram whose
  caption reports the number and share of minimum and maximum scores.
  Item parameters switch from CML to MML under sparse response
  categories, with a note, as in the targeting analysis.

  **This analysis replaces the former "Sum Score to Logit
  Transformation" analysis**, which is removed from the menu: the
  score-to-theta lookup table and its figure (the option previously
  titled just "Figure") are now the *Sum score to logit table* and
  *Sum score to logit figure* options here, sharing the WLE/EAP method
  choice and theta range with the person estimates, and computed by
  `easyRasch2::RMscoreSE()` exactly as before. Saved analyses from
  earlier module versions that used the old analysis need to be
  re-created. Note that with EAP, the lookup table's sum-score EAP
  (mirt, standard-normal prior) and the pattern-based EAP person
  estimates use different engines, so values for the same sum score
  can differ slightly (a table note explains this).
- **New analysis: Item Characteristic Curves (CICC)** — model-expected
  item score curves with observed class-interval means overlaid
  (Kreiner-style graphical item fit), via
  `easyRasch2::RMitemICCPlot()`. Grouping happens on the total-score
  scale (Buchardt, Christensen & Jensen, 2023): quantile groups
  (default), equal-width score intervals, or each total score
  separately. Observed means can carry confidence intervals,
  complemented by a model-based error band around the expected curve
  showing where the means should fall if the model holds; a minimum
  observations-per-interval filter is available. An optional DIF
  variable draws the observed means separately per group, with
  partial-gamma DIF annotations per item.
- **Conditional item infit and Q3: optional bootstrap p-values.** With
  simulation-based cutoffs enabled, a checkbox adds Monte-Carlo
  p-values comparing each observed statistic against its simulated
  distribution (two-sided per item for infit; one-sided per item pair
  for Q3, testing excess positive local dependence) and adjusted
  p-values, defaulting to the Westfall-Young step-down familywise
  correction (Ferreira, 2024) with Benjamini-Hochberg and
  Benjamini-Yekutieli FDR alternatives; the adjusted-p column title
  names the chosen method. Flagging then follows the adjusted p-value
  (< 0.05), with the expected range kept as the effect-size reference.
  Matches `easyRasch2::RMitemInfit()` / `RMlocdepQ3()` with
  `p_value = TRUE`; at least 1000 iterations are recommended when
  reporting p-values (noted in the output). The other simulation-based
  analyses will gain the same option in later releases; the asymptotic
  p-values in the item-restscore and partial-gamma analyses keep their
  fixed Benjamini-Hochberg correction (the cheap first screen).
- **Q3 heatmap** — a lower-triangle tile plot of the observed Q3
  matrix with a diverging fill centred on the mean off-diagonal Q3 and
  pairs above the global dynamic cutoff outlined in black (shown when
  the cutoff simulation is enabled).
- **Partial gamma local dependence: simulation-based expected ranges
  and a per-pair figure**, matching the machinery its DIF sibling
  already had: parametric bootstrap under the fitted model (no true LD
  by construction), Lower/Upper columns and above-range flags on both
  direction tables, and simulated per-pair distributions as dot clouds
  with the observed values overlaid.
- **CFA: standardized-loadings table and figure** — each item's
  observed one-factor loading against its simulated expected range,
  flagged above/below, pointing to the items driving
  multidimensionality.
- **Targeting: the threshold table defaults to a wide layout** — one
  row per item with a column per threshold plus the mean location, and
  no SE/CI. The optional **long layout** (one row per threshold) adds
  the SE and Wald CI columns at the chosen confidence level, matching
  the intervals drawn in the figure.

Other changes:

- **Reliability**: the "Empirical" row is replaced by **"Marginal"**
  reliability (Green, 1984) — the model-based marginal coefficient is
  complementary to the PSI (a large gap between them serves as an
  off-target diagnostic), whereas the previous
  `mirt::empirical_rxx()` estimate was largely redundant with it. The
  **PSI is now the native WLE-based person separation index** (values
  shift, noticeably for samples with many extreme scores), and the
  **bootstrap yields HDCIs for Cronbach's alpha, PSI, and Marginal
  reliability**, not just alpha.

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
