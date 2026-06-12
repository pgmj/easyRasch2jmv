# easyRasch2jmv 2.0.0 — complete change log (from 1.0.0)

This release is the result of a function-by-function consistency review
covering UX, analytical options, and documentation detail. This document
holds the full record; NEWS.md carries the summary.

## New analysis: Partial Gamma Local Dependence (`locdepgamma`)

Based on `easyRasch2::RMlocdepGamma()`. Tests each item pair for residual
association via `iarm::partgam_LD()`, controlling for the rest score in
both directions (two tables). Options:

- Number of pairs to display (top-N by |gamma|; 0 = show all).
- Show only significant pairs (BH-adjusted p < .05).
- Show only pairs with |gamma| at or above a threshold.
- Sort by |gamma| descending.
- Show SE and 95% CI (off by default).

The simulation-based cutoff part of the R-package function is not
included. Both tables carry a footnote explaining the BH abbreviation.
Requires at least 3 items (with 2 items the rest score reduces to the
other item and the statistic is undefined).

## Module-wide standardizations

- **Observed/Expected column-label convention**: short metric name in
  the column title, full name in the table title/footnotes. Conditional
  infit: "Observed infit"; Q3 pair table: "Observed Q3"; partial gamma
  DIF: "Observed gamma"; item-restscore: "Observed gamma" / "Expected
  gamma" (were "Observed value"/"Expected value"). Simulation interval
  columns are titled Lower/Upper under an **"Expected range"**
  super-title; asymptotic CI columns under a **"95% CI"** super-title.
- **Misfit/flag vocabulary settled**: "overfit"/"underfit" for the
  infit and restscore analyses (footnotes state each rule — note the
  direction is inverted between the two families: high restscore
  difference = overfit, but *low* infit = overfit); "above"/"below" for
  the Q3 pair table; "TRUE" for partial gamma DIF / CFA / LR-DIF flag
  columns.
- **HDCI width defaults unified at 99%** (previously 99.9 in both
  conditional infit analyses, 99.5 planned for Q3, 99 in partial gamma
  DIF). This changes default flagging in the infit analyses.
- **BH p-value adjustment hardcoded everywhere**: the adjustment-method
  dropdown in item-restscore was removed (Holm/Bonferroni no longer
  selectable); columns titled "Adj. p-value (BH)"; footnotes spell out
  the Benjamini-Hochberg false-discovery-rate method and the star
  coding (*** p < .001, ** p < .01, * p < .05, . p < .10, adjusted).
- **Random seed always applied** (default 42). Previously Q3 and the
  CFA cutoff treated `seed = 0` as "unseeded" (undocumented);
  simulation results are now reproducible by default everywhere.
- **Raw values + `zto` formatting**: no analysis pre-rounds values any
  more, so jamovi's Number format preferences apply; `zto` formatting
  added where missing.
- **Informative minimum-item messages**: every analysis now explains
  its minimum item requirement in the output instead of returning
  silently (3 items for the restscore/LD/infit/DIF families, with the
  statistical reason stated; 2 for targeting and the ICC plot; 4 for
  CFA per lavaan's identification requirement).
- **Sample-size / missing-data notes**: every analysis reports n and
  how missing responses are handled. Where two samples are involved the
  note says so precisely — e.g. iarm refits on complete cases for
  restscore/infit statistics while eRm CML retains partially missing
  rows for item locations.

## Simulation robustness guard (module-wide)

All simulation-based analyses (conditional infit, partial gamma DIF,
Q3, bootstrap item-restscore) now require **at least 20 successful
iterations and a 50% success rate**. Previously a single successful
iteration could silently produce degenerate cutoffs (HDCI lower =
upper, flagging every item) or percentages over tiny denominators.

- Conditional infit, partial gamma DIF, Q3 degrade gracefully: the
  observed statistics are still shown, and the note reports the
  dominant failure reason (e.g. "fewer than 8 positive responses in at
  least one item").
- Bootstrap item-restscore stops with an informative error (the
  bootstrap is the entire analysis; there is no observed-only
  fallback).

## Sparse-category warnings (module-wide)

Every item-based analysis checks whether any item has a response
category with fewer than 3 observations (the same trigger as the
targeting fallback and upstream easyRasch2) and shows a footnote naming
the affected items ("Estimates may be unstable"). The DIF analyses
(partial gamma DIF and LR-test DIF) check counts **within each DIF
group** — including categories a group never uses — and point users to
the response-distribution tileplot, so the check no longer depends on
the tileplot being displayed. Shared helpers live in `utils-sparse.R`.

## Per-analysis changes

### Item-restscore

- New **Misfit** column: "overfit" (observed > expected) / "underfit"
  (observed < expected) when adjusted p < .05 — same rule and labels as
  the bootstrap analysis.
- Removed the **Location** column; only **Rel. location** remains
  (footnote defines it).
- Sample-size note documents the dual-sample situation with missing
  data: restscore statistics use complete cases (iarm refits the
  model), item locations use all available responses (eRm CML). iarm's
  refit console message is suppressed.
- Footnotes for the Difference column (sign interpretation), BH/stars,
  and Rel. location.
- Requires >= 3 items with an explanatory message (with 2 items the
  restscore reduces to the other item; observed = expected by
  construction, p = 1).
- Raw values; `zto` on all numeric columns; `compilerMode: tame`;
  harmonized variable-type suggestions; `rejectInf`.

### Bootstrap item-restscore

- Removed the **% No misfit** column (redundant: always 100 minus the
  other two) and the **Cond. infit MSQ** column.
- The yes/blank **Flagged** column became a **Misfit** column
  ("overfit"/"underfit" when that classification share exceeds the
  display cutoff; if both exceed it, the more frequent one wins).
- New **Sort items by** dropdown: none (data order), % overfit, or
  % underfit (descending). The plot follows the table's item order.
- Plot: y-axis flipped to **observed − expected** (overfit positive,
  matching the asymptotic analysis); `coord_flip()` with items on the
  y-axis, first item on top (matching the infit plot); text size
  matches the infit plot (base_size 15); larger jitter points; caption
  documents the colour coding (blue = overfit, red = underfit, grey =
  no misfit) and the sign direction.
- HTML note documents that rows with missing responses can be drawn
  into bootstrap samples but are excluded by the within-iteration
  refit, shrinking the effective per-iteration n; iarm's per-iteration
  console message is suppressed.
- Footnotes for the classification rule, the Misfit rule, and Rel.
  location. Requires >= 3 items with an explanatory message.

### Conditional item infit

- **Flagged → Misfit** column: "overfit" = infit *below* the expected
  range (item more predictable than the model expects), "underfit" =
  above. The footnote notes the direction is inverted relative to the
  restscore analyses.
- Footnote documents the dual-sample situation (complete-case infit vs
  all-rows locations) and defines Rel. location.
- Requires >= 3 items (with 2 items conditional infit equals 1 by
  construction); sort option added to `clearWith`; plot x-axis breaks
  no longer hardcoded to 0.5–1.5; plot caption says "diamonds".

### Q3 residual correlations

- **New "Q3 by Item Pair" table**, shown automatically (with the
  per-pair figure above it) whenever the simulation-based cutoff is
  enabled: Observed Q3 per pair, an "Expected range" (per-pair HDCI of
  the values simulated under local independence; new **HDCI width**
  option), and a Flagged column marking pairs "above" (stronger
  residual association than the model predicts; positive LD) or
  "below" (weaker; can indicate multidimensionality). Pairs sort by
  deviation from the simulated per-pair median, matching the figure.
- The global cutoff summary table moved up to sit directly below the
  correlation matrix; both now state that they use a **global** cutoff
  derived from all item pairs, in contrast to the per-pair intervals
  in the pair table; footnote explains the percentiles refer to
  (max Q3 − mean Q3) across simulated datasets; the duplicate **p99**
  row was removed (identical to the suggested cutoff).
- **Observed Q3 precision increased from 2 to 4 decimals**
  (`mirt::residuals(digits = 4)`), matching the simulated values; the
  mean Q3, dynamic cutoff, and flags now use full precision. The same
  fix was applied upstream in easyRasch2 (`RMlocdepQ3()`).
- Sample-size note: N, complete cases, that MML retains partially
  missing rows, and the mean Q3. The per-pair plot caption no longer
  claims "complete responses".
- Requires >= 3 items (mirt cannot estimate the 2-item model).

### Partial gamma DIF

- **SE and 95% CI columns are opt-in** ("Show SE and 95% confidence
  interval", default off), matching partial gamma LD; CI columns under
  a "95% CI" super-title. The note states the CI construction (Wald,
  gamma ± 1.96 SE, verified against iarm).
- **Plot draws the observed gamma's 95% Wald CI** as an orange segment
  behind each diamond, documented in the caption; caption wording
  corrected to "diamonds".
- Always-visible note reporting n complete cases and exclusions
  (missing item responses or missing DIF value).
- Footnotes: BH + star coding; expected-range/Flagged rule including
  how the null is generated (DIF variable randomly reassigned with
  preserved group proportions).
- Requires >= 3 items (with 2 items `partgam_DIF` returns
  mirror-duplicate rows of a single test).
- Tileplot gained a "Note." caption stating group sizes and the
  highlight cutoff; `sortByGamma` added to `clearWith`; fixed a stray
  "no sink to remove" warning (double sink removal).

### Targeting plot

- **CML-to-MML fallback ported from `easyRasch2::RMtargeting()`**:
  when any item response category has fewer than 3 observations, item
  thresholds (and SEs) are estimated with MML (mirt), which is more
  numerically stable under sparse categories. Person locations always
  come from the eRm model. The HTML note states which method was used
  and why.
- The threshold table gained an **SE** column and `zto` formatting.
- Requires >= 2 items with an explanatory message; the plot caption's
  n counts respondents with at least one valid response.

### Item probability curves (ICC plot)

- **Rebuilt as a ggplot, porting `easyRasch2::RMitemCatProb()`**:
  faceted per item, viridis category colours, theta-range and legend
  options retained, "Note." caption with n and model. Curve values are
  numerically identical to the upstream function (verified).
- **Model is now auto-selected** (RM for dichotomous, PCM for
  polytomous); previously PCM was always used. The old base-graphics
  `eRm::plotICC(empICC = "raw")` call never actually drew its
  empirical overlay — eRm only supports `empICC` for dichotomous RM
  objects and silently warned every run — so no displayed output is
  lost by the rewrite.
- **Dichotomous data uses a joint-ICC design** (cf.
  `eRm::plotjointICC()`): all items in a single panel, one curve per
  item showing P(response = 1), item legend. The faceted category view
  is redundant for dichotomous items (two mirror-image curves per
  panel) and is now used for polytomous data only. Note this is a
  deliberate presentation divergence from `RMitemCatProb()`, which
  facets dichotomous data too.
- New always-on HTML note (n, model, CML, missing-data handling, with
  the sparse-category warning folded in); explanatory message when
  fewer than 2 items are selected.

### Sum score to logit transformation

- Corrected the sample description in the footnote, HTML note, and plot
  caption: item parameters are fitted on **all** rows (eRm CML / mirt
  MML retain partially missing responses), not on complete cases as
  previously claimed. Output now reports N respondents with the
  complete-case count in parentheses.
- When the EAP method fails because the MML model is not estimable
  (very few items), the error now suggests using WLE or selecting more
  items instead of only showing mirt's "too few degrees of freedom".
- Requires at least 2 items with an explanatory message; values are no
  longer pre-rounded; `zto` added to the SE column; `rejectInf` added;
  plot text size aligned at base_size 15.

### Reliability

- New sample-size note documents the mixed-samples situation with
  missing data: Cronbach's alpha is closed-form on complete cases,
  while PSI / Empirical / RMU retain partially missing rows (eRm CML /
  mirt MML). The previous note incorrectly stated that incomplete rows
  were "excluded from model fitting".
- New table footnote defining the HDCI and stating which estimates
  carry intervals (bootstrapped alpha and RMU) vs point estimates only
  (PSI, Empirical).
- **Bootstrap-CI robustness**: the alpha CI now requires at least 20
  usable resamples and a 50% success rate (same thresholds as the
  simulation analyses); otherwise the Notes column reports how many
  resamples were usable. Previously 2 usable resamples sufficed.
- Requires at least 3 items with an explanatory message (with 2
  dichotomous items the mirt model is not estimable; 2-item
  reliability is generally not informative).
- Estimates are no longer pre-rounded; `rejectInf` added to the
  variable filter.

### CFA cutoff

- Since `lavaan::cfa()` requires at least 4 items for an
  over-identified model, an explanatory message is shown when fewer
  are selected (with 3 items the model is just-identified and every
  fit index is perfect by construction).
- **Simulation robustness guard** (same thresholds as the other
  analyses): if fewer than 20 iterations or under 50% succeed, the
  observed fit indices are still shown (the lavaan fit is independent
  of the simulation) without cutoffs, and the note reports the
  dominant failure reason.
- The table's three index rows (CFI/RMSEA/SRMR) are pre-filled in
  `.init()` so the table renders meaningfully before computation
  finishes; new footnote states the flagging direction (CFI below the
  cutoff; RMSEA/SRMR above).
- Cleanup: removed an unreachable estimator check (the dropdown only
  offers WLSMV/ULSMV); error texts no longer mention DWLS; the cutoff
  percentile bounds in code now match the option's 50-99.9 range
  (previously exactly 50 was accepted by the option but rejected by
  the code).
- Random seed always applied (see module-wide changes).
- Sparse-category warning evaluated on the complete cases actually
  analysed.

### LR-test DIF

- Sparse-category warning within DIF groups (see module-wide changes),
  pointing to the response-distribution tileplot.
- Table footnotes added: MaxDiff definition (highest minus lowest group
  location, All column excluded) and the flagging rule.
- The Andersen LR-test note now reports how many rows were excluded
  (missing item responses or DIF value; eRm::LRtest requires complete
  cases).
- The panel figure dropped its duplicated in-plot titles (the jamovi
  element already carries the title); the "item locations are means of
  threshold locations" explanation moved into the caption (shown for
  polytomous data at item level).
- Tileplot gained the "Note." caption with group sizes and the
  highlight cutoff (same as partial gamma DIF); the partial gamma DIF
  tileplot adopted lrdif's refinements in return (wrapped facet labels
  for long group names, same text sizes) -- the two DIF tileplots are
  now identical in style, sharing the `er2_wrap_labels()` helper.
- Requires at least 2 items with an explanatory message.
