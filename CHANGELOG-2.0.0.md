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

All simulation/bootstrap analyses require **at least 20 successful
iterations** before cutoffs/percentages are produced. Previously a
single successful iteration could silently produce degenerate cutoffs
(HDCI lower = upper, flagging every item) or percentages over tiny
denominators.

- Conditional infit, partial gamma DIF, Q3, residual-PCA and the CFA
  cutoff degrade gracefully below the floor: the observed statistics
  are still shown and the note reports the dominant failure reason
  (e.g. "fewer than 8 positive responses in at least one item").
- Bootstrap item-restscore stops with an informative error (the
  bootstrap is the entire analysis; there is no observed-only
  fallback).
- **No success-rate requirement.** An earlier version also required a
  50% success rate, which blocked legitimate small-sample runs where
  many simulated datasets fail sparse-category validation even though
  the survivors form a usable null distribution. The rate condition was
  the same "validity heuristic blocking a computable result" pattern
  removed for the n < 30 gate, so it is gone. Instead, when fewer than
  100 iterations succeed the result is shown with a caveat
  (`low_iteration_caveat()`) noting the cutoffs rest on a small null
  distribution and may be unstable — informing rather than blocking.

## Sparse-category warnings (module-wide)

Every item-based analysis checks whether any item has a response
category with fewer than 3 observations (the same trigger as the
targeting fallback and upstream easyRasch2) and shows a footnote naming
the affected items ("Estimates may be unstable"). The DIF analyses
(partial gamma DIF and LR-test DIF) check counts **within each DIF
group** — including categories a group never uses — and point users to
the response-distribution tileplot, so the check no longer depends on
the tileplot being displayed. Shared helpers live in `utils-sparse.R`.

## Test suite

Added a `testthat` (edition 3) suite under `tests/` — the module had
none. Three layers: (1) unit tests for the pure validation/sparse/note
helpers; (2) smoke tests that run every analysis on the bundled example
datasets (polytomous and dichotomous, cutoffs on/off) and check results
are populated, not just error-free; (3) behavioural tests pinning
decisions from this release (too-few-items guard shows a note rather
than erroring and clears once enough items are added; a perfectly
correlated item pair yields a footnote rather than a halt; ICC category
probabilities sum to 1; CFA requires 4 items). The heavy engines are
Imports, so they are assumed present rather than skipped; the suite runs
in ~20 s via `devtools::test()`.

## Shared validation helper

All analyses now call a single `prepare_item_data()` helper
(`utils-validation.R`) instead of carrying a copy-pasted ~45-line
validation block: numeric conversion, all-NA column check,
sentinel-value check (> 20 looks like an unmarked missing-value code),
response validation, per-item variation check, and the
**identical-items check (correlation = 1), which previously existed in
only four analyses and now applies module-wide**. Complete-case
handling intentionally remains per-analysis (it legitimately varies:
retained by CML/MML, dropped for lavaan/prcomp/LRtest, joint with a
DIF variable).

## Iteration guidance

The simulation-based notes now state that **more iterations are
generally recommended for publication-ready results** — with the
documented exception for **conditional infit**, where with small
samples around 100 iterations can yield better detection power than a
larger number (Johansson, 2025); the infit notes say so explicitly.
The advice is shown only when running **at or below the analysis
default** (a user who lowered the count needs it even more; one who
raised it is spared the nag) and is generated by a shared
`iteration_note()` helper. The differing per-analysis defaults (Q3
100, infit 200, most others 250, MI-infit 500 total) are deliberate,
reflecting very different per-iteration costs.

## Analysis descriptions

All analyses now have a `description` block in their option
definitions (added for conditional infit, item probability curves,
targeting, partial gamma DIF, and partial gamma LD), summarising the
method, estimation choices, and — where relevant — the iteration
guidance above.

## Removed the small-sample (n < 30) gate

The "Warning: Only N complete cases found. Results may be unreliable."
check — present in 13 analyses via `jmvcore::reject()` — has been
removed. `reject()` is fatal (it halts the analysis and greys out the
results pane), so despite the caveat-style wording it produced **no**
output for any sample under 30, contradicting the message. More
fundamentally, n = 30 is a generic small-sample heuristic, not a
computational limit: the real adequacy floor depends on targeting and
analysis type (and is generally higher), so a fixed gate was both too
strict (blocking a computable 29-case analysis) and too lenient
(passing an underpowered 35-case DIF analysis). Sample-size adequacy
is left to the analyst. Genuinely fatal floors are unchanged
(`n_complete == 0`; the per-analysis minimum-item guards), and any
unworkably small n now surfaces the estimator's own error rather than a
blanket block. The cell-size (sparse response category) warning is
kept — it *is* tied to computational feasibility and is already a
non-fatal footnote. (jamovi library audit, 2.0.0.)

## Stale "requires at least N items" message cleared on re-run

The minimum-item guard messages set an HTML note and returned early.
When the user then added enough items the guard was skipped, and the
stale "requires at least N items" text lingered — because jamovi clears
HTML content neither via `clearWith` nor via `setContent("")` (an
empty-string update is a no-op in the frontend); the only thing that
replaces a note is a later write of *non-empty* content. Twelve of the
analyses already rewrite their note with a real sample-size/result
caption on every successful run, so they cleared correctly. Conditional
item infit did not, when simulation-based cutoffs were off (its note is
otherwise written only on the cutoff path), so it now always writes a
brief base note in that case. The message disappears as soon as enough
items are selected.

## Checkbox labels as noun phrases

The 25 checkbox / toggle labels across 12 analyses were rephrased from
action-verb instructions ("Show ...", "Compute ...", "Use ...", "Sort by
...", "Bootstrap ...") to noun phrases naming the feature or resulting
state ("Classification plot", "Simulation-based cutoffs", "Robust
statistics ...", "Sorted by ...", "Bootstrap CI (Cronbach's alpha)"),
following jamovi's UI convention that a checkbox names the thing being
toggled on. (jamovi library audit, 2.0.0.)

## Identical-items check no longer halts the analysis

When three or more items are selected and two are perfectly correlated
(r = 1, e.g. a coincidentally duplicated or re-entered item), the
analysis previously stopped via `jmvcore::reject()` despite the
caveat-style wording — so it produced no output. Since all estimators
(eRm, iarm, mirt) compute fine with such a pair, this is now a
**non-fatal footnote** (`duplicate_items_note()`, parallel to the
sparse-category footnote): the analysis runs and the table/note flags
the correlated pair. The genuinely-degenerate two-item case (where the
correlated pair leaves effectively one item) keeps its fatal `stop()`
in `prepare_item_data()`. (jamovi library audit, 2.0.0.)

## Output text: literal Unicode symbols

Math/Greek symbols in results notes now use literal Unicode characters
(≤, ≥, ±, ×, χ) instead of HTML named entities (`&le;`, `&ge;`,
`&plusmn;`, `&times;`, `&chi;`). The named entities only rendered by
side-effect of jamovi's current results renderer and already failed on
non-HTML export (showing literal `&plusmn;` etc.); a documented upcoming
renderer fix would have made them display literally everywhere. The
genuinely HTML-special `&lt;`/`&gt;` escapes are kept as-is. (jamovi
library audit, 2.0.0.)

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

### Conditional item infit (multiple imputation)

- **Flagged → Misfit** column ("overfit"/"underfit", same rule and
  inverted-direction footnote as the standard infit analysis).
- **Pooled-SE caveat footnote**: the within-imputation variance
  component is iarm's asymptotic infit SE, which Mueller (2020) showed
  can be unreliable; the footnote directs fit decisions to the
  simulation-based expected range instead.
- **Rel. location column added**: the pooled item location was computed
  in every run but never displayed.
- **Simulation robustness guard** on the stacked distribution (at
  least 20 successes and a 50% success rate across imputations), with
  graceful degradation: pooled estimates remain visible and the note
  explains the failure.
- Column retitled "Pooled infit" (Option B convention); requires at
  least 3 items with an explanatory message; values no longer
  pre-rounded; `zto` on all numeric columns; `sortByInfit` added to
  `clearWith`; plot caption says "diamonds"; plot x-axis breaks no
  longer hardcoded; `rejectInf` added.

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

### PCA of standardized Rasch residuals

- **Seed standardization completed**: this analysis still had the old
  `seed = 0` = "unseeded" behaviour, overlooked in the earlier sweep;
  the seed (default 42) is now always applied.
- **Simulation robustness guard** with graceful degradation: when fewer
  than 20 iterations or under 50% succeed, the eigenvalue table and
  variance partition are still shown and the note reports the dominant
  failure reason (previously only total failure was caught).
- Eigenvalue rows are pre-created in `.init()` (the row count --
  min(nComponents, number of items) -- is options-derivable).
- Requires at least 3 items with an explanatory message (with 2 items
  the residual PCA has a single possible contrast, carrying no
  information about multidimensionality).
- Values are no longer pre-rounded; `zto` formatting added to all
  numeric columns; `rejectInf` added.

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
