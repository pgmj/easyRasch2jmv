# easyRasch2jmv 3.2.0, detailed record

Companion to the NEWS entry, which carries the summary. Design record and
rationale in `dev/DESIGN-3.2.0.md`. Previous released version: 3.1.0. The
3.1.1 section was never released and is folded in here.

## Person Change (new analysis)

Files: `jamovi/personchange.{a,r,u}.yaml`, `R/personchange.b.R`,
`tests/testthat/test-personchange.R`, helper `change_data()` in
`tests/testthat/helper-er2.R`.

### Input handling

- **Two variable boxes paired by position.** `vars1` and `vars2` must be the
  same length; an unequal pair produces an explanatory note rather than an
  error, as with the other too-few-items guards.
- **Both occasions are prepared as one frame.** `prepare_item_data()` decides
  the 1-based recode from the minimum of every item it is given, so preparing
  the occasions separately could shift one and not the other and silently
  destroy comparability. Names are made unique with `make.unique()` first, so a
  variable dragged into both boxes does not collide, then split back.
- **Both frames are renamed to the Time 1 variable names before the package
  call.** `RMpersonChange()` requires the two occasions to carry identical item
  names in the same order, which two sets of jamovi variables never do. The
  pairing table records what each name stands for.
- **Rows with no responses at one occasion are dropped up front.** The bundled
  easyRasch2 does not survive them: its CML fit raises "subscript out of
  bounds" rather than returning `NA`. Reported upstream, see below. The note
  distinguishes rows that answered at one occasion only.
- Fewer than 10 respondents with data at both occasions is a hard stop naming
  the reason, since the single-respondent path is unreachable in jamovi.

### Statistics

- `critical = "exact"` is the only path and is not exposed as an option.
  Simulation and the normal approximation are both omitted: the normal cutoff
  is wrong in a way a dropdown cannot explain, and simulation buys nothing once
  enumeration is both fast and not approximate. 1.96 appears nowhere in the
  module.
- Consequently the analysis has no iteration count and no seed on its main
  path. The seed lives in the retest-SD section, the only part that simulates.
- `conditional_crit` is exposed, default `FALSE`. When on, the summary reports
  the median and range of the per-respondent values rather than one pair.

### Output

- Order: pairing table, figure, summary, conditional-critical-value figure,
  per-respondent table, retest SD table, note. **The pairing table sits above
  the figure**, deliberately breaking the figure-first layout of the Targeting
  Plot: nothing in Targeting can be silently mis-specified by the order of a
  variable box, and a mispaired analysis here draws a plausible scatter.
- Per-respondent table off by default, with a flagged-only sub-toggle.
- Seven output variables, following Person Parameters and Person Fit.
- The Html note carries the three things a jamovi user has no man page for:
  which null is in force and that a detected change is a necessary condition
  rather than evidence, that change is in logits while the index is unitless,
  and that respondents must not be ranked by the index.

### Caching

`self$results$changeCache` mirrors the `simCache` element in Person Fit. The
enumeration does not depend on the output-variable toggles or the table
switches, which jamovi nonetheless reruns `.run()` for, and a signature check
makes reuse self-validating rather than resting on `clearWith` alone. The
figure is cached alongside the result because `RMpersonChange()` has no output
mode returning both, so a cold run costs two enumerations and a warm one none.

Measured cost of one enumeration, sequential: n = 500 with 10 items, 1.2 s;
n = 1000 with 20 items, 13.4 s.

### Conditional critical values figure

Module-drawn rather than delegated, since the package has no such view yet.
Plots each respondent's critical value against their null location, which is
reconstructed from the public output as the precision-weighted mean of
`theta_t1` and `theta_t2`, matching `easyRasch2:::.pc_theta_null()`. No
internals are touched.

A curve rather than a histogram because the critical value is a deterministic
function of `(theta_null, key1, key2)`: with complete data every respondent
lies on one curve, and a histogram would show that curve marginalised over the
sample's targeting with no way to tell items from sample. Three details:

- the step line is grouped on the pair of answered-item sets, since joining
  across patterns draws interleaved step functions as one zigzag;
- `geom_step`, because the value is a quantile of a discrete distribution and
  slanted connectors draw transitions that do not exist;
- lines are dropped past three distinct pattern pairs, points only. Ragged
  missingness produced 44 pairs at n = 250 in testing.

## Reliability

- `.init()` row label `"Marginal"` becomes `"Marginal (curve mean)"`; the
  `context` note, the `a.yaml` description and the module description are
  rewritten to the new formula. A `marginalchange` note states that values
  moved in 3.2.0 and points at the curve for the superseded coefficient.
- `RMreliabilityCurve()` folded in as an optional figure plus summary. First
  put in its own collapsed section, following the
  `RMscoreSE()`-in-Person-Parameters precedent, then **promoted to the main
  panel**: the curve is a headline result rather than an estimation detail,
  and a collapsed box makes it easy to miss. `showCurve` now sits alongside
  the bootstrap checkbox with its sub-options nested under it, which is the
  same nesting `bootAlpha`/`bootIter` already use in this panel. Default
  still off.
- **The summary is computed with `boot = FALSE` regardless of the option.**
  None of the reported quantities depends on the bootstrap, and the figure is
  redrawn in the render function, so computing both with the band would run the
  only expensive part twice for one figure.
- `benchmark_range` is a data.frame of `xmin`/`xmax`, not a length-2 vector,
  and is `NULL` when nothing reaches the benchmark. The qualifying region can
  be disjoint, so every row is reported.
- The HDCI grouped heading names the width in force (`95% HDCI`), set with
  `getColumn()$setSuperTitle()` in `.init()` since the width is an option.
- The curve summary gained **Mean theta** and **SD theta**, from
  `RMpersonParameters()` on the same data, labelled as in Person Parameters.
  These describe where the sample sits, which is what the density behind the
  curve draws, and are deliberately kept separate from the latent SD. The
  spread of WLE estimates is wider than sigma because each estimate carries
  measurement error, and a test asserts that ordering.

  There is no latent mean to report because `.latent_sd()` does not estimate
  one: `.estimate_prior_sd()` optimises sigma alone with the prior mean held at
  0. **This is an assumption, not a consequence of centring the thresholds**,
  which fixes the item mean and says nothing about where the persons sit. An
  earlier draft of the row note claimed otherwise and was wrong.

  The consequence is reportable and is now stated in an `offtarget` note:
  sigma absorbs any distance between the sample and the item locations, so
  marginal reliability **rises** as targeting worsens while PSI falls.
  Measured on 8 items, n = 1500, true person SD 0.90: shifting the sample by
  +2 logits took sigma from 0.99 to 2.02, marginal from .810 to .902 and PSI
  from .761 to .600. Freeing the mean recovers the true SD at every shift
  (0.906, 0.910, 0.958). The existing PSI-versus-Marginal guidance survives
  and is stronger than stated, since the gap opens from both ends, but the
  marginal value itself is optimistic off-target. Recorded as an upstream
  issue.

  **Fixed upstream in easyRasch2 1.3.1 during this same work**, so the module
  requires that version and the caveat became unnecessary. Both the mean and
  the SD of the latent distribution are now fitted, marginal reliability falls
  with targeting as it should, and the curve summary reports the estimated
  **Latent mean (logits)** beside the latent SD. The note on the main
  reliability table now explains the PSI-Marginal gap rather than warning
  about a bias.

  Worth keeping in mind for the release notes: a user coming from 3.1.0 sees
  marginal reliability move for two reasons at once, the formula change (up)
  and the mistargeting fix (down, and only if their sample is off target).
- The curve summary also gained **Average test information**, computed as
  `1 / sem_average^2`. That is the package's own convention: `RMreliabilityCurve()`
  draws exactly this as the flat reference line on the information axis. It is
  the average SEM restated rather than a separate average over the grid, so the
  two rows cannot drift apart, and a test pins them together.
- The superseded Green row is a migration aid with a shelf life and should be
  removed once reconciling against 3.1.0 output stops being a live concern.
  Recorded in memory with the accompanying pieces to take with it.
- **Cell notes are kept to a few words.** The first draft explained each row
  in the Notes column, up to 496 characters on the latent mean. jamovi styles
  every table cell `white-space: nowrap` (`.jmv-results-table-cell` in the
  results view stylesheet), so a cell note is one unbreakable line and a long
  one drags the table far past the width of the results pane. There is no way
  to break a line inside a cell either: the client sanitises all displayed
  text and keeps only `<i>`, `<em>`, `<b>`, `<strong>`, `<sub>` and `<sup>`,
  so `<br>` arrives escaped and shows as literal text. Footnotes are a
  different matter, rendered one per row with no `nowrap`, so they wrap.
  The explanations therefore live in the footnotes and the cells carry a
  short tag each, matching the ≤27-character notes the package returns for
  the Reliability Estimates table. A test caps cell notes at 40 characters.
- Footnotes are split one topic per note (`context`, `views`, `extremes`,
  `marginal`, then the conditional `bench` and `notest`). Each renders as its
  own line under the table, which is the only line-breaking mechanism jamovi
  offers here.
- Row labels shortened to `Marginal (curve mean)` and
  `Marginal (Green, superseded)`, matching the wording in the Reliability
  Estimates table above rather than restating "reliability" in a column whose
  every row is one.
- **Two hidden cache elements, `relCache` and `curveCache`.** jamovi reruns
  `.run()` on every option change, so folding the curve into this analysis
  meant that ticking it on recomputed `RMreliability()` in full, mirt
  plausible values and bootstrap included, for a table none of the curve
  options can move. The cache splits on what each half depends on: `rel_sig`
  covers `estim`, `draws`, `rmuIter`, `confInt`, `bootAlpha`, `bootIter`,
  `seed` and the theta range; `curve_sig` covers the statistic, benchmark,
  reference line, density, `nNodes` and the theta range. The data frame is
  compared once for both. Changing the HDCI width therefore recomputes the
  estimates and reuses the curve; changing the benchmark does the reverse.
  When the curve is switched off the curve pieces are carried forward
  untouched, so switching it back on is free.

  Measured on the bundled 5-item polytomous data with the default 1000
  bootstrap iterations: 17.2 s cold, then 0.06 s to tick the curve on, 0.03 s
  to change its statistic or add a benchmark. Before the cache every one of
  those cost the full 17.2 s again.

  **`clearWith` is what keeps a cache alive, and omitting it is a silent
  bug.** The first version of this used one element declaring no `clearWith`,
  on the reasoning that the signature check in `.run()` already decides
  validity. It does, but it never gets the chance: `jmvcore::Html$new()`
  defaults `clearWith` to `"*"`, and `ResultsElement$fromProtoBuf()` returns
  early without restoring the state when any option has changed. The cache
  was therefore cold on every run, which is exactly the case it exists for.
  Reported from jamovi, where ticking the curve on visibly re-estimated the
  first table.

  Two elements rather than one, each declaring the options its own half
  depends on, so the split survives: a `clearWith` union on a single element
  would drop the curve whenever the HDCI width changed. Neither list names
  `showCurve`.

  `changeCache` in Person Change had the same omission and is fixed with it.
  The effect there was only speed, since `.run()` recomputes and the figure
  reads whatever is current, but the cache was equally dead.

  A test asserts the declarations directly, by reading `.clearWith` off the
  generated elements. The behavioural test could not catch this: `setState()`
  writes straight past `fromProtoBuf()`, so an in-process test sees a cache
  that works no matter what the yaml says.

  The WLE person locations behind Mean/SD theta are cached with the curve,
  not the estimates. They are a separate `RMpersonParameters()` fit and
  depend only on the data and the theta range.

## Targeting Plot

- New `panel` option, defaulting to the package's new `"categories"`. The
  previous dot-and-whisker view is `"thresholds"`.
- New `rowGap`, `viridisOption`, `viridisBegin`, `viridisEnd`, all passed
  straight through.
- New `categoryLabels`, which recovers response-category labels from the jamovi
  factor levels via `shared_category_labels()` in `R/utils-validation.R`. The
  helper is deliberately narrow: labels are taken only from factors whose
  levels cannot be parsed as numbers, which is the one case where
  `to_numeric_responses()` maps responses to 0-based factor positions and the
  level order therefore is the category order. `haven_labelled` vectors are
  keyed by value rather than position, so a shift applied by
  `prepare_item_data()` would misalign them, and numeric-string levels add
  nothing over the category scores. Every item must agree, and there must be
  one label per observed category. Anything else falls back to the scores and
  says so in the note.

  **Unverified end to end.** The `vars` option is `permitted: numeric`, and
  jmvcore rejects a text-labelled factor before the analysis runs, so the
  recovery path could not be exercised from R. The helper is unit-tested
  directly and the fallback is clean, but whether jamovi ever hands the module
  a variable this fires on needs checking in the app.

## Dynamic CFA Fit-Index Cutoffs

`run_observed_cfa_fit()` **stays**, contrary to the open question in the
design. The upstream 1.3.0 fix makes `RMdimCFACutoff()` and `RMdimCFA()`
tolerate non-syntactic item names, which covers the delegated main path, but
the module helper exists for a different reason: `RMdimCFA()` refuses to run
without a simulated reference distribution, so the observed-fit-only fallback
has no package equivalent. Its placeholder renaming is incidental.

## Figure captions

Captions were 10 pt when the module drew them and 9 pt when easyRasch2 did.
Both are now **10.5 pt**, from `er2_caption_size()` in `R/utils-theme.R`.

The catch: ggplot2 refuses to merge two theme elements of different classes,
and easyRasch2 draws captions with `ggtext::element_markdown()` when `ggtext`
is installed so that the italic "*Note.*" prefix can sit beside roman body
text. Handing a package figure a plain `element_text()` therefore fails with
"Only elements of the same class can be merged". `er2_caption_element()` runs
the same `requireNamespace()` test the package does and returns the matching
class.

Two consequences:

- **`ggtext` is added to `Imports:`**, and documented in `R/utils-imports.R`
  alongside the other bundling pins. It sits in easyRasch2's `Suggests:`, so
  without this declaration jmvtools would not bundle it and a jamovi user with
  no personal R library would get plain captions while everyone else got
  markdown ones. The appearance is now the same for everyone.
- **patchwork figures take the caption in two passes.** Panels that blank their
  caption cannot accept a text element, so `er2_bump_text()` applies the text
  bump per panel with `&` and the caption size to the assembled plot.

## Packaging

- `Version` 3.2.0, `Date` 2026-09-13, in `DESCRIPTION` and `jamovi/0000.yaml`.
- `Imports: easyRasch2 (>= 1.3.0)`.
- `Remotes:` still pins a GitHub commit: jamovi builds against a frozen CRAN
  snapshot rather than tracking CRAN, so easyRasch2 has to come from GitHub
  despite being on CRAN.
- **`grDevices` is not added.** It is new in easyRasch2 1.3.0 and called only
  by the new targeting band panel, but it is base-priority, so it ships with
  every R and is never bundled into a `.jmo`, and it sits in easyRasch2's
  `Imports:` rather than `Suggests:`, so the bundling rule documented in
  `R/utils-imports.R` does not apply. No module code calls it.
- No other new dependency. easyRasch2's `Suggests:` is unchanged between 1.2.0
  and 1.3.0, and the new or changed functions reach only for packages already
  in the module's `Imports:` (`scales` for the viridis palettes, `ggdist` for
  the density overlay, `mirai` only on the parallel path the module never
  takes).
- `ggtext` added to `Imports:`, see Figure captions above.
- Six references added to `jamovi/00refs.yaml`: `jacobsontruax1991`,
  `caronni2026`, `maassen2004`, `zumbo2026` for Person Change;
  `mcneishdumas2025`, `milanzi2015` for the reliability curve. 40 entries, all
  with links.

## Tests

347 passing. New: `test-personchange.R` (16 tests) and a 3.2.0 block in
`test-behavior.R` covering the marginal relabel, the curve fold and its
agreement with the table row, the benchmark row, the targeting panel default,
and `shared_category_labels()`.

Note for future tests: `Output` options cannot be set from R. They are driven
by the jamovi UI and are absent from the generated wrapper function, and
passing `TRUE` to the generated `Options$new()` leaves the value `FALSE`. The
row-alignment logic that output variables depend on is therefore asserted
through the default respondent IDs, which are built from the same index.

## Upstream observations for easyRasch2

Neither is fixed. (A third, the latent mean held at 0, was fixed in easyRasch2
1.3.1 as part of this work.)

1. **All-NA rows crash `RMpersonChange()`.** A respondent with no responses at
   one occasion makes the CML fit raise "CML estimation failed (subscript out
   of bounds)". This is the same class as the all-NA bug fixed in 1.0.0.9000
   for the other functions, which `.drop_empty_respondents()` addressed; the
   two-occasion entry points appear not to be covered.

2. **The pooled critical value is far slower than the conditional one.** At
   n = 1000 with 20 items, pooled takes 13.4 s against 0.8 s conditional, so
   the default path is about 17 times slower than the option. Pooling
   concatenates every respondent's grid and takes one weighted quantile over
   roughly n × (R+1)² values; with complete data all respondents share the same
   grid and only the weights differ, which is not currently exploited.
