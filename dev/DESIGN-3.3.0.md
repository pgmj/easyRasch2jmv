# Design: easyRasch2jmv 3.3.0

Status: **proposed 2026-09-30, no code written** beyond the item-restscore
calibration footnote and the version bump, which the user asked for directly.

Brings the module onto the easyRasch2 development version after 1.3.1, whose
main addition is `RMitemRestscoreCutoff()`: a parametric bootstrap null for
the item-restscore test, with Westfall-Young or FDR corrected p-values.
Background in `easyRasch2/dev/restscore-cutoff-design.md` and
`easyRasch2/dev/restscore_asymptotic_null.qmd`.

## Blocker: the package dependency

`DESCRIPTION` pins `Remotes: pgmj/easyRasch2@238b130...`, and
`RMitemRestscoreCutoff()` so far exists only as uncommitted code in the local
easyRasch2 repository. Nothing below can be built until the package is
committed and pushed and the pin is moved. Then:

- `Imports: easyRasch2 (>= 1.3.1.9000)` or the next release number
- `Remotes:` to the new commit
- `0000.yaml` description: "Built against easyRasch2 ..." updated

## What else the new pin brings in

Moving the pin brings every change in the package's development version, not
only the new function. From the package NEWS:

| package change | module effect |
|---|---|
| `RMlocdepQ3Cutoff()` resample DGP keeps missingness patterns | **Q3 cutoffs and p-values change (wider) for data with missing responses.** Complete data identical. Needs a module NEWS entry. |
| Cutoff sample-size check | None. The module always passes the same data to the cutoff and the consumer. |
| Reliability latent moments 2.5 times faster | Reliability analysis faster, identical values. |
| CFA 1.4 to 2 times faster | CFA analysis faster, identical values. |
| `RMitemRestscore()` kable caption | None, the module uses the dataframe output. |

## Added 2026-10-04: Partial Gamma DIF update (easyRasch2 1.3.1.9001)

Done, uncommitted. Background in `easyRasch2/dev/dif-impact-design.md`, Part A.

- **Required.** `RMdifGamma()` now defaults `p_value = NULL`, which resolves to
  TRUE with a full cutoff object and drops `padj_bh` / `Significance`. The
  module now passes `p_value` explicitly, gated on `computeCutoff`, as
  locdepgamma does. `Imports: easyRasch2 (>= 1.3.1.9001)`.
- **Inherited.** `RMdifGammaCutoff()` keeps each respondent's group and total
  score (conditional DGP, the package default; the module does not pass
  `dgp`). Footnote, comments and a.yaml description rewritten; the old text
  said the DIF variable is randomly reassigned. The description also said
  "rest score"; partial gamma DIF conditions on the total score.
- **Consistency.** `pValues` (default TRUE) and `correction` options with
  p-value and adjusted p-value columns, as in itemrestscore; `pAdjBoot`
  title set in `.init()` via `padjusted_title()`. Defaults 400 iterations,
  95% HDCI (max 99.9). Notes: `iteration_note(..., 400L, corrected =)`,
  `pvalue_iteration_caveat()` without a floor (the 400 floor is from item
  fit, as for Q3), `interval_flagging_note()` when p-values are off.
- **BH labels.** The package now applies a real Benjamini-Hochberg
  adjustment in `RMdifGamma()` and `RMlocdepGamma()` (iarm's column was
  Bonferroni), so the existing "(BH)" titles and footnotes in partgamdif and
  locdepgamma become correct without module code changes. NEWS entry.
- Permutation null not exposed.
- **Column visibility (2026-10-05).** jmvcore's `Options$eval()` only treats
  a `visible:` expression as a condition when it matches
  `^\([$A-Za-z].*\)$`. A negation such as `(!computeCutoff)` returns the
  string itself, which counts as visible, so the asymptotic columns showed
  empty in itemrestscore, partgamdif and locdepgamma. Rewritten as
  `(computeCutoff == FALSE ...)`, verified in all four option states and
  pinned by a behaviour test. Rule: never start a `visible:` expression
  with `!`.
- **Visibility after a failed simulation (2026-10-05).** `visible:` follows
  the options, so a failed simulation (error or < 20 iterations, which falls
  back to the asymptotic output) showed empty simulation columns and plots
  and hid filled asymptotic columns. New helpers `set_columns_visible()` and
  `set_elements_visible()` (utils-validation.R) set visibility from
  `cutoff_res` / `use_pvalues` once the fallback is known, in iteminfit,
  iteminfitmi, itemrestscore, locdepq3, locdepgamma, partgamdif and
  residualpca. The yaml expressions stay as the initial state. Behaviour
  test mocks the easyRasch2 cutoff functions to fail (all but MI infit).
- Tests: new behaviour test checks gamma, p, adjusted p and flags against
  easyRasch2 for fwer and fdr_bh, the interval path, and a stale pValues.
  Module suite 431 expectations, 1 failure (the known stale reliability-curve
  test), run against easyRasch2 1.3.1.9001 in a temporary library.
- `Remotes:` pin still points at 238b130 (CRAN 1.3.1); move it after the
  package is pushed.

## Proposal: simulation-based p-values in Item-Restscore Correlations

### Where

**Recommended: in the existing Item-Restscore Correlations analysis**, as an
option that mirrors the Conditional Item Infit analysis. The two are the
module's item-fit pair, and `RMitemRestscore(cutoff = )` is the package
counterpart of `RMitemInfit(cutoff = )`.

The alternative is the Bootstrap Item-Restscore analysis. That analysis wraps
`RMitemRestscoreBoot()`, which resamples the data and reports how often each
item is flagged by the asymptotic test. It answers a stability question with a
different mechanism, and putting a parametric null there would mix two
bootstraps in one analysis. **Decision A.**

### Options (mirroring iteminfit.a.yaml)

| name | title | type | default |
|---|---|---|---|
| `computeCutoff` | Simulation-based p-values | Bool | false |
| `hdciWidth` | HDCI width | Number | 95 (50 to 99.9) |
| `iterations` | Number of simulation iterations | Integer | 400 (50 to 5000) |
| `seed` | Random seed | Integer | 42 |
| `pValues` | Bootstrap p-values | Bool | true |
| `correction` | Multiple-comparison correction | List | fwer (fwer, fdr_bh, fdr_by) |

- **Checkbox label.** Infit uses "Simulation-based cutoffs". For item-restscore
  the thing turned on is better described by the p-values, since the interval
  is descriptive only. Following the house label rule (name what is turned on,
  no verb), "Simulation-based p-values" fits, but it differs from infit.
  **Decision B:** match infit's label exactly, or use the more accurate one.
- **`pValues` nested under `computeCutoff`.** Kept for parity with infit,
  where unchecking it flags on the interval. With restscore it is less useful,
  and could be dropped so the checkbox means one thing. **Decision C.**
- **`dgp`** is not exposed, as in infit. Resample is the package default.

### Table

When `computeCutoff` is on:

- `pAdjusted` (asymptotic BH) is hidden, as the package drops it.
- New columns **Lower** and **Upper** under superTitle "Expected range" for
  the difference (house convention "Option B"), then **p-value** and
  **Adj. p-value** with the correction named via `padjusted_title()`.
- `Flagged` from the corrected p-value. Its direction comes from the side of
  the simulated mean the item falls on, which the footnote states.

Footnotes: the asymptotic calibration note is shown only on the asymptotic
path. On the bootstrap path, the same set as infit: iterations and correction
label (`correction_label()`), `iteration_note()` and
`pvalue_iteration_caveat()` for runs below 400 and 1000, and attrition.

### Figure

`RMitemRestscorePlot(simfit, data)` as an Image, visible with
`computeCutoff`, following the infit plot:

- The cutoff object lives in the plot state with a signature check (iterations,
  seed, width) and a `clearWith` list, so option changes that do not touch the
  simulation (sorting, correction, p-values) reuse it. Per the module rule, a
  state without `clearWith` defaults to "*" and never carries.
- Stored as a built grob via `er2_plot_grob()`, as the other figures are.

### Guards

Same as infit: at least 20 successful iterations, else the observed table is
shown and the note explains why the simulation part is missing.

### Cost

About 0.06 to 0.1 s per iteration for 7 to 9 polytomous items, so 25 to 40 s
at the default 400. Infit is about half that. The note should say so, as the
infit note does.

## Timing

**Decision D:** ship 3.3.0 before or after the validation study
(`easyRasch2/dev/restscore_cutoff_validation.qmd`, running 2026-09-30). The
package function is built but its calibration is not yet verified. The
recommendation is to implement now and release after the validation shows
nominal family-wise error.

## Tests

Following `tests/testthat/test-smoke.R` and `test-behavior.R`:

- smoke: polytomous and dichotomous, cutoff on and off, populated table
- behaviour: calibration footnote only on the asymptotic path, `pAdjusted`
  hidden with the cutoff, the cache is reused when only `correction` changes
- identity: module output equals `RMitemRestscore(df, cutoff = RMitemRestscoreCutoff(df, ...))`

Pre-existing and unrelated: `test-behavior.R:830` expects the reliability
figure in `curvePlot$state`, but 3.2.2 moved it to `curveCache`. The test is
stale.
