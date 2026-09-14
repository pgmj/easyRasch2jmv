# easyRasch2jmv 3.2.0 design

Design record for bringing the module onto easyRasch2 1.3.0. Written
2026-09-13, decisions agreed in the same session.

**BUILT 2026-09-13.** This file is the design; `CHANGELOG-3.2.0.md` is the
as-built record and takes precedence where they differ. Four things the build
changed or settled are noted inline below.

Scope is the 1.2.0 to 1.3.0 delta in easyRasch2 NEWS.md, plus the jamovi
library audit items already drafted in the unreleased 3.1.1 section of the
module NEWS.

## Version and NEWS

3.1.1 is unreleased. 3.1.0 is what users have. Shipping new analyses and two
results-changing items inside a maintenance patch would misrepresent it, so
**this release is 3.2.0** and the 3.1.1 NEWS section is folded into it.

Per the house convention, an unreleased version's NEWS describes the final
state against the last released version and never records intra-development
churn. So the 3.2.0 section is written against 3.1.0 as a whole, with the audit
items and the 1.3.0 items in one narrative rather than two. The detailed record
goes in a `CHANGELOG-3.2.0.md` alongside the existing `CHANGELOG-2.0.0.md`.

3.2.0 confirmed 2026-09-13, over the argument for 4.0.0 on the strength of the
two results changes. Against 4.0.0: 3.1.0 carried a larger results change
(flagging moved to corrected bootstrap p-values across three analyses) and was
a minor bump.

`jamovi/0000.yaml` picks up analysis descriptions from the `a.yaml` files when
`jmvtools::prepare(".")` runs, but its module-level `version`, `date` and
`description` still need editing by hand.

## 1. New analysis: Person Change

`easyRasch2::RMpersonChange()`, new `personchange` analysis under the existing
**Person-level** menu subgroup, beside Person Parameters, Person Fit and
Reliability. All four are per-respondent statistics off the same CML/WLE
pipeline. A separate "Change over time" subgroup would hold one analysis and
fragment a menu that was only just grouped.

### Input, and the problem it creates

The package takes `data_t1` and `data_t2`: same items, same order, row *i* the
same person at both occasions. jamovi has one rectangular sheet, so the
analysis needs **two variable boxes, "Time 1 items" and "Time 2 items", paired
by position in the list**, plus an optional ID variable box.

Consequences to design around:

- **Pairing is by order, not by name.** Items dragged into the two boxes in
  different orders produce a silent wrong answer that looks plausible.
  Mitigation is a **pairing table** (`PHQ1_t1` with `PHQ1_t2`, one row per
  pair) placed **above** the figure, so the check precedes the result it
  validates. Settled 2026-09-13.

  **Found during the build:** `RMpersonChange()` requires the two occasions to
  carry identical item names in the same order, which two sets of jamovi
  variables never do. Both frames are renamed to the Time 1 names before the
  package call.
- **Wide format only.** Long format cannot be reshaped inside jamovi. State it
  in the `a.yaml` description and refuse cleanly when the boxes differ in
  length.
- **Validation runs over the pair**, not per box: equal item counts, matching
  category ranges within each pair, and the 1-based recode rule applied
  identically across both occasions rather than independently per box.
- **Missingness has two levels.** `prepare_item_data()` drops all-NA
  respondents, but a respondent can be all-NA at one occasion only, which kills
  the change without killing their contribution to the calibration. The caption
  needs a three-part count: rows total, rows contributing to the calibration,
  rows with a testable change.

  **Found during the build:** this is not a nicety. An all-NA row at one
  occasion makes the package's CML fit raise "subscript out of bounds" rather
  than returning `NA`, so such rows must be dropped module-side. Reported as an
  upstream observation in the changelog.

### Out of scope for jamovi

`item_params` has no UI, so the **single-respondent case cannot be reached**.
That is the headline use in the package docs (one clinician, one patient,
external calibration), and jamovi gets only the sample case via `anchor`. Say
so in the analysis description rather than letting users find it through an
`n < 10` error.

A later route would be a Variables box holding an item-parameter table produced
by the Targeting analysis. Separate project, not this release.

### Options (`personchange.a.yaml`)

| Option | Package argument | Notes |
|---|---|---|
| `vars1`, `vars2`, `id` | `data_t1`, `data_t2`, `id` | paired by order |
| `anchor` | `anchor` | stack (default) / t1 / t2 |
| `method` | `method` | WLE / EAP, wording copied from Person Parameters |
| (fixed CML) | `estimator` | no dropdown, matching the rest of the module |
| `nullType` | `null` | "Measurement error only" / "Measurement plus retest fluctuation" |
| `retestSd` | `retest_sd` | visible under `nullType == retest`, labelled **per occasion** |
| `alpha`, `direction` | same | direction labelled in theta terms, not clinical terms |
| `conditionalCrit` | `conditional_crit` | default FALSE, see section 2 |
| `estimateRetestSd` | (runs `RMretestSD()`) | see section 3 |
| `seed` | `seed` | lives inside the retest-SD section only |
| `thetaMin`, `thetaMax` | `theta_range` | house pattern |
| `showTable`, `flaggedOnly` | (display) | per-respondent table, default off |

**No `critical` control.** `critical = "exact"` is the only path. The normal
approximation is documented in the package as wrong in a direction that matters
most on short scales (1.62 at four items, 1.70 at six), and a user picking it
from a dropdown has no man page telling them so. Simulation earns nothing here
either: exact enumeration is fast, it is not approximate, and it is not a
fallback path, since `.pc_exact_null()` groups respondents by answered-item set
pair and degrades linearly rather than exploding. 1.96 is not mentioned
anywhere in the module. The caption states the critical value in force.

**No iterations, no seed on the main path.** Without simulation the analysis is
deterministic. Say so in the description, because every other simulation-based
analysis in the module has spent three versions teaching users the opposite.
The seed reappears only inside the retest-SD section, which does simulate.

### Results (`personchange.r.yaml`)

Figure-led, as in the Targeting Plot, with one deliberate deviation: the
pairing table sits above the figure rather than below it. The Targeting
precedent puts the image first, but nothing there can be silently mis-specified
by the order of a variable box, and a mispaired analysis here produces a
perfectly plausible-looking scatter. The check goes before the result.

1. **pairingTable**, which item went with which. Settled 2026-09-13.
2. **changePlot**, Image, always visible, no gate. `RMpersonChange(output =
   "ggplot")` restyled the usual way. The figure is the default output.
3. **summaryTable**, counts and percentages by `change_class`, the critical
   values in force, null, anchor, alpha, calibration n.
4. **critPlot**, Image, visible under `conditionalCrit` only. See section 2.
5. **changeTable**, one row per respondent, **off by default**
   (`visible: (showTable)`, `showTable` default false) with a flagged-only
   sub-toggle. Columns id, sum_t1, sum_t2, theta_t1, theta_t2, change, se_diff,
   rci, p_value, change_class, crit_lower, crit_upper, retest_sd_tip. Settled
   2026-09-13: the figure plus the saved output variables carry the common
   case, and a 300-row table on first run is noise. Single-subject readers turn
   it on.
6. **retestTable**, visible under `estimateRetestSd`.
7. **changeNote**, Html.

**Output variables into the dataset**, following Person Parameters and Person
Fit: theta_t1, se_t1, theta_t2, se_t2, change, se_diff, rci, p_value,
change_class.

### Captions and notes

`.pc_caption()` already emits the null in force with its `SE_diff` formula,
then the exact pooled critical values, then the calibration clause naming the
anchor and item count. Keep it rather than building a module caption. The
module adds the `.n_caption()` clause with the three-part count above.

Three things the Html note must carry, because jamovi users have no man page:

- which null is in force and that a flagged result under the measurement null
  is a necessary condition for change rather than evidence of it
- that `change` is logits, `rci` is unitless, and `change_class` is a threshold
  decision
- that respondents must not be ranked by `rci`. The package docs spend a page
  on this misreading and it is likelier in a point-and-click audience.

**Footnote under the figure, `conditional_crit` only:** the no-change band is
drawn at the median critical values, so a point can sit on the wrong side of
the band relative to its row in the table. Read classification from the table.

**Footnote under the figure, always:** with mixed missingness the band is drawn
for the most common pair of answered-item sets.

## 2. Conditional critical values, module-side figure

`conditional_crit` stays as an option, default FALSE. The per-respondent null
is the more defensible one for single-subject reading, since a respondent near
a boundary has a different null from one in the middle, and exact enumeration
makes it free. Default FALSE keeps the headline figure and caption exact.

When it is on, the module draws a second figure: **the per-respondent critical
value against that respondent's null location**, not a histogram. Prototype and
rendered examples in `dev/mockup_conditional_crit.R`.

Why the curve rather than a histogram of the same numbers: `crit` is a
deterministic function of `(theta_null, key1, key2)`, so with complete data
every respondent sits on one curve, and a histogram shows that curve
marginalised over the sample's targeting with no way for a reader to separate
item behaviour from sample spread. The curve also survives n = 25, which the
clinical use needs and a histogram does not, and it shares an x-axis with the
Reliability Curve going into the same release.

Measured on six items, four categories, n = 250: an arch peaking near 1.84
around theta -1 and falling to about 1.08 at both ends, against a pooled 1.73.

Three details the prototype settled:

1. **Group the step line on the pair of answered-item sets.** Joining across
   patterns draws two interleaved step functions as one zigzag that reads as
   noise.
2. **`geom_step`, not `geom_line`.** The critical value is a quantile of a
   discrete distribution, so it is a step function of the null location and
   slanted connectors draw transitions that do not exist.
3. **Drop the lines past three distinct pattern pairs**, points only. Ragged
   missingness produced 44 pairs at n = 250.

`theta_null` is reconstructed from the public output as the precision-weighted
mean of `theta_t1` and `theta_t2`, matching `.pc_theta_null()`, so no internals
are touched. The answered-item grouping needs the two response matrices, so the
image state has to carry both alongside the result.

**Summary table rows report median and range, not median and IQR.** On the
complete-data case the IQR was 1.65 to 1.78 while the range was 1.08 to 1.84.
The IQR says the critical values barely move and compresses away the tails,
where the whole effect lives.

**Upstream:** the same view belongs in `RMpersonChange()` in the next
easyRasch2 version, where it would also serve the single-subject case jamovi
cannot reach. Noted for easyRasch2 development. It needs a decision on how it
reaches the user, since `output = "ggplot"` already returns the occasion
scatter.

## 3. RMretestSD, folded not separate

`RMretestSD()` becomes a checkbox-gated section inside Person Change rather
than a nineteenth menu entry. It takes the same two datasets and the same
anchor and method arguments, and produces one number feeding `retest_sd` in the
same panel.

The label must say it is valid only when the two occasions are a stability
study with no expected change. Estimating it from trial data and feeding it
back is circular.

Cost is the catch: `sim_iter = 500` plus `boot_iter = 500`, both sequential
under jamovi, inside an analysis that is otherwise instant. The section carries
its own iteration and seed controls and its own run-time warning.

## 4. Reliability, results change (required)

Marginal reliability changes formula in 1.3.0 and **moves upward**, more so on
short scales. It is now the latent-density-weighted mean of
`sigma^2 / (sigma^2 + SEM(theta)^2)` rather than Green's subtractive
coefficient.

Four places name the old quantity and all four change to "Marginal (curve
mean)" and the new description:

- the hard-coded label in `.init()`, `R/reliability.b.R:24`
- the `context` table note, same file
- the `description.main` block in `jamovi/reliability.a.yaml`
- the analysis description in `jamovi/0000.yaml` (regenerated from the a.yaml)

`green1984` stays in the `refs` list only while the superseded coefficient is
still mentioned. NEWS must state that results move and roughly by how much
(phq9 .862 to .886).

The PSI-versus-Marginal gap sentence stays valid unchanged.

## 5. Reliability Curve, folded into Reliability

`RMreliabilityCurve()` draws the picture of the number the table now reports,
so it goes inside the Reliability analysis as an optional figure plus summary
rather than a new menu entry. Same move Person Parameters already makes with
`RMscoreSE()`.

It also gives somewhere to surface `marginal_green`, the superseded value, for
users comparing against results produced under 3.1.0.

New options: `statistic` (sem / information / reliability), `benchmark`,
`reference`, `show_density`, `boot` with `boot_iter`, `n_nodes`. That is six
more on a panel that already has nine, so it needs its own collapsible section.

Overlap to declare in the note: the score-to-logit table in Person Parameters
carries the same message on the score metric.

## 6. Targeting, breaking change

The bottom panel is now response-category bands, with thresholds as a coloured
error-bar row beneath. Categories never most likely, and threshold reversals,
are marked in red. Estimates are unchanged.

**Follow the package default** and add a "Bottom panel" dropdown offering the
previous dot-and-whisker view via `panel = "thresholds"`. Also expose
`row_gap` and the three viridis options.

**`category_labels` is the jamovi-specific win.** jamovi variables carry factor
level labels, so the module can populate the argument from the data instead of
leaving users to type them. The R package cannot do this for its own users.

## 7. CFA cutoff cleanup

**Answered during the build: the workaround stays.** The upstream fix covers
the delegated main path, but `run_observed_cfa_fit()` exists for a different
reason: `RMdimCFA()` refuses to run without a simulated reference distribution,
so the observed-fit-only fallback has no package equivalent. Its placeholder
renaming is incidental to that purpose.

## 8. Packaging

- `Imports: easyRasch2 (>= 1.3.0)`, done manually.
- **`Remotes: pgmj/easyRasch2@<commit>` stays.** jamovi builds against a
  fixed-in-time CRAN snapshot rather than tracking CRAN, so easyRasch2 has to
  come from GitHub despite being on CRAN.
- **`grDevices` does not need adding.** It is genuinely new in easyRasch2
  1.3.0 (1.2.0 Imports were graphics, knitr, mirt, psychotools, stats, utils,
  rlang) and is called in three places, all inside the new targeting
  category-band panel: `.ci_contrast()` and `.band_label_colour()` in
  `targeting_plot.R`. The module still does not declare it, for two independent
  reasons. It is base-priority, so it ships with every R and is never bundled
  into a `.jmo`. And it sits in easyRasch2's `Imports:` rather than
  `Suggests:`, so the rule written up in `R/utils-imports.R` does not apply,
  that rule existing for hard runtime requirements hiding in `Suggests:` where
  jmvtools will not see them.
- **No other new dependency.** easyRasch2's `Suggests:` is unchanged between
  1.2.0 and 1.3.0, and the three new or changed functions reach only for
  packages the module already Imports: the band panel uses
  `scales::viridis_pal()` (the viridis palettes come from scales, not the
  viridis package), `RMreliabilityCurve()` and `RMretestSD()` use `ggdist`, and
  `RMpersonChange()` uses `mirai` only on the parallel path the module does not
  take.

## 9. References to add to `jamovi/00refs.yaml`

`caronni2026`, `jacobsontruax1991`, `maassen2004`, `zumbo2026` for Person
Change. `mcneishdumas2025` and `milanzi2015` for the Reliability Curve. Every
entry needs a link, per the 3.1.1 audit item.

## 10. Cross-references between analyses

`RMreliabilityCurve(benchmark = )` reports the theta region reaching a given
reliability and the share of respondents inside it. That is the region where
change is detectable. The Person Change note points at it, and the Reliability
note points at Person Change.

## Open decisions

All closed. Version 3.2.0, per-respondent table off by default, pairing table
above the figure, CFA workaround kept.

**One thing shipped but not verified end to end:** `category_labels` in the
Targeting Plot. The `vars` option is `permitted: numeric` and jmvcore rejects a
text-labelled factor before the analysis runs, so the label-recovery path could
not be exercised from R. The helper is unit-tested directly and falls back
cleanly to category scores, but whether jamovi ever hands the module a variable
this fires on needs checking in the app.
