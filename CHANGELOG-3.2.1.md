# easyRasch2jmv 3.2.1 design notes

Audit responses. No computed value changes, only how often work runs and
when tables become readable.

## Person Change: two computations outside the cache

`changeCache` was written so that the option changes jamovi reruns `.run()`
for, the output-variable toggles and the table switches, cost nothing. Two
calls sat outside it and ran unconditionally.

**The pooled critical value.** `.pooledCrit()` is a second full
`RMpersonChange()` enumeration with `conditional_crit = FALSE`, producing one
number for the reference line on the conditional figure. Its comment said it
runs "only when the conditional figure is shown", which was true and beside
the point: the cost is how often, not whether.

Moved onto the cold path, computed alongside `res` and `plot` and stored as
`cached$pooled`. It depends on nothing outside `sig`, and `cond` is already in
`sig`, so ticking the conditional figure on misses the cache and computes it
then. It stays guarded on `cond` inside the cold path, so a user who never
ticks that figure still pays for two enumerations rather than three.

**The retest SD.** `RMretestSD()` runs `retestIter` simulations (500 by
default) plus an optional `retestBootIter` bootstrap (another 500). Given its
own element, `retestCache`, rather than a second signature inside
`changeCache`: it depends on a different set of options, and sharing one
element would mean each half discarded the other's work. Signature is
`anchor`, `method`, `retestIter`, `retestBoot`, `retestBootIter`, `confInt`,
`seed` and the theta range, with `clearWith` mirroring `retestTable`'s list.
`estimateRetestSd` appears in neither, which is the point.

A failed simulation is cached along with a successful one. The failure is
deterministic given the signature, so rerunning a doomed 500-iteration
simulation on every option change wastes as much as rerunning a good one.

Measured on the 150 x 6 test data with both boxes ticked: 6.4 s cold, 0.7 s
to toggle a table or an output variable. The retest simulation is nearly all
of it at this size (6.8 s with only that cache cold, 0.8 s with only the
enumeration cache cold), though the pooled enumeration scales with items and
respondents where the simulation scales with the iteration count.

## Person Change: two tables built after the work rather than before

`pairingTable` and `retestTable` were `rows: 0` with `addRow()` in `.run()`,
so both sat empty until the enumeration finished.

For the pairing table this undercut its own purpose. It is placed above the
figure because a mispaired analysis produces a plausible-looking result and
the check has to precede the result it validates. As built, it was the last
thing to appear.

Now `rows: (vars1)`, matching the six analyses that already declare rows from
a variable list, with `pair`, `item1` and `item2` filled in `.init()` along
with the pairing note. Only the response categories, which need the data, are
left to `.run()`, which now uses `setRow(rowNo = )`. The table is complete the
moment variables are dropped in.

`.init()` fills the pairing ahead of the two-items-per-occasion and
equal-length guards, rather than after them: a user in the middle of dragging
variables across is exactly who the note is for. When the lists are unequal
`item2` is left empty rather than mispaired.

`retestTable` gets its four fixed rows in `.init()` when the option is on, the
way `summaryTable` already did, and `.fillRetest()` switched to
`setRow(rowKey = )`.

`setRow()` sets only the columns named in `values` and leaves the rest alone
(`jmvcore::Table$setRow` loops over the columns present in `values`), so the
labels written in `.init()` survive the value fill. An earlier claim in
conversation that `setRow()` blanks unnamed columns was wrong: what looked
like a wiped note was an empty string, which `asDF` returns as `NA`.

## Reliability: the figure was recomputed on every redraw

`.curvePlot()` called `RMreliabilityCurve()` itself, with
`boot = isTRUE(self$options$curveBoot)`. Its comment said the figure was
"redrawn here rather than stored, so that only the figure pays for the
bootstrap band", which held for the first draw and no other: jamovi calls a
render function on every resize and every export.

Measured on the bundled 50 x 5 data: 0.05 s with no band, **2.6 s at the
default 200 iterations, 25 s at the 2000 maximum**, each time.

The audit suggested storing the curve data frame and adding an entry point in
easyRasch2 that turns a precomputed frame into the ggplot. That is the tidier
long-term shape, but it needs a package release, and 1.3.1 is on CRAN with no
development version open. The module can fix this alone by storing the
**ggplot object**, which is what Person Change already does with its own
figure, and needs no upstream change or duplicated plot code.

So the cold path now builds the figure alongside the summary and caches both,
and `.curvePlot()` only applies the theme bump and prints. A cold run costs
two curves, one cheap frame for the table and one bootstrapped figure, which
is the same total as before plus 0.05 s. Every redraw, export and cache hit
afterwards is free rather than 2.6 s.

`curve_sig` gained `boot`, `boot_iter`, `conf_int` and `seed`, all three of
the latter gated on the band being on: with it off none of them reaches the
figure or the summary, and an ungated `conf_int` threw the curve away
whenever the HDCI width changed. A test caught that, having been written for
3.2.0 to pin that the HDCI width moves the estimates and not the curve.

Storage was the thing worth checking before choosing this route. The stored
state is ~3.9 MB raw and ~0.36 MB gzipped, and it stays there: 300 x 10 and
1000 x 20 both come out at the same size, because the cost is the ggplot
object's fixed overhead rather than the data. Verified end to end through the
same serialise, gzip, `readRDS` path jmvcore uses for element state, and the
restored object renders.

## Results state: ggplot objects are much bigger than they look

Raised as a check rather than a diagnosis, correctly: it needed measuring.
jmvcore compresses element state and prints `WARNING: state object for <path>
is too large` above 500,000 bytes, confirmed in
`ResultsElement$asProtoBuf()`. Measured on 500 respondents and 20 items:

| element | before |
|---|---|
| `personchange` `changeCache` | **1181 KB** (the figure alone 1159) |
| `personfit` `simCache` | **1022 KB** (three figures, 998) |
| `reliability` `curveCache` | 374 KB (the figure 359) |
| `reliability` `curvePlot` | 359 KB (the same figure again) |

Two over the threshold, and the reliability pair only under it by storing one
figure in two elements, which the 3.2.1 redraw fix had just introduced.

The audit attributed the size to `plot_env` capturing the package function's
scope. That is not what it is here. Dissecting the person-change figure: no
component of the object exceeds 2 KB on its own, `plot_env` holds 15 objects
totalling 29 KB, and pointing `plot_env` at `globalenv()` saves nothing. The
cost is the ggplot2 object itself. **An empty ggplot in ggplot2 4.0.3
serialises to 102 KB compressed**, so every cached figure carries that floor
before any data.

That changes the remedy. Storing what the plot draws and rebuilding it, as
suggested, means reproducing easyRasch2's plot code in the module, which is
the thing the caching was there to avoid. Building the figure once and
storing the **grob** gets the same win with no duplicated code:

| figure | as ggplot | as grob |
|---|---|---|
| person change | 418 KB | 16 KB |
| person fit, infit | 287 KB | 14 KB |
| person fit, outfit | 287 KB | 14 KB |
| person fit, lz | 752 KB | 12 KB |
| reliability curve | 359 KB | 16 KB |

New `er2_plot_grob()` and `er2_draw_grob()` in `R/utils-theme.R`. The theme
bump moves from the render function into the build, since a built grob can no
longer take a theme, which is fine because `er2_bump_text()` uses fixed sizes
and does not depend on the render context. `grid` added to `Imports`.

After: `changeCache` 39 KB, `simCache` 62 KB, `curveCache` 32 KB, `curvePlot`
17 KB. All five figures were then rendered end to end through the module's own
render functions and checked to produce real output, not blank pages.

**Building a figure needs a graphics device, and `.run()` has none.**
Reported from jamovi: rendering the reliability curve opened a second jamovi
entry in the launch sidebar and an empty window titled "Quartz 2 [*]".
`ggplotGrob()` measures text, and with no device open ggplot2 opens the
platform default rather than failing, which on macOS is a visible Quartz
window; under `Rscript` the same bug shows up as a stray `pdf` device left
open after `.run()`. `er2_plot_grob()` now opens a file-less `pdf(NULL)` for
the build and closes it again, guarded on `dev.list()` being empty so a
render function's own device is never touched. `grDevices` added to
`Imports`. A test asserts that `dev.list()` is unchanged across a run of both
affected analyses.

Building against a stand-in device does freeze the layout against that
device's font metrics, which is worth knowing since jamovi draws to a
different one. Measured: the gtable's 13 layout column widths differ by at
most 0.18 pt between a pdf build and a png build, so it makes no visible
difference.

Worth keeping in view: most analyses in the module already store data and draw
in the render function, and `locdepq3`'s comment says so in as many words.
Only these three store figures, because their figure comes back from the same
expensive package call that produces the numbers.

## Images declared requiresData

`requiresData: true` makes jamovi read the whole dataset before calling a
render function, on every redraw and every export. All 27 `Image` elements
declared it. A sweep of the render functions in every `.b.R` found exactly one
that reads `self$data`: `.treePlot` in `R/diftree.b.R`, on its documented
fallback path. Removed from the other 26.

## Labels

Five checkboxes added in 3.2.0 opened with an action verb. jamovi's
convention, and the rest of this module after the 06-13 audit, is that a
checkbox names the thing it turns on.

| Was | Now |
|---|---|
| Draw the flat marginal summary as a reference line | Marginal reliability as a reference line |
| Show the respondent distribution behind the curve | Respondent distribution behind the curve |
| Shade the region reaching a reliability benchmark | Region reaching the reliability benchmark |
| Label categories with their value labels | Value labels as category labels |
| Estimate the retest SD from these data | Retest SD estimated from these data |

A sweep of every `Bool` title in every `.a.yaml` against a list of 24 leading
verbs found no others, confirming the audit's "everything else is already
right". Recorded as a convention in memory so it is applied when options are
written rather than when an audit catches them.

## Tests

Tampered-value probes for both caches, on the pattern used elsewhere: a
tampered number reaching the output proves the cached object was reused. Each
is paired with an assertion in the other direction, that an option the
computation genuinely depends on forces a real rerun. Plus `an$init()` called
on its own, without `run()`, to assert both tables are complete and labelled
before any computation, and that an unequal pairing shows the Time 1 side
rather than nothing.

`retestCache` added to the `clearWith` declaration test from 3.2.0.
