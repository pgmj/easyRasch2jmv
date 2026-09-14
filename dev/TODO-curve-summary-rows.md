# Deferred: rows to move between the curve summary table and the package

Status: **noted 2026-09-14, not implemented.** Counterpart to
`easyRasch2/dev/TODO-reliability-curve-annotations.md`, which carries the
full comparison.

The module's Conditional Precision Summary was built against the attributes
of `RMreliabilityCurve()` rather than against the package's own
`.curve_kable()`, so the two tables have drifted. Two rows the package has
and the module does not:

- **Minimum SEM (logits)**
- **Theta at minimum SEM**

Both come straight off `curve_df` (`which.min(curve_df$sem)`), so neither
needs anything the module does not already hold.

One row the package has that the module demotes:

- **Theta range reaching the benchmark.** The module carries this in the
  benchmark row's cell note (`theta -0.20 to 0.21`) rather than as a row of
  its own. A row reads better and would match the package. Note the cell
  cannot hold much: jamovi table cells do not wrap, so a disjoint range with
  several stretches will not fit as a note but would as a row.

Worth doing in the same pass as whichever easyRasch2 release adds the four
rows going the other way (latent mean, mean and SD theta, average test
information), so the two tables land in line rather than drifting again.
