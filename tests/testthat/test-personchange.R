# Person Change, new in 3.2.0. These tests encode the decisions that the
# analysis would otherwise lose silently: that items are paired by position,
# that the critical value is enumerated rather than taken from a normal
# distribution, and that rows without data at one occasion are dropped
# before they reach a package call that cannot survive them.

test_that("the analysis runs and reports an exact critical value below 1.96", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(), id = "pid"))

  sm <- r$summaryTable$asDF
  crit_hi <- sm$value[sm$statistic == "Critical value (upper)"]
  crit_lo <- sm$value[sm$statistic == "Critical value (lower)"]

  # The whole reason the module exposes no critical-value control: on a
  # six-item scale the enumerated value sits well below the normal one.
  expect_lt(crit_hi, 1.96)
  expect_gt(crit_hi, 1.4)
  # Two-sided with equal missingness is exchangeable, so the bounds agree.
  expect_equal(crit_hi, -crit_lo, tolerance = 1e-8)
})

test_that("the pairing table lists every pair in box order", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  pt <- r$pairingTable$asDF
  expect_equal(nrow(pt), 6L)
  expect_equal(pt$item1, change_t1())
  expect_equal(pt$item2, change_t2())
  expect_true("pairing" %in% names(r$pairingTable$notes))
})

test_that("reversing one box changes the pairing, and the results with it", {
  d <- change_data()
  straight <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  crossed <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = rev(change_t2())))

  # Pairing by position is the analysis's main hazard: a reversed box is
  # accepted and produces a different, plausible-looking answer. The
  # pairing table is the only thing that shows it, which is why it sits
  # above the figure.
  expect_equal(crossed$pairingTable$asDF$item2, rev(change_t2()))
  n_flagged <- function(x) {
    sm <- x$summaryTable$asDF
    sum(sm$value[sm$statistic %in% c("Increase in theta", "Decrease in theta")])
  }
  expect_false(identical(n_flagged(straight), n_flagged(crossed)))
})

test_that("unequal box lengths give a note rather than an error", {
  d <- change_data()
  expect_no_error(suppressWarnings(
    r <- er2$personchange(data = d, vars1 = change_t1(),
                          vars2 = change_t2()[1:5])))
  expect_match(r$changeNote$content, "same number of items")
})

test_that("fewer than two items per occasion gives the guard note", {
  d <- change_data()
  expect_no_error(suppressWarnings(
    r <- er2$personchange(data = d, vars1 = change_t1()[1],
                          vars2 = change_t2()[1])))
  expect_match(r$changeNote$content, "at least")
})

test_that("rows with no data at one occasion are dropped, not passed on", {
  d <- change_data()
  d[1:5, change_t2()] <- NA
  # The bundled easyRasch2 fails its CML fit on an all-NA row ("subscript
  # out of bounds"), so this has to be handled module-side.
  expect_no_error(suppressWarnings(
    r <- er2$personchange(data = d, vars1 = change_t1(),
                          vars2 = change_t2())))
  sm <- r$summaryTable$asDF
  expect_equal(sm$value[sm$statistic == "Respondents tested"], 145)
  expect_match(sm$notes[sm$statistic == "Respondents tested"], "of 150")
})

test_that("partial missingness at one occasion is retained", {
  d <- change_data()
  d[1:20, change_t2()[5:6]] <- NA
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  sm <- r$summaryTable$asDF
  expect_equal(sm$value[sm$statistic == "Respondents tested"], 150)
})

test_that("the class counts add up to the respondents tested", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  sm <- r$summaryTable$asDF
  n <- sm$value[sm$statistic == "Respondents tested"]
  parts <- sm$value[sm$statistic %in% c("Increase in theta",
                                        "No change detected",
                                        "Decrease in theta")]
  expect_equal(sum(parts), n)
})

test_that("a one-sided test leaves the untested bound empty", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(),
    direction = "increase"))
  sm <- r$summaryTable$asDF
  expect_true(is.na(sm$value[sm$statistic == "Critical value (lower)"]))
  expect_match(sm$notes[sm$statistic == "Critical value (lower)"], "not tested")
  # One tail carries the whole alpha, so the bound is lower than two-sided.
  two <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  two_hi <- two$summaryTable$asDF
  expect_lt(sm$value[sm$statistic == "Critical value (upper)"],
            two_hi$value[two_hi$statistic == "Critical value (upper)"])
})

test_that("the retest null widens the critical band and detects less change", {
  d <- change_data()
  meas <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  ret <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(),
    nullType = "retest", retestSd = 0.4))
  flagged <- function(x) {
    sm <- x$summaryTable$asDF
    sum(sm$value[sm$statistic %in% c("Increase in theta", "Decrease in theta")])
  }
  expect_lte(flagged(ret), flagged(meas))
  expect_match(ret$summaryTable$notes$null$note, "occasion")
})

test_that("per-respondent critical values vary and report a median", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(),
    conditionalCrit = TRUE, showTable = TRUE))
  ct <- r$changeTable$asDF
  expect_gt(length(unique(round(ct$critUpper, 6))), 1L)
  sm <- r$summaryTable$asDF
  expect_match(sm$notes[sm$statistic == "Critical value (upper)"], "Median")
})

test_that("the per-respondent table is off by default and filters when asked", {
  d <- change_data()
  off <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2()))
  expect_false(off$changeTable$visible)

  all_rows <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(), showTable = TRUE))
  only_flagged <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(),
    showTable = TRUE, flaggedOnly = TRUE))
  expect_equal(nrow(all_rows$changeTable$asDF), 150L)
  expect_lt(nrow(only_flagged$changeTable$asDF), 150L)
  expect_true(all(only_flagged$changeTable$asDF$changeClass != "none detected"))
})

test_that("the retest SD table reports the estimate without applying it", {
  d <- change_data()
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(),
    estimateRetestSd = TRUE, retestIter = 60, retestBootIter = 60))
  rt <- r$retestTable$asDF
  expect_equal(nrow(rt), 4L)
  expect_true(any(grepl("Retest SD", rt$quantity)))
  # Reported, never silently fed back into the test: the null in force
  # must still be the one shown in the options.
  expect_match(r$retestTable$notes$use$note, "not applied automatically")
  expect_match(r$retestTable$notes$valid$note, "stability study")
  expect_match(r$summaryTable$notes$null$note, "measurement error")
})

test_that("dropped rows keep their original row numbers", {
  # Output variables are written back to the spreadsheet by row number, so
  # the mapping from the analysed subset to the full data has to survive
  # the drop. jamovi drives Output options from the UI and they cannot be
  # set from R, so the mapping is asserted through the default IDs, which
  # are built from the same index.
  d <- change_data()
  d[1:5, change_t2()] <- NA
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(), showTable = TRUE))

  ids <- r$changeTable$asDF$id
  expect_equal(length(ids), 145L)
  expect_equal(ids[1], "6")
  expect_equal(ids[145], "150")
})

test_that("a supplied ID variable is carried through to the results", {
  d <- change_data()
  d[1:5, change_t2()] <- NA
  r <- suppressWarnings(er2$personchange(
    data = d, vars1 = change_t1(), vars2 = change_t2(), id = "pid",
    showTable = TRUE))
  ids <- r$changeTable$asDF$id
  expect_equal(ids[1], "P6")
  expect_equal(ids[145], "P150")
})
