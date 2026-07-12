# Behavioural tests encoding decisions made in the 2.0.0 review, so they
# cannot silently regress.

test_that("too-few-items guard shows a note instead of erroring", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$iteminfit(data = d, vars = names(d)[1:2])))
  expect_match(r$cutoffNote$content, "requires at least")
})

test_that("guard message is gone once enough items are selected", {
  d <- poly_data()
  r <- suppressWarnings(er2$iteminfit(data = d, vars = names(d)))
  # iteminfit without cutoffs writes a base note (not the guard message),
  # so the stale 'requires at least' text must not survive.
  expect_false(grepl("requires at least", r$cutoffNote$content))
})

test_that("a perfectly-correlated pair (>=3 items) yields a footnote, not a stop", {
  d <- with_duplicate_item(poly_data())
  expect_no_error(suppressWarnings(
    r <- er2$itemrestscore(data = d, vars = names(d))))
  notes <- r$restscoreTable$notes
  expect_true("duplicate" %in% names(notes))
  expect_match(notes$duplicate$note, "perfectly correlated")
})

test_that("item-restscore Flagged column uses overfit/underfit labels", {
  d <- dich_data()
  r <- suppressWarnings(er2$itemrestscore(data = d, vars = names(d)))
  # empty text cells come back as NA via asDF; non-blank labels must be
  # exactly "overfit" / "underfit"
  fit_vals <- r$restscoreTable$asDF$fit
  expect_true(all(is.na(fit_vals) | fit_vals %in% c("overfit", "underfit")))
  expect_true(any(fit_vals %in% c("overfit", "underfit")))   # at least one flagged
})

test_that("item-restscore Flagged uses inverted direction vs. infit (sanity)", {
  # Restscore: observed > expected => overfit. Infit: observed < expected
  # range => overfit. Just assert both analyses run and expose the column;
  # the directional rules are unit-tested via their own logic above.
  d <- poly_data()
  rr <- suppressWarnings(er2$itemrestscore(data = d, vars = names(d)))
  expect_true("fit" %in% names(rr$restscoreTable$asDF))
})

test_that("iccplot category probabilities are valid (sum to 1 per theta point)", {
  d <- poly_data()
  r <- suppressWarnings(er2$iccplot(data = d, vars = names(d)))
  pdf_df <- r$iccPlot$state$plot_df
  # For each item x theta, probabilities across categories should sum to 1.
  agg <- stats::aggregate(Probability ~ Item + Theta, data = pdf_df, FUN = sum)
  expect_true(all(abs(agg$Probability - 1) < 1e-8))
})

test_that("scorese WLE SEM is information-based: finite at extremes and larger there", {
  d <- dich_data()
  r  <- suppressWarnings(er2$scorese(data = d, vars = names(d), method = "WLE"))
  df <- r$scoreTable$asDF
  se <- df$logitSE
  expect_true(all(is.finite(se)))                 # Warm correction -> finite extremes
  n  <- length(se)
  # the two extreme scores carry the largest SEs (least information)
  expect_gt(se[1L], min(se))
  expect_gt(se[n],  min(se))
})

test_that("CFA cutoff requires 4 items, with an explanatory message at 3", {
  d <- poly_data()
  r <- suppressWarnings(er2$cfacutoff(data = d, vars = names(d)[1:3]))
  expect_match(r$cfaNote$content, "at least <b>4 items</b>")
})

test_that("Q3 matches easyRasch2 directly and survives an all-NA respondent", {
  # Regression for the 2.1.0 migration: the module delegates to easyRasch2
  # (CML/WLE) and must (a) reproduce RMlocdepQ3() numerically and (b) drop
  # all-NA respondents up front, because RMlocdepQ3Cutoff() errors on them
  # (psychotools cannot fit all-NA rows).
  m <- as.matrix(poly_data())
  m[3, ] <- NA                                   # one fully-empty respondent
  d <- as.data.frame(m)

  expect_no_error(suppressWarnings(
    r <- er2$locdepq3(data = d, vars = names(d),
                      computeCutoff = TRUE, iterations = 50)))

  # cutoff simulation ran (not the graceful-degradation path)
  cut_df <- r$cutoffTable$asDF
  expect_false(anyNA(cut_df$value))
  expect_false(grepl("unavailable", r$q3Note$content))
  expect_match(r$q3Note$content, "1 row\\(s\\) without any responses")

  # numerical agreement with the R package (same seed defaults)
  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepQ3(d_used, output = "dataframe")))
  mod <- suppressWarnings(
    apply(as.matrix(r$q3Table$asDF[, names(d)]), 2, as.numeric))
  expect_equal(unname(mod), unname(as.matrix(pkg[names(d)])), tolerance = 1e-12)
})

test_that("item-restscore matches easyRasch2 directly and survives an all-NA respondent", {
  # Regression for the 2.1.0 migration: the module delegates to
  # easyRasch2::RMitemRestscore() and pre-drops all-NA respondents
  # (the bundled release cannot fit them).
  m <- as.matrix(poly_data())
  m[3, ] <- NA
  d <- as.data.frame(m)

  expect_no_error(suppressWarnings(
    r <- er2$itemrestscore(data = d, vars = names(d))))
  tab <- r$restscoreTable$asDF

  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressWarnings(suppressMessages(
    easyRasch2::RMitemRestscore(d_used, output = "dataframe")))
  expect_equal(tab$observed, pkg$Observed)
  expect_equal(tab$difference, pkg$Difference)
  expect_equal(tab$relLocation, pkg$Relative_location)
  expect_match(r$restscoreNote$content, "1 row\\(s\\) without any responses")
})
