# Smoke tests: each analysis runs without error and populates results on
# the bundled example datasets. Heavy packages (eRm/mirt/iarm/lavaan/mice)
# are Imports, so they are assumed present rather than skipped. Iteration
# counts are kept small for speed; the fixed default seed (42) makes the
# simulations deterministic.
#
# Results are checked for "produced something" via rowCount / note content
# rather than specific values, which keeps the tests robust to estimator
# updates while still catching the common regression (an analysis erroring
# or returning empty after a refactor).

test_that("itemrestscore runs on polytomous and dichotomous data", {
  expect_no_error(suppressWarnings(
    r <- er2$itemrestscore(data = poly_data(), vars = names(poly_data()))))
  expect_gt(r$restscoreTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$itemrestscore(data = dich_data(), vars = names(dich_data()))))
})

test_that("bootstrap item-restscore runs", {
  d <- dich_data()
  expect_no_error(suppressWarnings(
    r <- er2$bootrestscore(data = d, vars = names(d),
                           iterations = 60, samplesize = 100)))
  expect_gt(r$bootstrapTable$rowCount, 0)
})

test_that("conditional infit runs with and without cutoffs", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    er2$iteminfit(data = d, vars = names(d))))                       # no cutoff
  expect_no_error(suppressWarnings(
    r <- er2$iteminfit(data = d, vars = names(d),
                       computeCutoff = TRUE, iterations = 60)))      # cutoff
  expect_gt(r$infitTable$rowCount, 0)
})

test_that("Q3 residual correlations run with and without cutoff", {
  d <- poly_data()
  expect_no_error(suppressWarnings(er2$locdepq3(data = d, vars = names(d))))
  expect_no_error(suppressWarnings(
    er2$locdepq3(data = d, vars = names(d),
                 computeCutoff = TRUE, iterations = 60)))
})

test_that("partial gamma LD runs (and with cutoff + SE columns)", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$locdepgamma(data = d, vars = names(d), showSE = TRUE)))
  expect_gt(r$dir1Table$rowCount, 0)
})

test_that("reliability runs, incl. bootstrap alpha CI", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$reliability(data = d, vars = names(d),
                         bootAlpha = TRUE, bootIter = 100,
                         draws = 100, rmuIter = 5)))
  expect_equal(r$relTable$rowCount, 4)
})

test_that("person parameters runs for WLE and EAP, incl. score-to-logit table", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$personparams(data = d, vars = names(d), method = "WLE",
                          showScoreTable = TRUE, showFigure = TRUE)))
  expect_gt(r$scoreTable$rowCount, 0)
  expect_gt(r$summaryTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$personparams(data = d, vars = names(d), method = "EAP",
                     showScoreTable = TRUE)))
})

test_that("residual PCA runs with and without cutoff", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$residualpca(data = d, vars = names(d))))
  expect_gt(r$pcaTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$residualpca(data = d, vars = names(d),
                    computeCutoff = TRUE, iterations = 60)))
})

test_that("targeting runs (CML path)", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$targeting(data = d, vars = names(d))))
  expect_gt(r$thresholdTable$rowCount, 0)
})

test_that("item probability curves run for polytomous and dichotomous data", {
  expect_no_error(suppressWarnings(
    er2$iccplot(data = poly_data(), vars = names(poly_data()))))
  expect_no_error(suppressWarnings(
    er2$iccplot(data = dich_data(), vars = names(dich_data()))))
})

test_that("partial gamma DIF runs (and with cutoff)", {
  d <- dif_data()
  expect_no_error(suppressWarnings(
    r <- er2$partgamdif(data = d, vars = dif_items(), difVar = "dif")))
  expect_gt(r$pgdifTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$partgamdif(data = d, vars = dif_items(), difVar = "dif",
                   computeCutoff = TRUE, iterations = 60)))
})

test_that("LR-test DIF runs at item and threshold level", {
  d <- dif_data()
  expect_no_error(suppressWarnings(
    er2$lrdif(data = d, vars = dif_items(), difVar = "dif", level = "item")))
  expect_no_error(suppressWarnings(
    er2$lrdif(data = d, vars = dif_items(), difVar = "dif",
              level = "threshold")))
})

test_that("MI conditional infit runs on data with missingness", {
  d <- with_missing(poly_data())
  expect_no_error(suppressWarnings(
    r <- er2$iteminfitmi(data = d, vars = names(poly_data()),
                         m = 3, maxit = 5)))
  expect_gt(r$infitTable$rowCount, 0)
})

test_that("CFA cutoff runs (>= 4 items)", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$cfacutoff(data = d, vars = names(d), iterations = 60)))
  expect_equal(r$cfaTable$rowCount, 3)
})

test_that("CICC runs on polytomous and dichotomous data, with and without DIF", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$cicc(data = d, vars = names(d))))
  expect_true(!is.null(r$ciccPlot$state))
  expect_match(r$ciccNote$content, "total-score group")
  # all three grouping methods run
  expect_no_error(suppressWarnings(
    er2$cicc(data = d, vars = names(d), method = "width")))
  r_sc <- suppressWarnings(er2$cicc(data = d, vars = names(d),
                                    method = "score"))
  expect_match(r_sc$ciccNote$content, "does not apply")
  expect_no_error(suppressWarnings(
    er2$cicc(data = dich_data(), vars = names(dich_data()))))
  dd <- dif_data()
  expect_no_error(suppressWarnings(
    r2 <- er2$cicc(data = dd, vars = setdiff(names(dd), "dif"),
                   difVar = "dif")))
  expect_match(r2$ciccNote$content, "partial-gamma")
})

test_that("partial gamma LD expected ranges + plot state (new in 2.1.0)", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$locdepgamma(data = d, vars = names(d),
                         computeCutoff = TRUE, iterations = 60)))
  t1 <- r$dir1Table$asDF
  expect_true(all(c("gammaLow", "gammaHigh") %in% names(t1)))
  expect_false(anyNA(t1$gammaLow))
  # the plot reads the cutoff object from the hidden simCache element
  # (which doubles as the simulation cache)
  expect_true(!is.null(r$simCache$state$cutoff_res))
})

test_that("person fit runs on polytomous and dichotomous data", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$personfit(data = d, vars = names(d), iterations = 100)))
  expect_gt(r$summaryTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$personfit(data = dich_data(), vars = names(dich_data()),
                  iterations = 100, statLz = FALSE)))
})

test_that("Martin-Loef test runs, incl. sequential stopping", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$martinlof(data = d, subscale1 = names(d)[1:3],
                       subscale2 = names(d)[4:5], iterations = 150)))
  expect_gt(r$summaryTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$martinlof(data = d, subscale1 = names(d)[1:3],
                  subscale2 = names(d)[4:5], iterations = 150,
                  sequential = TRUE)))
})

test_that("tree-based DIF runs on polytomous and dichotomous data", {
  d <- dif_data()
  expect_no_error(suppressWarnings(
    r <- er2$diftree(data = d, vars = dif_items(), covariates = "dif")))
  expect_true(nzchar(r$difNote$content))
  dd <- dich_data()
  dd$grp <- factor(rep(c("x", "y"), length.out = nrow(dd)))
  expect_no_error(suppressWarnings(
    er2$diftree(data = dd, vars = setdiff(names(dd), "grp"),
                covariates = "grp")))
})

test_that("cicc note explains the grouping rule without claiming a count", {
  d <- dich_data()
  strip <- function(x) gsub("<[^>]+>", "", x)

  q <- suppressWarnings(suppressMessages(
    er2$cicc(data = d, vars = names(d), method = "quantile",
             classIntervals = 4)))
  qt <- strip(q$ciccNote$content)
  expect_match(qt, "aiming for 4 groups")
  expect_match(qt, "the groups either side merge")
  expect_match(qt, "figure caption reports the grouping actually used")

  w <- suppressWarnings(suppressMessages(
    er2$cicc(data = d, vars = names(d), method = "width",
             classIntervals = 4)))
  wt <- strip(w$ciccNote$content)
  expect_match(wt, "4 equal-width intervals")
  expect_match(wt, "contributes no point")
  expect_match(wt, "figure caption reports the grouping actually used")

  # score-level grouping has nothing that can differ, so no pointer
  s <- suppressWarnings(suppressMessages(
    er2$cicc(data = d, vars = names(d), method = "score")))
  st <- strip(s$ciccNote$content)
  expect_match(st, "does not apply")
  expect_false(grepl("figure caption reports", st))
})
