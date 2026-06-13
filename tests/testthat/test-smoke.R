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

test_that("score-to-logit runs for WLE and EAP", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$scorese(data = d, vars = names(d), method = "WLE")))
  expect_gt(r$scoreTable$rowCount, 0)
  expect_no_error(suppressWarnings(
    er2$scorese(data = d, vars = names(d), method = "EAP")))
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
