# Unit tests for the iteration-guidance note helpers.

test_that("iteration_note is silent above the default and shown at/below it", {
  expect_identical(er2$iteration_note(300L, 250L), "")
  expect_match(er2$iteration_note(250L, 250L), "publication-ready")
  expect_match(er2$iteration_note(60L, 250L),  "publication-ready")
})

test_that("iteration_note corrected variant names the floor and the final range", {
  # `corrected` covers every analysis that flags on the Westfall-Young
  # p-value: conditional infit and, from module 3.1.0, both local
  # dependence analyses.
  msg <- er2$iteration_note(100L, 200L, corrected = TRUE)
  expect_match(msg, "400")
  expect_match(msg, "1000 to 2000")
  expect_match(msg, "Johansson, 2026")
  # the withdrawn small-sample exception is gone
  expect_false(grepl("detection power", msg))
  # the generic variant carries neither
  expect_false(grepl("400", er2$iteration_note(100L, 200L)))
})

test_that("pvalue_iteration_caveat is two-tier when given a floor", {
  # below the calibrated floor: the correction itself is off
  low <- er2$pvalue_iteration_caveat(200L, floor = 400L)
  expect_match(low, "below the calibrated floor of 400")
  expect_match(low, "liberal")
  # calibrated but not yet reproducible
  mid <- er2$pvalue_iteration_caveat(400L, floor = 400L)
  expect_match(mid, "calibrated")
  expect_match(mid, "seed-dependent")
  expect_false(grepl("liberal", mid))
  # silent once high enough, with or without a floor
  expect_identical(er2$pvalue_iteration_caveat(1000L, floor = 400L), "")
  expect_identical(er2$pvalue_iteration_caveat(2000L), "")
})

test_that("pvalue_iteration_caveat keeps the generic wording without a floor", {
  # Q3 and any other analysis: the 400 floor was measured on item fit
  # statistics, so it is not claimed for them until their own study lands.
  msg <- er2$pvalue_iteration_caveat(200L)
  expect_match(msg, "At least 1000 iterations")
  expect_false(grepl("calibrated floor", msg))
  expect_false(grepl("Johansson, 2026", msg))
})

test_that("iteration_attrition_note explains the shortfall and the remedy", {
  msg <- er2$iteration_attrition_note(391L, 400L)
  expect_match(msg, "9 of the 400 simulated datasets")
  expect_match(msg, "unused response category")
  expect_match(msg, "rest on 391 datasets")
  expect_match(msg, "Increase the number of simulation iterations")
  # no instability clause when the surviving count is still healthy
  expect_false(grepl("may be unstable", msg))
  # but there is one when it is not
  expect_match(er2$iteration_attrition_note(53L, 400L), "may be unstable")
  # silent when nothing was lost
  expect_identical(er2$iteration_attrition_note(400L, 400L), "")
  expect_identical(er2$iteration_attrition_note(400L, NULL), "")
})

test_that("interval_flagging_note reports the familywise rate it implies", {
  # 1 - .95^9 = 37%
  msg <- er2$interval_flagging_note(0.95, 9L)
  expect_match(msg, "familywise error rate of about 37%")
  expect_match(msg, "9 items")
  expect_match(msg, "Bootstrap p-values")
  # 1 - .99^20 = 18%
  expect_match(er2$interval_flagging_note(0.99, 20L), "about 18%")
  # unusable inputs stay silent rather than guessing
  expect_identical(er2$interval_flagging_note(NULL, 9L), "")
  expect_identical(er2$interval_flagging_note(1, 9L), "")
  expect_identical(er2$interval_flagging_note(NA_real_, 9L), "")
  expect_identical(er2$interval_flagging_note(0.95, 0L), "")
})

test_that("interval_flagging_note counts item pairs when asked to", {
  # Pairs grow quadratically in items, so the same width buys a far worse
  # rate: 1 - .95^36 = 84% over the 36 pairs of a nine-item scale.
  msg <- er2$interval_flagging_note(0.95, 36L, unit = "item pairs")
  expect_match(msg, "familywise error rate of about 84%")
  expect_match(msg, "36 item pairs")
  expect_match(msg, "^ Item pairs are flagged")
  # the default is still items, so existing callers are unaffected
  expect_match(er2$interval_flagging_note(0.95, 9L), "^ Items are flagged")
})

test_that("iteration_attrition_note also covers a short run that lost nothing", {
  # Absorbed from the former low_iteration_caveat(), so one helper now serves
  # every analysis and no call site has to branch between the two.
  expect_identical(er2$iteration_attrition_note(100L, 100L), "")
  expect_identical(er2$iteration_attrition_note(250L, 250L), "")
  msg <- er2$iteration_attrition_note(40L, 40L)
  expect_match(msg, "only 40 iterations")
  expect_match(msg, "may be unstable")
  # and without a requested count at all
  expect_match(er2$iteration_attrition_note(40L), "only 40 iterations")
  expect_identical(er2$iteration_attrition_note(250L), "")
  # a run that is both short and lossy reports the loss, with the caveat
  both <- er2$iteration_attrition_note(40L, 400L)
  expect_match(both, "360 of the 400 simulated datasets")
  expect_match(both, "may be unstable")
})

test_that("every reference key an analysis lists is defined in 00refs.yaml", {
  # A key with no entry is dropped silently at build time, so the citation in
  # the note has nothing to resolve against.
  root <- testthat::test_path("..", "..")
  if (!dir.exists(file.path(root, "jamovi"))) root <- "."
  skip_if_not(dir.exists(file.path(root, "jamovi")), "jamovi/ not found")

  defined <- names(yaml::read_yaml(file.path(root, "jamovi", "00refs.yaml"))$refs)
  for (f in list.files(file.path(root, "jamovi"), pattern = "[.]r[.]yaml$",
                       full.names = TRUE)) {
    y <- yaml::read_yaml(f)
    used <- unlist(lapply(y$items, function(it) it$refs))
    expect_true(all(used %in% defined),
                info = paste(basename(f), ":",
                             paste(setdiff(used, defined), collapse = ", ")))
  }
})

test_that("analyses that cite the 2026 preprint list it as a reference", {
  # The iteration and interval-flagging notes name Johansson (2026) in text,
  # so the analyses that can emit them have to carry the reference.
  root <- testthat::test_path("..", "..")
  if (!dir.exists(file.path(root, "jamovi"))) root <- "."
  skip_if_not(dir.exists(file.path(root, "jamovi")), "jamovi/ not found")

  for (a in c("iteminfit", "iteminfitmi", "locdepq3", "locdepgamma")) {
    y <- yaml::read_yaml(file.path(root, "jamovi", paste0(a, ".r.yaml")))
    used <- unlist(lapply(y$items, function(it) it$refs))
    expect_true("johansson2026_cutoffs" %in% used, info = a)
  }
})
