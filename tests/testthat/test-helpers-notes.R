# Unit tests for the iteration-guidance note helpers.

test_that("iteration_note is silent above the default and shown at/below it", {
  expect_identical(er2$iteration_note(300L, 250L), "")
  expect_match(er2$iteration_note(250L, 250L), "publication-ready")
  expect_match(er2$iteration_note(60L, 250L),  "publication-ready")
})

test_that("iteration_note infit variant adds the small-sample exception", {
  msg <- er2$iteration_note(100L, 200L, infit = TRUE)
  expect_match(msg, "detection power")
  expect_match(msg, "Johansson")
  # non-infit variant does not
  expect_false(grepl("detection power", er2$iteration_note(100L, 200L)))
})

test_that("low_iteration_caveat fires below 100 successful, silent at/above", {
  expect_identical(er2$low_iteration_caveat(100L), "")
  expect_identical(er2$low_iteration_caveat(250L), "")
  msg <- er2$low_iteration_caveat(40L)
  expect_match(msg, "only 40 iterations")
  expect_match(msg, "may be unstable")
})
