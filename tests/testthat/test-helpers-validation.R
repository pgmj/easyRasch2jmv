# Unit tests for the shared data-preparation / validation helpers.

test_that("prepare_item_data converts and returns a clean numeric data.frame", {
  raw <- data.frame(A = c(0, 1, 2, 1), B = c(1, 1, 0, 2))
  out <- er2$prepare_item_data(raw, c("A", "B"))
  expect_s3_class(out, "data.frame")
  expect_true(all(vapply(out, is.numeric, logical(1L))))
  expect_identical(names(out), c("A", "B"))
})

test_that("prepare_item_data converts SPSS-style factor responses", {
  raw <- data.frame(
    A = factor(c("0", "1", "2", "1")),
    B = factor(c("1", "1", "0", "2"))
  )
  out <- er2$prepare_item_data(raw, c("A", "B"))
  expect_true(all(vapply(out, is.numeric, logical(1L))))
  expect_equal(out$A, c(0, 1, 2, 1))
})

test_that("prepare_item_data errors on an all-NA item", {
  raw <- data.frame(A = c(0, 1, 2), B = c(NA, NA, NA))
  expect_error(er2$prepare_item_data(raw, c("A", "B")),
               "no valid numeric data")
})

test_that("prepare_item_data flags sentinel-like values (> 20)", {
  raw <- data.frame(A = c(0, 1, 2, 999), B = c(1, 0, 2, 1))
  expect_error(er2$prepare_item_data(raw, c("A", "B")),
               "missing-value codes")
})

test_that("prepare_item_data errors on an item with no variation", {
  raw <- data.frame(A = c(1, 1, 1, 1), B = c(0, 1, 2, 1))
  expect_error(er2$prepare_item_data(raw, c("A", "B")),
               "no variation")
})

test_that("prepare_item_data stops on exactly two identical items", {
  raw <- data.frame(A = c(0, 1, 2, 1, 0), B = c(0, 1, 2, 1, 0))
  expect_error(er2$prepare_item_data(raw, c("A", "B")),
               "identical")
})

test_that("prepare_item_data does NOT stop with >= 3 items incl. a duplicate", {
  raw <- data.frame(A = c(0, 1, 2, 1, 0, 2),
                    B = c(0, 1, 2, 1, 0, 2),   # identical to A
                    C = c(1, 0, 2, 2, 1, 0))
  expect_no_error(er2$prepare_item_data(raw, c("A", "B", "C")))
})

test_that("identical_item_pairs detects perfect correlation and nothing else", {
  dup  <- data.frame(A = c(0, 1, 2, 1), B = c(0, 1, 2, 1), C = c(2, 0, 1, 1))
  none <- data.frame(A = c(0, 1, 2, 1), B = c(1, 0, 2, 0), C = c(2, 0, 1, 1))
  expect_match(er2$identical_item_pairs(dup), "'A' and 'B'")
  expect_identical(er2$identical_item_pairs(none), character(0))
})
