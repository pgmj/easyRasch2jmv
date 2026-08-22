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

# ---------------------------------------------------------------------
# 1-based coding
# ---------------------------------------------------------------------
test_that("is_one_based only fires when every item's minimum is exactly 1", {
  d0 <- data.frame(a = c(0L, 1L, 2L), b = c(0L, 2L, 1L))
  d1 <- data.frame(a = c(1L, 2L, 3L), b = c(1L, 3L, 2L))

  expect_true(er2$is_one_based(d1))
  expect_false(er2$is_one_based(d0))
  # ragged: one item never uses 1, so the reading is not obvious
  expect_false(er2$is_one_based(data.frame(a = c(1L, 2L), b = c(2L, 3L))))
  # a minimum above 1 is left alone, it is more likely truncation
  expect_false(er2$is_one_based(data.frame(a = c(3L, 4L), b = c(3L, 5L))))
  # a gap in the middle does not block the shift, it is handled downstream
  expect_true(er2$is_one_based(data.frame(a = c(1L, 3L), b = c(1L, 2L))))
  # all-NA column cannot be judged
  expect_false(er2$is_one_based(
    data.frame(a = c(1L, 2L), b = c(NA_integer_, NA_integer_))))
})

test_that("prepare_item_data shifts 1-based data and leaves the rest alone", {
  d1 <- data.frame(a = c(1L, 2L, 3L), b = c(1L, 3L, 2L))
  out <- er2$prepare_item_data(d1, names(d1))
  expect_equal(out$a, c(0, 1, 2))
  expect_equal(out$b, c(0, 2, 1))

  # already 0-based: untouched
  d0 <- data.frame(a = c(0L, 1L, 2L), b = c(0L, 2L, 1L))
  expect_equal(er2$prepare_item_data(d0, names(d0))$a, c(0, 1, 2))

  # ragged minima still hit the existing error rather than being shifted
  ragged <- data.frame(a = c(1L, 2L, 3L), b = c(2L, 3L, 2L))
  expect_error(er2$prepare_item_data(ragged, names(ragged)),
               regexp = "scored starting at 0")
})

test_that("recode_note reads the raw data, not a flag on the prepared frame", {
  d1 <- data.frame(a = c(1L, 2L, 3L), b = c(1L, 3L, 2L))
  d0 <- data.frame(a = c(0L, 1L, 2L), b = c(0L, 2L, 1L))
  expect_match(er2$recode_note(d1, names(d1)), "lowest response was 1")
  expect_null(er2$recode_note(d0, names(d0)))

  # survives the frame rebuilds some analyses do straight after preparing
  prepared <- er2$prepare_item_data(d1, names(d1))
  rebuilt <- as.data.frame(matrix(as.numeric(as.matrix(prepared)),
                                  nrow = nrow(prepared),
                                  dimnames = list(NULL, names(prepared))))
  expect_null(attr(rebuilt, "recoded_from_one"))       # an attribute would be gone
  expect_match(er2$recode_note(d1, names(d1)), "lowest response was 1")
})

test_that("the 1-based shift does not change results", {
  d0 <- poly_data()
  d1 <- as.data.frame(lapply(d0, function(x) x + 1L))
  a <- suppressWarnings(suppressMessages(
    er2$iteminfit(data = d0, vars = names(d0))))
  b <- suppressWarnings(suppressMessages(
    er2$iteminfit(data = d1, vars = names(d1))))
  expect_equal(a$infitTable$asDF, b$infitTable$asDF)
})
