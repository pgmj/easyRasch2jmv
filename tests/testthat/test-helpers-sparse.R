# Unit tests for sparse-category and duplicate-item footnote helpers.

test_that("sparse_note is NULL when every category is well populated", {
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:2, 60 * 4, replace = TRUE),
                             ncol = 4, dimnames = list(NULL, paste0("I", 1:4))))
  expect_null(er2$sparse_note(df))
})

test_that("sparse_note names the item with a thin category", {
  df <- data.frame(
    I1 = c(rep(0, 20), rep(1, 20), 2, 2),   # category 2 has only 2 obs
    I2 = rep(c(0, 1, 2), length.out = 42)
  )
  msg <- er2$sparse_note(df)
  expect_type(msg, "character")
  expect_match(msg, "I1")
  expect_match(msg, "fewer than 3")
})

test_that("sparse_note_grouped flags a category sparse within one group", {
  # Whole sample uses categories 0/1/2; group B never reaches 2.
  df <- data.frame(I1 = c(rep(0:2, 10), 0, 1))
  grp <- factor(c(rep("A", 30), rep("B", 2)))
  msg <- er2$sparse_note_grouped(df, grp)
  expect_match(msg, "I1")
  expect_match(msg, "\\(B\\)")
})

test_that("duplicate_items_note fires only on perfect correlation", {
  dup  <- data.frame(A = c(0, 1, 2, 1, 0), B = c(0, 1, 2, 1, 0),
                     C = c(1, 0, 2, 2, 1))
  none <- data.frame(A = c(0, 1, 2, 1, 0), B = c(1, 0, 2, 0, 1),
                     C = c(1, 0, 2, 2, 1))
  expect_match(er2$duplicate_items_note(dup), "perfectly correlated")
  expect_null(er2$duplicate_items_note(none))
})
