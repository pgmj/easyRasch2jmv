# Unit tests for run_cutoff_simulation(), the shared cache / error /
# iteration-floor step of every simulation-based analysis.

fake_cutoff <- function(n) list(actual_iterations = n, tag = "fresh")

test_that("run_cutoff_simulation does nothing when not requested", {
  called <- FALSE
  out <- er2$run_cutoff_simulation(
    compute = FALSE, simulate = function() { called <<- TRUE },
    requested = 100)
  expect_null(out$cutoff_res)
  expect_null(out$fail_msg)
  expect_false(called)
})

test_that("run_cutoff_simulation reuses a cache only when signature and keys match", {
  df  <- data.frame(a = 1:3)
  sig <- list(iterations = 100, seed = 42L)
  cached <- list(cutoff_res = list(tag = "cached"), sig = sig, df = df,
                 dif = c("x", "y", "x"))
  run <- function(...) er2$run_cutoff_simulation(
    compute = TRUE, cached = cached, simulate = function() fake_cutoff(100),
    requested = 100, ...)

  expect_identical(run(sig = sig, keys = list(df = df))$cutoff_res$tag,
                   "cached")
  expect_identical(
    run(sig = sig, keys = list(df = df, dif = c("x", "y", "x")))$cutoff_res$tag,
    "cached")
  # any mismatch reruns the simulation
  expect_identical(
    run(sig = list(iterations = 100, seed = 1L),
        keys = list(df = df))$cutoff_res$tag, "fresh")
  expect_identical(
    run(sig = sig, keys = list(df = data.frame(a = 3:1)))$cutoff_res$tag,
    "fresh")
  expect_identical(
    run(sig = sig, keys = list(df = df, dif = c("y", "y", "x")))$cutoff_res$tag,
    "fresh")
})

test_that("run_cutoff_simulation captures errors and applies the 20-iteration floor", {
  err <- er2$run_cutoff_simulation(
    compute = TRUE, simulate = function() stop("boom"), requested = 100)
  expect_null(err$cutoff_res)
  expect_identical(err$fail_msg, "boom")

  low <- er2$run_cutoff_simulation(
    compute = TRUE, simulate = function() fake_cutoff(19L), requested = 100,
    hint = er2$sim_rule_hint("positive"))
  expect_null(low$cutoff_res)
  expect_match(low$fail_msg, "^Only 19 of 100 simulation iterations")
  expect_match(low$fail_msg, "at least 8 positive responses")

  ok <- er2$run_cutoff_simulation(
    compute = TRUE, simulate = function() fake_cutoff(20L), requested = 100)
  expect_identical(ok$cutoff_res$actual_iterations, 20L)
  expect_null(ok$fail_msg)

  # warnings and messages from the package are silenced
  expect_silent(er2$run_cutoff_simulation(
    compute = TRUE, requested = 100,
    simulate = function() { warning("w"); message("m"); fake_cutoff(50L) }))
})

test_that("sim_rule_hint names each package discard rule", {
  expect_match(er2$sim_rule_hint("positive"), "at least 8 positive responses")
  expect_match(er2$sim_rule_hint("both"), "8 positive and 8 negative")
  expect_match(er2$sim_rule_hint("cfa"), "at least 2 responses in each")
  expect_error(er2$sim_rule_hint("dif"))
})
