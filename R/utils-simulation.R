# Observed CFA fit helper, retained solely for the cfacutoff analysis'
# graceful-degradation path: when the cutoff simulation fails,
# easyRasch2::RMdimCFA() refuses to run (observed CFA indices are not
# interpretable without the simulated reference distribution), but the
# module still shows the observed indices with an explanation.
# Adapted from easyRasch2 (https://github.com/pgmj/easyRasch2), GPL >= 3.

#' Compute observed CFA fit on a complete-cases item data.frame
#'
#' Returns the same `(cfi, rmsea, srmr)` triple as `run_single_cfa_sim()`
#' so the populating code can format both uniformly. Returns a character
#' message on failure (typically a non-converging WLSMV fit).
#'
#' The items are refitted under placeholder names `V1...Vk`. jamovi variable
#' names may contain spaces or start with a digit (`Item 1`, `3 months`),
#' which lavaan's model syntax cannot express: it rejects such names as
#' undefined variables, and neither back-quoting nor double-quoting them
#' parses (verified with lavaan 0.6-21). Only the fit indices are returned,
#' so the names never surface; the fit is unchanged by the renaming.
#'
#' @noRd
run_observed_cfa_fit <- function(df, estimator) {
  safe <- paste0("V", seq_along(df))
  names(df) <- safe
  fmla <- paste0("F1 =~ ", paste(safe, collapse = " + "))
  tryCatch({
    fit <- suppressWarnings(suppressMessages(
      lavaan::cfa(
        model     = fmla,
        data      = df,
        ordered   = safe,
        estimator = estimator,
        warn      = FALSE,
        verbose   = FALSE
      )
    ))
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      return("convergence_failed: lavaan did not converge")
    }
    suppressWarnings(extract_cfa_fit(fit, estimator))
  }, error = function(e) as.character(conditionMessage(e)))
}

#' Pull (CFI, RMSEA, SRMR) from a fitted lavaan object
#'
#' Uses the Satorra-Bentler-scaled CFI / RMSEA (`cfi.scaled`,
#' `rmsea.scaled`) when available, falling back to the uncorrected
#' `cfi` / `rmsea` when not. Scaled variants are preferred over the
#' Yuan-Bentler `.robust` variants because the latter return `NA` for a
#' non-trivial fraction of small-n fits, which would produce holes in
#' the simulated null distribution. For percentile-based comparison the
#' binding requirement is internal consistency across iterations (same
#' metric on both observed and simulated data), which `.scaled` provides
#' more reliably. SRMR is unaffected by either correction.
#'
#' @noRd
extract_cfa_fit <- function(fit, estimator = NULL) {
  fm <- lavaan::fitMeasures(fit)
  pick_scaled <- function(name) {
    key_s <- paste0(name, ".scaled")
    if (key_s %in% names(fm) && is.finite(fm[[key_s]])) {
      return(as.numeric(fm[[key_s]]))
    }
    if (name %in% names(fm) && is.finite(fm[[name]])) {
      return(as.numeric(fm[[name]]))
    }
    NA_real_
  }
  c(pick_scaled("cfi"),
    pick_scaled("rmsea"),
    as.numeric(fm[["srmr"]]))
}

#' Run a cutoff simulation with caching, error capture and an iteration floor
#'
#' Shared by every analysis with *Simulation-based cutoffs*. jamovi reruns
#' `.run()` on every option change, including display options that do not
#' touch the simulation, so each analysis keeps the cutoff object in a
#' results element's state (a plot, or a hidden `simCache`) whose `clearWith`
#' lists the simulation-relevant options. A surviving state is reused only
#' when its signature and data keys match, so reuse does not rely on
#' `clearWith` alone. Saving the state stays with each analysis, since plot
#' states also carry plot-specific fields.
#'
#' On an error, or fewer than 20 successful iterations, the cutoff object is
#' dropped and `fail_msg` explains why. The caller falls back to its
#' asymptotic output and passes `fail_msg` to sim_failure_notice().
#'
#' @param compute Logical: was the simulation requested?
#' @param cached The state of the caching element, or `NULL` for none.
#' @param sig List of simulation settings (iterations, seed, HDCI width).
#' @param keys Named list of data the result depends on, compared with the
#'   same-named fields of `cached` (e.g. `list(df = df, dif = dif_vec)`).
#' @param simulate Zero-argument closure calling the `easyRasch2::RM*Cutoff()`
#'   function. Its warnings and messages are suppressed; the module states
#'   its own caveats.
#' @param requested Number of iterations requested, for the failure message.
#' @param hint Optional sentence naming the analysis' discard rule, from
#'   sim_rule_hint().
#' @return `list(cutoff_res, fail_msg)`, at most one of them non-`NULL`.
#' @noRd
run_cutoff_simulation <- function(compute, cached = NULL, sig = NULL,
                                  keys = list(), simulate, requested,
                                  hint = NULL) {
  out <- list(cutoff_res = NULL, fail_msg = NULL)
  if (!isTRUE(compute)) return(out)

  keys_match <- function() {
    all(vapply(names(keys), function(k) identical(cached[[k]], keys[[k]]),
               logical(1)))
  }
  if (!is.null(cached) && !is.null(cached$cutoff_res) &&
      identical(cached$sig, sig) && keys_match()) {
    out$cutoff_res <- cached$cutoff_res
    return(out)
  }

  res <- tryCatch(
    suppressWarnings(suppressMessages(simulate())),
    error = function(e) {
      out$fail_msg <<- conditionMessage(e)
      NULL
    }
  )
  # With very few successful iterations the intervals and percentiles
  # collapse onto a handful of values and every flag is spurious.
  if (!is.null(res) && res$actual_iterations < 20L) {
    out$fail_msg <- paste0(
      "Only ", res$actual_iterations, " of ", requested, " simulation ",
      "iterations succeeded, and at least 20 are needed for reliable ",
      "results. A simulated dataset is discarded when the model cannot be ",
      "refitted to it, which is most likely when items have very low or ",
      "very high endorsement rates relative to the sample size.",
      if (!is.null(hint)) paste0(" ", hint)
    )
    res <- NULL
  }
  out$cutoff_res <- res
  out
}

#' The discard rule of an easyRasch2 cutoff simulation, as a sentence
#'
#' The package validates each simulated dataset before refitting the model
#' and discards those that fail. The rules differ by analysis; these are
#' the ones in easyRasch2 1.4.0. The DIF simulation has no such rule (its
#' iterations fail only on errors) and takes no hint.
#'
#' @param rule `"positive"` (infit, item-restscore, Q3, partial gamma LD),
#'   `"both"` (residual PCA) or `"cfa"` (CFA fit indices).
#' @return A single sentence.
#' @noRd
sim_rule_hint <- function(rule = c("positive", "both", "cfa")) {
  rule <- match.arg(rule)
  dich <- switch(rule,
    positive = "at least 8 positive responses",
    both     = "at least 8 positive and 8 negative responses",
    cfa      = "at least 2 responses in each category"
  )
  paste0("Each simulated item needs ", dich, " for dichotomous data, and ",
         "every response category must occur for polytomous data.")
}
