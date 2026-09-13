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
