# Internal helpers for Warm's weighted-likelihood (WLE) person estimation.
#
# Ported from easyRasch2's utils-theta.R so the score-to-logit table uses
# the same solver and information-based SEM as easyRasch2::RMscoreSE() /
# RMpersonParameters(). They operate on a normalised item-parameter
# representation: a list `thr_list` with one element per item, each a
# numeric vector of Andrich thresholds (category-boundary locations on the
# logit difficulty scale). A dichotomous item is an item with a single
# threshold equal to its difficulty, so one set of routines handles both
# the Rasch and Partial Credit cases.

#' Category probabilities for one item at a given theta (PCM / Rasch)
#'
#' @param theta Numeric person location.
#' @param thr Numeric vector of Andrich thresholds (length = n_cat - 1).
#' @return Numeric vector of length `length(thr) + 1` (categories 0..m),
#'   computed with a log-sum-exp shift for numerical stability.
#' @noRd
.pcm_cat_probs <- function(theta, thr) {
  cum_eta <- c(0, cumsum(theta - thr))
  shift   <- max(cum_eta)
  exp(cum_eta - (shift + log(sum(exp(cum_eta - shift)))))
}

#' Centre item thresholds to grand mean zero
#'
#' Standard Rasch identification: the mean of all thresholds is set to
#' zero so person locations are on a common scale, matching the
#' convention in easyRasch2.
#'
#' @param thr_list List of threshold vectors, one per item.
#' @return The list with every threshold shifted by the grand mean.
#' @noRd
.center_thresholds <- function(thr_list) {
  shift <- mean(unlist(thr_list))
  lapply(thr_list, function(x) x - shift)
}

#' Warm's weighted likelihood estimate for a single response pattern
#'
#' Solves the Warm-corrected score equation
#' `l'(theta) + J(theta) / (2 I(theta)) = 0` (J = summed third central
#' moment of the item scores), which is finite even at the extreme
#' scores. The SEM is the information-based `1 / sqrt(I(theta))` evaluated
#' at the estimate -- matching catR / TAM and easyRasch2 (>= dev), and
#' replacing the iarm "expected SEM" used previously.
#'
#' @param resp Numeric response vector for one person (may contain `NA`).
#' @param thr_list List of threshold vectors, one per item, aligned with
#'   `resp`.
#' @param theta_range Length-2 numeric search/boundary range.
#' @param tol uniroot tolerance.
#' @return Named numeric `c(theta, sem)`; `sem` is `NA` when the root lies
#'   outside `theta_range` (the boundary is then returned for `theta`).
#' @noRd
.theta_wle <- function(resp, thr_list, theta_range, tol = 1e-6) {
  ok <- !is.na(resp)
  if (!any(ok)) return(c(theta = NA_real_, sem = NA_real_))

  resp     <- resp[ok]
  thr_list <- thr_list[ok]
  r        <- sum(resp)

  # Test information at theta (also the inverse-square of the SEM).
  info_at <- function(theta) {
    s <- 0
    for (i in seq_along(thr_list)) {
      cats <- 0:length(thr_list[[i]])
      p    <- .pcm_cat_probs(theta, thr_list[[i]])
      E    <- sum(cats * p)
      s    <- s + sum((cats - E)^2 * p)
    }
    s
  }

  # Warm's weighted-likelihood score function.
  score <- function(theta) {
    dll <- r; info <- 0; third <- 0
    for (i in seq_along(thr_list)) {
      cats <- 0:length(thr_list[[i]])
      p    <- .pcm_cat_probs(theta, thr_list[[i]])
      E    <- sum(cats * p)
      dll   <- dll  - E
      info  <- info + sum((cats - E)^2 * p)
      third <- third + sum((cats - E)^3 * p)
    }
    dll + third / (2 * info)
  }

  lo <- theta_range[1L]; hi <- theta_range[2L]
  g_lo <- score(lo); g_hi <- score(hi)

  if (is.finite(g_lo) && is.finite(g_hi) && g_lo * g_hi < 0) {
    theta <- stats::uniroot(score, c(lo, hi), tol = tol)$root
    info  <- info_at(theta)
    sem   <- if (info < 1e-12) NA_real_ else 1 / sqrt(info)
  } else {
    # Root outside the search range: clamp to the nearer boundary.
    theta <- if (r <= sum(vapply(thr_list, length, integer(1))) / 2) lo else hi
    sem   <- NA_real_
  }
  c(theta = theta, sem = sem)
}

#' Build a representative response pattern summing to a given raw score
#'
#' For complete data the raw total is the sufficient statistic, so any
#' pattern with that sum yields the same WLE; this greedily fills items up
#' to their category maximum.
#'
#' @param r Target raw score.
#' @param steps Integer vector of per-item threshold counts (category max).
#' @return Integer response vector summing to `r`.
#' @noRd
.score_pattern <- function(r, steps) {
  resp <- integer(length(steps))
  rem  <- r
  for (i in seq_along(steps)) {
    take    <- min(steps[i], rem)
    resp[i] <- take
    rem     <- rem - take
  }
  resp
}

#' Extract centred Andrich thresholds from a complete response data.frame
#'
#' Fits RM (dichotomous) or PCM (polytomous) via eRm and returns the
#' per-item threshold list, grand-mean-zero centred. Mirrors the relevant
#' part of easyRasch2:::.rasch_fit_cml(se = FALSE); sparse-category
#' reporting is handled separately in the analysis via sparse_note().
#'
#' @param df Numeric, complete-case response data.frame.
#' @return List of centred threshold vectors, one per item.
#' @noRd
.wle_thresholds <- function(df) {
  is_poly <- max(as.matrix(df), na.rm = TRUE) > 1L
  if (!is_poly) {
    fit      <- eRm::RM(df)
    thr_list <- as.list(-as.numeric(fit$betapar))
  } else {
    fit      <- eRm::PCM(df)
    thr_tab  <- eRm::thresholds(fit)$threshtable[[1L]]
    thr_cols <- grep("^Threshold", colnames(thr_tab))
    thr_mat  <- thr_tab[, thr_cols, drop = FALSE]
    thr_list <- lapply(seq_len(nrow(thr_mat)), function(i) {
      v <- as.numeric(thr_mat[i, ])
      v[!is.na(v)]
    })
  }
  .center_thresholds(thr_list)
}
