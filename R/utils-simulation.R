# Internal simulation helpers for the Q3 cutoff parametric bootstrap.
# Adapted from easyRasch2 (https://github.com/pgmj/easyRasch2), GPL >= 3.

#' Simulate responses for a single polytomous item
#'
#' @param deltas Numeric vector of threshold parameters.
#' @param thetas Numeric vector of person parameters.
#' @return Integer vector of simulated responses (0-indexed).
#' @noRd
sim_poly_item <- function(deltas, thetas) {
  k <- length(deltas) + 1L
  n <- length(thetas)
  Y <- outer(thetas, deltas, "-")
  cumsums <- t(rbind(rep(0, times = n), apply(X = Y, MARGIN = 1, FUN = cumsum)))
  expcumsums <- exp(cumsums)
  norms <- apply(X = expcumsums, MARGIN = 1, FUN = sum)
  z <- expcumsums / norms
  vapply(X = seq_len(n), FUN = function(x) {
    sample(x = 0L:(k - 1L), size = 1L, replace = TRUE, prob = z[x, ])
  }, FUN.VALUE = 1L)
}

#' Simulate partial score (polytomous) response matrix
#'
#' @param deltaslist List of threshold parameter vectors, one per item.
#' @param thetavec Numeric vector of person parameters.
#' @return Integer matrix of simulated responses (rows = persons, cols = items).
#' @noRd
sim_partial_score <- function(deltaslist, thetavec) {
  deltaslist <- lapply(X = deltaslist, FUN = unlist)
  sapply(X = deltaslist, FUN = sim_poly_item, thetas = thetavec)
}

#' Extract item threshold parameters from polytomous data via CML
#'
#' @param data A data.frame or matrix of polytomous item responses (min = 0).
#' @return Numeric matrix of threshold parameters.
#' @noRd
extract_item_thresholds <- function(data) {
  pcm_fit <- eRm::PCM(data)
  thresh <- eRm::thresholds(pcm_fit)
  thresh_mat <- thresh$threshtable$`1`
  if ("Location" %in% colnames(thresh_mat)) {
    thresh_mat <- thresh_mat[, colnames(thresh_mat) != "Location", drop = FALSE]
  }
  thresh_mat
}

#' Run a single Q3 simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside `.runCutoffSim()`. Should include
#'   `item_names` (character vector) so per-pair `pair_q3` rows are labelled
#'   with the user's item names rather than mirt's auto V1/V2/... .
#' @return A list with `mean`, `max`, and a long-format `pair_q3`
#'   data.frame (one row per upper-triangle pair: `Item1`, `Item2`, `Q3`),
#'   or a character string on failure.
#' @noRd
run_single_q3_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n, replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta  = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      # Fewer than 8 positive responses in an item causes numerical instability
      # during mirt model fitting (near-zero or zero variance in the item).
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df  <- as.data.frame(sim_mat)

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return(paste0("validation_failed: not all categories represented in item ", j))
        }
      }
    }

    # Preserve the user's item labels so per-pair Q3 rows use them
    # (rather than psychotools / mirt's auto-generated V1, V2, ...).
    if (!is.null(data_list$item_names) &&
        length(data_list$item_names) == ncol(sim_df)) {
      colnames(sim_df) <- data_list$item_names
    }

    mirt_fit <- mirt::mirt(
      sim_df,
      model    = 1,
      itemtype = "Rasch",
      verbose  = FALSE,
      accelerate = "squarem",
      quadpts = 29, # for about 20-25% speed improvement without noticeable
      TOL = 0.005   # precision loss for Q3 residuals
    )

    q3_mat <- mirt::residuals(mirt_fit, type = "Q3", digits = 4, verbose = FALSE)
    diag(q3_mat) <- NA

    mean_q3 <- mean(q3_mat, na.rm = TRUE)
    max_q3  <- max(q3_mat,  na.rm = TRUE)

    # Per-pair Q3 (upper triangle, long format) for the per-pair plot.
    item_names_q3 <- colnames(q3_mat)
    if (is.null(item_names_q3)) {
      item_names_q3 <- as.character(seq_len(ncol(q3_mat)))
    }
    upper_idx <- which(upper.tri(q3_mat), arr.ind = TRUE)
    pair_q3 <- data.frame(
      Item1 = item_names_q3[upper_idx[, "row"]],
      Item2 = item_names_q3[upper_idx[, "col"]],
      Q3    = q3_mat[upper_idx],
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    list(mean = mean_q3, max = max_q3, pair_q3 = pair_q3)
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

#' Run Q3 simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each iteration.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_q3_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                  verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_q3_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' Run a single infit simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside the cutoff simulation.
#' @return A data.frame with columns `Item`, `InfitMSQ`, `OutfitMSQ`, or a
#'   character string on failure.
#' @noRd
run_single_infit_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n, replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }

      model_fit <- eRm::RM(sim_df, se = FALSE)
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }

      model_fit <- psychotools::pcmodel(sim_df, hessian = FALSE)
    }

    cfit <- iarm::out_infit(model_fit)

    data.frame(
      Item      = data_list$item_names,
      InfitMSQ  = round(cfit$Infit,  3),
      OutfitMSQ = round(cfit$Outfit, 3),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

#' Run infit simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each iteration.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_infit_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                     verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_infit_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' Run a single partial gamma DIF simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside the cutoff simulation.
#' @return A data.frame with columns `Item` and `gamma`, or a character string
#'   on failure.
#' @noRd
run_single_partgam_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L,
                       integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    # Create a random DIF variable with the same group proportions
    # but no actual relationship to item responses (no true DIF)
    random_dif <- sample(
      data_list$dif_levels,
      size = data_list$sample_n,
      replace = TRUE,
      prob = data_list$dif_proportions
    )

    # Compute partial gamma DIF via iarm
    pgam <- iarm::partgam_DIF(sim_df, random_dif)
    pgam_df <- as.data.frame(pgam)

    data.frame(
      Item  = data_list$item_names,
      gamma = as.numeric(pgam_df[seq_along(data_list$item_names), "gamma"]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

#' Run partial gamma DIF simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each iteration.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_partgam_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                       verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_partgam_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' Run a single residual-PCA simulation iteration
#'
#' Simulates a dataset of size `data_list$sample_n` from the fitted Rasch
#' or partial-credit model, fits the same model to the simulated data,
#' and returns the largest eigenvalue of the unrotated PCA on the
#' standardised residuals (Chou & Wang, 2010).
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside `.runCutoffSim()`.
#' @return A numeric scalar (the first eigenvalue), or a character string
#'   on failure.
#' @noRd
run_single_pca_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta  = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
      neg_counts <- nrow(sim_df) - pos_counts
      if (any(neg_counts < 8L)) {
        return("validation_failed: fewer than 8 negative responses in at least one item")
      }

      model_fit <- eRm::RM(sim_df, se = FALSE)
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df  <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L,
                       integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }

      model_fit <- eRm::PCM(sim_df, se = FALSE)
    }

    pp        <- eRm::person.parameter(model_fit)
    ifit      <- eRm::itemfit(pp)
    st_resids <- ifit$st.res

    if (anyNA(st_resids)) {
      keep <- stats::complete.cases(st_resids)
      st_resids <- st_resids[keep, , drop = FALSE]
    }
    if (nrow(st_resids) < ncol(st_resids)) {
      return("validation_failed: too few rows in residual matrix for PCA")
    }

    pca_fit <- stats::prcomp(st_resids)
    as.numeric(pca_fit$sdev[1L]^2)
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

#' Run residual-PCA simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each iteration.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration); each element
#'   is either a numeric scalar (the first eigenvalue) or a character
#'   string describing the failure.
#' @noRd
run_pca_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                   verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_pca_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' Compute PCM expected-score matrix
#'
#' Used by the variance-partition computation in the residualpca module.
#' Returns an n_persons x n_items matrix of model-expected scores under
#' the partial-credit model.
#'
#' @param thetas Numeric vector of person locations.
#' @param thresh_mat Numeric matrix of thresholds (rows = items, cols =
#'   thresholds; pad shorter rows with NA).
#' @return Numeric matrix of expected scores.
#' @noRd
pcm_expected_scores <- function(thetas, thresh_mat) {
  n        <- length(thetas)
  n_items  <- nrow(thresh_mat)
  expected <- matrix(NA_real_, nrow = n, ncol = n_items)

  for (j in seq_len(n_items)) {
    taus <- thresh_mat[j, !is.na(thresh_mat[j, ])]
    K_j  <- length(taus)
    cumsum_taus <- c(0, cumsum(taus))
    log_num <- outer(thetas, 0:K_j) -
      matrix(cumsum_taus, nrow = n, ncol = K_j + 1L, byrow = TRUE)
    log_num <- log_num - apply(log_num, 1L, max)
    probs   <- exp(log_num)
    probs   <- probs / rowSums(probs)
    expected[, j] <- as.numeric(probs %*% (0:K_j))
  }
  expected
}

#' Run a single posterior-predictive CFA-fit-index simulation iteration
#'
#' Simulates a dataset under the supplied PCM/RM parameters, refits a
#' one-factor lavaan CFA on it, and returns the fit indices.
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List built inside `.runCutoffSim()`. Includes
#'   `estimator` so the appropriate robust/scaled fit-index variants can
#'   be returned.
#' @return A length-3 numeric vector `(cfi, rmsea, srmr)` on success, or
#'   a character string describing the failure.
#' @noRd
run_single_cfa_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(theta = thetas_res,
                                  beta  = data_list$item_params)
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      neg_counts <- nrow(sim_df) - pos_counts
      if (any(pos_counts < 2L) || any(neg_counts < 2L)) {
        return("validation_failed: an item has < 2 responses in one of the two categories")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df  <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist,
                       function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    fmla <- paste0("F1 =~ ",
                   paste(data_list$item_names, collapse = " + "))

    fit <- suppressWarnings(suppressMessages(
      lavaan::cfa(
        model     = fmla,
        data      = sim_df,
        ordered   = data_list$item_names,
        estimator = data_list$estimator,
        warn      = FALSE,
        verbose   = FALSE
      )
    ))

    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      return("convergence_failed: lavaan did not converge")
    }

    suppressWarnings(extract_cfa_fit(fit, data_list$estimator))
  }, error = function(e) as.character(conditionMessage(e)))
}

#' Run CFA simulations sequentially
#'
#' @noRd
run_cfa_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                   verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_cfa_sim(sim_seeds[sim], sim_data_list)
    if (verbose) utils::setTxtProgressBar(pb, sim)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' Compute observed CFA fit on a complete-cases item data.frame
#'
#' Returns the same `(cfi, rmsea, srmr)` triple as `run_single_cfa_sim()`
#' so the populating code can format both uniformly. Returns a character
#' message on failure (typically a non-converging WLSMV fit).
#'
#' @noRd
run_observed_cfa_fit <- function(df, estimator) {
  fmla <- paste0("F1 =~ ", paste(names(df), collapse = " + "))
  tryCatch({
    fit <- suppressWarnings(suppressMessages(
      lavaan::cfa(
        model     = fmla,
        data      = df,
        ordered   = names(df),
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
