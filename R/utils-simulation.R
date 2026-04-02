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
#' @param data_list List produced inside `.runCutoffSim()`.
#' @return A list with `mean` and `max` Q3, or a character string on failure.
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

    mirt_fit <- mirt::mirt(
      sim_df,
      model    = 1,
      itemtype = "Rasch",
      verbose  = FALSE,
      accelerate = "squarem"
    )

    q3_mat <- mirt::residuals(mirt_fit, type = "Q3", digits = 4, verbose = FALSE)
    diag(q3_mat) <- NA

    mean_q3 <- mean(q3_mat, na.rm = TRUE)
    max_q3  <- max(q3_mat,  na.rm = TRUE)

    list(mean = mean_q3, max = max_q3)
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
