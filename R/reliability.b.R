#' @export
reliabilityClass <- R6::R6Class(
  "reliabilityClass",
  inherit = reliabilityBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2)
        return()

      # 2. Required suggested package
      if (!requireNamespace("ggdist", quietly = TRUE))
        stop("Package 'ggdist' is required. Install with: install.packages(\"ggdist\")")

      # 3. Extract data and convert to numeric
      data <- self$data
      vars <- self$options$vars
      df   <- data[, vars, drop = FALSE]

      df <- to_numeric_responses_df(df)

      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste("The following variables contain no valid numeric data:",
                   paste(bad_vars, collapse = ", ")))
      }

      max_obs <- max(as.matrix(df), na.rm = TRUE)
      if (is.finite(max_obs) && max_obs > 20) {
        bad_cols <- names(df)[
          vapply(df, function(x) {
            mx <- suppressWarnings(max(x, na.rm = TRUE))
            is.finite(mx) && mx > 20
          }, logical(1L))
        ]
        stop(paste0(
          "Item(s) ", paste(bad_cols, collapse = ", "),
          " contain values > 20, which look like missing-value codes ",
          "(e.g., 999, 8888) rather than ordinal responses. ",
          "Mark these codes as missing in the data editor, or recode your data."
        ))
      }

      validate_response_data(df)

      n_complete <- sum(complete.cases(df))
      n_total    <- nrow(df)
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # 4. Read options
      estim       <- self$options$estim
      draws       <- self$options$draws
      rmu_iter    <- self$options$rmuIter
      conf_int    <- self$options$confInt / 100
      theta_range <- c(self$options$thetaMin, self$options$thetaMax)
      boot_alpha  <- isTRUE(self$options$bootAlpha)
      boot_iter   <- self$options$bootIter
      seed        <- self$options$seed

      if (theta_range[1L] >= theta_range[2L])
        stop("Theta lower bound must be less than upper bound.")

      tryCatch({
        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        data_mat      <- as.matrix(df)
        is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

        # 5. Full-sample fits
        mirt_fit <- suppressMessages(
          mirt::mirt(
            data       = df,
            model      = 1,
            itemtype   = "Rasch",
            verbose    = FALSE,
            accelerate = "squarem"
          )
        )
        erm_fit <- if (is_polytomous) eRm::PCM(df) else eRm::RM(df)

        # 6. Cronbach's alpha
        alpha <- .reliab_cronbach_alpha(df)

        # 7. PSI from eRm::SepRel (canonical)
        psi <- as.numeric(eRm::SepRel(eRm::person.parameter(erm_fit))$sep.rel)

        # 8. Empirical reliability (mirt)
        emp_rel <- as.numeric(
          mirt::empirical_rxx(
            mirt::fscores(mirt_fit,
                          method         = estim,
                          theta_lim      = theta_range,
                          full.scores.SE = TRUE,
                          verbose        = FALSE)
          )
        )

        # 9. RMU (mirt PVs + iterated split-half correlations)
        if (!is.null(seed)) set.seed(seed)
        pvs <- mirt::fscores(
          mirt_fit,
          method          = estim,
          theta_lim       = theta_range,
          plausible.draws = draws,
          plausible.type  = "MH",
          verbose         = FALSE
        )
        rmu_input <- do.call(cbind, lapply(pvs, as.numeric))

        rmu_iter_results <- do.call(
          rbind,
          lapply(seq_len(rmu_iter), function(i) {
            .reliab_rmu(rmu_input, level = conf_int)
          })
        )
        rmu_estimate <- mean(rmu_iter_results$rmu_estimate)
        rmu_lower    <- mean(rmu_iter_results$hdci_lowerbound)
        rmu_upper    <- mean(rmu_iter_results$hdci_upperbound)

        # 10. Bootstrap CI for Cronbach's alpha (sequential, closed-form)
        alpha_lower <- NA_real_
        alpha_upper <- NA_real_
        actual_boot <- NA_integer_
        if (isTRUE(boot_alpha)) {
          set.seed(seed + 1L)
          boot_seeds <- sample.int(.Machine$integer.max, boot_iter)
          alpha_vec <- vapply(seq_len(boot_iter), function(i) {
            set.seed(boot_seeds[i])
            idx <- sample.int(nrow(df), nrow(df), replace = TRUE)
            .reliab_cronbach_alpha(df[idx, , drop = FALSE])
          }, numeric(1L))
          alpha_vec <- alpha_vec[is.finite(alpha_vec)]
          actual_boot <- length(alpha_vec)
          if (actual_boot >= 2L) {
            alpha_int   <- ggdist::hdci(alpha_vec, .width = conf_int)
            alpha_lower <- alpha_int[1L, 1L]
            alpha_upper <- alpha_int[1L, 2L]
          }
        }

        # 11. Build result rows
        boot_note <- if (isTRUE(boot_alpha)) {
          if (!is.na(actual_boot) && actual_boot >= 2L) {
            paste0(actual_boot, " bootstrap resamples")
          } else {
            "bootstrap failed"
          }
        } else {
          "no bootstrap"
        }
        rmu_note <- paste0(draws, " PVs, ", rmu_iter, " RMU iterations")

        rows <- list(
          list(metric   = "Cronbach's alpha",
               estimate = round(alpha, 3),
               lower    = round(alpha_lower, 3),
               upper    = round(alpha_upper, 3),
               notes    = boot_note),
          list(metric   = "PSI",
               estimate = round(psi, 3),
               lower    = NA_real_,
               upper    = NA_real_,
               notes    = "eRm::SepRel"),
          list(metric   = paste0("Empirical (", estim, ")"),
               estimate = round(emp_rel, 3),
               lower    = NA_real_,
               upper    = NA_real_,
               notes    = "mirt::empirical_rxx"),
          list(metric   = paste0("RMU (", estim, ")"),
               estimate = round(rmu_estimate, 3),
               lower    = round(rmu_lower, 3),
               upper    = round(rmu_upper, 3),
               notes    = rmu_note)
        )

        # 12. Populate table
        table <- self$results$relTable
        for (i in seq_along(rows)) {
          table$addRow(rowKey = i, values = rows[[i]])
        }
        table$setNote(
          "context",
          paste0(
            "PSI uses eRm CML item parameters and excludes respondents with ",
            "min/max raw scores. Empirical reliability and RMU use MML item ",
            "parameters from mirt; theta estimator: ", estim, "."
          )
        )

        # 13. Caption
        n_excluded <- n_total - n_complete
        excluded_msg <- if (n_excluded > 0L) {
          paste0(" (", n_excluded, " of ", n_total,
                 " row(s) had a missing response and were excluded from ",
                 "model fitting)")
        } else {
          ""
        }
        self$results$relNote$setContent(
          paste0(
            "<p>Reliability based on n = ", n_complete,
            " complete responses across ", ncol(df), " items",
            excluded_msg, ". HDCI width = ",
            round(conf_int * 100, 1), "%.</p>"
          )
        )

      }, error = function(e) {
        stop(paste("Error in reliability analysis:", e$message))
      })
    }
  )
)

# ---------------------------------------------------------------------------
# Internal helpers (free functions, not on the R6 class)
# ---------------------------------------------------------------------------

#' Cronbach's alpha (closed-form, complete cases)
#'
#' @keywords internal
#' @noRd
.reliab_cronbach_alpha <- function(data) {
  d <- stats::na.omit(as.data.frame(data))
  k <- ncol(d)
  if (k < 2L || nrow(d) < 2L) return(NA_real_)
  total_var <- stats::var(rowSums(d))
  if (!is.finite(total_var) || total_var == 0) return(NA_real_)
  item_vars <- vapply(d, stats::var, numeric(1L))
  (k / (k - 1L)) * (1 - sum(item_vars) / total_var)
}

#' Relative Measurement Uncertainty from posterior / plausible-value draws
#'
#' Adapted from gbtoolbox::reliability() (GPL-2/3) and the easyRaschBayes
#' implementation. Random column split, paired Pearson correlations, summary
#' via ggdist::mean_hdci.
#'
#' @keywords internal
#' @noRd
.reliab_rmu <- function(input_draws, level = 0.95) {
  input_draws <- as.matrix(input_draws)
  if (ncol(input_draws) < 2L) {
    stop("`input_draws` must have at least 2 columns.", call. = FALSE)
  }

  col_select <- sample.int(ncol(input_draws), replace = FALSE)
  half       <- floor(length(col_select) / 2)
  cols_a     <- col_select[seq_len(half)]
  cols_b     <- col_select[(half + 1L):(2L * half)]

  draws_a <- input_draws[, cols_a, drop = FALSE]
  draws_b <- input_draws[, cols_b, drop = FALSE]

  rel_post <- vapply(seq_len(ncol(draws_a)), function(i) {
    x <- draws_a[, i]
    y <- draws_b[, i]
    if (stats::var(x, na.rm = TRUE) == 0 ||
        stats::var(y, na.rm = TRUE) == 0) {
      return(0)
    }
    stats::cor(x, y, method = "pearson", use = "complete.obs")
  }, numeric(1L))

  hdci <- ggdist::mean_hdci(rel_post, .width = level)
  colnames(hdci)[1:3] <- c("rmu_estimate", "hdci_lowerbound", "hdci_upperbound")
  hdci
}
