#' @export
reliabilityClass <- R6::R6Class(
  "reliabilityClass",
  inherit = reliabilityBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- build the fixed 4-row table structure up front so the table
    # renders immediately instead of flickering from a blank placeholder to a
    # populated table when .run() finishes. The row count (always 4) and their
    # metric labels are fully determined by options (estim), not by the data,
    # so they belong here. .run() fills in the computed estimates via setRow().
    # ---------------------------------------------------------------------
    .init = function() {
      if (is.null(self$options$vars) || length(self$options$vars) < 3)
        return()

      estim <- self$options$estim
      table <- self$results$relTable

      labels <- c(
        "Cronbach's alpha",
        "PSI",
        paste0("Empirical (", estim, ")"),
        paste0("RMU (", estim, ")")
      )
      for (i in seq_along(labels)) {
        table$addRow(rowKey = i, values = list(
          metric   = labels[i],
          estimate = NA_real_,
          lower    = NA_real_,
          upper    = NA_real_,
          notes    = ""
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met. With 2
      # dichotomous items the MML model (mirt) is not estimable (too few
      # degrees of freedom), and reliability estimates from 2-item
      # scales are generally not informative, so require 3 items.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$relNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With 2 ",
          "dichotomous items the latent model used for the Empirical and ",
          "RMU estimates cannot be estimated (too few degrees of ",
          "freedom), and reliability estimates from 2-item scales are ",
          "generally not informative. Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Required suggested package
      if (!requireNamespace("ggdist", quietly = TRUE))
        stop("Package 'ggdist' is required. Install with: install.packages(\"ggdist\")")

      # 3. Extract data and convert to numeric
      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$relTable$setNote("sparse", sparse_msg)

      n_complete <- sum(complete.cases(df))
      n_total    <- nrow(df)
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )


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
          # Same robustness thresholds as the simulation analyses: a
          # CI from a handful of usable resamples would be misleading.
          if (actual_boot >= 20L && actual_boot >= boot_iter / 2) {
            alpha_int   <- ggdist::hdci(alpha_vec, .width = conf_int)
            alpha_lower <- alpha_int[1L, 1L]
            alpha_upper <- alpha_int[1L, 2L]
          }
        }

        # 11. Build result rows
        boot_note <- if (isTRUE(boot_alpha)) {
          if (!is.na(actual_boot) &&
              actual_boot >= 20L && actual_boot >= boot_iter / 2) {
            paste0(actual_boot, " bootstrap resamples")
          } else {
            paste0("bootstrap failed (only ", actual_boot, " of ",
                   boot_iter, " resamples usable)")
          }
        } else {
          "no bootstrap"
        }
        rmu_note <- paste0(draws, " PVs, ", rmu_iter, " RMU iterations")

        # Raw values (no pre-rounding) so the jamovi frontend applies
        # the user's "Number format" preferences.
        rows <- list(
          list(metric   = "Cronbach's alpha",
               estimate = alpha,
               lower    = alpha_lower,
               upper    = alpha_upper,
               notes    = boot_note),
          list(metric   = "PSI",
               estimate = psi,
               lower    = NA_real_,
               upper    = NA_real_,
               notes    = "eRm::SepRel"),
          list(metric   = paste0("Empirical (", estim, ")"),
               estimate = emp_rel,
               lower    = NA_real_,
               upper    = NA_real_,
               notes    = "mirt::empirical_rxx"),
          list(metric   = paste0("RMU (", estim, ")"),
               estimate = rmu_estimate,
               lower    = rmu_lower,
               upper    = rmu_upper,
               notes    = rmu_note)
        )

        # 12. Populate table (rows created in .init(); fill values here)
        table <- self$results$relTable
        for (i in seq_along(rows)) {
          table$setRow(rowKey = i, values = rows[[i]])
        }
        table$setNote(
          "context",
          paste0(
            "PSI uses eRm CML item parameters and excludes respondents with ",
            "min/max raw scores. Empirical reliability and RMU use MML item ",
            "parameters from mirt; theta estimator: ", estim, "."
          )
        )
        table$setNote(
          "hdci",
          paste0(
            "HDCI = highest-density continuous interval (width set by the ",
            "HDCI width option; here ", round(conf_int * 100, 1),
            "%). Available for Cronbach's alpha (when bootstrapped) and ",
            "RMU; PSI and Empirical reliability are reported as point ",
            "estimates only."
          )
        )

        # 13. Caption. The estimates use different samples when data are
        # missing: Cronbach's alpha is closed-form on complete cases,
        # while the model-based estimates retain partially missing rows
        # (eRm CML / mirt MML).
        n_used <- sum(rowSums(!is.na(df)) > 0)
        missing_msg <- if (n_used > n_complete) {
          paste0(
            " Cronbach's alpha is computed from the ", n_complete,
            " complete cases; the model-based estimates (PSI, Empirical, ",
            "RMU) use all ", n_used, " rows with at least one response ",
            "(eRm's CML and mirt's MML estimation accommodate partially ",
            "missing responses)."
          )
        } else {
          ""
        }
        self$results$relNote$setContent(
          paste0(
            "<p>Reliability based on N = ", n_used,
            " respondents across ", ncol(df), " items",
            if (n_used > n_complete)
              paste0(" (", n_complete, " with complete responses)")
            else "",
            ".", missing_msg, "</p>"
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
