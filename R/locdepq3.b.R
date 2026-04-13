#' @export
locdepq3Class <- R6::R6Class(
  "locdepq3Class",
  inherit = locdepq3Base,
  private = list(

    # .init() runs immediately when options change — sets up table structure
    # so column headers appear instantly and don't flicker on re-run.
    .init = function() {
      vars <- self$options$vars
      if (is.null(vars) || length(vars) == 0)
        return()

      if (length(self$options$vars) < 2)
        return()

      q3_table <- self$results$q3Table

      q3_table$addColumn(name = "item", title = "", type = "text")
      for (v in vars) {
        q3_table$addColumn(name = v, title = v, type = "number")
      }
      if (isTRUE(self$options$computeCutoff)) {
        q3_table$addColumn(name = "above_cutoff", title = "Above cutoff", type = "text")
      }

      # Pre-populate rows with em dashes on the diagonal and blanks above
      for (i in seq_along(vars)) {
        row_vals <- list(item = vars[i])
        for (j in seq_along(vars)) {
          if (i == j) {
            row_vals[[ vars[j] ]] <- "\u2014"
          } else if (j > i) {
            row_vals[[ vars[j] ]] <- ""
          }
        }
        q3_table$addRow(rowKey = i, values = row_vals)
      }
    },

    .run = function() {
      # Return early if not enough variables
      if (is.null(self$options$vars) || length(self$options$vars) < 2) {
        return()
      }

      data            <- self$data
      vars            <- self$options$vars
      compute_cutoff  <- self$options$computeCutoff
      iterations      <- self$options$iterations
      seed_val        <- self$options$seed

      # Select only the specified variables
      df <- data[, vars, drop = FALSE]

      # Convert to numeric — handle factors properly
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          df[[col]] <- as.numeric(as.character(df[[col]]))
        } else {
          df[[col]] <- as.numeric(df[[col]])
        }
      }

      # Check for all-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste(
          "The following variables contain no valid numeric data:",
          paste(bad_vars, collapse = ", ")
        ))
      }

      # Validate data (non-negative integers, min = 0)
      validate_response_data(df)

      # Check for sufficient complete cases
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. Each row must have responses for all selected items.")
      }
      if (n_complete < 30) {
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable with small samples.",
          n = n_complete
        )
      }

      tryCatch(
        {
          # --- Step 1: Fit mirt Rasch model and compute Q3 matrix ----------------
          mirt_model <- mirt::mirt(
            df,
            model      = 1,
            itemtype   = "Rasch",
            verbose    = FALSE,
            accelerate = "squarem"
          )

          resid_mat <- mirt::residuals(mirt_model, type = "Q3", digits = 2, verbose = FALSE)
          diag(resid_mat) <- NA
          mean_resid <- mean(resid_mat, na.rm = TRUE)

          # --- Step 2: Simulation-based cutoff (optional) -------------------------
          cutoff_val <- NULL
          dyn_cutoff <- NULL

          if (compute_cutoff) {
            actual_seed <- if (seed_val == 0) NULL else as.integer(seed_val)
            cutoff_res  <- private$.runCutoffSim(df, iterations, actual_seed)
            cutoff_val  <- cutoff_res$suggested_cutoff
            dyn_cutoff  <- mean_resid + cutoff_val

            # Populate cutoff summary table
            ct <- self$results$cutoffTable
            ct$addRow(rowKey = "sug",  values = list(parameter = "Suggested cutoff (p99)", value = round(cutoff_res$suggested_cutoff, 4)))
            ct$addRow(rowKey = "p95",  values = list(parameter = "p95",  value = round(cutoff_res$p95,  4)))
            ct$addRow(rowKey = "p99",  values = list(parameter = "p99",  value = round(cutoff_res$p99,  4)))
            ct$addRow(rowKey = "p995", values = list(parameter = "p99.5", value = round(cutoff_res$p995, 4)))
            ct$addRow(rowKey = "p999", values = list(parameter = "p99.9", value = round(cutoff_res$p999, 4)))
            ct$addRow(rowKey = "iter", values = list(parameter = "Actual iterations", value = cutoff_res$actual_iterations))
            ct$addRow(rowKey = "n",    values = list(parameter = "Sample N",          value = cutoff_res$sample_n))
          }

          # --- Step 3: Determine above_cutoff flags using the lower triangle only
          # (matching what is displayed in the table)
          above_flags <- rep(FALSE, length(vars))
          if (!is.null(dyn_cutoff)) {
            for (i in seq_along(vars)) {
              if (i > 1) {
                above_flags[i] <- any(resid_mat[i, 1:(i - 1)] > dyn_cutoff, na.rm = TRUE)
              }
              # Row 1 has no lower-triangle cells, so it stays FALSE
            }
          }

          # --- Step 4: Populate the Q3 table (structure set up in .init()) --------
          q3_table <- self$results$q3Table

          # Round for display
          resid_rounded <- round(resid_mat, 2)

          # Add note when cutoff is applied
          if (!is.null(dyn_cutoff)) {
            q3_table$setNote(
              "cutoff",
              paste0(
                "Dynamic cut-off: ", round(dyn_cutoff, 3),
                " (mean Q3 = ", round(mean_resid, 3),
                " + ", round(cutoff_val, 3),
                " (p99 from simulation)). ",
                "Rows marked * contain at least one value above the cut-off."
              )
            )
          }

          # Populate rows — diagonal gets em dash, upper triangle gets blank,
          # lower triangle gets the Q3 value
          for (i in seq_along(vars)) {
            row_vals <- list(item = vars[i])
            for (j in seq_along(vars)) {
              if (i == j) {
                row_vals[[ vars[j] ]] <- "\u2014"
              } else if (j > i) {
                row_vals[[ vars[j] ]] <- ""
              } else {
                row_vals[[ vars[j] ]] <- resid_rounded[i, j]
              }
            }
            if (compute_cutoff) {
              row_vals[["above_cutoff"]] <- if (isTRUE(above_flags[i])) "*" else ""
            }
            q3_table$setRow(rowNo = i, values = row_vals)
          }
        },
        error = function(e) {
          stop(paste("Error in Q3 analysis:", e$message))
        }
      )
    },

    .runCutoffSim = function(df, iterations, seed) {
      if (!is.null(seed)) set.seed(seed)

      sim_seeds <- sample.int(.Machine$integer.max, iterations)
      data_mat  <- as.matrix(df)
      sample_n  <- nrow(data_mat)
      is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

      if (is_polytomous) {
        pcm_fit     <- eRm::PCM(data_mat)
        pp          <- eRm::person.parameter(pcm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores  <- rowSums(data_mat, na.rm = TRUE)
        thetas      <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        thresh_mat  <- extract_item_thresholds(data_mat)
        deltaslist  <- lapply(seq_len(nrow(thresh_mat)), function(i) {
          as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
        })
        sim_data_list <- list(
          type      = "polytomous",
          thetas    = thetas,
          deltaslist = deltaslist,
          n_items   = ncol(data_mat),
          sample_n  = sample_n
        )
      } else {
        rm_fit      <- eRm::RM(data_mat)
        pp          <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores  <- rowSums(data_mat, na.rm = TRUE)
        thetas      <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        sim_data_list <- list(
          type        = "dichotomous",
          thetas      = thetas,
          item_params = item_params,
          n_items     = ncol(data_mat),
          sample_n    = sample_n
        )
      }

      results_raw <- run_q3_sim_sequential(iterations, sim_seeds, sim_data_list,
                                           verbose = FALSE)

      ok         <- vapply(results_raw, function(x) is.list(x), logical(1L))
      successful <- results_raw[ok]

      if (length(successful) == 0L) {
        stop(
          "All simulation iterations failed. Check your data: items must have ",
          "sufficient response variation (at least 8 positive responses per item ",
          "for dichotomous data; all response categories represented for ",
          "polytomous data), and the sample must be large enough for stable ",
          "Rasch model estimation.",
          call. = FALSE
        )
      }

      actual_iterations <- length(successful)
      mean_q3 <- vapply(successful, function(x) x$mean, numeric(1L))
      max_q3  <- vapply(successful, function(x) x$max,  numeric(1L))
      diff_q3 <- max_q3 - mean_q3

      list(
        actual_iterations = actual_iterations,
        sample_n          = sample_n,
        max_diff          = max(diff_q3),
        sd_diff           = stats::sd(diff_q3),
        p95               = stats::quantile(diff_q3, 0.95),
        p99               = stats::quantile(diff_q3, 0.99),
        p995              = stats::quantile(diff_q3, 0.995),
        p999              = stats::quantile(diff_q3, 0.999),
        suggested_cutoff  = stats::quantile(diff_q3, 0.99)
      )
    }
  )
)
