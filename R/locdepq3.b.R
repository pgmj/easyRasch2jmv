#' @export
locdepq3Class <- R6::R6Class(
  "locdepq3Class",
  inherit = locdepq3Base,
  private = list(

    # .init() runs immediately when options change — sets up table structure
    # so column headers appear instantly and don't flicker on re-run.
    .init = function() {
      vars <- self$options$vars
      if (is.null(vars) || length(vars) < 3)
        return()

      q3_table <- self$results$q3Table

      q3_table$addColumn(name = "item", title = "", type = "text")
      for (v in vars) {
        q3_table$addColumn(name = v, title = v, type = "number", format = "zto")
      }
      if (isTRUE(self$options$computeCutoff)) {
        q3_table$addColumn(name = "above_cutoff", title = "Above cutoff", type = "text")
      }

      # Pre-populate rows with em dashes on the diagonal and blanks above
      for (i in seq_along(vars)) {
        row_vals <- list(item = vars[i])
        for (j in seq_along(vars)) {
          if (i == j) {
            row_vals[[ vars[j] ]] <- "—"
          } else if (j > i) {
            row_vals[[ vars[j] ]] <- ""
          }
        }
        q3_table$addRow(rowKey = i, values = row_vals)
      }

      # Pre-create the item-pair table rows: one per pair, count known
      # from the options (choose(k, 2)). The pair table is only shown
      # when the simulation-based cutoff is computed.
      if (isTRUE(self$options$computeCutoff)) {
        k <- length(vars)
        blank <- list(item1 = "", item2 = "", q3 = NA_real_,
                      q3Low = NA_real_, q3High = NA_real_, flagged = "")
        pt <- self$results$pairTable
        for (i in seq_len(k * (k - 1L) / 2L)) {
          pt$addRow(rowKey = i, values = blank)
        }
      }

      # The cutoff summary table has a fixed 6-row structure determined by the
      # design (not by data), so build it here with NA placeholders. .run()
      # fills in the bootstrap results via setRow() once the simulation
      # finishes -- this avoids the blank-then-populated UI jump.
      if (isTRUE(self$options$computeCutoff)) {
        ct <- self$results$cutoffTable
        ct$addRow(rowKey = "sug",  values = list(parameter = "Suggested cutoff (p99)", value = NA_real_))
        ct$addRow(rowKey = "p95",  values = list(parameter = "p95",                    value = NA_real_))
        ct$addRow(rowKey = "p995", values = list(parameter = "p99.5",                  value = NA_real_))
        ct$addRow(rowKey = "p999", values = list(parameter = "p99.9",                  value = NA_real_))
        ct$addRow(rowKey = "iter", values = list(parameter = "Actual iterations",      value = NA_real_))
        ct$addRow(rowKey = "n",    values = list(parameter = "Sample N",               value = NA_real_))
      }
    },

    .run = function() {
      # Return early / explain if requirements not met. With 2 items the
      # unidimensional mirt model is not estimable (too few degrees of
      # freedom), so at least 3 items are required.
      if (is.null(self$options$vars) || length(self$options$vars) == 0) {
        return()
      }
      if (length(self$options$vars) < 3) {
        self$results$q3Note$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the unidimensional model cannot be estimated (too few ",
          "degrees of freedom), so no Q3 residual correlations can be ",
          "computed. Select at least 3 items.</p>"
        ))
        return()
      }

      data            <- self$data
      vars            <- self$options$vars
      compute_cutoff  <- self$options$computeCutoff
      iterations      <- self$options$iterations
      seed_val        <- self$options$seed

      # Select only the specified variables
      df <- data[, vars, drop = FALSE]

      # Robust conversion: handles factors with text labels (SPSS),
      # haven_labelled vectors, and numerics.
      df <- to_numeric_responses_df(df)

      # Check for all-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste(
          "The following variables contain no valid numeric data:",
          paste(bad_vars, collapse = ", ")
        ))
      }

      # Sentinel-value sanity check (e.g., 999, 8888 unmarked as missing)
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

      # Validate data (non-negative integers, min = 0)
      validate_response_data(df)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$q3Table$setNote("sparse", sparse_msg)

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

          # digits = 4 matches the precision used for the simulated Q3
          # values (utils-simulation.R); the previous digits = 2 limited
          # the table display and the cutoff computation to 2 decimals.
          resid_mat <- mirt::residuals(mirt_model, type = "Q3", digits = 4, verbose = FALSE)
          diag(resid_mat) <- NA
          mean_resid <- mean(resid_mat, na.rm = TRUE)

          # --- Step 2: Simulation-based cutoff (optional) -------------------------
          cutoff_val   <- NULL
          dyn_cutoff   <- NULL
          cutoff_res   <- NULL
          sim_fail_msg <- NULL

          if (compute_cutoff) {
            # Seed is always applied (default 42) so results are
            # reproducible by default, consistent with the other
            # simulation-based analyses in this module. If the simulation
            # cannot deliver reliable cutoffs (e.g. too few successful
            # iterations), degrade gracefully: the Q3 matrix is still
            # shown and the note below explains why the cutoff, pair
            # table, and figure are unavailable.
            cutoff_res <- tryCatch(
              private$.runCutoffSim(df, iterations, as.integer(seed_val)),
              error = function(e) {
                sim_fail_msg <<- e$message
                NULL
              }
            )
          }

          if (!is.null(cutoff_res)) {
            cutoff_val  <- cutoff_res$suggested_cutoff
            dyn_cutoff  <- mean_resid + cutoff_val

            # Fill cutoff summary table (rows created in .init())
            ct <- self$results$cutoffTable
            ct$setRow(rowKey = "sug",  values = list(parameter = "Suggested cutoff (p99)", value = cutoff_res$suggested_cutoff))
            ct$setRow(rowKey = "p95",  values = list(parameter = "p95",  value = cutoff_res$p95))
            ct$setRow(rowKey = "p995", values = list(parameter = "p99.5", value = cutoff_res$p995))
            ct$setRow(rowKey = "p999", values = list(parameter = "p99.9", value = cutoff_res$p999))
            ct$setRow(rowKey = "iter", values = list(parameter = "Actual iterations", value = cutoff_res$actual_iterations))
            ct$setRow(rowKey = "n",    values = list(parameter = "Sample N",          value = cutoff_res$sample_n))
            ct$setNote("pctl", paste0(
              "Global cutoff based on all item pairs: percentiles of ",
              "(max Q3 - mean Q3), where the max and mean are taken over ",
              "all pairs within each simulated dataset. The suggested ",
              "cutoff (99th percentile) is added to the observed mean Q3 ",
              "to give the dynamic cut-off applied in the correlation ",
              "matrix above. See the item-pair table below for per-pair ",
              "intervals."
            ))
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

          # --- Step 3b: Save state for the per-pair Q3 plot ----------------------
          if (!is.null(cutoff_res)) {
            # Build observed_pair_df from the upper-triangle of the observed
            # mirt-fitted Q3 matrix, matching the pair labels used in the
            # simulated pair_results.
            obs_upper <- which(upper.tri(resid_mat), arr.ind = TRUE)
            observed_pair_df <- data.frame(
              Item1       = colnames(resid_mat)[obs_upper[, "row"]],
              Item2       = colnames(resid_mat)[obs_upper[, "col"]],
              observed_Q3 = resid_mat[obs_upper],
              stringsAsFactors = FALSE
            )
            observed_pair_df$Pair <- paste(observed_pair_df$Item1, "-",
                                           observed_pair_df$Item2)

            self$results$q3Plot$setState(list(
              pair_results      = cutoff_res$pair_results,
              actual_iterations = cutoff_res$actual_iterations,
              sample_n          = cutoff_res$sample_n,
              n_complete        = n_complete,
              observed_pair_df  = observed_pair_df,
              n_pairs           = self$options$nPairs
            ))
          }

          # --- Step 3c: Item-pair table (simulation-based) -------------------
          if (!is.null(cutoff_res)) {
            pairs_idx <- which(upper.tri(resid_mat), arr.ind = TRUE)
            pair_df <- data.frame(
              Item1 = colnames(resid_mat)[pairs_idx[, "row"]],
              Item2 = colnames(resid_mat)[pairs_idx[, "col"]],
              q3    = resid_mat[pairs_idx],
              stringsAsFactors = FALSE
            )
            pair_df$expected <- NA_real_
            pair_df$q3Low    <- NA_real_
            pair_df$q3High   <- NA_real_
            pair_df$flagged  <- ""

            # Per-pair credible intervals from the simulated null
            # distributions (same approach as easyRasch2's pair_cutoffs:
            # HDCI over the per-pair simulated Q3 values).
            hdci_width <- self$options$hdciWidth / 100
            pr  <- cutoff_res$pair_results
            prk <- paste(pr$Item1, pr$Item2, sep = "___")
            pdk <- paste(pair_df$Item1, pair_df$Item2, sep = "___")
            for (i in seq_len(nrow(pair_df))) {
              sims <- pr$Q3[prk == pdk[i]]
              sims <- sims[is.finite(sims)]
              if (length(sims) >= 2L) {
                iv <- ggdist::hdci(sims, .width = hdci_width)
                pair_df$expected[i] <- stats::median(sims)
                pair_df$q3Low[i]    <- iv[1L, 1L]
                pair_df$q3High[i]   <- iv[1L, 2L]
                if (!is.na(pair_df$q3[i])) {
                  if (pair_df$q3[i] > pair_df$q3High[i]) {
                    pair_df$flagged[i] <- "above"
                  } else if (pair_df$q3[i] < pair_df$q3Low[i]) {
                    pair_df$flagged[i] <- "below"
                  }
                }
              }
            }
            # Sort by deviation from the simulated null, matching the
            # ranking used by the per-pair plot.
            pair_df <- pair_df[order(-abs(pair_df$q3 - pair_df$expected)), ,
                               drop = FALSE]
            rownames(pair_df) <- NULL

            pt <- self$results$pairTable
            for (i in seq_len(nrow(pair_df))) {
              pt$setRow(rowNo = i, values = list(
                item1    = pair_df$Item1[i],
                item2    = pair_df$Item2[i],
                q3       = pair_df$q3[i],
                q3Low    = pair_df$q3Low[i],
                q3High   = pair_df$q3High[i],
                flagged  = pair_df$flagged[i]
              ))
            }
            pt$setNote("flag", paste0(
              "Expected range = ", self$options$hdciWidth, "% HDCI of the ",
              "per-pair Q3 values simulated under local independence. ",
              "Flagged: 'above' = stronger residual association than the ",
              "model predicts (positive local dependence); 'below' = ",
              "weaker than predicted (can indicate multidimensionality). ",
              "These per-pair intervals complement the global cutoff used ",
              "by the tables above. Pairs are sorted by deviation from ",
              "the simulated per-pair median (the black dots in the ",
              "figure), descending."
            ))
          }

          # --- Step 3d: Sample-size note ------------------------------------
          n_total <- sum(rowSums(!is.na(df)) > 0)
          missing_clause <- if (n_total > n_complete) {
            paste0(", of whom ", n_complete, " had complete responses; ",
                   "rows with partially missing responses are retained by ",
                   "MML estimation")
          } else ""
          fail_clause <- if (!is.null(sim_fail_msg)) {
            paste0(" <b>Simulation-based cutoffs unavailable:</b> ",
                   sim_fail_msg)
          } else ""
          self$results$q3Note$setContent(paste0(
            "<p>Q3 residual correlations from a unidimensional ",
            if (max(as.matrix(df), na.rm = TRUE) > 1L) "partial credit"
            else "Rasch",
            " model estimated with MML (mirt) on N = ", n_total,
            " respondents", missing_clause, ". Mean Q3 = ",
            round(mean_resid, 3), ".", fail_clause, "</p>"
          ))

          # --- Step 4: Populate the Q3 table (structure set up in .init()) --------
          # Raw (unrounded) values: the jamovi frontend applies the user's
          # "Number format" preferences.
          q3_table <- self$results$q3Table

          # Add note when cutoff is applied
          if (!is.null(dyn_cutoff)) {
            q3_table$setNote(
              "cutoff",
              paste0(
                "Dynamic cut-off: ", round(dyn_cutoff, 3),
                " (mean Q3 = ", round(mean_resid, 3),
                " + ", round(cutoff_val, 3),
                " (p99 from simulation)). ",
                "Rows marked * contain at least one value above the ",
                "cut-off. This is a global cut-off derived from all item ",
                "pairs; the item-pair table below applies per-pair ",
                "intervals instead."
              )
            )
          }

          # Populate rows — diagonal gets em dash, upper triangle gets blank,
          # lower triangle gets the Q3 value
          for (i in seq_along(vars)) {
            row_vals <- list(item = vars[i])
            for (j in seq_along(vars)) {
              if (i == j) {
                row_vals[[ vars[j] ]] <- "—"
              } else if (j > i) {
                row_vals[[ vars[j] ]] <- ""
              } else {
                row_vals[[ vars[j] ]] <- resid_mat[i, j]
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

      item_names_vec <- colnames(data_mat)
      if (is.null(item_names_vec)) {
        item_names_vec <- paste0("V", seq_len(ncol(data_mat)))
      }

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
          type       = "polytomous",
          thetas     = thetas,
          deltaslist = deltaslist,
          n_items    = ncol(data_mat),
          sample_n   = sample_n,
          item_names = item_names_vec
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
          sample_n    = sample_n,
          item_names  = item_names_vec
        )
      }

      results_raw <- run_q3_sim_sequential(iterations, sim_seeds, sim_data_list,
                                           verbose = FALSE)

      ok         <- vapply(results_raw, function(x) is.list(x), logical(1L))
      successful <- results_raw[ok]

      # Guard against degenerate cutoffs: with very few successful
      # iterations the global percentiles and the per-pair HDCIs collapse
      # and flagging becomes meaningless. Require at least 20 successes
      # and a 50% success rate; otherwise report the dominant failure
      # reason.
      n_ok <- length(successful)
      if (n_ok < 20L || n_ok < iterations / 2) {
        fail_msgs <- unlist(results_raw[!ok])
        top_reason <- if (length(fail_msgs) > 0L) {
          names(sort(table(fail_msgs), decreasing = TRUE))[1L]
        } else NULL
        stop(paste0(
          "Only ", n_ok, " of ", iterations, " simulation iterations ",
          "succeeded -- too few to estimate reliable cutoffs.",
          if (!is.null(top_reason))
            paste0(" Most common failure: ", top_reason, ".") else "",
          " Check your data: items must have sufficient response ",
          "variation (at least 8 positive responses per item for ",
          "dichotomous data; all response categories represented for ",
          "polytomous data), and the sample must be large enough for ",
          "stable Rasch model estimation."
        ), call. = FALSE)
      }

      actual_iterations <- length(successful)
      mean_q3 <- vapply(successful, function(x) x$mean, numeric(1L))
      max_q3  <- vapply(successful, function(x) x$max,  numeric(1L))
      diff_q3 <- max_q3 - mean_q3

      # Per-pair aggregation: stack pair_q3 frames + add iteration index
      pair_iter_dfs <- lapply(seq_along(successful), function(i) {
        d <- successful[[i]]$pair_q3
        d$iteration <- i
        d
      })
      pair_results <- do.call(rbind, pair_iter_dfs)
      rownames(pair_results) <- NULL

      list(
        actual_iterations = actual_iterations,
        sample_n          = sample_n,
        item_names        = item_names_vec,
        pair_results      = pair_results,
        max_diff          = max(diff_q3),
        sd_diff           = stats::sd(diff_q3),
        p95               = stats::quantile(diff_q3, 0.95),
        p99               = stats::quantile(diff_q3, 0.99),
        p995              = stats::quantile(diff_q3, 0.995),
        p999              = stats::quantile(diff_q3, 0.999),
        suggested_cutoff  = stats::quantile(diff_q3, 0.99)
      )
    },

    # ---------------------------------------------------------------------
    # Per-pair Q3 simulation plot — mirrors iteminfit's plot design
    # (ggdist::stat_dots dot cloud + black per-pair median + orange diamonds
    # for the observed Q3 from the mirt fit).
    # ---------------------------------------------------------------------
    .q3Plot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
      if (!requireNamespace("ggdist", quietly = TRUE))  return(FALSE)
      if (!requireNamespace("scales", quietly = TRUE))  return(FALSE)

      state <- image$state
      pair_results      <- state$pair_results
      actual_iterations <- state$actual_iterations
      sample_n          <- state$sample_n
      n_complete        <- state$n_complete
      observed_pair_df  <- state$observed_pair_df   # data.frame: Pair, observed_Q3
      n_pairs           <- state$n_pairs            # NULL or positive integer

      pair_results$Pair <- paste(pair_results$Item1, "-", pair_results$Item2)

      # --- Top-N filter by |observed Q3 - median(simulated Q3 per pair)| -------
      # Mirrors RMlocdepGammaPlot's ranking: rank pairs by how far the observed
      # value lies from the simulation null, keep the top N, and order the
      # y-axis so the most-deviant pair sits at the top of the plot.
      if (is.numeric(n_pairs) && length(n_pairs) == 1L && n_pairs >= 1L) {
        n_pairs <- as.integer(n_pairs)
        pair_names_all <- unique(pair_results$Pair)
        med_sim <- vapply(pair_names_all, function(pp) {
          stats::median(pair_results$Q3[pair_results$Pair == pp], na.rm = TRUE)
        }, numeric(1L))
        obs_lookup <- stats::setNames(observed_pair_df$observed_Q3,
                                      observed_pair_df$Pair)
        deviation  <- abs(obs_lookup[pair_names_all] - med_sim[pair_names_all])
        ord        <- order(deviation, decreasing = TRUE)
        keep_n     <- min(n_pairs, length(pair_names_all))
        keep_pairs <- pair_names_all[ord[seq_len(keep_n)]]

        pair_results     <- pair_results[pair_results$Pair %in% keep_pairs, ,
                                         drop = FALSE]
        observed_pair_df <- observed_pair_df[observed_pair_df$Pair %in%
                                             keep_pairs, , drop = FALSE]
        pair_levels      <- rev(keep_pairs)  # largest deviation at the top
      } else {
        pair_levels <- rev(unique(pair_results$Pair))
      }

      lo_hi <- do.call(rbind, lapply(unique(pair_results$Pair), function(pp) {
        sub <- pair_results[pair_results$Pair == pp, ]
        data.frame(
          Pair      = pp,
          min_q3    = stats::quantile(sub$Q3, 0.001, na.rm = TRUE),
          max_q3    = stats::quantile(sub$Q3, 0.999, na.rm = TRUE),
          p66lo_q3  = stats::quantile(sub$Q3, 0.167, na.rm = TRUE),
          p66hi_q3  = stats::quantile(sub$Q3, 0.833, na.rm = TRUE),
          median_q3 = stats::median  (sub$Q3,         na.rm = TRUE),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }))
      lo_hi$Pair_f <- factor(lo_hi$Pair, levels = pair_levels)

      q3_sim <- data.frame(
        Pair  = pair_results$Pair,
        Value = pair_results$Q3,
        stringsAsFactors = FALSE
      )
      q3_sim <- merge(q3_sim, observed_pair_df[, c("Pair", "observed_Q3")],
                      by = "Pair", sort = FALSE)
      q3_sim$Pair <- factor(q3_sim$Pair, levels = pair_levels)

      # sample_n counts all rows used by the MML fit (rows with partially
      # missing responses are retained), so it is the accurate n for both
      # the observed and the simulated Q3 values.
      caption_text <- er2_caption(paste0(
        "Observed and simulated Q3 are based on n = ", sample_n,
        " respondents.\n",
        "Simulated distributions: ", actual_iterations,
        " parametric-bootstrap datasets using the same n.\n",
        "Orange diamonds: observed Q3. Black dots: simulation median."
      ))

      p <- ggplot2::ggplot(q3_sim, ggplot2::aes(x = .data$Value, y = .data$Pair)) +
        ggdist::stat_dots(
          ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
          quantiles  = actual_iterations,
          layout     = "weave",
          slab_color = NA,
          .width     = c(0.666, 0.999)
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$min_q3, xend = .data$max_q3,
                       y = .data$Pair_f, yend = .data$Pair_f),
          color = "black", linewidth = 0.7
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$p66lo_q3, xend = .data$p66hi_q3,
                       y = .data$Pair_f, yend = .data$Pair_f),
          color = "black", linewidth = 1.2
        ) +
        ggplot2::geom_point(
          data = lo_hi,
          ggplot2::aes(x = .data$median_q3, y = .data$Pair_f),
          size = 3.6
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = .data$observed_Q3),
          color = "sienna2", shape = 18,
          position = ggplot2::position_nudge(y = -0.1),
          size = 7
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype   = "dashed",
          color      = "grey50",
          linewidth  = 0.4
        ) +
        ggplot2::labs(x = "Q3 residual correlation", y = "Item pair",
                      caption = caption_text) +
        ggplot2::scale_color_manual(
          values     = scales::brewer_pal()(3)[-1],
          aesthetics = "slab_fill", guide = "none"
        ) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(
          panel.spacing = ggplot2::unit(0.7, "cm")
        ) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
