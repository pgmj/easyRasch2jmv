#' @export
locdepq3Class <- R6::R6Class(
  "locdepq3Class",
  inherit = locdepq3Base,
  private = list(

    # .init() runs immediately when options change — sets up table structure
    # so column headers appear instantly and don't flicker on re-run.
    .init = function() {
      # The adjusted-p column names its correction method, extending the
      # "Adj. p-value (BH)" title convention of the asymptotic p columns.
      self$results$pairTable$getColumn("pAdjusted")$setTitle(
        padjusted_title(self$options$correction)
      )

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
      # unidimensional Rasch model has too few degrees of freedom for a
      # meaningful residual analysis, so at least 3 items are required.
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
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      # Respondents with no responses on any selected item are dropped up
      # front: easyRasch2::RMlocdepQ3() would drop them internally anyway,
      # but RMlocdepQ3Cutoff() currently errors on all-NA rows (psychotools
      # cannot fit them), so both calls must see the same reduced sample.
      n_total <- nrow(df)
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
      n_used <- nrow(df)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$q3Table$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$q3Table$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$q3Table$setNote("duplicate", dup_msg)

      tryCatch(
        {
          # --- Step 1: Simulation-based cutoff (optional) ------------------
          # All computation is delegated to the easyRasch2 package (CML item
          # estimation via psychotools, WLE person locations), so results are
          # numerically identical to RMlocdepQ3Cutoff() / RMlocdepQ3() with
          # the same seed and iterations. Seed is always applied (default 42)
          # so results are reproducible by default, consistent with the other
          # simulation-based analyses in this module. If the simulation cannot
          # deliver reliable cutoffs (e.g. too few successful iterations),
          # degrade gracefully: the Q3 matrix is still shown and the note
          # below explains why the cutoff, pair table, and figures are
          # unavailable. Package warnings (e.g. sparse categories) are
          # suppressed -- the module surfaces its own sparse-category
          # footnote above.
          cutoff_res   <- NULL
          sim_fail_msg <- NULL

          # The simulation is the expensive part, and jamovi reruns .run on
          # every option change -- including changes (pValues, correction,
          # nPairs) that do not affect the simulation. The heatmap state
          # carries the cutoff object and jamovi clears it exactly when a
          # simulation-relevant option changes (its clearWith list; the
          # per-pair plot is unsuitable as the cache because its state is
          # also cleared by nPairs). The signature check makes reuse
          # self-validating rather than relying on clearWith alone -- in
          # particular, hdciWidth changes the cutoff object but is not in
          # the heatmap's clearWith, and the signature catches it.
          sim_sig <- list(
            iterations = iterations,
            seed       = as.integer(seed_val),
            hdci_width = self$options$hdciWidth / 100
          )
          cached <- self$results$matrixPlot$state
          if (compute_cutoff &&
              !is.null(cached) && !is.null(cached$cutoff_res) &&
              identical(cached$sig, sim_sig) && identical(cached$df, df)) {
            cutoff_res <- cached$cutoff_res
          } else if (compute_cutoff) {
            cutoff_res <- tryCatch(
              suppressWarnings(suppressMessages(
                easyRasch2::RMlocdepQ3Cutoff(
                  df,
                  iterations = iterations,
                  parallel   = FALSE,
                  seed       = as.integer(seed_val),
                  hdci_width = self$options$hdciWidth / 100
                )
              )),
              error = function(e) {
                sim_fail_msg <<- e$message
                NULL
              }
            )
            # Guard against degenerate cutoffs: with very few successful
            # iterations the global percentiles and the per-pair HDCIs
            # collapse and flagging becomes meaningless.
            if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
              sim_fail_msg <- paste0(
                "Only ", cutoff_res$actual_iterations, " of ", iterations,
                " simulation iterations succeeded -- too few to estimate ",
                "reliable cutoffs. Check your data: items must have ",
                "sufficient response variation (at least 8 positive ",
                "responses per item for dichotomous data; all response ",
                "categories represented for polytomous data), and the ",
                "sample must be large enough for stable Rasch model ",
                "estimation."
              )
              cutoff_res <- NULL
            }
          }

          # --- Step 2: Observed Q3 matrix (and pair table with cutoff) -----
          # Bootstrap p-values need the full cutoff object (it carries the
          # simulated per-pair distributions), so they are only computed
          # when the simulation succeeded. The pValues option can hold a
          # stale TRUE while greyed out (jamovi disables but does not reset
          # nested options), hence the explicit computeCutoff gate. The
          # package's below-1000-iterations warning is suppressed with the
          # rest; the module states the same caveat in the note below.
          use_pvalues <- compute_cutoff &&
            isTRUE(self$options$pValues) &&
            !is.null(cutoff_res)

          if (!is.null(cutoff_res)) {
            res <- suppressWarnings(suppressMessages(easyRasch2::RMlocdepQ3(
              df, cutoff = cutoff_res, output = "dataframe",
              p_value = use_pvalues, correction = self$options$correction
            )))
            matrix_df <- res$matrix
            pairs_df  <- res$pairs
          } else {
            matrix_df <- suppressWarnings(suppressMessages(easyRasch2::RMlocdepQ3(
              df, output = "dataframe"
            )))
            pairs_df  <- NULL
          }

          # Lower-triangle matrix (unrounded) and its symmetric completion
          # for per-pair lookups. Mean of the lower triangle equals the mean
          # off-diagonal Q3 (the matrix is symmetric).
          resid_mat <- as.matrix(matrix_df[vars])
          full_mat  <- resid_mat
          full_mat[upper.tri(full_mat)] <- t(full_mat)[upper.tri(full_mat)]
          mean_resid <- mean(resid_mat, na.rm = TRUE)

          cutoff_val <- NULL
          dyn_cutoff <- NULL
          if (!is.null(cutoff_res)) {
            cutoff_val <- as.numeric(cutoff_res$suggested_cutoff)
            dyn_cutoff <- mean_resid + cutoff_val
          }

          # --- Step 3: Cutoff summary table (rows created in .init()) ------
          if (!is.null(cutoff_res)) {
            ct <- self$results$cutoffTable
            ct$setRow(rowKey = "sug",  values = list(parameter = "Suggested cutoff (p99)", value = as.numeric(cutoff_res$suggested_cutoff)))
            ct$setRow(rowKey = "p95",  values = list(parameter = "p95",   value = as.numeric(cutoff_res$p95)))
            ct$setRow(rowKey = "p995", values = list(parameter = "p99.5", value = as.numeric(cutoff_res$p995)))
            ct$setRow(rowKey = "p999", values = list(parameter = "p99.9", value = as.numeric(cutoff_res$p999)))
            ct$setRow(rowKey = "iter", values = list(parameter = "Actual iterations", value = cutoff_res$actual_iterations))
            ct$setRow(rowKey = "n",    values = list(parameter = "Sample N",          value = cutoff_res$sample_n))
            ct$setNote("pctl", paste0(
              "Global cutoff based on all item pairs: percentiles of ",
              "(max Q3 - mean Q3), where the max and mean are taken over ",
              "all pairs within each simulated dataset. The suggested ",
              "cutoff (99th percentile) is added to the observed mean Q3 ",
              "to give the dynamic cut-off applied in the correlation ",
              "matrix above. See the item-pair table below for per-pair ",
              "intervals.",
              iteration_note(iterations, 400L, corrected = TRUE),
              iteration_attrition_note(cutoff_res$actual_iterations,
                                       iterations)
            ))
          }

          # --- Step 3b: Save state for the plots ---------------------------
          # Both figures are drawn by easyRasch2::RMlocdepQ3Plot() inside the
          # render functions; the observed CML/WLE refit it performs there is
          # cheap. Storing the data + cutoff object (rather than ggplot
          # objects) keeps the saved analysis small.
          if (!is.null(cutoff_res)) {
            # sig makes the heatmap state double as the simulation cache
            # consulted at the top of this function.
            plot_state <- list(df = df, cutoff_res = cutoff_res,
                               sig = sim_sig)
            self$results$q3Plot$setState(plot_state)
            self$results$matrixPlot$setState(plot_state)
          }

          # --- Step 3c: Item-pair table (simulation-based) -----------------
          # Row order and flags come from RMlocdepQ3()$pairs (sorted by
          # departure from the per-pair simulated median); the displayed
          # values are joined back unrounded from the cutoff object and the
          # observed matrix, per the module's raw-values + format "zto"
          # convention ($pairs pre-rounds to 3 decimals).
          if (!is.null(pairs_df)) {
            pc  <- cutoff_res$pair_cutoffs
            pck <- paste(pc$Item1, pc$Item2, sep = "___")
            pdk <- paste(pairs_df$Item1, pairs_df$Item2, sep = "___")
            idx <- match(pdk, pck)

            pt <- self$results$pairTable
            for (i in seq_len(nrow(pairs_df))) {
              vals <- list(
                item1   = pairs_df$Item1[i],
                item2   = pairs_df$Item2[i],
                q3      = full_mat[pairs_df$Item1[i], pairs_df$Item2[i]],
                q3Low   = pc$Q3_low[idx[i]],
                q3High  = pc$Q3_high[idx[i]],
                flagged = pairs_df$Flagged[i]
              )
              if (use_pvalues) {
                vals$pValue    <- pairs_df$p_q3[i]
                vals$pAdjusted <- pairs_df$padj_q3[i]
              }
              pt$setRow(rowNo = i, values = vals)
            }
            if (use_pvalues) {
              # With bootstrap p-values, flagging follows the adjusted
              # p-value (upstream behavior; one-sided test for excess
              # positive local dependence, so only 'above' pairs are
              # flagged); the expected-range columns stay as the
              # effect-size reference.
              pt$setNote("flag", paste0(
                "Flagged: adjusted p-value < 0.05, indicating stronger ",
                "residual association than the model predicts (positive ",
                "local dependence; the test is one-sided, so only 'above' ",
                "pairs are flagged). The expected-range columns remain as ",
                "the effect-size reference. Pairs are sorted by deviation ",
                "from the simulated per-pair median (the black dots in ",
                "the figure), descending."
              ))
              pt$setNote("pvalues", paste0(
                "p-value: probability of a Q3 at least as large as ",
                "observed under local independence, computed from the ",
                cutoff_res$actual_iterations, " simulated datasets ",
                "(Monte-Carlo, one-sided). Adj. p-value: corrected for ",
                "multiple comparisons across the ", nrow(pairs_df),
                " item pairs using ",
                correction_label(self$options$correction), "."
              ))
            } else {
              pt$setNote("flag", paste0(
                "Expected range = ", self$options$hdciWidth, "% HDCI of the ",
                "per-pair Q3 values simulated under local independence. ",
                "Flagged: 'above' = stronger residual association than the ",
                "model predicts (positive local dependence); 'below' = ",
                "weaker than predicted (can indicate multidimensionality). ",
                "These per-pair intervals complement the global cutoff used ",
                "by the tables above. Pairs are sorted by deviation from ",
                "the simulated per-pair median (the black dots in the ",
                "figure), descending.",
                # Flagging falls back to the expected range, whose width sets
                # a familywise rate over pairs that the user has not chosen.
                interval_flagging_note(cutoff_res$hdci_width, nrow(pairs_df),
                                       unit = "item pairs")
              ))
            }
          }

          # --- Step 3d: Sample-size note ------------------------------------
          # Respondents with no responses were dropped above (n_total vs
          # n_used); incomplete response patterns are retained (CML/WLE
          # handles missingness directly).
          n_complete <- sum(stats::complete.cases(df))
          missing_clause <- if (n_complete < n_used) {
            paste0(", of whom ", n_complete, " had complete responses; ",
                   "rows with partially missing responses are retained by ",
                   "CML/WLE estimation")
          } else ""
          drop_clause <- if (n_used < n_total) {
            paste0(" (", n_total - n_used, " row(s) without any responses ",
                   "on the selected items excluded)")
          } else ""
          fail_clause <- if (!is.null(sim_fail_msg)) {
            paste0(" <b>Simulation-based cutoffs unavailable:</b> ",
                   sim_fail_msg)
          } else ""
          pvalue_clause <- if (use_pvalues) {
            pvalue_iteration_caveat(cutoff_res$actual_iterations)
          } else ""
          self$results$q3Note$setContent(paste0(
            "<p>Q3 residual correlations from a unidimensional ",
            if (max(as.matrix(df), na.rm = TRUE) > 1L) "partial credit"
            else "Rasch",
            " model, computed with the easyRasch2 R package: conditional ",
            "maximum likelihood (CML) item estimation with weighted ",
            "likelihood (WLE) person estimates, on n = ", n_used,
            " respondents", drop_clause, missing_clause, ". Mean Q3 = ",
            round(mean_resid, 3), ".", pvalue_clause, fail_clause, "</p>"
          ))

          # --- Step 4: Populate the Q3 table (structure set up in .init()) --
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
          # lower triangle gets the Q3 value. The above-cutoff row flag comes
          # from RMlocdepQ3()$matrix (lower triangle vs the dynamic cut-off).
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
              row_vals[["above_cutoff"]] <-
                if (isTRUE(matrix_df$above_cutoff[i])) "*" else ""
            }
            q3_table$setRow(rowNo = i, values = row_vals)
          }
        },
        error = function(e) {
          stop(paste("Error in Q3 analysis:", e$message))
        }
      )
    },

    # ---------------------------------------------------------------------
    # Per-pair Q3 simulation plot — easyRasch2::RMlocdepQ3Plot()$pairs
    # (ggdist dot cloud + black per-pair median + orange diamonds for the
    # observed Q3), restyled to the module's plot conventions (base size 15).
    # ---------------------------------------------------------------------
    .q3Plot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      plots <- suppressWarnings(suppressMessages(
        easyRasch2::RMlocdepQ3Plot(
          image$state$cutoff_res,
          data    = image$state$df,
          n_pairs = self$options$nPairs
        )
      ))
      if (is.null(plots$pairs)) return(FALSE)

      # Re-apply the module's theme on top of the package theme: base size
      # 15 for jamovi's 500x600 canvas, then restore the theme pieces that
      # theme_minimal() resets.
      p <- plots$pairs +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    },

    # ---------------------------------------------------------------------
    # Lower-triangle Q3 heatmap — easyRasch2::RMlocdepQ3Plot()$matrix
    # (diverging fill centred on the mean off-diagonal Q3; pairs above the
    # global dynamic cut-off outlined in black).
    # ---------------------------------------------------------------------
    .matrixPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      plots <- suppressWarnings(suppressMessages(
        easyRasch2::RMlocdepQ3Plot(
          image$state$cutoff_res,
          data = image$state$df
        )
      ))
      if (is.null(plots$matrix)) return(FALSE)

      print(er2_bump_text(plots$matrix))
      TRUE
    }
  )
)
