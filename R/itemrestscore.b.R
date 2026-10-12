#' @export
itemrestscoreClass <- R6::R6Class(
  "itemrestscoreClass",
  inherit = itemrestscoreBase,
  private = list(
    .init = function() {
      # The bootstrap adjusted-p column names its correction method, as in
      # the infit analysis. The asymptotic column keeps its fixed BH title.
      self$results$restscoreTable$getColumn("pAdjBoot")$setTitle(
        padjusted_title(self$options$correction)
      )
    },

    .run = function() {
      # Return early / explain if requirements not met.
      # With 2 items the restscore degenerates to the other item of the
      # pair and observed = expected exactly (p = 1) for both items, so
      # at least 3 items are required for a meaningful analysis.
      if (is.null(self$options$vars) || length(self$options$vars) == 0) {
        return()
      }
      if (length(self$options$vars) < 3) {
        self$results$restscoreNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the restscore (total score minus the item) reduces to the ",
          "other item, making observed and expected correlations identical ",
          "by construction. Select at least 3 items.</p>"
        ))
        return()
      }

      data <- self$data
      vars <- self$options$vars
      sort_by_diff <- self$options$sortByDiff

      # Select the specified variables (drop = FALSE keeps it as data.frame)
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$restscoreTable$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$restscoreTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$restscoreTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item are dropped up
      # front: the bundled easyRasch2 release cannot fit all-NA rows
      # (psychotools errors on polytomous and crashes on dichotomous data).
      n_total <- nrow(df)
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]

      # Sufficient complete cases?
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. The item-restscore statistics require at least one row with responses to all selected items.")
      }

      tryCatch(
        {
          # Simulation-based cutoffs. All computation is delegated to the
          # easyRasch2 package, so results are numerically identical to
          # RMitemRestscore(df, cutoff = RMitemRestscoreCutoff(df, ...))
          # with the same seed and iterations. If the simulation cannot
          # deliver usable results, the observed table is still shown and
          # a notice (sim_failure_notice()) says why the simulation part is
          # missing.
          #
          # Cached in the plot state (see run_cutoff_simulation()), so changing
          # pValues, correction or sortByDiff does not rerun the simulation.
          sim_sig <- list(
            iterations = self$options$iterations,
            seed       = as.integer(self$options$seed),
            hdci_width = self$options$hdciWidth / 100
          )
          sim <- run_cutoff_simulation(
            compute   = self$options$computeCutoff,
            cached    = self$results$restscorePlot$state,
            sig       = sim_sig,
            keys      = list(df = df),
            simulate  = function() easyRasch2::RMitemRestscoreCutoff(
              df,
              iterations = self$options$iterations,
              parallel   = FALSE,
              seed       = as.integer(self$options$seed),
              hdci_width = self$options$hdciWidth / 100
            ),
            requested = self$options$iterations,
            hint      = sim_rule_hint("positive")
          )
          cutoff_res   <- sim$cutoff_res
          sim_fail_msg <- sim$fail_msg

          # The pValues option can hold a stale TRUE while greyed out
          # (jamovi disables but does not reset nested options), hence the
          # explicit computeCutoff gate.
          use_pvalues <- isTRUE(self$options$computeCutoff) &&
            isTRUE(self$options$pValues) &&
            !is.null(cutoff_res)

          # Visibility follows what was computed, not the options: if the
          # simulation failed, the simulation columns would otherwise show
          # empty (see set_columns_visible()).
          sim_ok <- !is.null(cutoff_res)
          tbl <- self$results$restscoreTable
          set_columns_visible(tbl, c("diffLow", "diffHigh"), sim_ok)
          set_columns_visible(tbl, c("pValue", "pAdjBoot"), use_pvalues)
          # With a cutoff the package drops the asymptotic p-value.
          set_columns_visible(tbl, "pAdjusted", !sim_ok)
          set_elements_visible(list(self$results$restscorePlot), sim_ok)

          # Without a cutoff this is the asymptotic test with hardcoded BH
          # adjustment (module-wide convention, also the package default).
          # Package warnings (e.g. sparse categories) are suppressed -- the
          # module surfaces its own sparse-category footnote above.
          results <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemRestscore(
              df,
              cutoff     = cutoff_res,
              p_value    = if (is.null(cutoff_res)) NULL else use_pvalues,
              correction = self$options$correction,
              output     = "dataframe"
            )
          ))

          # Sort by absolute magnitude when requested, so both over- and
          # underfit items rise to the top while the signed value remains
          # visible in the table.
          if (isTRUE(sort_by_diff)) {
            results <- results[order(abs(results$Difference), decreasing = TRUE), ]
            rownames(results) <- NULL
          }

          private$.fillTable(results, cutoff_res, use_pvalues)
          private$.setTableNotes(results, cutoff_res, use_pvalues)

          if (!is.null(sim_fail_msg)) {
            sim_failure_notice(self, sim_fail_msg,
                               "the asymptotic p-values instead")
          }

          private$.writeNote(df, n_total, n_complete, nrow(results),
                             cutoff_res, use_pvalues)

          # Save state for the plot. The figure is drawn by
          # easyRasch2::RMitemRestscorePlot() in the render function.
          # Storing the data + cutoff object keeps the saved analysis small
          # and doubles as the simulation cache above.
          if (!is.null(cutoff_res)) {
            self$results$restscorePlot$setState(list(
              df         = df,
              cutoff_res = cutoff_res,
              sig        = sim_sig
            ))
          }
        },
        error = function(e) {
          stop(paste("Error in item-restscore analysis:", e$message))
        }
      )
    },

    # ---------------------------------------------------------------------
    # Table rows. Asymptotic p-values without a cutoff object; expected
    # range with one, plus bootstrap p-values when requested.
    # ---------------------------------------------------------------------
    .fillTable = function(results, cutoff_res, use_pvalues) {
      table <- self$results$restscoreTable
      for (i in seq_len(nrow(results))) {
        vals <- list(
          item        = results$Item[i],
          observed    = results$Observed[i],
          expected    = results$Expected[i],
          difference  = results$Difference[i],
          fit         = results$Flagged[i],
          relLocation = results$Relative_location[i]
        )
        if (is.null(cutoff_res)) {
          vals$pAdjusted <- results$p_adjusted[i]
        } else {
          vals$diffLow  <- results$Diff_low[i]
          vals$diffHigh <- results$Diff_high[i]
        }
        if (use_pvalues) {
          vals$pValue   <- results$p_restscore[i]
          vals$pAdjBoot <- results$padj_restscore[i]
        }
        table$setRow(rowNo = i, values = vals)
      }
    },

    # ---------------------------------------------------------------------
    # Table footnotes. Three paths: bootstrap p-values, expected range, or
    # the asymptotic test (which also carries the calibration caveat).
    # ---------------------------------------------------------------------
    .setTableNotes = function(results, cutoff_res, use_pvalues) {
      table <- self$results$restscoreTable
      diff_txt <- paste0(
        "Difference = observed - expected gamma. Positive: ",
        "over-discrimination (overfit, often local dependence); ",
        "negative: under-discrimination (underfit, often ",
        "multidimensionality or noise)."
      )
      if (use_pvalues) {
        table$setNote("diff", paste0(
          diff_txt, " Items are labelled in the Flagged column when the ",
          "adjusted p-value < .05, with the direction taken from the ",
          "side of the simulated mean the difference falls on. In small ",
          "samples that mean is slightly above 0, so an underfitting ",
          "item can occasionally show a small positive difference."
        ))
        table$setNote("pvalues", paste0(
          "p-value: probability of a difference at least as far from ",
          "its simulated mean as observed if the model fits, computed ",
          "from the ", cutoff_res$actual_iterations, " simulated ",
          "datasets (Monte-Carlo). Adj. p-value: corrected for multiple ",
          "comparisons across the ", nrow(results), " items using ",
          correction_label(self$options$correction), "."
        ))
      } else if (!is.null(cutoff_res)) {
        table$setNote("diff", paste0(
          diff_txt, " Items are labelled in the Flagged column when the ",
          "difference falls outside the expected range: above = ",
          "overfit, below = underfit."
        ))
      } else {
        table$setNote("diff", paste0(
          diff_txt, " Items are labelled in the Flagged column when the ",
          "adjusted p-value < .05."
        ))
        table$setNote("sig", paste0(
          "P-values adjusted with the Benjamini-Hochberg (BH) ",
          "false-discovery-rate method."
        ))
        # The asymptotic test divides by the SE of the observed gamma
        # alone, which is too large for the observed minus expected
        # difference, and the difference is biased upwards in small
        # samples. See easyRasch2::RMitemRestscore() and
        # easyRasch2/dev/restscore_asymptotic_null.qmd.
        table$setNote("calibration", paste0(
          "The p-values are asymptotic and miscalibrated under the ",
          "Rasch model. In small samples too many items are flagged as ",
          "overfit, and at any sample size too few are flagged as ",
          "underfit. In simulation, 20 dichotomous items at N = 150 ",
          "gave at least one overfit flag in 13% of datasets that ",
          "fitted the model perfectly, while underfitting items need to ",
          "misfit clearly before they are flagged.",
          # A failed simulation is reported by the notice instead, so
          # the advice is only for users who have not asked for it.
          if (!isTRUE(self$options$computeCutoff)) {
            paste0(" Enable Simulation-based cutoffs for p-values that ",
                   "are calibrated to the data.")
          }
        ))
      }
      table$setNote("loc", paste0(
        "Rel. location = mean item (threshold) location relative to ",
        "the mean person location (weighted likelihood estimates, ",
        "WLE), in logits."
      ))
      table$setNote("bidirectional", bidirectional_flag_note(results$Flagged))
    },

    # ---------------------------------------------------------------------
    # Html note below the table: sample size, missing-data handling and the
    # simulation basis.
    # ---------------------------------------------------------------------
    .writeNote = function(df, n_total, n_complete, n_rows, cutoff_res,
                          use_pvalues) {
      # Sample-size / missing-data note. The two parts of the table
      # use different samples when responses are missing: iarm refits
      # the model on complete cases for the restscore statistics,
      # while the Location columns come from the CML fit on all
      # available responses (partial rows retained). This mirrors
      # easyRasch2::RMitemRestscore(). Rows with no valid responses at
      # all were dropped in .run() and are reported separately.
      n_used <- nrow(df)
      drop_clause <- if (n_used < n_total) {
        paste0(" (", n_total - n_used, " row(s) without any responses ",
               "on the selected items excluded)")
      } else ""
      missing_clause <- if (n_used > n_complete) {
        paste0(", of whom ", n_complete, " had complete responses on ",
               "all ", ncol(df), " items. Restscore correlations and ",
               "p-values are computed from the ", n_complete,
               " complete cases (the model is refitted on complete ",
               "cases for this purpose); item locations use all ",
               "available responses via conditional maximum ",
               "likelihood estimation.")
      } else {
        paste0(", all with complete responses on all ", ncol(df),
               " items.")
      }

      # Always rewrite the note on the success path: jamovi does not
      # clear set HTML content, so a stale guard message would survive.
      sim_html <- if (!is.null(cutoff_res)) {
        paste0(
          "<p>Expected ranges and p-values based on ",
          cutoff_res$actual_iterations, " simulated datasets (",
          cutoff_res$hdci_width * 100, "% HDCI) drawn from the same n = ",
          n_complete, " complete cases, with the model refitted in each.",
          # With p-values the caveat carries the iteration advice, so
          # the general recommendation would repeat it.
          if (!use_pvalues) {
            iteration_note(self$options$iterations, 400L)
          },
          iteration_attrition_note(
            cutoff_res$actual_iterations,
            self$options$iterations
          ),
          if (use_pvalues) {
            pvalue_iteration_caveat(cutoff_res$actual_iterations,
                                    floor = 400L, count_stated = TRUE)
          } else {
            # Flagging falls back to the expected range, whose width
            # sets a familywise error rate the user has not chosen.
            interval_flagging_note(cutoff_res$hdci_width, n_rows)
          },
          "</p>"
        )
      } else ""
      self$results$restscoreNote$setContent(paste0(
        "<p>Analysis based on N = ", n_used, " respondents",
        drop_clause, missing_clause, "</p>", sim_html
      ))
    },

    # ---------------------------------------------------------------------
    # Simulated difference plot -- easyRasch2::RMitemRestscorePlot() (ggdist
    # dot cloud + black per-item median + orange diamonds for the observed
    # difference), restyled to the module's plot conventions (base size 15).
    # ---------------------------------------------------------------------
    .restscorePlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMitemRestscorePlot(
          image$state$cutoff_res,
          data = image$state$df
        )
      ))
      if (is.null(p)) return(FALSE)

      p <- p +
        ggplot2::theme_minimal(base_size = 15) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
