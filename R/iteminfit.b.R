#' @export
iteminfitClass <- R6::R6Class(
  "iteminfitClass",
  inherit = iteminfitBase,
  private = list(
    .init = function() {
      # The adjusted-p column names its correction method, extending the
      # "Adj. p-value (BH)" title convention of the asymptotic p columns.
      self$results$infitTable$getColumn("pAdjusted")$setTitle(
        padjusted_title(self$options$correction)
      )
    },

    .run = function() {
      # 1. Return early / explain if requirements not met. With 2 items
      # the conditional infit is ~1 for both items by construction (no
      # degrees of freedom left for misfit once the total score is
      # conditioned on), so at least 3 items are required.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$cutoffNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the conditional infit equals 1 for both items by ",
          "construction and carries no information about item fit. ",
          "Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Extract and validate data
      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$infitTable$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$infitTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$infitTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item are dropped up
      # front: the bundled easyRasch2 release cannot fit all-NA rows
      # (psychotools errors on polytomous and crashes on dichotomous data).
      n_total <- nrow(df)
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]

      n_complete <- sum(complete.cases(df))
      n_excluded <- n_total - n_complete

      if (n_complete == 0)
        stop("No complete cases found in the data. Conditional infit requires at least one row with responses to all selected items.")

      # 3. Compute infit (and optional cutoffs) via easyRasch2
      tryCatch({
        # All computation is delegated to the easyRasch2 package (CML item
        # estimation via psychotools, WLE person locations), so results are
        # numerically identical to RMitemInfit() / RMitemInfitCutoff() with
        # the same seed and iterations. If the simulation cannot deliver
        # reliable cutoffs (e.g. too few successful iterations), degrade
        # gracefully: the observed infit table is still shown and the note
        # below explains why the expected range and figure are unavailable.
        # Package warnings (e.g. sparse categories) are suppressed -- the
        # module surfaces its own sparse-category footnote above.
        cutoff_res   <- NULL
        sim_fail_msg <- NULL

        # The simulation is the expensive part, and jamovi reruns .run on
        # every option change -- including changes (pValues, correction,
        # sortByInfit) that do not affect the simulation. The plot state
        # already carries the cutoff object and jamovi clears it exactly
        # when a simulation-relevant option changes (its clearWith list),
        # so a surviving state is a valid cache. The signature check makes
        # reuse self-validating rather than relying on clearWith alone.
        sim_sig <- list(
          iterations = self$options$iterations,
          seed       = as.integer(self$options$seed),
          hdci_width = self$options$hdciWidth / 100
        )
        cached <- self$results$infitPlot$state
        if (isTRUE(self$options$computeCutoff) &&
            !is.null(cached) && !is.null(cached$cutoff_res) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df)) {
          cutoff_res <- cached$cutoff_res
        } else if (isTRUE(self$options$computeCutoff)) {
          cutoff_res <- tryCatch(
            suppressWarnings(suppressMessages(
              easyRasch2::RMitemInfitCutoff(
                df,
                iterations = self$options$iterations,
                parallel   = FALSE,
                seed       = as.integer(self$options$seed),
                hdci_width = self$options$hdciWidth / 100
              )
            )),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )
          # Guard against degenerate cutoffs: with very few successful
          # iterations the HDCI collapses (lower = upper) and every item
          # is spuriously flagged.
          if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
            sim_fail_msg <- paste0(
              "Only ", cutoff_res$actual_iterations, " of ",
              self$options$iterations, " simulation iterations succeeded ",
              "-- too few to estimate reliable cutoff intervals. This ",
              "typically happens when items have very low or very high ",
              "endorsement rates relative to the sample size."
            )
            cutoff_res <- NULL
          }
        }

        # Bootstrap p-values need the full cutoff object (it carries the
        # simulated distributions), so they are only computed when the
        # simulation succeeded. The pValues option can hold a stale TRUE
        # while greyed out (jamovi disables but does not reset nested
        # options), hence the explicit computeCutoff gate.
        use_pvalues <- isTRUE(self$options$computeCutoff) &&
          isTRUE(self$options$pValues) &&
          !is.null(cutoff_res)

        # Observed conditional infit (+ expected range and flags when the
        # cutoff simulation succeeded, + bootstrap p-values when requested).
        # The package's below-1000-iterations warning is suppressed with the
        # rest; the module states the same caveat in the note below.
        results <- suppressWarnings(suppressMessages(
          easyRasch2::RMitemInfit(df, cutoff = cutoff_res,
                                  p_value    = use_pvalues,
                                  correction = self$options$correction,
                                  output     = "dataframe")
        ))

        # 4. Sort if requested
        if (isTRUE(self$options$sortByInfit)) {
          results <- results[order(results$Infit_MSQ, decreasing = TRUE), ]
          rownames(results) <- NULL
        }

        # 5. Populate table
        table <- self$results$infitTable
        for (i in seq_len(nrow(results))) {
          vals <- list(
            item = results$Item[i],
            infitMSQ = results$Infit_MSQ[i],
            relLocation = results$Relative_location[i]
          )
          if (!is.null(cutoff_res)) {
            vals$infitLow <- results$Infit_low[i]
            vals$infitHigh <- results$Infit_high[i]
            vals$misfit <- results$Flagged[i]
          }
          if (use_pvalues) {
            vals$pValue <- results$p_infit[i]
            vals$pAdjusted <- results$padj_infit[i]
          }
          table$setRow(rowNo = i, values = vals)
        }

        # 6. Always-visible footnotes: complete-case basis + column docs.
        # When responses are missing, the infit statistics come from
        # iarm's complete-case refit while Rel. location comes from the
        # CML fit on all available responses (partial rows retained).
        excluded_clause <- if (n_excluded > 0L) {
          paste0(
            " ", n_excluded, " of ", n_total,
            " row(s) had a missing response on at least one selected item ",
            "and were excluded from the infit computation; item locations ",
            "use all available responses via conditional maximum ",
            "likelihood estimation."
          )
        } else {
          ""
        }
        table$setNote(
          "ncomplete",
          paste0(
            "Conditional infit MSQ is computed from complete responses only ",
            "(n = ", n_complete, " row(s) with no missing values across the ",
            "selected items).", excluded_clause
          )
        )
        table$setNote(
          "loc",
          paste0(
            "Rel. location = mean item (threshold) location relative to ",
            "the mean person location (weighted likelihood estimates, WLE), ",
            "in logits."
          )
        )
        if (use_pvalues) {
          # With bootstrap p-values, flagging follows the adjusted p-value
          # (upstream behavior); the expected-range columns stay as the
          # effect-size reference.
          table$setNote(
            "misfit",
            paste0(
              "Flagged: adjusted p-value < 0.05; observed infit below 1 = ",
              "overfit, above 1 = underfit. Note the direction is inverted ",
              "relative to the item-restscore analyses."
            )
          )
          table$setNote(
            "pvalues",
            paste0(
              "p-value: probability of an infit at least as extreme as ",
              "observed if the model fits, computed from the ",
              cutoff_res$actual_iterations, " simulated datasets ",
              "(Monte-Carlo). Adj. p-value: corrected for multiple ",
              "comparisons across the ", nrow(results), " items using ",
              correction_label(self$options$correction), "."
            )
          )
        } else if (!is.null(cutoff_res)) {
          table$setNote(
            "misfit",
            paste0(
              "Flagged: infit below the expected range = overfit (item is ",
              "more predictable than the model expects); above = underfit ",
              "(noisier than expected). Note the direction is inverted ",
              "relative to the item-restscore analyses."
            )
          )
        }

        # 7. Set cutoff note (HTML element below the table). Always set
        # non-empty content on the success path so a stale "requires at
        # least 3 items" guard message is overwritten when enough items
        # are selected -- jamovi does not clear set HTML content via
        # setContent("") / clearWith, so the note must be actively
        # rewritten with real text.
        if (!is.null(cutoff_res)) {
          method_label <- paste0(cutoff_res$hdci_width * 100, "% HDCI")

          note_html <- paste0(
            "<p>Cutoff values based on ",
            cutoff_res$actual_iterations, " simulation iterations (",
            method_label, ") drawn from the same n = ", n_complete,
            " complete cases.",
            iteration_note(self$options$iterations, 400L, corrected = TRUE),
            iteration_attrition_note(
              cutoff_res$actual_iterations,
              self$options$iterations
            ),
            if (use_pvalues) {
              pvalue_iteration_caveat(cutoff_res$actual_iterations, floor = 400L)
            } else {
              # Flagging falls back to the expected range, whose width sets
              # a familywise error rate the user has not chosen explicitly.
              interval_flagging_note(cutoff_res$hdci_width, nrow(results))
            },
            "</p>"
          )
          self$results$cutoffNote$setContent(note_html)
        } else if (!is.null(sim_fail_msg)) {
          self$results$cutoffNote$setContent(paste0(
            "<p><b>Simulation-based cutoffs unavailable:</b> ",
            sim_fail_msg, "</p>"
          ))
        } else {
          self$results$cutoffNote$setContent(paste0(
            "<p>Conditional infit MSQ for n = ", n_complete,
            " complete cases. Enable <i>Simulation-based cutoffs</i> to ",
            "flag each item against a per-item expected range simulated ",
            "under the fitted model.</p>"
          ))
        }

        # 8. Save state for plot. The figure is drawn by
        # easyRasch2::RMitemInfitPlot() inside the render function; the
        # observed CML refit it performs there is cheap. Storing the data +
        # cutoff object (rather than a ggplot) keeps the saved analysis
        # small.
        if (!is.null(cutoff_res)) {
          self$results$infitPlot$setState(list(
            df         = df,
            cutoff_res = cutoff_res,
            n_complete = n_complete,
            sig        = sim_sig
          ))
        }

      }, error = function(e) {
        stop(paste("Error in conditional infit analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Simulated-infit dot plot — easyRasch2::RMitemInfitPlot() (ggdist dot
    # cloud + black per-item median + orange diamonds for the observed
    # conditional infit), restyled to the module's plot conventions
    # (base size 15).
    # ---------------------------------------------------------------------
    .infitPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMitemInfitPlot(
          image$state$cutoff_res,
          data = image$state$df
        )
      ))
      if (is.null(p)) return(FALSE)

      # Re-apply the module's theme on top of the package theme: base size
      # 15 for jamovi's 500x600 canvas, then restore the theme pieces that
      # theme_minimal() resets.
      p <- p +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
