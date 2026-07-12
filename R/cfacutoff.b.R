#' @export
cfacutoffClass <- R6::R6Class(
  "cfacutoffClass",
  inherit = cfacutoffBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- the three rows exist via r.yaml (rows: 3); pre-fill the
    # design-fixed index labels so the table renders meaningfully before
    # .run() finishes.
    # ---------------------------------------------------------------------
    .init = function() {
      if (is.null(self$options$vars) || length(self$options$vars) < 4)
        return()
      table <- self$results$cfaTable
      idx_names <- c("CFI", "RMSEA", "SRMR")
      for (i in seq_len(3L)) {
        table$setRow(rowNo = i, values = list(
          index    = idx_names[i],
          observed = NA_real_,
          cutoff   = NA_real_,
          flagged  = ""
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early if requirements not met.
      # A one-factor CFA needs at least 4 indicators to be over-identified
      # (df > 0). With 3 items the model is just-identified (df = 0) and every
      # fit index is degenerate (CFI = 1, RMSEA = 0, SRMR = 0) for both the
      # observed data and every simulated dataset, so the cutoff distribution
      # carries no information; with <= 2 items the model is not identified and
      # lavaan cannot invert the information matrix. Verified with lavaan
      # 0.6-21 (ordered = TRUE, WLSMV).
      vars <- self$options$vars
      if (is.null(vars) || length(vars) == 0)
        return()
      if (length(vars) < 4) {
        self$results$cfaNote$setContent(paste0(
          "<p>This analysis requires at least <b>4 items</b>. A one-factor ",
          "CFA with fewer indicators is just-identified (3 items, 0 degrees ",
          "of freedom) or not identified (≤ 2 items), so the model fit is ",
          "perfect by construction and the simulated fit-index cutoffs carry ",
          "no information. Select at least 4 items.</p>"
        ))
        return()
      }

      # 2. Extract data + standard validations
      data <- self$data
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      # Complete-case counts for the notes (easyRasch2 drops incomplete
      # rows internally -- lavaan + the simulation need complete cases)
      n_total     <- nrow(df)
      df_complete <- stats::na.omit(df)
      n_complete  <- nrow(df_complete)
      n_excluded  <- n_total - n_complete

      if (n_complete == 0L)
        stop("No complete cases found. CFA requires at least one row with responses to all selected items.")

      for (col in names(df_complete)) {
        if (length(unique(df_complete[[col]])) < 2L)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # Sparse-category warning on the complete cases actually analysed
      sparse_msg <- sparse_note(df_complete)
      if (!is.null(sparse_msg))
        self$results$cfaTable$setNote("sparse", sparse_msg)

      dup_msg <- duplicate_items_note(df_complete)
      if (!is.null(dup_msg))
        self$results$cfaTable$setNote("duplicate", dup_msg)

      # 3. Read options
      estimator  <- toupper(self$options$estimator)
      percentile <- self$options$percentile
      iterations <- self$options$iterations
      seed_val   <- self$options$seed

      if (!is.numeric(percentile) || percentile < 50 || percentile > 99.9) {
        stop("Cutoff percentile must be between 50 and 99.9. Common choices: 95, 99, 99.5.")
      }

      # rgl workaround
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      # 4. Simulate cutoffs + compare observed fit via easyRasch2.
      # Results are numerically identical to RMdimCFACutoff() / RMdimCFA()
      # with the same seed, iterations, percentile, and estimator. Package
      # warnings are suppressed -- the module surfaces its own footnotes.
      tryCatch({
        sim_fail_msg <- NULL
        cutoff_res <- tryCatch(
          suppressWarnings(suppressMessages(
            easyRasch2::RMdimCFACutoff(
              df,
              iterations = iterations,
              percentile = percentile,
              output     = "list",
              parallel   = FALSE,
              seed       = as.integer(seed_val),
              estimator  = estimator
            )
          )),
          error = function(e) {
            sim_fail_msg <<- e$message
            NULL
          }
        )
        # Guard against degenerate cutoffs: with very few successful
        # iterations the percentile cutoffs collapse onto a handful of
        # values.
        if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
          sim_fail_msg <- paste0(
            "Only ", cutoff_res$actual_iterations, " of ", iterations,
            " simulation iterations succeeded -- too few to estimate ",
            "reliable cutoffs. Often this indicates very sparse items; ",
            "inspect the per-item response distribution."
          )
          cutoff_res <- NULL
        }

        # Graceful degradation: the observed lavaan fit is independent of
        # the simulation, so it remains valid and is shown without
        # cutoffs. easyRasch2::RMdimCFA() deliberately refuses to run
        # without the simulated reference distribution, so this fallback
        # uses the module's retained observed-fit helper.
        if (is.null(cutoff_res)) {
          observed <- run_observed_cfa_fit(df_complete, estimator)
          if (!is.numeric(observed)) {
            stop(paste0(
              "Observed CFA fit failed (", observed, "). ",
              "lavaan WLSMV / ULSMV typically fails when items have ",
              "very rare or empty response categories. Inspect the ",
              "response distribution per item before running this ",
              "analysis."
            ))
          }
          names(observed) <- c("cfi", "rmsea", "srmr")
          table <- self$results$cfaTable
          idx_names <- c("CFI", "RMSEA", "SRMR")
          for (i in seq_len(3L)) {
            k <- c("cfi", "rmsea", "srmr")[i]
            table$setRow(rowNo = i, values = list(
              index    = idx_names[i],
              observed = observed[[k]],
              cutoff   = NA_real_,
              flagged  = ""
            ))
          }
          self$results$cfaNote$setContent(paste0(
            "<p><b>Simulation-based cutoffs unavailable:</b> ",
            sim_fail_msg,
            " The observed fit indices are shown without cutoffs.</p>"
          ))
          return()
        }

        # Observed fit + loadings against the simulated reference
        res <- suppressWarnings(suppressMessages(
          easyRasch2::RMdimCFA(df, cutoff = cutoff_res,
                               output = "dataframe")
        ))
        fit_df  <- res$fit       # Index, Observed, Cutoff, Direction, Flagged
        load_df <- res$loadings  # Item, Observed, Expected_low/high, Flagged

        is_polytomous <- max(as.matrix(df), na.rm = TRUE) > 1L
        actual_iterations <- cutoff_res$actual_iterations

        # 5. Populate the fit-index table (3 fixed rows, set in r.yaml)
        table <- self$results$cfaTable
        for (i in seq_len(3L)) {
          table$setRow(rowNo = i, values = list(
            index    = fit_df$Index[i],
            observed = fit_df$Observed[i],
            cutoff   = fit_df$Cutoff[i],
            flagged  = fit_df$Flagged[i]
          ))
        }
        table$setNote("flag", paste0(
          "Flagged = TRUE when the observed value lies beyond the cutoff ",
          "in the unfavourable direction (CFI below the cutoff; RMSEA ",
          "and SRMR above)."
        ))

        # 5b. Loadings table: observed standardized loadings vs the
        # per-item expected range from the same simulation.
        lt <- self$results$loadingsTable
        for (i in seq_len(nrow(load_df))) {
          lt$addRow(rowKey = i, values = list(
            item     = load_df$Item[i],
            observed = load_df$Observed[i],
            low      = load_df$Expected_low[i],
            high     = load_df$Expected_high[i],
            flagged  = load_df$Flagged[i]
          ))
        }
        lt$setNote("flag", paste0(
          "Expected range = central ", percentile, "% interval of each ",
          "item's simulated standardized loadings (tails of ",
          round((100 - percentile) / 2, 2), "% each). Flagged: 'below' = ",
          "the item loads weaker on the common factor than ",
          "unidimensionality predicts; 'above' = stronger. Deviating ",
          "loadings point to the items driving multidimensionality."
        ))

        # 6. Caption note (HTML below the table)
        excluded_clause <- if (n_excluded > 0L) {
          paste0(" (", n_excluded, " of ", n_total,
                 " row(s) excluded due to missing responses)")
        } else ""

        success_clause <- if (actual_iterations < iterations) {
          paste0(" Note: ", iterations - actual_iterations,
                 " of ", iterations, " iterations failed and were dropped; ",
                 actual_iterations, " contributed to the cutoffs.")
        } else ""

        note_html <- paste0(
          "<p><b>Posterior-predictive CFA fit-index check.</b> Observed ",
          "one-factor CFA fit (lavaan ", estimator, ", ordered = TRUE) ",
          "compared to a parametric-bootstrap null distribution simulated ",
          "under ", if (is_polytomous) "PCM" else "RM",
          " unidimensionality at n = ", n_complete,
          " complete cases", excluded_clause,
          ", with ", actual_iterations,
          " successful iterations. Cutoffs are one-sided at the ",
          percentile,
          "th percentile of the simulated distribution: CFI is flagged ",
          "when below the (", round(100 - percentile, 1),
          "th) lower-tail cutoff; RMSEA / SRMR are flagged when above ",
          "the upper-tail cutoff. Results are identical to ",
          "easyRasch2::RMdimCFACutoff() and RMdimCFA() with the same ",
          "seed.",
          success_clause,
          iteration_note(iterations, 250L),
          low_iteration_caveat(actual_iterations), "</p>"
        )
        self$results$cfaNote$setContent(note_html)

        # 7. Save state for the plots. Both figures are drawn by
        # easyRasch2::RMdimCFAPlot() inside the render functions (the
        # observed lavaan refit it performs there is a single fit).
        plot_state <- list(df = df, cutoff_res = cutoff_res)
        self$results$cfaPlot$setState(plot_state)
        self$results$cfaLoadingsPlot$setState(plot_state)
      }, error = function(e) {
        stop(paste("Error in CFA-cutoff analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # .cfaPlot -- observed fit indices vs the simulated null distributions
    # (easyRasch2::RMdimCFAPlot()$fit)
    # ---------------------------------------------------------------------
    .cfaPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      plots <- suppressWarnings(suppressMessages(
        easyRasch2::RMdimCFAPlot(
          image$state$cutoff_res,
          data = image$state$df
        )
      ))
      if (is.null(plots$fit)) return(FALSE)
      print(er2_bump_text(plots$fit))
      TRUE
    },

    # ---------------------------------------------------------------------
    # .cfaLoadingsPlot -- observed standardized loadings vs their simulated
    # expected ranges (easyRasch2::RMdimCFAPlot()$loadings)
    # ---------------------------------------------------------------------
    .cfaLoadingsPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      plots <- suppressWarnings(suppressMessages(
        easyRasch2::RMdimCFAPlot(
          image$state$cutoff_res,
          data = image$state$df
        )
      ))
      if (is.null(plots$loadings)) return(FALSE)
      print(er2_bump_text(plots$loadings))
      TRUE
    }
  )
)
