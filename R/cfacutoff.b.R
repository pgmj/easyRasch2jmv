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
          "of freedom) or not identified (&le; 2 items), so the model fit is ",
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

      # Drop incomplete rows -- lavaan + simulation both need complete cases
      n_total     <- nrow(df)
      df_complete <- stats::na.omit(df)
      n_complete  <- nrow(df_complete)
      n_excluded  <- n_total - n_complete

      if (n_complete == 0L)
        stop("No complete cases found. CFA requires at least one row with responses to all selected items.")
      if (n_complete < 30L)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )

      for (col in names(df_complete)) {
        if (length(unique(df_complete[[col]])) < 2L)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # Sparse-category warning on the complete cases actually analysed
      sparse_msg <- sparse_note(df_complete)
      if (!is.null(sparse_msg))
        self$results$cfaTable$setNote("sparse", sparse_msg)

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

      # 4. Fit observed and run simulation
      tryCatch({
        sim_data_list <- private$.buildSimDataList(df_complete, estimator)

        observed <- run_observed_cfa_fit(df_complete, estimator)
        if (!is.numeric(observed)) {
          stop(paste0(
            "Observed CFA fit failed (", observed, "). ",
            "lavaan WLSMV / ULSMV typically fails when items have ",
            "very rare or empty response categories. Inspect the response ",
            "distribution per item (e.g., the descriptives module) before ",
            "running this analysis."
          ))
        }
        names(observed) <- c("cfi", "rmsea", "srmr")

        # Seed is always applied (default 42) so results are reproducible
        # by default, consistent with the other simulation-based analyses.
        set.seed(as.integer(seed_val))
        sim_seeds <- sample.int(.Machine$integer.max, iterations)

        # Sequential only -- jamovi runs single-threaded
        results_raw <- run_cfa_sim_sequential(
          iterations    = iterations,
          sim_seeds     = sim_seeds,
          sim_data_list = sim_data_list,
          verbose       = FALSE
        )

        ok         <- vapply(results_raw, is.numeric, logical(1L))
        successful <- results_raw[ok]

        # Guard against degenerate cutoffs: with very few successful
        # iterations the percentile cutoffs collapse onto a handful of
        # values. Require at least 20 successes and a 50% success rate;
        # otherwise show the observed fit indices without cutoffs and
        # explain the dominant failure reason (the observed lavaan fit
        # is independent of the simulation, so it remains valid).
        n_ok <- length(successful)
        if (n_ok < 20L || n_ok < iterations / 2) {
          fail_msgs <- unlist(results_raw[!ok])
          top_reason <- if (length(fail_msgs) > 0L) {
            names(sort(table(fail_msgs), decreasing = TRUE))[1L]
          } else NULL

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
            "<p><b>Simulation-based cutoffs unavailable:</b> Only ", n_ok,
            " of ", iterations, " simulation iterations succeeded -- too ",
            "few to estimate reliable cutoffs.",
            if (!is.null(top_reason))
              paste0(" Most common failure: ", top_reason, ".") else "",
            " Often this indicates very sparse items; inspect the ",
            "per-item response distribution. The observed fit indices ",
            "are shown without cutoffs.</p>"
          ))
          return()
        }

        actual_iterations <- length(successful)
        sim_mat <- do.call(rbind, successful)
        colnames(sim_mat) <- c("cfi", "rmsea", "srmr")
        simulated_df <- data.frame(
          iteration = seq_len(actual_iterations),
          cfi       = as.numeric(sim_mat[, "cfi"]),
          rmsea     = as.numeric(sim_mat[, "rmsea"]),
          srmr      = as.numeric(sim_mat[, "srmr"]),
          stringsAsFactors = FALSE
        )

        cutoffs <- private$.computeCfaCutoffs(simulated_df, percentile)
        flagged <- private$.computeCfaFlagged(observed, cutoffs)

        is_polytomous <- sim_data_list$type == "polytomous"

        # 5. Populate the table (3 fixed rows, set in r.yaml)
        table <- self$results$cfaTable

        idx_names <- c("CFI", "RMSEA", "SRMR")
        for (i in seq_len(3L)) {
          k <- c("cfi", "rmsea", "srmr")[i]
          table$setRow(rowNo = i, values = list(
            index    = idx_names[i],
            observed = observed[[k]],
            cutoff   = cutoffs[[k]],
            flagged  = if (isTRUE(flagged[[k]])) "TRUE" else ""
          ))
        }
        table$setNote("flag", paste0(
          "Flagged = TRUE when the observed value lies beyond the cutoff ",
          "in the unfavourable direction (CFI below the cutoff; RMSEA ",
          "and SRMR above)."
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
          "the upper-tail cutoff. An item flagged 'TRUE' lies in the ",
          "worst ", round(100 - percentile, 1),
          "% of the simulated distribution in the unfavourable direction.",
          success_clause,
          iteration_note(iterations, 250L), "</p>"
        )
        self$results$cfaNote$setContent(note_html)

        # 7. Save state for the plot
        self$results$cfaPlot$setState(list(
          simulated         = simulated_df,
          observed          = observed,
          cutoffs           = cutoffs,
          flagged           = flagged,
          percentile        = percentile,
          actual_iterations = actual_iterations,
          n_complete        = n_complete,
          is_polytomous     = is_polytomous,
          estimator         = estimator
        ))
      }, error = function(e) {
        stop(paste("Error in CFA-cutoff analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # .buildSimDataList -- inlined from easyRasch2::RMdimCFACutoff
    # ---------------------------------------------------------------------
    .buildSimDataList = function(df, estimator) {
      data_mat       <- as.matrix(df)
      sample_n       <- nrow(data_mat)
      is_polytomous  <- max(data_mat, na.rm = TRUE) > 1L
      item_names_vec <- colnames(data_mat)
      if (is.null(item_names_vec))
        item_names_vec <- paste0("V", seq_len(ncol(data_mat)))

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
        list(type       = "polytomous",
             thetas     = thetas,
             deltaslist = deltaslist,
             n_items    = ncol(data_mat),
             sample_n   = sample_n,
             item_names = item_names_vec,
             estimator  = estimator)
      } else {
        rm_fit      <- eRm::RM(data_mat)
        pp          <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores  <- rowSums(data_mat, na.rm = TRUE)
        thetas      <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        list(type        = "dichotomous",
             thetas      = thetas,
             item_params = item_params,
             n_items     = ncol(data_mat),
             sample_n    = sample_n,
             item_names  = item_names_vec,
             estimator   = estimator)
      }
    },

    # ---------------------------------------------------------------------
    # Cutoff + flag computation -- one-sided, directional
    # ---------------------------------------------------------------------
    .computeCfaCutoffs = function(simulated_df, percentile) {
      pct <- percentile / 100
      # Filter to finite values -- lavaan can return Inf for RMSEA when
      # chi-square is exactly 0 (degenerate near-perfect fit), which
      # would otherwise propagate into the cutoff.
      cfi_v   <- simulated_df$cfi[is.finite(simulated_df$cfi)]
      rmsea_v <- simulated_df$rmsea[is.finite(simulated_df$rmsea)]
      srmr_v  <- simulated_df$srmr[is.finite(simulated_df$srmr)]
      c(
        cfi   = as.numeric(stats::quantile(cfi_v,   1 - pct, na.rm = TRUE)),
        rmsea = as.numeric(stats::quantile(rmsea_v, pct,     na.rm = TRUE)),
        srmr  = as.numeric(stats::quantile(srmr_v,  pct,     na.rm = TRUE))
      )
    },

    .computeCfaFlagged = function(observed, cutoffs) {
      c(
        cfi   = !is.na(observed[["cfi"]])   && observed[["cfi"]]   < cutoffs[["cfi"]],
        rmsea = !is.na(observed[["rmsea"]]) && observed[["rmsea"]] > cutoffs[["rmsea"]],
        srmr  = !is.na(observed[["srmr"]])  && observed[["srmr"]]  > cutoffs[["srmr"]]
      )
    },

    # ---------------------------------------------------------------------
    # .cfaPlot -- faceted histogram with diamond marker and cutoff line
    # ---------------------------------------------------------------------
    .cfaPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state <- image$state

      # Long-format data for faceting
      sim_long <- data.frame(
        Index = factor(rep(c("CFI", "RMSEA", "SRMR"),
                           each = nrow(state$simulated)),
                       levels = c("CFI", "RMSEA", "SRMR")),
        Value = c(state$simulated$cfi,
                  state$simulated$rmsea,
                  state$simulated$srmr),
        stringsAsFactors = FALSE
      )
      sim_long <- sim_long[is.finite(sim_long$Value), , drop = FALSE]

      obs_df <- data.frame(
        Index    = factor(c("CFI", "RMSEA", "SRMR"),
                          levels = c("CFI", "RMSEA", "SRMR")),
        Observed = c(state$observed[["cfi"]],
                     state$observed[["rmsea"]],
                     state$observed[["srmr"]]),
        Cutoff   = c(state$cutoffs[["cfi"]],
                     state$cutoffs[["rmsea"]],
                     state$cutoffs[["srmr"]]),
        Flagged  = c(state$flagged[["cfi"]],
                     state$flagged[["rmsea"]],
                     state$flagged[["srmr"]]),
        stringsAsFactors = FALSE
      )
      obs_df$Color <- ifelse(obs_df$Flagged, "red", "black")

      cfi_pct_lbl <- 100 - state$percentile
      caption <- er2_caption(paste0(
        "Histograms: ", state$actual_iterations,
        " parametric-bootstrap datasets simulated under ",
        if (state$is_polytomous) "PCM" else "RM",
        " unidimensionality at n = ", state$n_complete, ",\n",
        "refitted with lavaan::cfa(ordered = TRUE, estimator = \"",
        state$estimator, "\").\n",
        "Diamond: observed value (red = flagged at the ", state$percentile,
        "th percentile in the unfavourable direction).\n",
        "Dashed line: cutoff (CFI: ", round(cfi_pct_lbl, 1),
        "th pct; RMSEA / SRMR: ", state$percentile, "th pct)."
      ))

      p <- ggplot2::ggplot(sim_long,
                           ggplot2::aes(x = .data$Value)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey80",
                                colour = "white") +
        ggplot2::geom_vline(
          data = obs_df,
          ggplot2::aes(xintercept = .data$Cutoff),
          linetype = "dashed", colour = "grey40"
        ) +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = .data$Observed, colour = .data$Color),
          y = 0, size = 6, shape = 18
        ) +
        ggplot2::scale_colour_identity() +
        ggplot2::facet_wrap(~ Index, scales = "free", nrow = 1) +
        ggplot2::labs(
          x       = "Fit index value",
          y       = "Count of simulated datasets",
          caption = caption
        ) +
        ggplot2::theme_bw(base_size = 13) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
