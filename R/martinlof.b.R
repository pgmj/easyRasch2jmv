#' @export
martinlofClass <- R6::R6Class(
  "martinlofClass",
  inherit = martinlofBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- pre-create the fixed summary rows.
    # ---------------------------------------------------------------------
    .init = function() {
      st <- self$results$summaryTable
      rows <- list(
        t    = "Observed Martin-Löf statistic (T)",
        p    = "Monte Carlo p-value",
        iter = "Successful MC iterations",
        n    = "Complete cases used",
        k1   = "Items in subscale 1",
        k2   = "Items in subscale 2"
      )
      for (key in names(rows)) {
        st$addRow(rowKey = key,
                  values = list(statistic = rows[[key]], value = NA_real_))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {
      s1 <- self$options$subscale1
      s2 <- self$options$subscale2
      if ((is.null(s1) || length(s1) == 0) &&
          (is.null(s2) || length(s2) == 0))
        return()
      if (is.null(s1) || length(s1) < 2 ||
          is.null(s2) || length(s2) < 2) {
        self$results$mlNote$setContent(paste0(
          "<p>The Martin-Löf test requires <b>at least 2 items in each ",
          "subscale</b>. Assign the items of the two hypothesised ",
          "subscales to the two boxes. The partition must be specified ",
          "a priori (from theory or instrument design), not derived from ",
          "this data.</p>"
        ))
        return()
      }

      data <- self$data
      vars <- c(s1, s2)
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$summaryTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$summaryTable$setNote("duplicate", dup_msg)

      # Complete cases only (the conditional test statistic is defined on
      # the joint subscore distribution); upstream requires >= 30.
      n_total    <- nrow(df)
      n_complete <- sum(stats::complete.cases(df))
      if (n_complete < 30) {
        self$results$mlNote$setContent(paste0(
          "<p>The Martin-Löf test requires at least <b>30 complete ",
          "cases</b> (rows with no missing responses on the selected ",
          "items); the data contain ", n_complete, ".</p>"
        ))
        return()
      }

      partition  <- list(s1, s2)
      iterations <- self$options$iterations
      stopping   <- if (isTRUE(self$options$sequential)) "sequential"
                    else "none"

      tryCatch({
        # The Monte Carlo simulation is the expensive part; the hidden
        # simCache element carries the result, with a signature check so
        # reuse is self-validating (e.g. the heatmap's minExpected option
        # re-renders without rerunning the simulation).
        sim_sig <- list(
          subscale1  = s1,
          subscale2  = s2,
          iterations = iterations,
          stopping   = stopping,
          seed       = as.integer(self$options$seed)
        )
        cached <- self$results$simCache$state
        if (!is.null(cached) && !is.null(cached$res) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df)) {
          res <- cached$res
        } else {
          # All computation is delegated to the easyRasch2 package, so
          # results are identical to RMdimMartinLof() with the same seed
          # and iterations.
          res <- suppressWarnings(suppressMessages(
            easyRasch2::RMdimMartinLof(
              df,
              partition  = partition,
              iterations = iterations,
              stopping   = stopping,
              parallel   = FALSE,
              seed       = as.integer(self$options$seed)
            )
          ))
          self$results$simCache$setState(list(
            res = res, sig = sim_sig, df = df
          ))
        }

        # --- Summary table ------------------------------------------------
        st <- self$results$summaryTable
        st$setRow(rowKey = "t",    values = list(value = res$T_obs))
        st$setRow(rowKey = "p",    values = list(value = res$p_value))
        st$setRow(rowKey = "iter", values = list(value = res$actual_iterations))
        st$setRow(rowKey = "n",    values = list(value = res$sample_n))
        st$setRow(rowKey = "k1",   values = list(value = length(s1)))
        st$setRow(rowKey = "k2",   values = list(value = length(s2)))
        # Attainable p-values are k / (iterations + 1); a p-value at the
        # floor means no simulated statistic reached the observed one and
        # only bounds the true p-value from above.
        p_floor <- if (!is.null(res$p_value_floor)) res$p_value_floor
                   else 1 / (res$actual_iterations + 1)
        floor_clause <- if (isTRUE(all.equal(res$p_value, p_floor))) {
          paste0(
            " The reported p-value equals the smallest value attainable ",
            "with ", res$actual_iterations, " iterations (1/",
            res$actual_iterations + 1, "): no simulated statistic ",
            "reached the observed one, so read it as p < ",
            signif(p_floor, 2), " -- the true p-value may be much ",
            "smaller. Increase the iterations for finer resolution."
          )
        } else ""
        st$setNote("p", paste0(
          "Monte Carlo p-value with (exceedances + 1) / (iterations + 1) ",
          "correction: the proportion of datasets simulated under the ",
          "unidimensional null whose statistic is at least as large as ",
          "the observed T. A small p-value indicates that the joint ",
          "distribution of the two subscale scores departs from what a ",
          "single dimension predicts.", floor_clause,
          iteration_note(iterations, 250L),
          # Under sequential stopping, ending early is the intended
          # behaviour (the p-value stays valid), so the few-successes
          # caveat applies only to full runs.
          if (stopping == "none")
            iteration_attrition_note(res$actual_iterations, iterations)
          else ""
        ))

        # --- Subscale WLE correlation (effect size) -------------------------
        wc <- res$wle_correlation
        self$results$corrTable$setRow(rowNo = 1, values = list(
          r      = wc$r[1],
          ciLow  = wc$ci_lower[1],
          ciHigh = wc$ci_upper[1],
          n      = as.integer(wc$n[1])
        ))
        self$results$corrTable$setNote("es", paste0(
          "Pearson correlation between the two subscales' person ",
          "estimates (Warm's WLE from a separate CML fit per subscale); ",
          "N = persons with finite estimates on both subscales (persons ",
          "at a subscale's minimum or maximum score are excluded ",
          "pairwise). This is the effect size to read alongside the ",
          "p-value: a rejected test with a correlation near 1 indicates ",
          "a statistically detectable but trivial departure from ",
          "unidimensionality; a correlation clearly below 1 indicates ",
          "substantive multidimensionality."
        ))

        # --- Note -----------------------------------------------------------
        drop_clause <- if (res$sample_n < n_total) {
          paste0(" (", n_total - res$sample_n, " row(s) with missing ",
                 "responses excluded; the test is defined on complete ",
                 "response patterns)")
        } else ""
        seq_clause <- if (stopping == "sequential") {
          paste0(" Sequential stopping (Besag & Clifford, 1991) was used: ",
                 "the simulation stops once 50 simulated statistics ",
                 "exceed the observed one, which shortens computation ",
                 "when the null is compatible with the data while ",
                 "keeping the p-value valid; ", res$actual_iterations,
                 " of at most ", iterations, " iterations were run.")
        } else ""
        self$results$mlNote$setContent(paste0(
          "<p>Martin-Löf test of unidimensionality against the a priori ",
          "two-subscale partition, for n = ", res$sample_n,
          " complete cases", drop_clause, ", using a ",
          if (isTRUE(res$is_polytomous)) "partial credit" else "Rasch",
          " model (CML via psychotools). The p-value is obtained by ",
          "Monte Carlo simulation under the unidimensional null model ",
          "because the asymptotic chi-square approximation is biased ",
          "toward conservatism at realistic sample sizes.", seq_clause,
          " <b>The partition must be specified a priori</b> (from theory ",
          "or instrument design): testing a partition suggested by the ",
          "same data -- for example the split indicated by the first ",
          "residual PCA contrast -- inflates the Type-I error rate and ",
          "invalidates the p-value. Results are identical to ",
          "easyRasch2::RMdimMartinLof().</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in Martin-Löf test:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Monte Carlo null distribution: histogram of the simulated statistics
    # with the observed T marked. Pure presentation of the cached result.
    # ---------------------------------------------------------------------
    .nullPlot = function(image, ggtheme, theme, ...) {
      state <- self$results$simCache$state
      if (is.null(state) || is.null(state$res)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      res <- state$res
      caption_text <- er2_caption(paste0(
        "Martin-Löf statistics from ", res$actual_iterations,
        " datasets simulated under the unidimensional null ",
        "(n = ", res$sample_n, " per dataset). Dashed line: observed ",
        "T = ", round(res$T_obs, 2), "; Monte Carlo p = ",
        signif(res$p_value, 3), ".",
        iteration_note(self$options$iterations, 250L)
      ))

      p <- ggplot2::ggplot(
        data.frame(T_rep = res$T_rep),
        ggplot2::aes(x = .data$T_rep)
      ) +
        ggplot2::geom_histogram(bins = 30, fill = "#78b0a0",
                                colour = "white") +
        ggplot2::geom_vline(xintercept = res$T_obs, linetype = "dashed",
                            linewidth = 0.8) +
        ggplot2::labs(x = "Martin-Löf statistic (simulated under H0)",
                      y = "Simulated datasets", caption = caption_text) +
        ggplot2::theme_minimal(base_size = 15) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    },

    # ---------------------------------------------------------------------
    # Observed vs expected score cross-table residuals -- drawn by
    # easyRasch2::RMdimMartinLofResiduals() (no Monte Carlo involved; data
    # and partition read from the cache element).
    # ---------------------------------------------------------------------
    .residualPlot = function(image, ggtheme, theme, ...) {
      state <- self$results$simCache$state
      if (is.null(state) || is.null(state$df)) return(FALSE)

      min_exp <- self$options$minExpected
      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMdimMartinLofResiduals(
          state$df,
          partition    = list(state$sig$subscale1, state$sig$subscale2),
          output       = "ggplot",
          min_expected = if (min_exp > 0) min_exp else NULL
        )
      ))
      if (is.null(p)) return(FALSE)

      print(er2_bump_text(p))
      TRUE
    }
  )
)
