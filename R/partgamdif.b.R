#' @export
partgamdifClass <- R6::R6Class(
  "partgamdifClass",
  inherit = partgamdifBase,
  private = list(
    .init = function() {
      # The bootstrap adjusted-p column names its correction method, as in
      # the infit and item-restscore analyses. The asymptotic column keeps
      # its fixed BH title.
      self$results$pgdifTable$getColumn("pAdjBoot")$setTitle(
        padjusted_title(self$options$correction)
      )
    },

    .run = function() {
      # 1. Return early / explain if requirements not met. Partial gamma
      # conditions on the total score, and with 2 items the total score
      # fixes one item's response given the other, so partgam_DIF returns
      # mirror-duplicate rows (gamma_1 = -gamma_2, identical SE and p),
      # i.e. a single test shown twice -- so at least 3 items are required.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (is.null(self$options$difVar))
        return()
      if (length(self$options$vars) < 3) {
        self$results$cutoffNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the total score fixes the response to one item given the ",
          "other, so the two rows of the table become mirror images of a ",
          "single test. Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Extract and validate data
      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      # Extract DIF variable (kept in original class -- typically a factor)
      dif_raw <- data[[self$options$difVar]]

      # Check DIF variable levels
      dif_levels <- if (is.factor(dif_raw)) levels(dif_raw) else
        unique(stats::na.omit(dif_raw))
      if (length(dif_levels) < 2)
        stop("The DIF variable must have at least 2 distinct levels.")

      # Handle complete cases including DIF variable
      complete_mask <- complete.cases(df) & !is.na(dif_raw)
      n_complete <- sum(complete_mask)
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      df <- df[complete_mask, , drop = FALSE]
      dif_vec <- dif_raw[complete_mask]

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # Sparse-category warning per DIF group -- the relevant split for
      # this analysis. Shown even when the tileplot is off, pointing
      # users to it for visual inspection.
      sparse_msg <- sparse_note_grouped(df, dif_vec)
      if (!is.null(sparse_msg)) {
        self$results$pgdifTable$setNote("sparse", paste0(
          sparse_msg, " Enable 'Show response distribution by DIF group' ",
          "to inspect the counts."
        ))
      }

      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$pgdifTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$pgdifTable$setNote("duplicate", dup_msg)

      # rgl workaround
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      tryCatch({
        # All computation is delegated to the easyRasch2 package
        # (iarm::partgam_DIF() for the observed statistics; a conditional
        # parametric bootstrap that keeps each respondent's group and total
        # score for the expected ranges and p-values). Results are
        # numerically identical to RMdifGamma() / RMdifGammaCutoff() with
        # the same seed and iterations. Package warnings are suppressed --
        # the module surfaces its own footnotes.

        # Optionally compute cutoffs first (they feed RMdifGamma). If the
        # simulation cannot deliver reliable cutoffs, degrade gracefully:
        # show the observed gammas without the expected range and say why
        # in a notice (sim_failure_notice()).
        #
        # Cached in the plot state (see run_cutoff_simulation()), so changing
        # the tileplot options, sortByGamma or showSE does not rerun the
        # simulation. The DIF variable is part of the key.
        sim_sig <- list(
          iterations = self$options$iterations,
          seed       = as.integer(self$options$seed),
          hdci_width = self$options$hdciWidth / 100
        )
        sim <- run_cutoff_simulation(
          compute   = self$options$computeCutoff,
          cached    = self$results$pgdifPlot$state,
          sig       = sim_sig,
          keys      = list(df = df, dif = dif_vec),
          simulate  = function() easyRasch2::RMdifGammaCutoff(
            df,
            dif_var    = dif_vec,
            iterations = self$options$iterations,
            parallel   = FALSE,
            seed       = as.integer(self$options$seed),
            hdci_width = self$options$hdciWidth / 100
          ),
          requested = self$options$iterations
        )
        cutoff_res   <- sim$cutoff_res
        sim_fail_msg <- sim$fail_msg

        # `p_value` is passed explicitly rather than left to the package
        # default, which is NULL from easyRasch2 1.3.1.9001 and resolves on
        # the presence of a cutoff object. The pValues option can hold a
        # stale TRUE while greyed out (jamovi disables but does not reset
        # nested options), hence the explicit computeCutoff gate, as in the
        # other simulation analyses. The two branches return different
        # columns: with p-values the asymptotic `padj_bh` is dropped in
        # favour of `p_gamma` and `padj_gamma`.
        use_pvalues <- isTRUE(self$options$computeCutoff) &&
          isTRUE(self$options$pValues) &&
          !is.null(cutoff_res)

        # Visibility follows what was computed, not the options: if the
        # simulation failed, the simulation columns would otherwise show
        # empty (see set_columns_visible()).
        sim_ok <- !is.null(cutoff_res)
        tbl <- self$results$pgdifTable
        set_columns_visible(tbl, c("gammaLow", "gammaHigh"), sim_ok)
        set_columns_visible(tbl, c("pValue", "pAdjBoot"), use_pvalues)
        set_columns_visible(tbl, "padjBH", !use_pvalues)
        set_elements_visible(list(self$results$pgdifPlot), sim_ok)

        pgam_df <- suppressWarnings(suppressMessages(
          easyRasch2::RMdifGamma(
            df,
            dif_var    = dif_vec,
            cutoff     = cutoff_res,
            p_value    = use_pvalues,
            correction = self$options$correction,
            output     = "dataframe"
          )
        ))

        n_complete_used <- n_complete

        # Sort if requested
        if (isTRUE(self$options$sortByGamma)) {
          pgam_df <- pgam_df[order(abs(pgam_df$gamma), decreasing = TRUE), ]
          rownames(pgam_df) <- NULL
        }

        # Populate table
        table <- self$results$pgdifTable
        for (i in seq_len(nrow(pgam_df))) {
          vals <- list(
            item   = pgam_df$Item[i],
            gamma  = pgam_df$gamma[i],
            se     = pgam_df$se[i],
            lower  = pgam_df$lower[i],
            upper  = pgam_df$upper[i]
          )
          if (use_pvalues) {
            vals$pValue   <- pgam_df$p_gamma[i]
            vals$pAdjBoot <- pgam_df$padj_gamma[i]
          } else {
            vals$padjBH <- pgam_df$padj_bh[i]
          }
          # Without a cutoff object the package returns no flag, so the
          # asymptotic rule (BH-adjusted p < .05) is applied here, as in
          # the item-restscore analysis.
          flag <- if (!is.null(cutoff_res)) {
            isTRUE(pgam_df$flagged[i])
          } else {
            isTRUE(pgam_df$padj_bh[i] < 0.05)
          }
          vals$flagged <- if (flag) "beyond cutoff" else ""
          if (!is.null(cutoff_res)) {
            vals$gammaLow  <- pgam_df$gamma_low[i]
            vals$gammaHigh <- pgam_df$gamma_high[i]
          }
          table$setRow(rowNo = i, values = vals)
        }

        # Table footnotes
        if (!use_pvalues) {
          table$setNote("sig", paste0(
            "Asymptotic p-values adjusted across items with the ",
            "Benjamini-Hochberg (BH) false-discovery-rate method.",
            if (is.null(cutoff_res)) {
              " Flagged 'beyond cutoff' when the adjusted p-value < .05."
            }
          ))
        }
        if (!is.null(cutoff_res)) {
          range_txt <- paste0(
            "Expected range = ", cutoff_res$hdci_width * 100, "% HDCI of ",
            "partial gamma values simulated under no DIF (each simulated ",
            "dataset keeps every respondent's group and total score and ",
            "redraws the responses from the fitted model given that score)."
          )
          if (use_pvalues) {
            table$setNote("flag", paste0(
              range_txt, " Flagged 'beyond cutoff' when the adjusted ",
              "p-value < .05; the expected range is shown as description."
            ))
            table$setNote("pvalues", paste0(
              "p-value: two-sided probability of a partial gamma at least ",
              "as far from its simulated mean as observed if there is no ",
              "DIF, computed from the ", cutoff_res$actual_iterations,
              " simulated datasets (Monte-Carlo). Adj. p-value: corrected ",
              "for multiple comparisons across the ", nrow(pgam_df),
              " items using ", correction_label(self$options$correction), "."
            ))
          } else {
            table$setNote("flag", paste0(
              range_txt, " Flagged 'beyond cutoff' when the observed ",
              "gamma falls outside the expected range."
            ))
          }
        }

        # HTML note: n reporting (always), CI level when shown, cutoff
        # basis or simulation-failure explanation.
        n_total    <- length(complete_mask)
        n_excluded <- n_total - n_complete_used
        excluded_clause <- if (n_excluded > 0L) {
          paste0(" (", n_excluded, " of ", n_total, " row(s) excluded ",
                 "due to missing item responses or a missing DIF value)")
        } else ""
        se_clause <- if (isTRUE(self$options$showSE)) {
          paste0(" Confidence intervals are 95% Wald intervals ",
                 "(gamma ± 1.96 × SE).")
        } else ""
        cutoff_clause <- if (!is.null(cutoff_res)) {
          paste0(" Expected ranges", if (use_pvalues) " and p-values",
                 " based on ", cutoff_res$actual_iterations,
                 " simulation iterations (", cutoff_res$hdci_width * 100,
                 "% HDCI). Results are identical to ",
                 "easyRasch2::RMdifGamma() and RMdifGammaCutoff() with ",
                 "the same seed.",
                 # With p-values the caveat carries the iteration advice,
                 # so the general recommendation would repeat it.
                 if (!use_pvalues) {
                   iteration_note(self$options$iterations, 400L)
                 },
                 iteration_attrition_note(cutoff_res$actual_iterations,
                                          self$options$iterations),
                 if (use_pvalues) {
                   pvalue_iteration_caveat(cutoff_res$actual_iterations,
                                           count_stated = TRUE)
                 } else {
                   # Flagging falls back to the expected range, whose
                   # width sets a familywise error rate the user has not
                   # chosen.
                   interval_flagging_note(cutoff_res$hdci_width,
                                          nrow(pgam_df))
                 })
        } else ""
        if (!is.null(sim_fail_msg)) {
          sim_failure_notice(self, sim_fail_msg,
                             "the asymptotic p-values instead")
        }
        self$results$cutoffNote$setContent(paste0(
          "<p>Partial gamma DIF analysis based on n = ", n_complete_used,
          " complete cases", excluded_clause, ".", se_clause,
          cutoff_clause, "</p>"
        ))

        if (!is.null(cutoff_res)) {
          # Save state for plot (incl. the 95% Wald CI of the observed
          # gamma, drawn as a segment in the same colour as the diamond).
          # The full cutoff object plus sig/df/dif also serve as the
          # simulation cache consulted above.
          observed_gamma_vec <- pgam_df$gamma
          observed_lower_vec <- pgam_df$lower
          observed_upper_vec <- pgam_df$upper
          names(observed_gamma_vec) <- pgam_df$Item
          names(observed_lower_vec) <- pgam_df$Item
          names(observed_upper_vec) <- pgam_df$Item

          self$results$pgdifPlot$setState(list(
            cutoff_res        = cutoff_res,
            sig               = sim_sig,
            df                = df,
            dif               = dif_vec,
            observed_gamma    = observed_gamma_vec,
            observed_lower    = observed_lower_vec,
            observed_upper    = observed_upper_vec,
            item_names_data   = names(df)
          ))
        }

        # Tileplot: per-item × category × DIF-group response counts,
        # drawn by easyRasch2::RMplotTile() in the render function.
        if (isTRUE(self$options$showTileplot)) {
          self$results$tileplot$setState(list(df = df, dif = dif_vec))
        }

      }, error = function(e) {
        stop(paste("Error in partial gamma DIF analysis:", e$message))
      })
    },

    # ------------------------------------------------------------------
    # Simulated-gamma dot plot -- module-drawn: easyRasch2's plot shows
    # the observed diamonds but not the observed gamma's 95% Wald CI
    # segment, which is a module convention worth keeping. The simulated
    # distributions come from RMdifGammaCutoff()$results and the observed
    # values from RMdifGamma(), so the numbers match the package exactly.
    # ------------------------------------------------------------------
    .pgDIFplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
      if (!requireNamespace("ggdist",  quietly = TRUE)) return(FALSE)

      state             <- image$state
      results_df        <- state$cutoff_res$results
      item_names        <- state$cutoff_res$item_names
      actual_iterations <- state$cutoff_res$actual_iterations
      sample_n          <- state$cutoff_res$sample_n
      observed_gamma    <- state$observed_gamma
      item_names_data   <- state$item_names_data

      item_levels <- rev(item_names)

      # Compute per-item summary intervals
      lo_hi <- do.call(rbind, lapply(item_names, function(item) {
        sub <- results_df[results_df$Item == item, ]
        data.frame(
          Item            = item,
          min_gamma       = stats::quantile(sub$gamma, 0.005, na.rm = TRUE),
          max_gamma       = stats::quantile(sub$gamma, 0.995, na.rm = TRUE),
          p66lo_gamma     = stats::quantile(sub$gamma, 0.167, na.rm = TRUE),
          p66hi_gamma     = stats::quantile(sub$gamma, 0.833, na.rm = TRUE),
          median_gamma    = stats::median(sub$gamma, na.rm = TRUE),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }))
      rownames(lo_hi) <- NULL

      observed_df <- data.frame(
        Item           = item_names_data,
        observed_gamma = as.numeric(observed_gamma[item_names_data]),
        observed_lower = as.numeric(state$observed_lower[item_names_data]),
        observed_upper = as.numeric(state$observed_upper[item_names_data]),
        stringsAsFactors = FALSE
      )
      observed_df$Item_f <- factor(observed_df$Item, levels = item_levels)

      gamma_sim <- data.frame(
        Item  = results_df$Item,
        Value = results_df$gamma,
        stringsAsFactors = FALSE
      )
      gamma_sim <- merge(gamma_sim, observed_df[, c("Item", "observed_gamma")], by = "Item", sort = FALSE)
      gamma_sim$Item <- factor(gamma_sim$Item, levels = item_levels)
      lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

      caption_text <- er2_caption(paste0(
        "Results from ", actual_iterations,
        " simulated datasets with ", sample_n, " respondents.\n",
        "Orange diamonds indicate observed partial gamma, with the ",
        "orange line showing its 95% Wald CI.\n",
        "Black dots indicate median gamma from simulations."
      ))

      p <- ggplot2::ggplot(gamma_sim, ggplot2::aes(x = .data$Value, y = .data$Item)) +
        ggdist::stat_dots(
          ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
          quantiles = actual_iterations,
          layout = "weave",
          slab_color = NA,
          .width = c(0.666, 0.99)
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$min_gamma, xend = .data$max_gamma,
                       y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 0.7
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$p66lo_gamma, xend = .data$p66hi_gamma,
                       y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 1.2
        ) +
        ggplot2::geom_point(
          data = lo_hi,
          ggplot2::aes(x = .data$median_gamma, y = .data$Item_f),
          size = 3.6
        ) +
        ggplot2::geom_segment(
          data = observed_df,
          ggplot2::aes(x = .data$observed_lower, xend = .data$observed_upper,
                       y = .data$Item_f, yend = .data$Item_f),
          color = "sienna2", linewidth = 0.8,
          position = ggplot2::position_nudge(y = -0.1)
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = .data$observed_gamma),
          color = "sienna2", shape = 18,
          position = ggplot2::position_nudge(y = -0.1),
          size = 7
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::labs(x = "Partial gamma", y = "Item", caption = caption_text) +
        ggplot2::scale_color_manual(
          values = scales::brewer_pal()(3)[-1],
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
    },

    # ------------------------------------------------------------------
    # Faceted tileplot of item × category response counts by DIF group,
    # drawn by easyRasch2::RMplotTile() (the function the module's
    # previous hand-rolled version was ported from).
    # ------------------------------------------------------------------
    .tileplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMplotTile(
          image$state$df,
          group   = image$state$dif,
          cutoff  = self$options$tileCutoff,
          percent = isTRUE(self$options$tilePercent)
        )
      ))
      p <- er2_bump_text(p)

      print(p)
      TRUE
    }
  )
)
