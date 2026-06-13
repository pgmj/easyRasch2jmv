#' @export
partgamdifClass <- R6::R6Class(
  "partgamdifClass",
  inherit = partgamdifBase,
  private = list(
    .run = function() {
      # 1. Return early / explain if requirements not met. With 2 items
      # the rest score reduces to the other item of the pair and
      # partgam_DIF returns mirror-duplicate rows (gamma_1 = -gamma_2,
      # identical SE and p), i.e. a single coarsely-conditioned test
      # shown twice -- so at least 3 items are required.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (is.null(self$options$difVar))
        return()
      if (length(self$options$vars) < 3) {
        self$results$cutoffNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the rest score (total score minus the item) reduces to ",
          "the other item, and the two rows of the table become mirror ",
          "duplicates of a single test. Select at least 3 items.</p>"
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

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$pgdifTable$setNote("duplicate", dup_msg)

      # rgl workaround
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      tryCatch({
        # Compute partial gamma DIF (suppress iarm console output; the
        # finally handler guarantees the sink is removed exactly once)
        sink(nullfile())
        pgam_raw <- tryCatch(
          iarm::partgam_DIF(as.data.frame(df), dif_vec),
          finally = sink()
        )

        pgam_df <- data.frame(
          Item    = as.character(pgam_raw$Item),
          gamma   = as.numeric(pgam_raw$gamma),
          se      = as.numeric(pgam_raw$se),
          pvalue  = as.numeric(pgam_raw$pvalue),
          # Column 6 is the BH-adjusted p-value; iarm does not always give it
          # a stable name, so positional access matches easyRasch2 source behaviour.
          padj_bh = as.numeric(pgam_raw[[6]]),
          sig     = as.character(pgam_raw$sig),
          lower   = as.numeric(pgam_raw$lower),
          upper   = as.numeric(pgam_raw$upper),
          stringsAsFactors = FALSE
        )

        n_complete_used <- n_complete

        # Optionally compute cutoffs. If the simulation cannot deliver
        # reliable cutoffs (e.g. too few successful iterations), degrade
        # gracefully: show the observed gammas without the expected range
        # and explain why in the note below the table.
        cutoff_res <- NULL
        sim_fail_msg <- NULL
        if (isTRUE(self$options$computeCutoff)) {
          cutoff_res <- tryCatch(
            private$.runCutoffSim(df, dif_vec),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )

          if (!is.null(cutoff_res)) {
            cutoff_df <- cutoff_res$item_cutoffs
            data_items <- pgam_df$Item
            cutoff_sub <- cutoff_df[, c("Item", "gamma_low", "gamma_high")]
            pgam_df <- merge(pgam_df, cutoff_sub, by = "Item", sort = FALSE)
            pgam_df <- pgam_df[match(data_items, pgam_df$Item), ]
            rownames(pgam_df) <- NULL
            pgam_df$Flagged <- pgam_df$gamma < pgam_df$gamma_low | pgam_df$gamma > pgam_df$gamma_high
          }
        }

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
            upper  = pgam_df$upper[i],
            padjBH = pgam_df$padj_bh[i],
            # iarm's sig column is padded with leading spaces; trim for display.
            sig    = trimws(pgam_df$sig[i])
          )
          if (!is.null(cutoff_res)) {
            vals$gammaLow  <- pgam_df$gamma_low[i]
            vals$gammaHigh <- pgam_df$gamma_high[i]
            vals$flagged   <- ifelse(pgam_df$Flagged[i], "TRUE", "")
          }
          table$setRow(rowNo = i, values = vals)
        }

        # Table footnotes
        table$setNote("sig", paste0(
          "P-values adjusted with the Benjamini-Hochberg (BH) ",
          "false-discovery-rate method. ",
          "*** p < .001, ** p < .01, * p < .05, . p < .10 (adjusted)."
        ))
        if (!is.null(cutoff_res)) {
          table$setNote("flag", paste0(
            "Expected range = ", cutoff_res$hdci_width * 100, "% HDCI of ",
            "partial gamma values simulated under no DIF (the DIF ",
            "variable is randomly reassigned, preserving group ",
            "proportions). Flagged = TRUE when the observed gamma falls ",
            "outside the expected range."
          ))
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
          paste0(" Cutoff values based on ", cutoff_res$actual_iterations,
                 " simulation iterations (", cutoff_res$hdci_width * 100,
                 "% HDCI).",
                 iteration_note(self$options$iterations, 250L),
                 low_iteration_caveat(cutoff_res$actual_iterations))
        } else if (!is.null(sim_fail_msg)) {
          paste0(" <b>Simulation-based cutoffs unavailable:</b> ",
                 sim_fail_msg)
        } else ""
        self$results$cutoffNote$setContent(paste0(
          "<p>Partial gamma DIF analysis based on n = ", n_complete_used,
          " complete cases", excluded_clause, ".", se_clause,
          cutoff_clause, "</p>"
        ))

        if (!is.null(cutoff_res)) {
          # Save state for plot (incl. the 95% Wald CI of the observed
          # gamma, drawn as a segment in the same colour as the diamond)
          observed_gamma_vec <- pgam_df$gamma
          observed_lower_vec <- pgam_df$lower
          observed_upper_vec <- pgam_df$upper
          names(observed_gamma_vec) <- pgam_df$Item
          names(observed_lower_vec) <- pgam_df$Item
          names(observed_upper_vec) <- pgam_df$Item

          self$results$pgdifPlot$setState(list(
            results_df        = cutoff_res$results,
            item_names        = cutoff_res$item_names,
            actual_iterations = cutoff_res$actual_iterations,
            sample_n          = cutoff_res$sample_n,
            observed_gamma    = observed_gamma_vec,
            observed_lower    = observed_lower_vec,
            observed_upper    = observed_upper_vec,
            item_names_data   = names(df)
          ))
        }

        # Tileplot: per-item × category × DIF-group response counts
        if (isTRUE(self$options$showTileplot)) {
          tile_state <- private$.computeTileCounts(df, dif_vec)
          tile_state$cutoff   <- self$options$tileCutoff
          tile_state$percent  <- isTRUE(self$options$tilePercent)
          self$results$tileplot$setState(tile_state)
        }

      }, error = function(e) {
        stop(paste("Error in partial gamma DIF analysis:", e$message))
      })
    },

    # ------------------------------------------------------------------
    # Build per-(item × category × group) counts for the tileplot.
    # `df` is numeric, complete-cases item data; `dif_vec` is the
    # corresponding DIF variable values (length = nrow(df)).
    # ------------------------------------------------------------------
    .computeTileCounts = function(df, dif_vec) {
      item_names <- names(df)

      all_vals       <- unlist(df, use.names = FALSE)
      all_vals       <- all_vals[!is.na(all_vals)]
      min_val        <- min(all_vals)
      max_val        <- max(all_vals)
      all_categories <- seq(min_val, max_val)

      dif_factor <- if (is.factor(dif_vec)) {
        droplevels(dif_vec)
      } else {
        as.factor(dif_vec)
      }

      parts <- lapply(levels(dif_factor), function(g) {
        mask <- !is.na(dif_factor) & dif_factor == g
        sub  <- df[mask, , drop = FALSE]
        rows <- lapply(item_names, function(it) {
          vals <- sub[[it]]
          vals <- vals[!is.na(vals)]
          tab  <- table(factor(vals, levels = all_categories))
          data.frame(
            item     = it,
            category = as.integer(names(tab)),
            n        = as.integer(tab),
            group    = g,
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, rows)
      })
      count_df <- do.call(rbind, parts)
      rownames(count_df) <- NULL
      count_df$group <- factor(count_df$group, levels = levels(dif_factor))

      totals <- stats::aggregate(
        count_df$n,
        by  = count_df[, c("item", "group"), drop = FALSE],
        FUN = sum
      )
      colnames(totals)[ncol(totals)] <- "total"
      count_df <- merge(count_df, totals, by = c("item", "group"), sort = FALSE)
      count_df$percentage <- round(count_df$n / count_df$total * 100, 1)

      # Item ordering: top of y-axis = first column of df
      count_df$item_label <- factor(count_df$item, levels = rev(item_names))

      group_sizes <- table(dif_factor)

      list(
        count_df       = count_df,
        all_categories = all_categories,
        item_names     = item_names,
        group_sizes    = group_sizes
      )
    },

    .runCutoffSim = function(df, dif_vec) {
      # Implements RMdifGammaCutoff() logic (sequential only)
      hdci_width  <- self$options$hdciWidth / 100
      iterations  <- self$options$iterations
      seed        <- self$options$seed

      if (!requireNamespace("ggdist", quietly = TRUE)) {
        stop("Package 'ggdist' is required for HDCI cutoff method. Install with: install.packages(\"ggdist\")")
      }

      data_complete <- stats::na.omit(df)
      if (nrow(data_complete) == 0L)
        stop("No complete cases for simulation.")

      set.seed(seed)
      sim_seeds <- sample.int(.Machine$integer.max, iterations)

      data_mat   <- as.matrix(data_complete)
      sample_n   <- nrow(data_mat)
      is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
      item_names_vec <- colnames(data_mat)

      # DIF group structure
      dif_levels_sim      <- sort(unique(dif_vec))
      dif_table           <- table(dif_vec)
      dif_proportions_sim <- as.numeric(dif_table[as.character(dif_levels_sim)]) / length(dif_vec)

      if (is_polytomous) {
        pcm_fit    <- eRm::PCM(data_mat)
        pp         <- eRm::person.parameter(pcm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas     <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        thresh_mat <- extract_item_thresholds(data_mat)
        deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
          as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
        })
        sim_data_list <- list(
          type = "polytomous", thetas = thetas, deltaslist = deltaslist,
          n_items = ncol(data_mat), sample_n = sample_n, item_names = item_names_vec,
          dif_levels = dif_levels_sim, dif_proportions = dif_proportions_sim
        )
      } else {
        rm_fit  <- eRm::RM(data_mat)
        pp      <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas     <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        sim_data_list <- list(
          type = "dichotomous", thetas = thetas, item_params = item_params,
          n_items = ncol(data_mat), sample_n = sample_n, item_names = item_names_vec,
          dif_levels = dif_levels_sim, dif_proportions = dif_proportions_sim
        )
      }

      results_raw <- run_partgam_sim_sequential(iterations, sim_seeds, sim_data_list, verbose = FALSE)

      ok         <- vapply(results_raw, is.data.frame, logical(1L))
      successful <- results_raw[ok]

      # Guard against degenerate cutoffs: with very few successful
      # iterations the HDCI collapses and every item is spuriously
      # flagged. Require at least 20 successes and a 50% success rate;
      # otherwise report the dominant failure reason.
      n_ok <- length(successful)
      if (n_ok < 20L) {
        fail_msgs <- unlist(results_raw[!ok])
        top_reason <- if (length(fail_msgs) > 0L) {
          names(sort(table(fail_msgs), decreasing = TRUE))[1L]
        } else NULL
        stop(paste0(
          "Only ", n_ok, " of ", iterations, " simulation iterations ",
          "succeeded -- too few to estimate reliable cutoff intervals.",
          if (!is.null(top_reason))
            paste0(" Most common failure: ", top_reason, ".") else "",
          " This typically happens when items have very low or very ",
          "high endorsement rates relative to the sample size."
        ), call. = FALSE)
      }

      actual_iterations <- length(successful)
      iter_dfs <- lapply(seq_along(successful), function(i) {
        d <- successful[[i]]
        d$iteration <- i
        d
      })
      results_df <- do.call(rbind, iter_dfs)
      rownames(results_df) <- NULL

      item_names <- unique(results_df$Item)
      item_cutoffs <- do.call(rbind, lapply(item_names, function(item) {
        sub <- results_df[results_df$Item == item, ]
        gamma_interval <- ggdist::hdci(sub$gamma, .width = hdci_width)
        data.frame(
          Item       = item,
          gamma_low  = gamma_interval[1L, 1L],
          gamma_high = gamma_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }))
      rownames(item_cutoffs) <- NULL

      list(
        results           = results_df,
        item_cutoffs      = item_cutoffs,
        actual_iterations = actual_iterations,
        sample_n          = sample_n,
        item_names        = item_names_vec,
        hdci_width        = hdci_width
      )
    },

    .pgDIFplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
      if (!requireNamespace("ggdist",  quietly = TRUE)) return(FALSE)

      state             <- image$state
      results_df        <- state$results_df
      item_names        <- state$item_names
      actual_iterations <- state$actual_iterations
      sample_n          <- state$sample_n
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
    # Faceted tileplot of item × category response counts, faceted by
    # the DIF grouping variable. Diagnostic to inspect subgroup
    # response distributions before / alongside the DIF analysis.
    # Mirrors easyRasch2::RMplotTile logic.
    # ------------------------------------------------------------------
    .tileplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state          <- image$state
      count_df       <- state$count_df
      all_categories <- state$all_categories
      cutoff         <- state$cutoff
      use_pct        <- isTRUE(state$percent)

      p <- ggplot2::ggplot(
        count_df,
        ggplot2::aes(x = .data$category,
                     y = .data$item_label,
                     fill = .data$n)
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(
          expression(italic(n)),
          limits = c(0, NA)
        ) +
        ggplot2::scale_x_continuous(
          "Response category",
          expand = c(0, 0),
          breaks = all_categories
        )

      # Cell labels (counts or percentages) with low-count highlighting
      label_aes <- if (use_pct) {
        ggplot2::aes(label = paste0(.data$percentage, "%"),
                     color = ifelse(.data$n < cutoff, "red", "orange"))
      } else {
        ggplot2::aes(label = .data$n,
                     color = ifelse(.data$n < cutoff, "red", "orange"))
      }
      p <- p +
        ggplot2::geom_text(label_aes) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()

      caption_text <- if (!is.null(state$group_sizes)) {
        gs <- state$group_sizes
        er2_caption(paste0(
          "Response category counts by DIF group (",
          paste0(names(gs), ": n = ", as.integer(gs), collapse = "; "),
          "). Counts below ", cutoff, " are highlighted in red."
        ))
      } else NULL

      p <- p +
        ggplot2::facet_wrap(~ group,
                            labeller = ggplot2::labeller(
                              group = er2_wrap_labels
                            )) +
        ggplot2::labs(y = "Items", caption = caption_text) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(
          axis.text.x   = ggplot2::element_text(size = 10),
          panel.grid    = ggplot2::element_blank(),
          panel.spacing = ggplot2::unit(0.7, "cm")
        ) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
