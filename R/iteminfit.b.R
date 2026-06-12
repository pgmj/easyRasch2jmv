#' @export
iteminfitClass <- R6::R6Class(
  "iteminfitClass",
  inherit = iteminfitBase,
  private = list(
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

      n_complete <- sum(complete.cases(df))
      n_total    <- nrow(df)
      n_excluded <- n_total - n_complete

      if (n_complete == 0)
        stop("No complete cases found in the data. Conditional infit requires at least one row with responses to all selected items.")
      if (n_complete < 30)
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)


      # 3. Fit Rasch model and compute infit
      tryCatch({
        data_mat <- as.matrix(df)
        n_items <- ncol(df)

        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        if (max(data_mat, na.rm = TRUE) == 1L) {
          erm_out <- eRm::RM(df)
          item_avg_locations <- stats::coef(erm_out, "beta") * -1
          pp <- eRm::person.parameter(erm_out)
          person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
        } else {
          erm_out <- eRm::PCM(df)
          thresh_obj <- eRm::thresholds(erm_out)
          thresh_table <- thresh_obj$threshtable[[1]]
          if ("Location" %in% colnames(thresh_table)) {
            item_avg_locations <- thresh_table[, "Location"]
          } else {
            item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
          }
          pp <- eRm::person.parameter(erm_out)
          person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
        }

        relative_item_avg_locations <- item_avg_locations - person_avg_location

        # iarm refits the model on complete cases (na.omit) when any
        # responses are missing -- suppress its console message; the
        # behaviour is documented in the table footnote.
        cfit <- suppressMessages(iarm::out_infit(erm_out))

        # Build results data.frame. Raw values (no pre-rounding) so the
        # jamovi frontend applies the user's "Number format" preferences.
        results <- data.frame(
          Item = names(df),
          Infit_MSQ = cfit$Infit,
          Relative_location = relative_item_avg_locations,
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        # 4. Optionally compute cutoffs. If the simulation cannot deliver
        # reliable cutoffs (e.g. too few successful iterations), degrade
        # gracefully: show the observed infit table without the expected
        # range and explain why in the note below the table.
        cutoff_res <- NULL
        if (isTRUE(self$options$computeCutoff)) {
          cutoff_res <- tryCatch(
            private$.runCutoffSim(df),
            error = function(e) {
              self$results$cutoffNote$setContent(paste0(
                "<p><b>Simulation-based cutoffs unavailable:</b> ",
                e$message, "</p>"
              ))
              NULL
            }
          )

          if (!is.null(cutoff_res)) {
            cutoff_df <- cutoff_res$item_cutoffs
            data_items <- results$Item
            cutoff_sub <- cutoff_df[, c("Item", "infit_low", "infit_high")]
            results <- merge(results, cutoff_sub, by = "Item", sort = FALSE)
            results <- results[match(data_items, results$Item), ]
            rownames(results) <- NULL
            results$Infit_low <- results$infit_low
            results$Infit_high <- results$infit_high
            results$infit_low <- NULL
            results$infit_high <- NULL
            # Misfit direction is inverted relative to the restscore
            # analyses: infit BELOW the expected range = overfit (too
            # predictable), ABOVE = underfit (noisy).
            results$Misfit <- ifelse(
              results$Infit_MSQ < results$Infit_low, "overfit",
              ifelse(results$Infit_MSQ > results$Infit_high, "underfit", "")
            )
            results <- results[, c("Item", "Infit_MSQ", "Infit_low", "Infit_high", "Misfit", "Relative_location")]
          }
        }

        # 5. Sort if requested
        if (isTRUE(self$options$sortByInfit)) {
          results <- results[order(results$Infit_MSQ, decreasing = TRUE), ]
          rownames(results) <- NULL
        }

        # 6. Populate table
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
            vals$misfit <- results$Misfit[i]
          }
          table$setRow(rowNo = i, values = vals)
        }

        # 7. Always-visible footnotes: complete-case basis + column docs.
        # When responses are missing, the infit statistics come from
        # iarm's complete-case refit while Rel. location comes from the
        # eRm fit on all available responses (CML retains partial rows).
        excluded_clause <- if (n_excluded > 0L) {
          paste0(
            " ", n_excluded, " of ", n_total,
            " row(s) had a missing response on at least one selected item ",
            "and were excluded from the infit computation; item locations ",
            "use all available responses via eRm's conditional maximum ",
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
            "the mean person location, in logits."
          )
        )
        if (!is.null(cutoff_res)) {
          table$setNote(
            "misfit",
            paste0(
              "Misfit: infit below the expected range = overfit (item is ",
              "more predictable than the model expects); above = underfit ",
              "(noisier than expected). Note the direction is inverted ",
              "relative to the item-restscore analyses."
            )
          )
        }

        # 8. Set cutoff note (HTML element below the table)
        if (!is.null(cutoff_res)) {
          method_label <- paste0(cutoff_res$hdci_width * 100, "% HDCI")

          note_html <- paste0(
            "<p>Cutoff values based on ",
            cutoff_res$actual_iterations, " simulation iterations (",
            method_label, ") drawn from the same n = ", n_complete,
            " complete cases.",
            iteration_note(self$options$iterations, 200L, infit = TRUE),
            "</p>"
          )
          self$results$cutoffNote$setContent(note_html)
        }

        # 9. Save state for plot
        if (!is.null(cutoff_res)) {
          self$results$infitPlot$setState(list(
            results_df = cutoff_res$results,
            item_names = cutoff_res$item_names,
            actual_iterations = cutoff_res$actual_iterations,
            sample_n = cutoff_res$sample_n,
            n_complete = n_complete,
            observed_infit = cfit$Infit,
            item_names_data = names(df)
          ))
        }

      }, error = function(e) {
        stop(paste("Error in conditional infit analysis:", e$message))
      })
    },

    .runCutoffSim = function(df) {
      # Implements RMitemInfitCutoff() logic (sequential only)
      hdci_width <- self$options$hdciWidth / 100
      iterations <- self$options$iterations
      seed <- self$options$seed

      if (!requireNamespace("ggdist", quietly = TRUE)) {
        stop("Package 'ggdist' is required for HDCI cutoff method. Install with: install.packages(\"ggdist\")")
      }

      data_complete <- stats::na.omit(df)
      if (nrow(data_complete) == 0L)
        stop("No complete cases for simulation.")

      set.seed(seed)
      sim_seeds <- sample.int(.Machine$integer.max, iterations)

      data_mat <- as.matrix(data_complete)
      sample_n <- nrow(data_mat)
      is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
      item_names_vec <- colnames(data_mat)

      if (is_polytomous) {
        pcm_fit <- eRm::PCM(data_mat)
        pp <- eRm::person.parameter(pcm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        thresh_mat <- extract_item_thresholds(data_mat)
        deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
          as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
        })
        sim_data_list <- list(
          type = "polytomous", thetas = thetas, deltaslist = deltaslist,
          n_items = ncol(data_mat), sample_n = sample_n, item_names = item_names_vec
        )
      } else {
        rm_fit <- eRm::RM(data_mat)
        pp <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        sim_data_list <- list(
          type = "dichotomous", thetas = thetas, item_params = item_params,
          n_items = ncol(data_mat), sample_n = sample_n, item_names = item_names_vec
        )
      }

      results_raw <- run_infit_sim_sequential(iterations, sim_seeds, sim_data_list, verbose = FALSE)

      ok <- vapply(results_raw, is.list, logical(1L))
      successful <- results_raw[ok]

      # Guard against degenerate cutoffs: with very few successful
      # iterations the HDCI collapses (lower = upper) and every item is
      # spuriously flagged. Require at least 20 successes and a 50%
      # success rate; otherwise report the dominant failure reason.
      n_ok <- length(successful)
      if (n_ok < 20L || n_ok < iterations / 2) {
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

          infit_interval <- ggdist::hdci(sub$InfitMSQ, .width = hdci_width)
          outfit_interval <- ggdist::hdci(sub$OutfitMSQ, .width = hdci_width)
          data.frame(Item = item, infit_low = infit_interval[1L, 1L], infit_high = infit_interval[1L, 2L],
                     outfit_low = outfit_interval[1L, 1L], outfit_high = outfit_interval[1L, 2L],
                     stringsAsFactors = FALSE, row.names = NULL)
      }))
      rownames(item_cutoffs) <- NULL

      list(
        results = results_df, item_cutoffs = item_cutoffs,
        actual_iterations = actual_iterations, sample_n = sample_n,
        sample_summary = summary(thetas), item_names = item_names_vec,
        hdci_width = hdci_width
      )
    },

    .infitPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
      if (!requireNamespace("ggdist", quietly = TRUE)) return(FALSE)

      state <- image$state
      results_df <- state$results_df
      item_names <- state$item_names
      actual_iterations <- state$actual_iterations
      sample_n <- state$sample_n
      n_complete <- state$n_complete
      observed_infit <- state$observed_infit
      item_names_data <- state$item_names_data

      item_levels <- rev(item_names)

      lo_hi <- do.call(rbind, lapply(item_names, function(item) {
        sub <- results_df[results_df$Item == item, ]
        data.frame(
          Item = item,
          min_infit_msq = stats::quantile(sub$InfitMSQ, 0.001, na.rm = TRUE),
          max_infit_msq = stats::quantile(sub$InfitMSQ, 0.999, na.rm = TRUE),
          p66lo_infit_msq = stats::quantile(sub$InfitMSQ, 0.167, na.rm = TRUE),
          p66hi_infit_msq = stats::quantile(sub$InfitMSQ, 0.833, na.rm = TRUE),
          median_infit = stats::median(sub$InfitMSQ, na.rm = TRUE),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }))
      rownames(lo_hi) <- NULL

      observed_df <- data.frame(
        Item = item_names_data,
        observed_infit = observed_infit,
        stringsAsFactors = FALSE
      )

      infit_sim <- data.frame(
        Item = results_df$Item,
        Value = results_df$InfitMSQ,
        stringsAsFactors = FALSE
      )
      infit_sim <- merge(infit_sim, observed_df[, c("Item", "observed_infit")], by = "Item", sort = FALSE)
      infit_sim$Item <- factor(infit_sim$Item, levels = item_levels)
      lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

      caption_text <- er2_caption(paste0(
        "Observed and simulated infit are based on n = ",
        if (!is.null(n_complete)) n_complete else sample_n,
        " complete responses.\n",
        "Simulated distributions: ", actual_iterations,
        " parametric-bootstrap datasets using the same n.\n",
        "Orange diamonds: observed conditional infit. Black dots: simulation median."
      ))

      p <- ggplot2::ggplot(infit_sim, ggplot2::aes(x = .data$Value, y = .data$Item)) +
        ggdist::stat_dots(
          ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
          quantiles = actual_iterations,
          layout = "weave",
          slab_color = NA,
          .width = c(0.666, 0.999)
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$min_infit_msq, xend = .data$max_infit_msq, y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 0.7
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$p66lo_infit_msq, xend = .data$p66hi_infit_msq, y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 1.2
        ) +
        ggplot2::geom_point(
          data = lo_hi,
          ggplot2::aes(x = .data$median_infit, y = .data$Item_f),
          size = 3.6
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = .data$observed_infit),
          color = "sienna2", shape = 18,
          position = ggplot2::position_nudge(y = -0.1),
          size = 7
        ) +
        ggplot2::labs(x = "Conditional Infit MSQ", y = "Item", caption = caption_text) +
        ggplot2::scale_color_manual(
          values = scales::brewer_pal()(3)[-1],
          aesthetics = "slab_fill", guide = "none"
        ) +
        ggplot2::scale_x_continuous(minor_breaks = NULL) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
