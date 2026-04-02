#' @export
iteminfitClass <- R6::R6Class(
  "iteminfitClass",
  inherit = iteminfitBase,
  private = list(
    .run = function() {
      # 1. Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2)
        stop("You need at least two variables to run an analysis.")

      # 2. Extract and validate data (same pattern as itemrestscore.b.R)
      data <- self$data
      vars <- self$options$vars
      df <- data[, vars, drop = FALSE]

      # Convert factors to numeric
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          df[[col]] <- as.numeric(as.character(df[[col]]))
        } else {
          df[[col]] <- as.numeric(df[[col]])
        }
      }

      # Check for identical items (correlation = 1)
      n_vars <- ncol(df)
      identical_pairs <- list()
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          cor_val <- cor(df[[i]], df[[j]], use = "complete.obs")
          if (!is.na(cor_val) && cor_val == 1) {
            identical_pairs <- append(identical_pairs, list(c(names(df)[i], names(df)[j])))
          }
        }
      }
      if (length(identical_pairs) > 0) {
        pair_strings <- sapply(identical_pairs, function(p) paste0("'", p[1], "' and '", p[2], "'"))
        pair_msg <- paste(pair_strings, collapse = "; ")
        if (ncol(df) == 2) {
          stop(paste("The two selected items are identical:", pair_msg, "- please select different items."))
        } else {
          jmvcore::reject("Warning: Some items appear to be identical ({pairs}). This may affect results.", pairs = pair_msg)
        }
      }

      # Check for all-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste("The following variables contain no valid numeric data:", paste(bad_vars, collapse = ", ")))
      }

      validate_response_data(df)

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

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
        cfit <- iarm::out_infit(erm_out)
        n_complete_used <- nrow(stats::na.omit(df))

        # Build results data.frame
        results <- data.frame(
          Item = names(df),
          Infit_MSQ = round(cfit$Infit, 3),
          Relative_location = round(relative_item_avg_locations, 2),
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        # 4. Optionally compute cutoffs
        cutoff_res <- NULL
        if (isTRUE(self$options$computeCutoff)) {
          cutoff_res <- private$.runCutoffSim(df)

          if (!is.null(cutoff_res)) {
            cutoff_df <- cutoff_res$item_cutoffs
            data_items <- results$Item
            cutoff_sub <- cutoff_df[, c("Item", "infit_low", "infit_high")]
            results <- merge(results, cutoff_sub, by = "Item", sort = FALSE)
            results <- results[match(data_items, results$Item), ]
            rownames(results) <- NULL
            results$Infit_low <- round(results$infit_low, 3)
            results$Infit_high <- round(results$infit_high, 3)
            results$infit_low <- NULL
            results$infit_high <- NULL
            results$Flagged <- results$Infit_MSQ < results$Infit_low | results$Infit_MSQ > results$Infit_high
            results <- results[, c("Item", "Infit_MSQ", "Infit_low", "Infit_high", "Flagged", "Relative_location")]
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
            vals$flagged <- ifelse(results$Flagged[i], "TRUE", "")
          }
          table$addRow(rowKey = i, values = vals)
        }

        # 7. Set cutoff note
          method_label <- paste0(cutoff_res$hdci_width * 100, "% HDCI")

          note_html <- paste0(
            "<p>MSQ values based on conditional estimation (n = ", n_complete_used,
            " complete cases). Cutoff values based on ",
            cutoff_res$actual_iterations, " simulation iterations (",
            method_label, ").</p>"
          )
          self$results$cutoffNote$setContent(note_html)


        # 8. Save state for plot
        if (!is.null(cutoff_res)) {
          self$results$infitPlot$setState(list(
            results_df = cutoff_res$results,
            item_names = cutoff_res$item_names,
            actual_iterations = cutoff_res$actual_iterations,
            sample_n = cutoff_res$sample_n,
            observed_infit = cfit$Infit,
            item_names_data = names(df)
          ))
        }

      }, error = function(e) {
        stop(paste("Error in conditional infit analysis:", e$message))
      })
    },

    .runCutoffSim = function(df) {
      # Implements RMinfitcutoff() logic (sequential only)
      hdci_width <- self$options$hdciWidth
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
      if (length(successful) == 0L)
        stop("All simulation iterations failed.")

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
      observed_infit <- state$observed_infit
      item_names_data <- state$item_names_data

      item_levels <- rev(item_names)

      # Compute per-item summary intervals
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

      caption_text <- paste0(
        "Note: Results from ", actual_iterations,
        " simulated datasets with ", sample_n, " respondents.\n",
        "Orange dots indicate observed conditional item fit. ",
        "Black dots indicate median fit from simulations."
      )

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
          # Request 3 colors from Brewer palette, drop the first (lightest) one
          values = scales::brewer_pal()(3)[-1],
          aesthetics = "slab_fill", guide = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = seq(0.5, 1.5, 0.1), minor_breaks = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

      print(p)
      TRUE
    }
  )
)
