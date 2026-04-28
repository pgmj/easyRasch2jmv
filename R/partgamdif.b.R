#' @export
partgamdifClass <- R6::R6Class(
  "partgamdifClass",
  inherit = partgamdifBase,
  private = list(
    .run = function() {
      # 1. Return early if no variables selected or difVar not set
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (is.null(self$options$difVar))
        return()
      if (length(self$options$vars) < 2)
        return()

      # 2. Extract and validate data
      data <- self$data
      vars <- self$options$vars
      df <- data[, vars, drop = FALSE]

      # Robust conversion: handles factors with text labels (SPSS),
      # haven_labelled vectors, and numerics.
      df <- to_numeric_responses_df(df)

      # Extract DIF variable (kept in original class — typically a factor)
      dif_col_name <- self$options$difVar
      dif_raw <- data[[dif_col_name]]

      # Check DIF variable levels
      dif_levels <- if (is.factor(dif_raw)) levels(dif_raw) else unique(stats::na.omit(dif_raw))
      if (length(dif_levels) < 2)
        stop("The DIF variable must have at least 2 distinct levels.")

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

      # Sentinel-value sanity check (e.g., 999, 8888 unmarked as missing)
      max_obs <- max(as.matrix(df), na.rm = TRUE)
      if (is.finite(max_obs) && max_obs > 20) {
        bad_cols <- names(df)[
          vapply(df, function(x) {
            mx <- suppressWarnings(max(x, na.rm = TRUE))
            is.finite(mx) && mx > 20
          }, logical(1L))
        ]
        stop(paste0(
          "Item(s) ", paste(bad_cols, collapse = ", "),
          " contain values > 20, which look like missing-value codes ",
          "(e.g., 999, 8888) rather than ordinal responses. ",
          "Mark these codes as missing in the data editor, or recode your data."
        ))
      }

      validate_response_data(df)

      # Handle complete cases including DIF variable
      complete_mask <- complete.cases(df) & !is.na(dif_raw)
      n_complete <- sum(complete_mask)
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)

      df <- df[complete_mask, , drop = FALSE]
      dif_vec <- dif_raw[complete_mask]

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # rgl workaround
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      tryCatch({
        # Compute partial gamma DIF
        nullfile_path <- nullfile()
        sink(nullfile_path)
        on.exit(sink(), add = TRUE)
        pgam_raw <- iarm::partgam_DIF(as.data.frame(df), dif_vec)
        sink()

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

        # Optionally compute cutoffs
        cutoff_res <- NULL
        if (isTRUE(self$options$computeCutoff)) {
          cutoff_res <- private$.runCutoffSim(df, dif_vec)

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
            padjBH = pgam_df$padj_bh[i]
          )
          if (!is.null(cutoff_res)) {
            vals$gammaLow  <- pgam_df$gamma_low[i]
            vals$gammaHigh <- pgam_df$gamma_high[i]
            vals$flagged   <- ifelse(pgam_df$Flagged[i], "TRUE", "")
          }
          table$setRow(rowNo = i, values = vals)
        }

        # Set cutoff note
        if (!is.null(cutoff_res)) {
          method_label <- paste0(cutoff_res$hdci_width * 100, "% HDCI")
          note_html <- paste0(
            "<p>Partial gamma DIF analysis (n = ", n_complete_used,
            " complete cases). Cutoff values based on ",
            cutoff_res$actual_iterations, " simulation iterations (",
            method_label, ").</p>"
          )
          self$results$cutoffNote$setContent(note_html)

          # Save state for plot
          observed_gamma_vec <- pgam_df$gamma
          names(observed_gamma_vec) <- pgam_df$Item

          self$results$pgdifPlot$setState(list(
            results_df        = cutoff_res$results,
            item_names        = cutoff_res$item_names,
            actual_iterations = cutoff_res$actual_iterations,
            sample_n          = cutoff_res$sample_n,
            observed_gamma    = observed_gamma_vec,
            item_names_data   = names(df)
          ))
        }

      }, error = function(e) {
        stop(paste("Error in partial gamma DIF analysis:", e$message))
      })
    },

    .runCutoffSim = function(df, dif_vec) {
      # Implements RMpgDIFcutoff() logic (sequential only)
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
        stringsAsFactors = FALSE
      )

      gamma_sim <- data.frame(
        Item  = results_df$Item,
        Value = results_df$gamma,
        stringsAsFactors = FALSE
      )
      gamma_sim <- merge(gamma_sim, observed_df[, c("Item", "observed_gamma")], by = "Item", sort = FALSE)
      gamma_sim$Item <- factor(gamma_sim$Item, levels = item_levels)
      lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

      caption_text <- paste0(
        "Note: Results from ", actual_iterations,
        " simulated datasets with ", sample_n, " respondents.\n",
        "Orange dots indicate observed partial gamma.\n",
        "Black dots indicate median gamma from simulations."
      )

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
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"),
                       plot.caption = ggplot2::element_text(size = 11))

      print(p)
      TRUE
    }
  )
)
