#' @export
bootrestscoreClass <- R6::R6Class(
  "bootrestscoreClass",
  inherit = bootrestscoreBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2)
        return()

      # 2. Extract data and convert to numeric
      data <- self$data
      vars <- self$options$vars
      df <- data[, vars, drop = FALSE]

      # Robust conversion: handles factors with text labels (SPSS),
      # haven_labelled vectors, and numerics.
      df <- to_numeric_responses_df(df)

      # Identical-item check (same pattern as itemrestscore.b.R)
      n_vars <- ncol(df)
      identical_pairs <- list()
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          cor_val <- cor(df[[i]], df[[j]], use = "complete.obs")
          if (!is.na(cor_val) && cor_val == 1) {
            identical_pairs <- append(identical_pairs,
                                      list(c(names(df)[i], names(df)[j])))
          }
        }
      }
      if (length(identical_pairs) > 0) {
        pair_strings <- sapply(identical_pairs,
                               function(p) paste0("'", p[1], "' and '", p[2], "'"))
        pair_msg <- paste(pair_strings, collapse = "; ")
        if (ncol(df) == 2) {
          stop(paste("The two selected items are identical:", pair_msg,
                     "- please select different items."))
        } else {
          jmvcore::reject(
            "Warning: Some items appear to be identical ({pairs}). This may affect results.",
            pairs = pair_msg
          )
        }
      }

      # All-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste(
          "The following variables contain no valid numeric data:",
          paste(bad_vars, collapse = ", ")
        ))
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

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Bootstrap results may be unreliable.",
          n = n_complete
        )

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # 3. Read options
      iterations <- self$options$iterations
      samplesize <- self$options$samplesize
      cutoff     <- self$options$cutoff
      seed       <- self$options$seed

      # Clamp samplesize to nrow(df) instead of failing (Jamovi-friendly)
      samplesize_used <- min(samplesize, nrow(df))
      samplesize_clamped <- samplesize_used < samplesize

      # 4. Run analysis
      tryCatch({
        data_mat      <- as.matrix(df)
        item_names    <- colnames(data_mat)
        n_items       <- ncol(data_mat)
        is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        # Full-sample model: locations + conditional infit
        if (is_polytomous) {
          erm_full <- eRm::PCM(df)
          thresh_table <- eRm::thresholds(erm_full)$threshtable[[1L]]
          if ("Location" %in% colnames(thresh_table)) {
            item_avg_locations <- thresh_table[, "Location"]
          } else {
            item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
          }
        } else {
          erm_full <- eRm::RM(df)
          item_avg_locations <- stats::coef(erm_full, "beta") * -1
        }
        pp <- eRm::person.parameter(erm_full)
        person_avg_location <- mean(pp$theta.table[["Person Parameter"]],
                                    na.rm = TRUE)
        relative_item_avg_locations <- item_avg_locations - person_avg_location

        cfit <- iarm::out_infit(erm_full)
        n_complete_used <- nrow(stats::na.omit(df))

        # Per-iteration seeds for reproducibility
        set.seed(seed)
        boot_seeds <- sample.int(.Machine$integer.max, iterations)

        boot_data_list <- list(
          data          = df,
          samplesize    = samplesize_used,
          is_polytomous = is_polytomous,
          item_names    = item_names
        )

        results_raw <- run_boot_restscore_sequential(
          iterations, boot_seeds, boot_data_list, verbose = FALSE
        )

        ok <- vapply(results_raw, is.data.frame, logical(1L))
        successful <- results_raw[ok]
        if (length(successful) == 0L)
          stop("All bootstrap iterations failed. Check your data.")
        actual_iterations <- length(successful)

        # Stack and tag with iteration index
        iter_dfs <- lapply(seq_along(successful), function(i) {
          d <- successful[[i]]
          d$iteration <- i
          d[, c("iteration", "Item", "item_restscore", "diff", "diff_abs")]
        })
        fit_all <- do.call(rbind, iter_dfs)
        rownames(fit_all) <- NULL

        # Per-item classification counts
        classes <- c("overfit", "underfit", "no misfit")
        counts <- as.data.frame(
          table(Item = factor(fit_all$Item, levels = item_names),
                item_restscore = factor(fit_all$item_restscore, levels = classes)),
          responseName = "n",
          stringsAsFactors = FALSE
        )
        counts$Item           <- as.character(counts$Item)
        counts$item_restscore <- as.character(counts$item_restscore)
        per_item_total <- tapply(counts$n, counts$Item, sum)
        counts$percent <- round(counts$n * 100 / per_item_total[counts$Item], 1)

        # Wide per-item summary
        wide <- data.frame(
          Item = item_names,
          stringsAsFactors = FALSE
        )
        get_pct <- function(item, cls) {
          v <- counts$percent[counts$Item == item & counts$item_restscore == cls]
          if (length(v) == 0L) 0 else v
        }
        wide$pctOverfit   <- vapply(item_names, get_pct, numeric(1L), cls = "overfit")
        wide$pctUnderfit  <- vapply(item_names, get_pct, numeric(1L), cls = "underfit")
        wide$pctNoMisfit  <- vapply(item_names, get_pct, numeric(1L), cls = "no misfit")
        wide$infitMSQ     <- round(cfit$Infit, 2)
        wide$relLocation  <- round(relative_item_avg_locations, 2)
        wide$flagged      <- ifelse(
          wide$pctOverfit > cutoff | wide$pctUnderfit > cutoff,
          "yes", ""
        )

        # 5. Populate the results table
        table <- self$results$bootstrapTable
        for (i in seq_len(nrow(wide))) {
          table$setRow(rowNo = i, values = list(
            item        = wide$Item[i],
            pctOverfit  = wide$pctOverfit[i],
            pctUnderfit = wide$pctUnderfit[i],
            pctNoMisfit = wide$pctNoMisfit[i],
            infitMSQ    = wide$infitMSQ[i],
            relLocation = wide$relLocation[i],
            flagged     = wide$flagged[i]
          ))
        }

        # 6. Caption note
        clamp_msg <- if (samplesize_clamped) {
          paste0(" Requested sample size (", samplesize,
                 ") exceeded the number of available rows; clamped to ",
                 samplesize_used, ".")
        } else {
          ""
        }
        note_html <- paste0(
          "<p>Results based on ", actual_iterations,
          " successful bootstrap iterations with n = ", samplesize_used,
          " and ", n_items, " items.",
          clamp_msg,
          " Conditional MSQ infit based on complete responders only (n = ",
          n_complete_used,
          "). Items where % overfit or % underfit exceed the display cutoff (",
          cutoff, "%) are marked as flagged.</p>"
        )
        self$results$bootstrapNote$setContent(note_html)

        # 7. Save state for plot (raw per-iteration data)
        if (isTRUE(self$options$showPlot)) {
          self$results$bootstrapPlot$setState(list(
            fit_all           = fit_all,
            item_names        = item_names,
            actual_iterations = actual_iterations,
            samplesize_used   = samplesize_used
          ))
        }

      }, error = function(e) {
        stop(paste("Error in bootstrap item-restscore analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Plot: per-item violin + jitter of (expected - observed) across
    # bootstrap iterations, coloured by per-iteration classification
    # ---------------------------------------------------------------------
    .bootstrapPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state <- image$state
      d <- state$fit_all
      d$Item <- factor(d$Item, levels = state$item_names)
      d$item_restscore <- factor(d$item_restscore,
                                 levels = c("overfit", "underfit", "no misfit"))

      caption_text <- paste0(
        "Note: ", state$actual_iterations,
        " bootstrap iterations with n = ", state$samplesize_used, " per draw.\n",
        "Each point is one iteration; colour = per-iteration classification (BH-adj. p < .05)."
      )

      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(x = .data$Item, y = .data$diff)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey50") +
        ggplot2::geom_violin(fill = "grey90", colour = NA) +
        ggplot2::geom_jitter(
          ggplot2::aes(colour = .data$item_restscore),
          width = 0.15, alpha = 0.5, size = 0.9
        ) +
        ggplot2::scale_colour_manual(
          values = c("overfit"   = "#377eb8",
                     "underfit"  = "#e41a1c",
                     "no misfit" = "grey60"),
          name = NULL,
          drop = FALSE
        ) +
        ggplot2::labs(
          x = NULL,
          y = "Expected − observed restscore correlation",
          caption = caption_text
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
          plot.caption = ggplot2::element_text(size = 10)
        )

      print(p)
      TRUE
    }
  )
)

# ---------------------------------------------------------------------------
# Internal helpers (mirror easyRasch2::RMbootRestscore internals)
# ---------------------------------------------------------------------------

#' Run a single item-restscore bootstrap iteration
#'
#' @keywords internal
#' @noRd
run_single_boot_restscore <- function(seed, data_list) {
  set.seed(seed)
  idx <- sample.int(nrow(data_list$data), data_list$samplesize, replace = TRUE)
  d   <- data_list$data[idx, , drop = FALSE]

  tryCatch({
    if (data_list$is_polytomous) {
      model_fit <- psychotools::pcmodel(d, hessian = FALSE)
    } else {
      model_fit <- eRm::RM(d, se = FALSE)
    }

    i1 <- as.data.frame(iarm::item_restscore(model_fit))
    res_mat <- i1[[1L]]
    n_items <- length(data_list$item_names)

    observed <- as.numeric(res_mat[seq_len(n_items), 1L])
    expected <- as.numeric(res_mat[seq_len(n_items), 2L])
    p_adj    <- as.numeric(res_mat[seq_len(n_items), 5L])
    diff_val <- expected - observed

    cls <- ifelse(p_adj < 0.05 & diff_val < 0, "overfit",
           ifelse(p_adj < 0.05 & diff_val > 0, "underfit", "no misfit"))

    data.frame(
      Item           = data_list$item_names,
      item_restscore = cls,
      diff           = diff_val,
      diff_abs       = abs(diff_val),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) as.character(conditionMessage(e)))
}

#' Run item-restscore bootstrap iterations sequentially
#'
#' @keywords internal
#' @noRd
run_boot_restscore_sequential <- function(iterations, boot_seeds, boot_data_list,
                                          verbose = FALSE) {
  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    results[[i]] <- run_single_boot_restscore(boot_seeds[i], boot_data_list)
  }
  results
}
