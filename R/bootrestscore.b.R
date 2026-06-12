#' @export
bootrestscoreClass <- R6::R6Class(
  "bootrestscoreClass",
  inherit = bootrestscoreBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met.
      # Same minimum as the asymptotic item-restscore analysis: with 2
      # items the restscore reduces to the other item of the pair and
      # observed = expected by construction.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$bootstrapNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the restscore (total score minus the item) reduces to the ",
          "other item, making observed and expected correlations identical ",
          "by construction. Select at least 3 items.</p>"
        ))
        return()
      }

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

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$bootstrapTable$setNote("sparse", sparse_msg)

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

        # Guard against misleading percentages: with few successful
        # iterations the classification shares are based on tiny
        # denominators. Require at least 20 successes and a 50% success
        # rate; otherwise stop with the dominant failure reason (there
        # is no observed-only fallback -- the bootstrap is the analysis).
        n_ok <- length(successful)
        if (n_ok < 20L || n_ok < iterations / 2) {
          fail_msgs <- unlist(results_raw[!ok])
          top_reason <- if (length(fail_msgs) > 0L) {
            names(sort(table(fail_msgs), decreasing = TRUE))[1L]
          } else NULL
          stop(paste0(
            "Only ", n_ok, " of ", iterations, " bootstrap iterations ",
            "succeeded -- too few for trustworthy classification ",
            "percentages.",
            if (!is.null(top_reason))
              paste0(" Most common failure: ", top_reason, ".") else "",
            " This typically happens when items have very low or very ",
            "high endorsement rates, so that resampled datasets often ",
            "contain items without response variation. Consider a larger ",
            "bootstrap sample size."
          ), call. = FALSE)
        }
        actual_iterations <- n_ok

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

        # Wide per-item summary. Pass raw numerics (no pre-rounding) so
        # the jamovi frontend applies the user's "Number format"
        # preferences.
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
        wide$relLocation  <- relative_item_avg_locations
        # Misfit labels mirror the asymptotic item-restscore analysis:
        # an item is labelled when its classification share exceeds the
        # display cutoff. If both shares exceed it (possible only with a
        # cutoff below 50%), the more frequent classification wins.
        over  <- wide$pctOverfit > cutoff
        under <- wide$pctUnderfit > cutoff
        wide$misfit <- ""
        wide$misfit[over  & (!under | wide$pctOverfit >= wide$pctUnderfit)] <- "overfit"
        wide$misfit[under & (!over  | wide$pctUnderfit > wide$pctOverfit)]  <- "underfit"

        # Sort if requested (table and plot share the resulting order)
        sort_by <- self$options$sortBy
        if (sort_by == "overfit") {
          wide <- wide[order(-wide$pctOverfit), , drop = FALSE]
        } else if (sort_by == "underfit") {
          wide <- wide[order(-wide$pctUnderfit), , drop = FALSE]
        }
        rownames(wide) <- NULL

        # 5. Populate the results table
        table <- self$results$bootstrapTable
        for (i in seq_len(nrow(wide))) {
          table$setRow(rowNo = i, values = list(
            item        = wide$Item[i],
            pctOverfit  = wide$pctOverfit[i],
            pctUnderfit = wide$pctUnderfit[i],
            relLocation = wide$relLocation[i],
            misfit      = wide$misfit[i]
          ))
        }

        # Footnotes explaining classification, flagging, and location
        table$setNote("cls", paste0(
          "Items are classified per iteration as overfit (observed > ",
          "expected restscore correlation) or underfit (observed < ",
          "expected) when the BH-adjusted p-value < .05; the % columns ",
          "give the percentage of bootstrap iterations with each ",
          "classification."
        ))
        table$setNote("flag", paste0(
          "Misfit = the item was classified as overfit (or underfit) in ",
          "more than ", cutoff, "% of the bootstrap iterations."
        ))
        table$setNote("loc", paste0(
          "Rel. location = item location relative to the mean person ",
          "location (full sample)."
        ))

        # 6. Caption note
        clamp_msg <- if (samplesize_clamped) {
          paste0(" Requested sample size (", samplesize,
                 ") exceeded the number of available rows; clamped to ",
                 samplesize_used, ".")
        } else {
          ""
        }
        missing_msg <- if (n_complete < nrow(df)) {
          paste0(
            " Note: ", nrow(df) - n_complete, " row(s) have missing ",
            "responses; such rows can be drawn into bootstrap samples but ",
            "are excluded when the model is refitted within each ",
            "iteration, so the effective per-iteration n is smaller than ",
            "the bootstrap sample size."
          )
        } else {
          ""
        }
        note_html <- paste0(
          "<p>Results based on ", actual_iterations,
          " successful bootstrap iterations with n = ", samplesize_used,
          " and ", n_items, " items.",
          clamp_msg, missing_msg, "</p>"
        )
        self$results$bootstrapNote$setContent(note_html)

        # 7. Save state for plot (raw per-iteration data; item order
        # follows the table sort)
        if (isTRUE(self$options$showPlot)) {
          self$results$bootstrapPlot$setState(list(
            fit_all           = fit_all,
            item_names        = wide$Item,
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
      # coord_flip() below puts items on the y-axis like the conditional
      # infit plot; reversed levels place the first item at the top.
      d$Item <- factor(d$Item, levels = rev(state$item_names))
      d$item_restscore <- factor(d$item_restscore,
                                 levels = c("overfit", "underfit", "no misfit"))

      caption_text <- er2_caption(paste0(
        state$actual_iterations,
        " bootstrap iterations with n = ", state$samplesize_used, " per draw.\n",
        "Each point is one iteration, classified by BH-adjusted p < .05 ",
        "and sign: blue = overfit, red = underfit, grey = no misfit.\n",
        "Positive values indicate over-discrimination (overfit), negative ",
        "values under-discrimination (underfit)."
      ))

      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(x = .data$Item, y = .data$diff)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey50") +
        ggplot2::geom_violin(fill = "grey90", colour = NA) +
        ggplot2::geom_jitter(
          ggplot2::aes(colour = .data$item_restscore),
          width = 0.15, alpha = 0.5, size = 1.6
        ) +
        ggplot2::scale_colour_manual(
          values = c("overfit"   = "#377eb8",
                     "underfit"  = "#e41a1c",
                     "no misfit" = "grey60"),
          name = NULL,
          drop = FALSE
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          x = NULL,
          y = "Observed − expected restscore correlation",
          caption = caption_text
        ) +
        ggplot2::theme_minimal(base_size = 15) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)

# ---------------------------------------------------------------------------
# Internal helpers (mirror easyRasch2::RMitemRestscoreBoot internals)
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

    # iarm refits on complete cases when the resampled rows contain
    # missing responses -- suppress its per-iteration console message;
    # the behaviour is documented in the HTML note below the table.
    i1 <- suppressMessages(as.data.frame(iarm::item_restscore(model_fit)))
    res_mat <- i1[[1L]]
    n_items <- length(data_list$item_names)

    observed <- as.numeric(res_mat[seq_len(n_items), 1L])
    expected <- as.numeric(res_mat[seq_len(n_items), 2L])
    p_adj    <- as.numeric(res_mat[seq_len(n_items), 5L])
    # Signed as observed - expected so that overfit is positive,
    # matching the Difference column in the item-restscore analysis.
    diff_val <- observed - expected

    cls <- ifelse(p_adj < 0.05 & diff_val > 0, "overfit",
           ifelse(p_adj < 0.05 & diff_val < 0, "underfit", "no misfit"))

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
