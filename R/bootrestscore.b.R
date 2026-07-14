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
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$bootstrapTable$setNote("sparse", sparse_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$bootstrapTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item are dropped up
      # front: the bundled easyRasch2 release cannot fit all-NA rows
      # (psychotools errors on polytomous and crashes on dichotomous data),
      # and such rows would poison the bootstrap resampling pool.
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      # 3. Read options
      iterations <- self$options$iterations
      samplesize <- self$options$samplesize
      cutoff     <- self$options$cutoff
      seed       <- self$options$seed

      # Clamp samplesize to nrow(df) instead of failing (Jamovi-friendly;
      # easyRasch2::RMitemRestscoreBoot errors when samplesize > nrow)
      samplesize_used <- min(samplesize, nrow(df))
      samplesize_clamped <- samplesize_used < samplesize

      # 4. Run analysis
      tryCatch({
        item_names <- colnames(df)
        n_items    <- ncol(df)

        # --- Bootstrap via easyRasch2 ---------------------------------------
        # output = "raw" returns the per-iteration long data (iteration,
        # Item, item_restscore, diff, diff_abs) that both the percentage
        # table and the violin plot are built from -- one bootstrap run,
        # numerically identical to RMitemRestscoreBoot() with the same seed.
        # (The package's summary output cannot feed the plot, so the
        # counts-to-percentages aggregation stays module-side; percentages
        # are computed unrounded from the raw counts.) Package warnings are
        # suppressed -- the module surfaces its own footnotes.
        # The bootstrap is the expensive part, and jamovi reruns .run on
        # every option change -- including changes (cutoff, sortBy,
        # showPlot) that do not affect the resampling. The hidden simCache
        # element carries the raw per-iteration data and the full-sample
        # restscore fit; jamovi clears its state exactly when a
        # bootstrap-relevant option changes (its clearWith list), and the
        # signature check makes reuse self-validating rather than relying
        # on clearWith alone. samplesize enters via its clamped value.
        sim_sig <- list(
          iterations = iterations,
          samplesize = samplesize_used,
          seed       = as.integer(seed)
        )
        cached <- self$results$simCache$state
        if (!is.null(cached) && !is.null(cached$fit_all) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df)) {
          fit_all           <- cached$fit_all
          obs_df            <- cached$obs_df
          actual_iterations <- cached$actual_iterations
        } else {
          fit_all <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemRestscoreBoot(
              df,
              iterations = iterations,
              samplesize = samplesize_used,
              parallel   = FALSE,
              seed       = as.integer(seed),
              output     = "raw"
            )
          ))
          # Upstream signs diff as (expected - observed); the module's
          # Difference convention (matching the asymptotic item-restscore
          # analysis) is (observed - expected), so flip the sign here.
          fit_all$diff <- -fit_all$diff

          # Guard against misleading percentages: with few successful
          # iterations the classification shares are based on tiny
          # denominators (there is no observed-only fallback -- the
          # bootstrap is the analysis). Failed iterations are discarded
          # inside the package, so count what came back.
          actual_iterations <- length(unique(fit_all$iteration))
          if (actual_iterations < 20L) {
            stop(paste0(
              "Only ", actual_iterations, " of ", iterations, " bootstrap ",
              "iterations succeeded -- too few for trustworthy ",
              "classification percentages. This typically happens when ",
              "items have very low or very high endorsement rates, so that ",
              "resampled datasets often contain items without response ",
              "variation. Consider a larger bootstrap sample size."
            ), call. = FALSE)
          }

          # --- Full-sample relative locations --------------------------------
          # Same full-sample CML/WLE fit the package's own summary output
          # reports (numerically identical to RMitemRestscoreBoot()'s
          # Relative_location column).
          obs_df <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemRestscore(df, output = "dataframe")
          ))

          self$results$simCache$setState(list(
            fit_all           = fit_all,
            obs_df            = obs_df,
            actual_iterations = actual_iterations,
            sig               = sim_sig,
            df                = df
          ))
        }
        relative_item_avg_locations <-
          obs_df$Relative_location[match(item_names, obs_df$Item)]

        # --- Per-item classification counts ----------------------------------
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
        counts$percent <- counts$n * 100 / per_item_total[counts$Item]

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
          "Flagged = the item was classified as overfit (or underfit) in ",
          "more than ", cutoff, "% of the bootstrap iterations."
        ))
        table$setNote("loc", paste0(
          "Rel. location = item location relative to the mean person ",
          "location (weighted likelihood estimates, WLE; full sample)."
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
          clamp_msg, missing_msg,
          iteration_note(iterations, 250L),
          low_iteration_caveat(actual_iterations), "</p>"
        )
        self$results$bootstrapNote$setContent(note_html)

        # 7. Save state for plot (item order follows the table sort; the
        # raw per-iteration data is read from the simCache element inside
        # the render function -- single storage)
        if (isTRUE(self$options$showPlot)) {
          self$results$bootstrapPlot$setState(list(
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
    # Plot: per-item violin + jitter of (observed - expected) across
    # bootstrap iterations, coloured by per-iteration classification
    # ---------------------------------------------------------------------
    .bootstrapPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state <- image$state
      # Raw per-iteration data lives in the hidden simCache element
      # (single storage; also serves as the bootstrap cache for .run).
      d <- self$results$simCache$state$fit_all
      if (is.null(d)) return(FALSE)
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
        er2_plot_caption() +
        ggplot2::theme(legend.position = "none")

      print(p)
      TRUE
    }
  )
)
