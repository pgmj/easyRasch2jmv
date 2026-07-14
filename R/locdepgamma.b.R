#' @export
locdepgammaClass <- R6::R6Class(
  "locdepgammaClass",
  inherit = locdepgammaBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- pre-create the table rows so both tables render their full
    # structure immediately instead of flickering from a blank placeholder
    # once .run() finishes. This is only possible when no result-dependent
    # filter is active: without filters the row count is fully determined
    # by options (choose(k, 2) pairs per direction, capped by nPairs when
    # > 0). When the significance filter or the |gamma| threshold is on,
    # the row count depends on the computed results and rows must be added
    # in .run() instead (defensible Level 3 case -- cannot be moved here).
    # ---------------------------------------------------------------------
    .init = function() {
      vars <- self$options$vars
      if (is.null(vars) || length(vars) < 3)
        return()
      if (isTRUE(self$options$sigOnly) || self$options$gammaThreshold > 0)
        return()

      k           <- length(vars)
      total_pairs <- k * (k - 1L) / 2L
      n_pairs     <- self$options$nPairs
      n_rows      <- if (n_pairs > 0L) min(n_pairs, total_pairs) else total_pairs

      blank <- list(item1 = "", item2 = "", gamma = NA_real_,
                    se = NA_real_, lower = NA_real_, upper = NA_real_,
                    padjBH = NA_real_, sig = "")
      for (tbl in list(self$results$dir1Table, self$results$dir2Table)) {
        for (i in seq_len(n_rows)) {
          tbl$addRow(rowKey = i, values = blank)
        }
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / reject if requirements not met.
      # Partial gamma LD conditions each item pair on the rest score (total
      # minus one of the items in the pair). With only 2 items the rest
      # score degenerates to the other item in the pair and iarm returns
      # all-NaN results, so at least 3 items are required.
      vars <- self$options$vars
      if (is.null(vars) || length(vars) == 0)
        return()
      if (length(vars) < 3) {
        self$results$ldNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. Partial gamma ",
          "conditions each item pair on the rest score (total score minus ",
          "one of the items in the pair); with only 2 items the rest score ",
          "reduces to the other item in the pair and the statistic is ",
          "undefined. Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Extract and validate data
      data <- self$data
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg)) {
        self$results$dir1Table$setNote("sparse", sparse_msg)
        self$results$dir2Table$setNote("sparse", sparse_msg)
      }

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg)) {
        self$results$dir1Table$setNote("duplicate", dup_msg)
        self$results$dir2Table$setNote("duplicate", dup_msg)
      }

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      # rgl workaround (iarm depends on vcdExtra -> rgl)
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      tryCatch({
        # 3. Optional simulation-based expected ranges (computed first --
        # they feed RMlocdepGamma()). If the simulation cannot deliver
        # reliable ranges, degrade gracefully: the tables are shown
        # without the expected range / flags and the note explains why.
        compute_cutoff <- isTRUE(self$options$computeCutoff)
        cutoff_res   <- NULL
        sim_fail_msg <- NULL

        # The simulation is the expensive part, and jamovi reruns .run on
        # every option change -- including changes (filters, sorting,
        # plotPairs, showSE) that do not affect the simulation. The hidden
        # simCache element carries the cutoff object; jamovi clears its
        # state exactly when a simulation-relevant option changes (its
        # clearWith list), and the signature check makes reuse
        # self-validating rather than relying on clearWith alone. (The
        # plot state cannot serve as the cache here: its clearWith also
        # includes plotPairs, which must keep forcing a re-render.)
        sim_sig <- list(
          iterations = self$options$iterations,
          seed       = as.integer(self$options$seed),
          hdci_width = self$options$hdciWidth / 100
        )
        cached <- self$results$simCache$state
        if (compute_cutoff &&
            !is.null(cached) && !is.null(cached$cutoff_res) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df)) {
          cutoff_res <- cached$cutoff_res
        } else if (compute_cutoff) {
          cutoff_res <- tryCatch(
            suppressWarnings(suppressMessages(
              easyRasch2::RMlocdepGammaCutoff(
                df,
                iterations = self$options$iterations,
                parallel   = FALSE,
                seed       = as.integer(self$options$seed),
                hdci_width = self$options$hdciWidth / 100
              )
            )),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )
          if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
            sim_fail_msg <- paste0(
              "Only ", cutoff_res$actual_iterations, " of ",
              self$options$iterations, " simulation iterations succeeded ",
              "-- too few to estimate reliable expected ranges. Check ",
              "your data: items must have sufficient response variation, ",
              "and the sample must be large enough for stable Rasch ",
              "model estimation."
            )
            cutoff_res <- NULL
          }
        }

        # 4. Partial gamma LD via easyRasch2 (numerically identical to
        # RMlocdepGamma(); a wrapper around iarm::partgam_LD()). Missing
        # values are allowed -- iarm handles NA itself.
        result_list <- suppressWarnings(suppressMessages(
          easyRasch2::RMlocdepGamma(
            df,
            cutoff = cutoff_res,
            output = "dataframe"
          )
        ))
        result_list <- list(result_list$direction1, result_list$direction2)

        # The bundled easyRasch2 release omits the se / lower / upper (95%
        # Wald CI) columns from the dataframe output (added in a later
        # version); when absent, source them from the same underlying
        # iarm::partgam_LD() call and join by pair, so the values are
        # identical by construction. iarm's console output is sunk; the
        # finally handler guarantees the sink is removed exactly once.
        if (!"se" %in% names(result_list[[1L]])) {
          sink(nullfile())
          pgam_raw <- tryCatch(
            iarm::partgam_LD(as.data.frame(df)),
            finally = sink()
          )
          for (idx in seq_along(result_list)) {
            d   <- result_list[[idx]]
            raw <- pgam_raw[[idx]]
            key_d   <- paste(d$Item1, d$Item2, sep = "___")
            key_raw <- paste(as.character(raw$Item1),
                             as.character(raw$Item2), sep = "___")
            m <- match(key_d, key_raw)
            d$se    <- as.numeric(raw$se)[m]
            d$lower <- as.numeric(raw$lower)[m]
            d$upper <- as.numeric(raw$upper)[m]
            result_list[[idx]] <- d
          }
        }

        # Non-finite gammas (degenerate pairs) display as empty cells
        for (idx in seq_along(result_list)) {
          d <- result_list[[idx]]
          d$gamma <- ifelse(is.finite(d$gamma), d$gamma, NA_real_)
          result_list[[idx]] <- d
        }

        # 5. Filter / sort pipeline, applied per direction:
        #    significance filter -> |gamma| threshold -> top-N by |gamma|
        #    (which sorts) -> otherwise sort by |gamma| if requested.
        # With filters active the two directions can retain different
        # numbers of rows (each direction has its own gammas / p-values).
        sig_only    <- isTRUE(self$options$sigOnly)
        gamma_thr   <- self$options$gammaThreshold
        n_pairs     <- self$options$nPairs
        sort_gamma  <- isTRUE(self$options$sortByGamma)
        total_pairs <- nrow(result_list[[1]])

        filter_applied <- FALSE
        for (idx in seq_along(result_list)) {
          d <- result_list[[idx]]
          if (sig_only) {
            d <- d[!is.na(d$padj_bh) & d$padj_bh < 0.05, , drop = FALSE]
          }
          if (gamma_thr > 0) {
            d <- d[!is.na(d$gamma) & abs(d$gamma) >= gamma_thr, , drop = FALSE]
          }
          if (n_pairs > 0L && n_pairs < nrow(d)) {
            filter_applied <- TRUE
            ord <- order(abs(d$gamma), decreasing = TRUE)
            d   <- d[ord[seq_len(n_pairs)], , drop = FALSE]
          } else if (sort_gamma) {
            d <- d[order(abs(d$gamma), decreasing = TRUE), , drop = FALSE]
          }
          rownames(d) <- NULL
          result_list[[idx]] <- d
        }

        # 6. Populate tables. Rows were pre-created in .init() only in the
        # unfiltered case (same condition as there); with filters active
        # the row count is result-dependent and rows are added here.
        rows_pre_created <- !sig_only && gamma_thr == 0
        tables <- list(self$results$dir1Table, self$results$dir2Table)
        for (idx in seq_along(tables)) {
          d <- result_list[[idx]]
          for (i in seq_len(nrow(d))) {
            vals <- list(
              item1  = d$Item1[i],
              item2  = d$Item2[i],
              gamma  = d$gamma[i],
              se     = d$se[i],
              lower  = d$lower[i],
              upper  = d$upper[i],
              padjBH = d$padj_bh[i],
              sig    = d$Significance[i]
            )
            if (!is.null(cutoff_res)) {
              vals$gammaLow  <- d$gamma_low[i]
              vals$gammaHigh <- d$gamma_high[i]
              vals$flagged   <- ifelse(isTRUE(d$flagged[i]), "TRUE", "")
            }
            if (rows_pre_created) {
              tables[[idx]]$setRow(rowNo = i, values = vals)
            } else {
              tables[[idx]]$addRow(rowKey = i, values = vals)
            }
          }
          tables[[idx]]$setNote(
            "bh",
            paste0(
              "BH = Benjamini-Hochberg false-discovery-rate correction ",
              "for multiple testing."
            )
          )
          tables[[idx]]$setNote(
            "direction",
            paste0(
              "Partial gamma between Item 1 and Item 2, controlling for ",
              "the rest score (total score minus Item 2). Each item pair ",
              "appears in both tables with the two items swapped, so the ",
              "two tables together test both rest-score directions for ",
              "every pair."
            )
          )
          if (!is.null(cutoff_res)) {
            tables[[idx]]$setNote("flag", paste0(
              "Expected range = ", cutoff_res$hdci_width * 100, "% HDCI ",
              "of partial gamma values simulated under the fitted ",
              "unidimensional model (no true local dependence). ",
              "Flagged = TRUE when the observed gamma falls outside the ",
              "expected range."
            ))
          }
          if (nrow(d) == 0L) {
            tables[[idx]]$setNote(
              "empty",
              "No item pairs met the filter criteria."
            )
          }
        }

        # 7. Save the simulation cache. The per-pair plot reads the cutoff
        # object and data from here too (via .ldPlot), so nothing is
        # stored twice. On simulation failure the cache is cleared so the
        # plot cannot render from a stale cutoff object.
        if (!is.null(cutoff_res)) {
          self$results$simCache$setState(list(
            cutoff_res = cutoff_res, sig = sim_sig, df = df
          ))
        } else if (compute_cutoff) {
          self$results$simCache$setState(NULL)
        }

        # 8. Caption note
        filter_clauses <- character(0)
        if (sig_only) {
          filter_clauses <- c(filter_clauses, paste0(
            "Showing only pairs with BH-adjusted p &lt; .05; note that ",
            "statistical significance depends on sample size (large samples ",
            "flag trivially small gammas, small samples may miss ",
            "substantial LD)."
          ))
        }
        if (gamma_thr > 0) {
          filter_clauses <- c(filter_clauses, paste0(
            "Showing only pairs with |gamma| ≥ ", gamma_thr, "."
          ))
        }
        if (filter_applied) {
          filter_clauses <- c(filter_clauses, paste0(
            "Showing the top ", n_pairs, " of ", total_pairs,
            " pairs by |gamma| per direction."
          ))
        }
        se_clause <- if (isTRUE(self$options$showSE)) {
          " Confidence intervals are 95% Wald intervals (gamma ± 1.96 × SE)."
        } else ""
        cutoff_clause <- if (!is.null(cutoff_res)) {
          paste0(
            " Expected ranges based on ", cutoff_res$actual_iterations,
            " simulation iterations (", cutoff_res$hdci_width * 100,
            "% HDCI); results are identical to easyRasch2::RMlocdepGamma() ",
            "and RMlocdepGammaCutoff() with the same seed.",
            iteration_note(self$options$iterations, 250L),
            low_iteration_caveat(cutoff_res$actual_iterations)
          )
        } else if (!is.null(sim_fail_msg)) {
          paste0(" <b>Simulation-based expected ranges unavailable:</b> ",
                 sim_fail_msg)
        } else ""
        self$results$ldNote$setContent(paste0(
          "<p>Partial gamma LD analysis (n = ", n_complete,
          " complete cases). Values near 0 indicate no local dependence; ",
          "large positive values suggest positive LD (items share variance ",
          "beyond the latent trait), large negative values suggest negative ",
          "LD. Because it matters which item of a pair is subtracted from ",
          "the total score, each pair is tested in both rest-score ",
          "directions.", se_clause, cutoff_clause,
          if (length(filter_clauses) > 0)
            paste0(" ", paste(filter_clauses, collapse = " "))
          else "",
          "</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in partial gamma LD analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Per-pair simulated-gamma dot plot — easyRasch2::RMlocdepGammaPlot()
    # (ggdist dot cloud + black per-pair median + orange diamonds for the
    # observed partial gamma), restyled to the module's plot conventions
    # (base size 15).
    # ---------------------------------------------------------------------
    .ldPlot = function(image, ggtheme, theme, ...) {
      # The cutoff object and data live in the hidden simCache element
      # (single storage; also serves as the simulation cache for .run).
      state <- self$results$simCache$state
      if (is.null(state) || is.null(state$cutoff_res)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMlocdepGammaPlot(
          state$cutoff_res,
          data    = state$df,
          n_pairs = self$options$plotPairs
        )
      ))
      if (is.null(p)) return(FALSE)

      # Re-apply the module's theme on top of the package theme: base size
      # 15 for jamovi's 500x600 canvas, then restore the theme pieces that
      # theme_minimal() resets.
      p <- p +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
