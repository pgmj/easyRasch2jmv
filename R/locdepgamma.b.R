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
        # 3. Compute partial gamma LD via iarm (suppress its console output).
        # Missing values are allowed: iarm::partgam_LD() handles NA itself,
        # matching easyRasch2::RMlocdepGamma() which passes data as-is.
        sink(nullfile())
        pgam_raw <- tryCatch(
          iarm::partgam_LD(as.data.frame(df)),
          finally = sink()
        )

        # pgam_raw is a list of two data.frames, one per rest-score
        # direction. Columns from iarm::partgam_LD():
        #   1 Item1, 2 Item2, 3 gamma, 4 se, 5 pvalue, 6 padj.BH, 7 sig,
        #   8 lower, 9 upper. Column 6 (BH-adjusted p) and 7 (a star string
        #   padded with spaces) are accessed positionally because iarm does
        #   not always give them stable names -- matches easyRasch2 source.
        # lower/upper are the 95% CI (gamma +/- 1.96 * se).
        process_pgam_df <- function(raw_df) {
          g <- as.numeric(raw_df$gamma)
          data.frame(
            Item1   = as.character(raw_df$Item1),
            Item2   = as.character(raw_df$Item2),
            gamma   = ifelse(is.finite(g), g, NA_real_),
            se      = as.numeric(raw_df$se),
            lower   = as.numeric(raw_df$lower),
            upper   = as.numeric(raw_df$upper),
            padj_bh = as.numeric(raw_df[[6]]),
            sig     = trimws(as.character(raw_df[[7]])),
            stringsAsFactors = FALSE
          )
        }

        result_list <- list(
          process_pgam_df(pgam_raw[[1]]),
          process_pgam_df(pgam_raw[[2]])
        )

        # 4. Filter / sort pipeline, applied per direction:
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

        # 5. Populate tables. Rows were pre-created in .init() only in the
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
              sig    = d$sig[i]
            )
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
          if (nrow(d) == 0L) {
            tables[[idx]]$setNote(
              "empty",
              "No item pairs met the filter criteria."
            )
          }
        }

        # 6. Caption note
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
        self$results$ldNote$setContent(paste0(
          "<p>Partial gamma LD analysis (n = ", n_complete,
          " complete cases). Values near 0 indicate no local dependence; ",
          "large positive values suggest positive LD (items share variance ",
          "beyond the latent trait), large negative values suggest negative ",
          "LD. Because it matters which item of a pair is subtracted from ",
          "the total score, each pair is tested in both rest-score ",
          "directions.", se_clause,
          if (length(filter_clauses) > 0)
            paste0(" ", paste(filter_clauses, collapse = " "))
          else "",
          "</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in partial gamma LD analysis:", e$message))
      })
    }
  )
)
