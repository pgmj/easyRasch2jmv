#' @export
itemrestscoreClass <- R6::R6Class(
  "itemrestscoreClass",
  inherit = itemrestscoreBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met.
      # With 2 items the restscore degenerates to the other item of the
      # pair and observed = expected exactly (p = 1) for both items, so
      # at least 3 items are required for a meaningful analysis.
      if (is.null(self$options$vars) || length(self$options$vars) == 0) {
        return()
      }
      if (length(self$options$vars) < 3) {
        self$results$restscoreNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the restscore (total score minus the item) reduces to the ",
          "other item, making observed and expected correlations identical ",
          "by construction. Select at least 3 items.</p>"
        ))
        return()
      }

      data <- self$data
      vars <- self$options$vars
      sort_by_diff <- self$options$sortByDiff

      # Select the specified variables (drop = FALSE keeps it as data.frame)
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$restscoreTable$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$restscoreTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$restscoreTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item are dropped up
      # front: the bundled easyRasch2 release cannot fit all-NA rows
      # (psychotools errors on polytomous and crashes on dichotomous data).
      n_total <- nrow(df)
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]

      # Sufficient complete cases?
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. The item-restscore statistics require at least one row with responses to all selected items.")
      }

      tryCatch(
        {
          n_items <- ncol(df)

          # All computation is delegated to the easyRasch2 package: CML item
          # estimation via psychotools with WLE person estimates for the
          # relative locations, and iarm::item_restscore() with hardcoded BH
          # adjustment (module-wide convention, also the package default).
          # Results are numerically identical to RMitemRestscore(). Values
          # are as reported by the package's dataframe output. Package
          # warnings (e.g. sparse categories) are suppressed -- the module
          # surfaces its own sparse-category footnote above.
          results <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemRestscore(df, output = "dataframe")
          ))

          # Sort by absolute magnitude when requested, so both over- and
          # underfit items rise to the top while the signed value remains
          # visible in the table.
          if (isTRUE(sort_by_diff)) {
            results <- results[order(abs(results$Difference), decreasing = TRUE), ]
            rownames(results) <- NULL
          }

          # Populate the results table
          table <- self$results$restscoreTable
          for (i in seq_len(nrow(results))) {
            table$setRow(rowNo = i, values = list(
              item        = results$Item[i],
              observed    = results$Observed[i],
              expected    = results$Expected[i],
              difference  = results$Difference[i],
              pAdjusted   = results$p_adjusted[i],
              fit         = results$Flagged[i],
              relLocation = results$Relative_location[i]
            ))
          }

          # Footnotes explaining columns and symbols
          table$setNote("diff", paste0(
            "Difference = observed - expected gamma. Positive: ",
            "over-discrimination (overfit, often local dependence); ",
            "negative: under-discrimination (underfit, often ",
            "multidimensionality or noise). Items are labelled in the ",
            "Flagged column when the adjusted p-value < .05."
          ))
          table$setNote("sig", paste0(
            "P-values adjusted with the Benjamini-Hochberg (BH) ",
            "false-discovery-rate method."
          ))
          table$setNote("loc", paste0(
            "Rel. location = mean item (threshold) location relative to ",
            "the mean person location (weighted likelihood estimates, ",
            "WLE), in logits."
          ))

          # Sample-size / missing-data note. The two parts of the table
          # use different samples when responses are missing: iarm refits
          # the model on complete cases for the restscore statistics,
          # while the Location columns come from the CML fit on all
          # available responses (partial rows retained). This mirrors
          # easyRasch2::RMitemRestscore(). Rows with no valid responses at
          # all were dropped above and are reported separately.
          n_used <- nrow(df)
          drop_clause <- if (n_used < n_total) {
            paste0(" (", n_total - n_used, " row(s) without any responses ",
                   "on the selected items excluded)")
          } else ""
          missing_clause <- if (n_used > n_complete) {
            paste0(", of whom ", n_complete, " had complete responses on ",
                   "all ", n_items, " items. Restscore correlations and ",
                   "p-values are computed from the ", n_complete,
                   " complete cases (the model is refitted on complete ",
                   "cases for this purpose); item locations use all ",
                   "available responses via conditional maximum ",
                   "likelihood estimation.")
          } else {
            paste0(", all with complete responses on all ", n_items,
                   " items.")
          }
          self$results$restscoreNote$setContent(paste0(
            "<p>Analysis based on N = ", n_used, " respondents",
            drop_clause, missing_clause, "</p>"
          ))
        },
        error = function(e) {
          stop(paste("Error in item-restscore analysis:", e$message))
        }
      )
    }
  )
)
