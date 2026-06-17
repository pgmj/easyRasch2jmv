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

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$restscoreTable$setNote("duplicate", dup_msg)

      # Sufficient complete cases?
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. Each row must have responses for all selected items.")
      }

      # Run analysis (logic inlined from easyRasch2::RMitemRestscore)
      tryCatch(
        {
          data_mat <- as.matrix(df)
          n_items <- ncol(df)

          # Fit Rasch model and compute item/person locations
          if (max(data_mat, na.rm = TRUE) == 1L) {
            # Dichotomous: Rasch model
            erm_out <- eRm::RM(df)
            item_avg_locations <- stats::coef(erm_out, "beta") * -1
            pp <- eRm::person.parameter(erm_out)
            person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
          } else {
            # Polytomous: Partial Credit Model
            erm_out <- eRm::PCM(df)
            thresh_obj <- eRm::thresholds(erm_out)
            thresh_table <- thresh_obj$threshtable[[1]]
            item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
            pp <- eRm::person.parameter(erm_out)
            person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
          }

          relative_item_avg_locations <- item_avg_locations - person_avg_location

          # Temporarily set rgl.useNULL to avoid rgl device issues during iarm fitting
          old_rgl <- getOption("rgl.useNULL")
          options(rgl.useNULL = TRUE)
          on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

          # iarm refits the model on complete cases (na.omit) when any
          # responses are missing -- suppress its console message; the
          # behaviour is documented in the HTML note below the table.
          # BH adjustment is hardcoded module-wide.
          i1 <- suppressMessages(iarm::item_restscore(erm_out, p.adj = "BH"))
          i1 <- as.data.frame(i1)

          # Pass raw numerics (no pre-rounding) so the jamovi frontend
          # applies the user's "Number format" preferences -- matches the
          # convention used by the newer analyses in this module.
          res_mat <- i1[[1]]
          observed <- as.numeric(res_mat[seq_len(n_items), 1L])
          expected <- as.numeric(res_mat[seq_len(n_items), 2L])
          p_adjusted <- as.numeric(res_mat[seq_len(n_items), 5L])

          # Assemble result data.frame. `Difference` is signed (observed -
          # expected): positive = item over-discriminates (often LD),
          # negative = item under-discriminates (often noise / multi-dim).
          # `Fit` classifies items as overfit (observed > expected) or
          # underfit (observed < expected) when adj. p < .05 -- the same
          # rule and labels as the bootstrap item-restscore analysis.
          difference <- observed - expected
          fit_class <- ifelse(
            !is.na(p_adjusted) & p_adjusted < 0.05 & difference > 0, "overfit",
            ifelse(!is.na(p_adjusted) & p_adjusted < 0.05 & difference < 0,
                   "underfit", "")
          )
          results <- data.frame(
            Item              = names(df),
            Observed          = observed,
            Expected          = expected,
            Difference        = difference,
            p_adjusted        = p_adjusted,
            Fit               = fit_class,
            Relative_location = relative_item_avg_locations,
            stringsAsFactors  = FALSE,
            row.names         = NULL
          )

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
              fit         = results$Fit[i],
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
            "the mean person location, in logits."
          ))

          # Sample-size / missing-data note. The two parts of the table
          # use different samples when responses are missing: iarm refits
          # the model on complete cases for the restscore statistics,
          # while the Location columns come from the eRm fit on all
          # available responses (CML accommodates partial missingness).
          # This mirrors easyRasch2::RMitemRestscore(). Rows with no valid
          # responses at all contribute nothing and are not counted in N.
          n_total <- sum(rowSums(!is.na(df)) > 0)
          missing_clause <- if (n_total > n_complete) {
            paste0(", of whom ", n_complete, " had complete responses on ",
                   "all ", n_items, " items. Restscore correlations and ",
                   "p-values are computed from the ", n_complete,
                   " complete cases (the model is refitted on complete ",
                   "cases for this purpose); item locations use all ",
                   "available responses via eRm's conditional maximum ",
                   "likelihood estimation.")
          } else {
            paste0(", all with complete responses on all ", n_items,
                   " items.")
          }
          self$results$restscoreNote$setContent(paste0(
            "<p>Analysis based on N = ", n_total, " respondents",
            missing_clause, "</p>"
          ))
        },
        error = function(e) {
          stop(paste("Error in item-restscore analysis:", e$message))
        }
      )
    }
  )
)
