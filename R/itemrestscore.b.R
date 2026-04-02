#' @export
itemrestscoreClass <- R6::R6Class(
  "itemrestscoreClass",
  inherit = itemrestscoreBase,
  private = list(
    .run = function() {
      # Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0) {
        return()
      }

      # Check minimum number of variables
      if (length(self$options$vars) < 2) {
        stop("You need at least two variables to run an analysis.")
      }

      data <- self$data
      vars <- self$options$vars
      p_adj <- self$options$pAdj
      sort_by_diff <- self$options$sortByDiff

      # Select only the specified variables (drop = FALSE keeps it as dataframe)
      df <- data[, vars, drop = FALSE]

      # Convert to numeric - handle factors properly
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          df[[col]] <- as.numeric(as.character(df[[col]]))
        } else {
          df[[col]] <- as.numeric(df[[col]])
        }
      }

      # Check for duplicate/identical variables using correlation
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
        pair_strings <- sapply(identical_pairs, function(p) paste0("'", p[1], "' and '", p[2], "'"))
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

      # Check for all-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste(
          "The following variables contain no valid numeric data:",
          paste(bad_vars, collapse = ", ")
        ))
      }

      # Validate data using shared helper (checks non-negative integers, starts at 0)
      validate_response_data(df)

      # Check for sufficient complete cases
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. Each row must have responses for all selected items.")
      }
      if (n_complete < 30) {
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable with small samples.",
          n = n_complete
        )
      }

      # Check for sufficient response variation per item
      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2) {
          stop(paste0(
            "Item '", col, "' has no variation in responses. ",
            "Each item needs at least two different response values."
          ))
        }
      }

      # Run analysis (logic inlined from easyRasch2::RMitemrestscore)
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

          i1 <- iarm::item_restscore(erm_out, p.adj = p_adj)
          i1 <- as.data.frame(i1)

          res_mat <- i1[[1]]
          observed <- round(as.numeric(res_mat[seq_len(n_items), 1L]), 2)
          expected <- round(as.numeric(res_mat[seq_len(n_items), 2L]), 2)
          p_adjusted <- round(as.numeric(res_mat[seq_len(n_items), 5L]), 3)
          significance <- as.character(res_mat[seq_len(n_items), 6L])

          # Assemble result data.frame
          results <- data.frame(
            Item                = names(df),
            Observed            = observed,
            Expected            = expected,
            Absolute_difference = round(abs(expected - observed), 3),
            p_adjusted          = p_adjusted,
            Significance        = significance,
            Location            = round(item_avg_locations, 2),
            Relative_location   = round(relative_item_avg_locations, 2),
            stringsAsFactors    = FALSE,
            row.names           = NULL
          )

          # Sort if requested
          if (isTRUE(sort_by_diff)) {
            results <- results[order(results$Absolute_difference, decreasing = TRUE), ]
            rownames(results) <- NULL
          }

          # Populate the results table
          table <- self$results$restscoreTable
          for (i in seq_len(nrow(results))) {
            table$setRow(rowNo = i, values = list(
              item        = results$Item[i],
              observed    = results$Observed[i],
              expected    = results$Expected[i],
              absDiff     = results$Absolute_difference[i],
              pAdjusted   = results$p_adjusted[i],
              significance = results$Significance[i],
              location    = results$Location[i],
              relLocation = results$Relative_location[i]
            ))
          }
        },
        error = function(e) {
          stop(paste("Error in item-restscore analysis:", e$message))
        }
      )
    }
  )
)
