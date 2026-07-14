#' @export
targetingClass <- R6::R6Class(
  "targetingClass",
  inherit = targetingBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met (the model needs
      # at least 2 items).
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$targetingNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model. Select at least 2 items.</p>"
        ))
        return()
      }

      # Extract and validate data
      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      # Strip jmvcore S4 column wrappers
      col_names <- names(df)
      df <- as.data.frame(
        matrix(as.numeric(as.matrix(df)),
               nrow = nrow(df),
               ncol = ncol(df),
               dimnames = list(NULL, col_names))
      )

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      tryCatch({
        item_names <- names(df)
        ci_level   <- self$options$ciLevel / 100

        # Estimation-method selection mirrors easyRasch2::RMtargeting():
        # CML (psychotools) when all response categories have at least 3
        # observations; MML (mirt) otherwise, which is more numerically
        # stable under sparse-category conditions. The same estimator is
        # passed to RMitemParameters() below so the threshold table and
        # the figure agree.
        sparse_items <- sparse_category_items(df, min_n = 3L)
        use_mml      <- length(sparse_items) > 0L
        estimator    <- if (use_mml) "MML" else "CML"

        # Threshold table via easyRasch2: grand-mean-centred threshold
        # locations, numerically identical to RMitemParameters() with the
        # same estimator. The default is the wide layout (one row per
        # item, one column per threshold, no SE/CI); the long layout (one
        # row per threshold) adds the delta-method SE and Wald CI.
        # Package messages (e.g. the empty-respondent drop) are
        # suppressed; the module note reports the sample.
        # NOTE: columns and rows are built here rather than in .init()
        # because the number of thresholds depends on the observed
        # response categories per item -- it cannot be derived from the
        # options alone (defensible Level 3 case).
        long_format <- isTRUE(self$options$longFormat)
        table <- self$results$thresholdTable

        if (long_format) {
          thr <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemParameters(
              df,
              estimator = estimator,
              format    = "long",
              se        = TRUE,
              ci_level  = ci_level,
              output    = "dataframe"
            )
          ))

          table$addColumn(name = "item", title = "Item", type = "text",
                          combineBelow = TRUE)
          table$addColumn(name = "threshold", title = "Threshold",
                          type = "text")
          table$addColumn(name = "location", title = "Location",
                          type = "number", format = "zto")
          table$addColumn(name = "se", title = "SE",
                          type = "number", format = "zto")
          table$addColumn(name = "ciLow", title = "Lower",
                          type = "number", format = "zto",
                          superTitle = "CI")
          table$addColumn(name = "ciHigh", title = "Upper",
                          type = "number", format = "zto",
                          superTitle = "CI")

          for (i in seq_len(nrow(thr))) {
            table$addRow(rowKey = i, values = list(
              item      = thr$item[i],
              threshold = paste0("T", thr$threshold[i]),
              location  = thr$location[i],
              se        = thr$se[i],
              ciLow     = thr$ci_lower[i],
              ciHigh    = thr$ci_upper[i]
            ))
          }
          table$setNote("ci", paste0(
            "Wald confidence intervals (", self$options$ciLevel,
            "%): location ± z × SE."
          ))
        } else {
          thr_w <- suppressWarnings(suppressMessages(
            easyRasch2::RMitemParameters(
              df,
              estimator = estimator,
              format    = "wide",
              se        = FALSE,
              output    = "dataframe"
            )
          ))
          t_cols <- setdiff(names(thr_w), c("item", "location"))

          table$addColumn(name = "item", title = "Item", type = "text")
          for (tc in t_cols) {
            table$addColumn(
              name       = tc,
              title      = toupper(tc),
              type       = "number",
              format     = "zto",
              superTitle = "Threshold location"
            )
          }
          table$addColumn(name = "location", title = "Location",
                          type = "number", format = "zto")

          for (i in seq_len(nrow(thr_w))) {
            vals <- list(item = thr_w$item[i],
                         location = thr_w$location[i])
            for (tc in t_cols) vals[[tc]] <- thr_w[[tc]][i]
            table$addRow(rowKey = i, values = vals)
          }
          table$setNote("wide", paste0(
            if (length(t_cols) > 0) {
              "Location = mean of the item's threshold locations. "
            } else {
              "Location = item difficulty (dichotomous items have a single threshold). "
            },
            "Enable 'Long format threshold table' for one row per ",
            "threshold with SE and confidence interval."
          ))
        }

        # Rows with no valid responses contribute nothing to estimation
        n_total <- sum(rowSums(!is.na(df)) > 0)

        # Save state for the plot: the figure is drawn by
        # easyRasch2::RMtargeting() inside the render function (the CML/WLE
        # refit it performs there is cheap); storing the data + options
        # keeps the saved analysis small.
        self$results$targetingPlot$setState(list(df = df))

        # Estimation-method / sample-size note
        method_clause <- if (use_mml) {
          paste0(
            "Item thresholds were estimated with <b>MML (mirt)</b> ",
            "because item(s) ", paste(sparse_items, collapse = ", "),
            " have at least one response category with fewer than 3 ",
            "observations; MML is more numerically stable under sparse ",
            "categories. Person locations are weighted likelihood ",
            "estimates (WLE), which are finite at extreme scores."
          )
        } else {
          paste0(
            "Item thresholds estimated with CML (psychotools); all ",
            "response categories have at least 3 observations. Person ",
            "locations are weighted likelihood estimates (WLE), which ",
            "are finite at extreme scores."
          )
        }
        self$results$targetingNote$setContent(paste0(
          "<p>Analysis based on N = ", n_total, " respondents (rows with ",
          "partially missing responses are retained by the estimation). ",
          method_clause,
          " Results are identical to easyRasch2::RMtargeting() and ",
          "RMitemParameters().</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in targeting plot analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Person-item targeting plot (Wright map) — drawn by
    # easyRasch2::RMtargeting(): back-to-back person / threshold histograms
    # plus a per-item threshold dot plot with optional CIs. Every module
    # option maps directly; xlim auto-expands inside the package so
    # nothing is clipped.
    # ---------------------------------------------------------------------
    .targetingPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required but is not installed.")
      }

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMtargeting(
          image$state$df,
          robust     = isTRUE(self$options$robust),
          sort_items = self$options$sortItems,
          bins       = self$options$bins,
          xlim       = c(self$options$xlimLow, self$options$xlimHigh),
          ci_level   = if (isTRUE(self$options$showCi))
                         self$options$ciLevel / 100
                       else NULL,
          output     = "patchwork"
        )
      ))
      p <- er2_bump_text(p)

      print(p)
      TRUE
    }
  )
)
