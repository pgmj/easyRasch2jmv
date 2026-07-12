#' @export
reliabilityClass <- R6::R6Class(
  "reliabilityClass",
  inherit = reliabilityBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- build the fixed 4-row table structure up front so the table
    # renders immediately instead of flickering from a blank placeholder to a
    # populated table when .run() finishes. The row count (always 4) and their
    # metric labels are fully determined by options (estim), not by the data,
    # so they belong here. .run() fills in the computed estimates via setRow().
    # ---------------------------------------------------------------------
    .init = function() {
      if (is.null(self$options$vars) || length(self$options$vars) < 3)
        return()

      estim <- self$options$estim
      table <- self$results$relTable

      labels <- c(
        "Cronbach's alpha",
        "PSI",
        "Marginal",
        paste0("RMU (", estim, ")")
      )
      for (i in seq_along(labels)) {
        table$addRow(rowKey = i, values = list(
          metric   = labels[i],
          estimate = NA_real_,
          lower    = NA_real_,
          upper    = NA_real_,
          notes    = ""
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met. With 2
      # dichotomous items the MML model (mirt) behind the RMU estimate is
      # not estimable (too few degrees of freedom), and reliability
      # estimates from 2-item scales are generally not informative, so
      # require 3 items.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$relNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With 2 ",
          "dichotomous items the latent model used for the RMU estimate ",
          "cannot be estimated (too few degrees of freedom), and ",
          "reliability estimates from 2-item scales are generally not ",
          "informative. Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Required suggested package
      if (!requireNamespace("ggdist", quietly = TRUE))
        stop("Package 'ggdist' is required. Install with: install.packages(\"ggdist\")")

      # 3. Extract data and convert to numeric
      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$relTable$setNote("sparse", sparse_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$relTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item are dropped up
      # front: the bundled easyRasch2 release cannot fit all-NA rows
      # (psychotools errors on polytomous and crashes on dichotomous data).
      n_total <- nrow(df)
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      # 4. Read options
      estim       <- self$options$estim
      draws       <- self$options$draws
      rmu_iter    <- self$options$rmuIter
      conf_int    <- self$options$confInt / 100
      theta_range <- c(self$options$thetaMin, self$options$thetaMax)
      boot_cis    <- isTRUE(self$options$bootAlpha)
      boot_iter   <- self$options$bootIter
      seed        <- self$options$seed

      if (theta_range[1L] >= theta_range[2L])
        stop("Theta lower bound must be less than upper bound.")

      tryCatch({
        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        # 5. All computation is delegated to the easyRasch2 package:
        # Cronbach's alpha (closed-form, complete cases), the WLE-based
        # PSI (native CML/WLE engine, min/max scorers excluded), the
        # native marginal reliability (Green, 1984; CML test information
        # integrated over the estimated normal latent density), and RMU
        # from mirt plausible values. When the bootstrap is enabled,
        # respondents are resampled and alpha / PSI / Marginal are
        # recomputed natively per resample for HDCIs. Results are
        # numerically identical to RMreliability() with the same seed.
        # Package warnings are suppressed -- the module surfaces its own
        # footnotes.
        results <- suppressWarnings(suppressMessages(
          easyRasch2::RMreliability(
            df,
            conf_int    = conf_int,
            draws       = draws,
            rmu_iter    = rmu_iter,
            estim       = estim,
            boot        = boot_cis,
            boot_iter   = boot_iter,
            parallel    = FALSE,
            seed        = as.integer(seed),
            theta_range = theta_range,
            output      = "dataframe"
          )
        ))

        # 6. Populate table (rows created in .init(); fill values here).
        # Row order from the package is fixed: alpha, PSI, Marginal, RMU.
        table <- self$results$relTable
        for (i in seq_len(nrow(results))) {
          table$setRow(rowKey = i, values = list(
            metric   = results$metric[i],
            estimate = results$estimate[i],
            lower    = results$lower[i],
            upper    = results$upper[i],
            notes    = results$notes[i]
          ))
        }
        table$setNote(
          "context",
          paste0(
            "PSI is the WLE-based person-separation reliability (CML item ",
            "parameters via psychotools) and excludes respondents with ",
            "min/max raw scores. Marginal is the model-based marginal ",
            "reliability (Green, 1984): CML test information integrated ",
            "over the estimated latent distribution -- a large gap between ",
            "PSI and Marginal suggests the sample is off-target relative ",
            "to the scale. RMU uses plausible values from an MML model ",
            "(mirt); theta estimator for the RMU draws: ", estim, "."
          )
        )
        table$setNote(
          "hdci",
          paste0(
            "HDCI = highest-density continuous interval (width set by the ",
            "HDCI width option; here ", round(conf_int * 100, 1),
            "%). Available for Cronbach's alpha, PSI, and Marginal when ",
            "the bootstrap is enabled, and always for RMU. RMU = Relative ",
            "Measurement Uncertainty."
          )
        )

        # 7. Caption. The estimates use different samples when data are
        # missing: Cronbach's alpha is closed-form on complete cases,
        # while the model-based estimates retain partially missing rows
        # (CML / MML estimation).
        n_used <- nrow(df)
        drop_clause <- if (n_used < n_total) {
          paste0(" (", n_total - n_used, " row(s) without any responses ",
                 "on the selected items excluded)")
        } else ""
        missing_msg <- if (n_used > n_complete) {
          paste0(
            " Cronbach's alpha is computed from the ", n_complete,
            " complete cases; the model-based estimates (PSI, Marginal, ",
            "RMU) use all ", n_used, " rows with at least one response ",
            "(CML and MML estimation accommodate partially missing ",
            "responses)."
          )
        } else {
          ""
        }
        self$results$relNote$setContent(
          paste0(
            "<p>Reliability based on N = ", n_used,
            " respondents across ", ncol(df), " items",
            drop_clause,
            if (n_used > n_complete)
              paste0(" (", n_complete, " with complete responses)")
            else "",
            ".", missing_msg, "</p>"
          )
        )

      }, error = function(e) {
        stop(paste("Error in reliability analysis:", e$message))
      })
    }
  )
)
