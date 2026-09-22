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

      # The interval width is an option, so the grouped heading names it
      # rather than leaving a bare "HDCI" that could mean any width.
      hdci_title <- paste0(format(self$options$confInt, trim = TRUE),
                           "% HDCI")
      table$getColumn("lower")$setSuperTitle(hdci_title)
      table$getColumn("upper")$setSuperTitle(hdci_title)

      labels <- c(
        "Cronbach's alpha",
        "PSI",
        "Marginal (curve mean)",
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

      # Conditional-precision summary rows. Fully determined by options, so
      # they are built here; the benchmark rows only exist when asked for.
      if (isTRUE(self$options$showCurve)) {
        ct <- self$results$curveTable
        rows <- list(
          thetamean = "Mean theta",
          thetasd   = "SD theta",
          latmean   = "Latent mean (logits)",
          sigma     = "Latent SD (logits)",
          margin    = "Marginal (curve mean)",
          green     = "Marginal (Green, superseded)",
          sem       = "Average SEM (logits)",
          info      = "Average test information"
        )
        if (isTRUE(self$options$useBenchmark))
          rows$bench <- "Respondents reaching the benchmark (%)"
        for (key in names(rows)) {
          ct$addRow(rowKey = key, values = list(
            metric = rows[[key]], value = NA_real_, notes = ""
          ))
        }
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
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$relTable$setNote("recode", recode_msg)

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
        # native marginal reliability (the latent-density-weighted mean of
        # the conditional reliability, which the curve below averages), and RMU
        # from mirt plausible values. When the bootstrap is enabled,
        # respondents are resampled and alpha / PSI / Marginal are
        # recomputed natively per resample for HDCIs. Results are
        # numerically identical to RMreliability() with the same seed.
        # Package warnings are suppressed -- the module surfaces its own
        # footnotes.
        #
        # Cached on the options it actually depends on. None of the curve
        # options can move any of these four estimates, so ticking the
        # curve on and off, changing its statistic or its benchmark must
        # not pay for the MML fit and the bootstrap again.
        rel_cached <- self$results$relCache$state
        if (!is.null(rel_cached) && !identical(rel_cached$df, df))
          rel_cached <- NULL
        rel_sig <- list(
          estim = estim, draws = draws, rmu_iter = rmu_iter,
          conf_int = conf_int, boot = boot_cis,
          # Only counts when the bootstrap is on, so editing the iteration
          # count with the box unticked does not throw the cache away.
          boot_iter = if (boot_cis) boot_iter else NA_integer_,
          seed = seed, theta = theta_range
        )
        if (!is.null(rel_cached) && identical(rel_cached$rel_sig, rel_sig)) {
          results <- rel_cached$results
        } else {
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
        }
        self$results$relCache$setState(
          list(df = df, rel_sig = rel_sig, results = results))

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
            "min/max raw scores. Marginal is the latent-density-weighted ",
            "mean of the conditional reliability, ",
            "sigma² / (sigma² + SEM(theta)²), over the ",
            "estimated latent distribution -- a large gap between PSI and ",
            "Marginal suggests the sample is off-target relative to the ",
            "scale. RMU uses plausible values from an MML model (mirt); ",
            "theta estimator for the RMU draws: ", estim, "."
          )
        )
        table$setNote(
          "offtarget",
          paste0(
            "Marginal reliability averages the conditional reliability over ",
            "the estimated latent distribution, whose mean and SD are both ",
            "fitted to these data. PSI divides instead by the observed ",
            "spread of the person estimates, so the two approach the same ",
            "coefficient by different routes and a wide gap between them ",
            "flags an off-target or non-normal sample."
          )
        )
        table$setNote(
          "marginalchange",
          paste0(
            "Marginal reliability changed formula in module version 3.2.0 ",
            "and its values are higher than those reported by earlier ",
            "versions, more so on short scales. Enable <i>Conditional ",
            "precision curve</i> to see the superseded Green (1984) ",
            "coefficient alongside it."
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

        # --- Conditional precision curve (easyRasch2::RMreliabilityCurve)
        # The figure is redrawn in the render function, so the summary
        # quantities here are computed with boot = FALSE: none of them
        # depends on the bootstrap, and running it twice for one band
        # would double the only expensive part of this option.
        if (isTRUE(self$options$showCurve)) {
          benchmark <- if (isTRUE(self$options$useBenchmark))
            self$options$benchmark else NULL

          # Cached on its own signature, separately from the reliability
          # estimates: the benchmark and the statistic move the curve
          # summary but nothing above it, and the HDCI width moves the
          # table but nothing here.
          curve_sig <- list(
            statistic = self$options$curveStatistic,
            benchmark = benchmark,
            reference = isTRUE(self$options$curveReference),
            density   = isTRUE(self$options$showDensity),
            n_nodes   = self$options$nNodes,
            theta     = theta_range,
            # The figure is cached alongside the summary, so the options
            # that move only the bootstrap band belong in this signature
            # too. They cannot touch curve_df, which is always computed
            # with boot = FALSE, so a band change recomputes a 0.05 s
            # data frame it did not have to. That is the cheaper mistake.
            # All three count only while the band is drawn. With it off,
            # neither the interval width nor the seed can move the figure
            # or the summary, so editing them must not cost a refit.
            boot      = isTRUE(self$options$curveBoot),
            boot_iter = if (isTRUE(self$options$curveBoot))
                          self$options$curveBootIter else NA_integer_,
            conf_int  = if (isTRUE(self$options$curveBoot))
                          conf_int else NA_real_,
            seed      = if (isTRUE(self$options$curveBoot))
                          seed else NA_integer_
          )
          curve_cached <- self$results$curveCache$state
          if (!is.null(curve_cached) && !identical(curve_cached$df, df))
            curve_cached <- NULL
          if (!is.null(curve_cached) &&
              identical(curve_cached$curve_sig, curve_sig)) {
            curve_df   <- curve_cached$curve_df
            loc        <- curve_cached$loc
            curve_plot <- curve_cached$curve_plot
          } else {
            curve_df <- suppressWarnings(suppressMessages(
              easyRasch2::RMreliabilityCurve(
                df,
                statistic   = self$options$curveStatistic,
                benchmark   = benchmark,
                reference   = if (isTRUE(self$options$curveReference))
                                "marginal" else "none",
                boot        = FALSE,
                show_density = isTRUE(self$options$showDensity),
                n_nodes     = self$options$nNodes,
                output      = "dataframe"
              )
            ))
            # Where this sample actually sits, which is what the density
            # behind the curve draws. Distinct from the latent SD below:
            # these are WLE point estimates, whose spread is inflated by
            # measurement error, and the labels match Person Parameters.
            loc <- suppressWarnings(suppressMessages(
              easyRasch2::RMpersonParameters(
                df, method = "WLE", estimator = "CML",
                theta_range = theta_range, output = "dataframe"
              )
            ))
            # The figure is built here and stored, not redrawn in the render
            # function. jamovi calls that function on every resize and every
            # export, so redrawing there made each one pay for the bootstrap
            # band again: 2.6 s at the default 200 iterations on 50
            # respondents and 5 items, 25 s at the 2000-iteration maximum.
            # RMreliabilityCurve() has no output mode returning both the
            # summary and the figure, so a cold run costs two curves and a
            # warm one none, which is how Person Change caches its own
            # figure. Stored built: 359 KB as a ggplot, 16 KB as a grob.
            # See er2_plot_grob().
            curve_plot <- er2_plot_grob(suppressWarnings(suppressMessages(
              easyRasch2::RMreliabilityCurve(
                df,
                statistic    = self$options$curveStatistic,
                benchmark    = benchmark,
                reference    = if (isTRUE(self$options$curveReference))
                                 "marginal" else "none",
                boot         = isTRUE(self$options$curveBoot),
                boot_iter    = self$options$curveBootIter,
                conf_int     = conf_int,
                parallel     = FALSE,
                seed         = as.integer(seed),
                show_density = isTRUE(self$options$showDensity),
                n_nodes      = self$options$nNodes,
                output       = "ggplot"
              )
            )))
          }
          self$results$curveCache$setState(list(
            df = df, curve_sig = curve_sig, curve_df = curve_df, loc = loc,
            curve_plot = curve_plot))

          get <- function(nm) attr(curve_df, nm, exact = TRUE)
          ct <- self$results$curveTable
          # jamovi renders every table cell with white-space: nowrap, so a
          # cell note is one unbreakable line and a long one widens the
          # table past the results pane, taking the column headers out of
          # view with it. Cell notes are therefore kept to a few words, in
          # line with the Reliability Estimates table above, and anything
          # that needs a sentence goes into the footnotes, which do wrap.
          ct$setRow(rowKey = "thetamean", values = list(
            value = mean(loc$theta, na.rm = TRUE),
            notes = paste0("WLE, n = ", sum(!is.na(loc$theta)),
                           " of ", nrow(loc))
          ))
          ct$setRow(rowKey = "thetasd", values = list(
            value = stats::sd(loc$theta, na.rm = TRUE),
            notes = "includes measurement error"
          ))

          ct$setRow(rowKey = "latmean", values = list(
            value = get("latent_mean"),
            notes = "0 = centre of the item scale"
          ))
          ct$setRow(rowKey = "sigma", values = list(
            value = get("sigma"),
            notes = "measurement error removed"
          ))
          ct$setRow(rowKey = "margin", values = list(
            value = get("marginal_ratio"),
            notes = "as reported above"
          ))
          ct$setRow(rowKey = "green",  values = list(
            value = get("marginal_green"),
            notes = "superseded in 3.2.0"
          ))
          ct$setRow(rowKey = "sem", values = list(
            value = get("sem_average"),
            notes = "root mean error variance"
          ))
          # The information counterpart of the average SEM, and the same
          # quantity the package draws as the flat reference line on the
          # information axis. Derived from the SEM rather than averaged
          # over the grid, so the two rows describe one summary.
          ct$setRow(rowKey = "info", values = list(
            value = 1 / get("sem_average")^2,
            notes = "1 / (average SEM)²"
          ))
          if (!is.null(benchmark)) {
            # benchmark_range is a data.frame of xmin/xmax, one row per
            # qualifying stretch of the scale, and NULL when none reaches
            # the benchmark. The qualifying region can be disjoint, so all
            # rows are reported rather than just the outermost bounds.
            rng <- get("benchmark_range")
            span <- if (is.null(rng) || !is.data.frame(rng) || nrow(rng) == 0L) {
              NULL
            } else {
              paste(
                sprintf("%.2f to %.2f", rng$xmin, rng$xmax),
                collapse = ", and "
              )
            }
            ct$setRow(rowKey = "bench", values = list(
              value = get("benchmark_percent"),
              notes = if (is.null(span)) "not reached" else paste0("theta ", span)
            ))
            bench_note <- if (is.null(span)) {
              paste0("No part of the scale reaches a conditional ",
                     "reliability of ", benchmark, ".")
            } else {
              paste0("Conditional reliability reaches ", benchmark,
                     " over the range of theta shown, in logits. The ",
                     "percentage is the share of respondents located ",
                     "inside that region.")
            }
          }

          # One topic per footnote. jamovi gives each note its own row
          # under the table, so short notes read as separate lines rather
          # than as one block of prose.
          ct$setNote("context", paste0(
            "Computed by easyRasch2::RMreliabilityCurve(). The bootstrap ",
            "band, when enabled, affects the figure only; none of these ",
            "quantities depends on it."
          ))
          ct$setNote("views", paste0(
            "The first four rows are two views of the same distribution. ",
            "Mean and SD theta summarise the person estimates; the latent ",
            "mean and SD describe the distribution fitted to them without ",
            "estimating anyone's location, which is the distribution ",
            "marginal reliability is averaged over. The latent SD is the ",
            "smaller of the two SDs, measurement error having been removed. ",
            "The curve itself is a property of the items and does not ",
            "depend on either."
          ))
          ct$setNote("extremes", paste0(
            "The two means usually agree. They separate when many ",
            "respondents score at the minimum or maximum: those locations ",
            "are extrapolations that cannot run past a bound, so they pull ",
            "the mean of the estimates inward, while the latent mean is ",
            "fitted without assigning anyone a location. Its distance from ",
            "0, the origin of the item scale, is how far off target the ",
            "sample is."
          ))
          ct$setNote("marginal", paste0(
            "Marginal (curve mean) is the same quantity as the Marginal ",
            "row in the table above. Green's subtractive coefficient is ",
            "what that row reported before module version 3.2.0, shown ",
            "here so earlier results can be reconciled."
          ))
          if (!is.null(benchmark)) ct$setNote("bench", bench_note)

          n_ne <- get("n_not_estimable")
          if (!is.null(n_ne) && is.finite(n_ne) && n_ne > 0) {
            ct$setNote("notest", paste0(
              n_ne, " of the ", self$options$nNodes, " grid points had no ",
              "estimable information and are omitted from the curve."
            ))
          }
        }
        # curveCache is left untouched when the curve is off. showCurve is
        # absent from its clearWith, so the state survives and turning the
        # curve back on is free.

      }, error = function(e) {
        stop(paste("Error in reliability analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Conditional precision across the latent scale, drawn by
    # easyRasch2::RMreliabilityCurve(). The counterpart to the single
    # coefficients above: the marginal reliability in the table is the
    # latent-density-weighted mean of this curve's reliability axis.
    # Built in .run() and read from curveCache here, because jamovi calls
    # this function on every resize and export and the bootstrap band is
    # not something to pay for twice. The theme bump moved into .run() with
    # it: a built grob can no longer take a theme. Reading the cache rather
    # than a copy on this element matches .renderMap() in Person Fit and
    # .changePlot() in Person Change, and keeps one grob in state instead
    # of two. Every option that moves the figure is in curvePlot's
    # clearWith, so the redraw still happens without state of its own.
    # ---------------------------------------------------------------------
    .curvePlot = function(image, ggtheme, theme, ...) {
      cached <- self$results$curveCache$state
      if (is.null(cached) || is.null(cached$curve_plot)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      er2_draw_grob(cached$curve_plot)
    }
  )
)
