#' @export
residualpcaClass <- R6::R6Class(
  "residualpcaClass",
  inherit = residualpcaBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- pre-create the eigenvalue rows: the row count is
    # min(nComponents, number of items), fully determined by options.
    # ---------------------------------------------------------------------
    .init = function() {
      vars <- self$options$vars
      if (is.null(vars) || length(vars) < 3)
        return()
      k_show <- min(self$options$nComponents, length(vars))
      table <- self$results$pcaTable
      for (i in seq_len(k_show)) {
        table$addRow(rowKey = i, values = list(
          component  = paste0("PC", i),
          eigenvalue = NA_real_,
          propVar    = NA_real_,
          cutoff     = NA_real_,
          flagged    = ""
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met. With 2 items
      # the residual PCA has a single possible contrast (the two items
      # against each other), which carries no information about
      # multidimensionality.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$pcaNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only ",
          "2 items the residual PCA has a single possible contrast (the ",
          "two items against each other), which carries no information ",
          "about multidimensionality. Select at least 3 items.</p>"
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
        self$results$pcaTable$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$pcaTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$pcaTable$setNote("duplicate", dup_msg)

      # Complete-case counts for the notes (easyRasch2 drops incomplete
      # rows internally -- prcomp does not handle NA)
      n_total     <- nrow(df)
      n_complete  <- sum(complete.cases(df))
      n_excluded  <- n_total - n_complete

      if (n_complete == 0L)
        stop("No complete cases found. PCA on residuals requires at least one row with responses to all selected items.")

      # 3. Read options
      n_components   <- self$options$nComponents
      coord_flip     <- isTRUE(self$options$coordFlip)
      compute_cutoff <- isTRUE(self$options$computeCutoff)

      # rgl workaround
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      # 4. Run analysis via easyRasch2 (WLE-based standardized residuals
      # + unrotated PCA; results are numerically identical to
      # RMdimResidualPCA() / RMdimResidualPCACutoff() with the same seed
      # and iterations). Package warnings are suppressed -- the module
      # surfaces its own footnotes.
      tryCatch({

        # 5. Optional simulation-based cutoff. Seed is always applied
        # (default 42) so results are reproducible by default. If the
        # simulation cannot deliver a reliable cutoff, degrade
        # gracefully: eigenvalues are shown without the cutoff/flag and
        # the note explains why.
        cutoff_res   <- NULL
        sim_fail_msg <- NULL

        # The simulation is the expensive part, and jamovi reruns .run on
        # every option change -- including changes (nComponents, coordFlip)
        # that do not affect the simulation. The hidden simCache element
        # carries the cutoff object; jamovi clears its state exactly when
        # a simulation-relevant option changes (its clearWith list), and
        # the signature check makes reuse self-validating rather than
        # relying on clearWith alone.
        sim_sig <- list(
          iterations = self$options$iterations,
          seed       = as.integer(self$options$seed)
        )
        cached <- self$results$simCache$state
        if (compute_cutoff &&
            !is.null(cached) && !is.null(cached$cutoff_res) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df)) {
          cutoff_res <- cached$cutoff_res
        } else if (compute_cutoff) {
          cutoff_res <- tryCatch(
            suppressWarnings(suppressMessages(
              easyRasch2::RMdimResidualPCACutoff(
                df,
                iterations = self$options$iterations,
                parallel   = FALSE,
                seed       = as.integer(self$options$seed)
              )
            )),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )
          # Guard against a degenerate cutoff: with very few successful
          # iterations the 99th percentile collapses onto a handful of
          # values.
          if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
            sim_fail_msg <- paste0(
              "Only ", cutoff_res$actual_iterations, " of ",
              self$options$iterations, " simulation iterations succeeded ",
              "-- too few to estimate a reliable cutoff. Check your data: ",
              "items must have sufficient response variation, and the ",
              "sample must be large enough for stable Rasch model ",
              "estimation."
            )
            cutoff_res <- NULL
          }
        }

        # Save (or clear, on failure) the simulation cache.
        if (!is.null(cutoff_res)) {
          self$results$simCache$setState(list(
            cutoff_res = cutoff_res, sig = sim_sig, df = df
          ))
        } else if (compute_cutoff) {
          self$results$simCache$setState(NULL)
        }

        result_df <- suppressWarnings(suppressMessages(
          easyRasch2::RMdimResidualPCA(
            df,
            cutoff       = cutoff_res,
            n_components = n_components,
            output       = "dataframe"
          )
        ))
        vp <- attr(result_df, "variance_partition")

        cutoff_value <- if (!is.null(cutoff_res)) {
          as.numeric(cutoff_res$suggested_cutoff)
        } else NULL

        # 6. Populate eigenvalue table (rows created in .init(); raw
        # values so jamovi's Number format preferences apply)
        table <- self$results$pcaTable
        for (i in seq_len(nrow(result_df))) {
          vals <- list(
            component  = result_df$Component[i],
            eigenvalue = result_df$Eigenvalue[i],
            propVar    = result_df$Proportion_of_variance[i]
          )
          if (!is.null(cutoff_value)) {
            vals$cutoff  <- cutoff_value
            vals$flagged <- if (isTRUE(result_df$Flagged[i])) "TRUE" else ""
          }
          table$setRow(rowNo = i, values = vals)
        }

        table$setNote(
          "ncomplete",
          paste0(
            "Eigenvalues are unrotated; expressed as a proportion of total ",
            "unexplained (residual) variance. PCA performed on n = ",
            n_complete, " complete responses",
            if (n_excluded > 0L)
              paste0(" (", n_excluded, " of ", n_total,
                     " row(s) excluded due to missing values)")
            else "",
            "."
          )
        )

        # 7. Variance partition + cutoff note (HTML)
        partition_avail <- !is.null(vp) && is.finite(vp$pct_explained)
        var_text <- if (partition_avail) {
          paste0(
            "<p><b>Variance partition.</b> ",
            round(vp$pct_explained * 100, 1),
            "% of the total observed variance is explained by the fitted ",
            "Rasch model; ",
            round(vp$pct_unexplained * 100, 1),
            "% is unexplained (residual) and is what the PCA above ",
            "decomposes (n = ", vp$n_persons, " respondents; weighted ",
            "likelihood person estimates retain extreme scorers).</p>"
          )
        } else {
          paste0(
            "<p><b>Variance partition.</b> Unavailable -- too few persons ",
            "with finite theta estimates to compute the partition.</p>"
          )
        }

        cutoff_text <- if (!is.null(cutoff_value)) {
          paste0(
            "<p><b>Simulation-based cutoff.</b> ",
            cutoff_res$actual_iterations,
            " parametric-bootstrap datasets drawn from the fitted ",
            "unidimensional model at the same n. Suggested cutoff is the ",
            "99th percentile of the simulated first-contrast eigenvalues ",
            "(= ", round(cutoff_value, 3), ").",
            iteration_note(self$options$iterations, 250L),
            iteration_attrition_note(cutoff_res$actual_iterations,
                                     self$options$iterations), "</p>"
          )
        } else if (!is.null(sim_fail_msg)) {
          paste0(
            "<p><b>Simulation-based cutoff unavailable:</b> ",
            sim_fail_msg, "</p>"
          )
        } else {
          ""
        }

        self$results$pcaNote$setContent(paste0(var_text, cutoff_text))

        # 8. Save state for the loadings plot. The PC1 loadings and item
        # locations are not part of the dataframe output, so they are
        # taken from the data underlying the package's own loadings plot
        # (RMdimResidualPCA(output = "ggplot")); the module keeps its own
        # rendering (coloured labels + optional coord_flip).
        loadings_df <- suppressWarnings(suppressMessages(
          easyRasch2::RMdimResidualPCA(df, output = "ggplot")
        ))$data

        variance_text <- if (partition_avail) {
          paste0(
            "Total observed variance: ",
            round(vp$pct_explained * 100, 1), "% explained by measures, ",
            round(vp$pct_unexplained * 100, 1),
            "% unexplained\n(basis for PCA; n = ", vp$n_persons,
            " respondents, WLE)."
          )
        } else {
          "Variance partition unavailable."
        }

        self$results$pcaPlot$setState(list(
          loadings      = loadings_df,
          variance_text = variance_text,
          coord_flip    = coord_flip
        ))
      }, error = function(e) {
        stop(paste("Error in residual PCA:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # .pcaPlot -- loadings plot, optionally with coord_flip. The loadings
    # and item locations come from easyRasch2::RMdimResidualPCA(); the
    # rendering (coloured, repelled labels; optional flip) is
    # module-specific.
    # ---------------------------------------------------------------------
    .pcaPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state    <- image$state
      loadings <- state$loadings

      p <- ggplot2::ggplot(
        loadings,
        ggplot2::aes(x = .data$PC1, y = .data$Location, label = .data$Item)
      ) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
        ggplot2::geom_point(ggplot2::aes(color = .data$Item),
                            size = 4) +
        ggplot2::scale_x_continuous(limits = c(-1, 1))

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(color = .data$Item),
                                          size = 4.5, max.overlaps = Inf)
      } else {
        p <- p + ggplot2::geom_text(ggplot2::aes(color = .data$Item),
                                    nudge_y = 0.05, size = 4.5)
      }

      p <- p +
        ggplot2::labs(
          x       = "Loading on first residual contrast (PC1)",
          y       = "Item location (logit scale)",
          caption = er2_caption(state$variance_text)
        ) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(
          legend.position = "none"
        ) +
        er2_axis_margins() +
        er2_plot_caption() +
        ggplot2::scale_color_viridis_d()

      if (isTRUE(state$coord_flip)) {
        p <- p + ggplot2::coord_flip()
      }

      print(p)
      TRUE
    }
  )
)
