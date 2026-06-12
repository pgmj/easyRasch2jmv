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

      # Drop incomplete rows (prcomp does not handle NA)
      n_total     <- nrow(df)
      df_complete <- stats::na.omit(df)
      n_complete  <- nrow(df_complete)
      n_excluded  <- n_total - n_complete

      if (n_complete == 0L)
        stop("No complete cases found. PCA on residuals requires at least one row with responses to all selected items.")
      if (n_complete < 30L)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )

      for (col in names(df_complete)) {
        if (length(unique(df_complete[[col]])) < 2L)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # 3. Read options
      n_components   <- self$options$nComponents
      coord_flip     <- isTRUE(self$options$coordFlip)
      compute_cutoff <- isTRUE(self$options$computeCutoff)

      # rgl workaround for any eRm dependency chains
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      # 4. Run analysis (logic inlined from easyRasch2::RMdimResidualPCA)
      tryCatch({
        pca_res <- private$.runResidualPCA(df_complete, n_components)

        # 5. Optional simulation-based cutoff. Seed is always applied
        # (default 42) so results are reproducible by default. If the
        # simulation cannot deliver a reliable cutoff, degrade
        # gracefully: eigenvalues are shown without the cutoff/flag and
        # the note explains why.
        cutoff_value <- NULL
        cutoff_iters <- NULL
        sim_fail_msg <- NULL
        if (compute_cutoff) {
          cutoff_res <- tryCatch(
            private$.runCutoffSim(df_complete,
                                  self$options$iterations,
                                  as.integer(self$options$seed)),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )
          if (!is.null(cutoff_res)) {
            cutoff_value <- cutoff_res$suggested_cutoff
            cutoff_iters <- cutoff_res$actual_iterations
            pca_res$result_df$Flagged <-
              pca_res$result_df$Eigenvalue > cutoff_value
          }
        }

        # 6. Populate eigenvalue table (rows created in .init(); raw
        # values so jamovi's Number format preferences apply)
        table <- self$results$pcaTable
        for (i in seq_len(nrow(pca_res$result_df))) {
          vals <- list(
            component  = pca_res$result_df$Component[i],
            eigenvalue = pca_res$result_df$Eigenvalue[i],
            propVar    = pca_res$result_df$Proportion_of_variance[i]
          )
          if (!is.null(cutoff_value)) {
            vals$cutoff  <- cutoff_value
            vals$flagged <- if (isTRUE(pca_res$result_df$Flagged[i])) "TRUE" else ""
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
        vp <- pca_res$variance_partition
        var_text <- if (isTRUE(vp$available)) {
          paste0(
            "<p><b>Variance partition.</b> ",
            round(vp$pct_explained * 100, 1),
            "% of the total observed variance is explained by the fitted ",
            "Rasch model; ",
            round(vp$pct_unexplained * 100, 1),
            "% is unexplained (residual) and is what the PCA above ",
            "decomposes (n = ", vp$n_persons, " non-extreme cases).</p>"
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
            cutoff_iters,
            " parametric-bootstrap datasets drawn from the fitted ",
            "unidimensional model at the same n. Suggested cutoff is the ",
            "99th percentile of the simulated first-contrast eigenvalues ",
            "(= ", round(cutoff_value, 3), ").",
            iteration_note(self$options$iterations, 250L), "</p>"
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

        # 8. Save state for the loadings plot
        self$results$pcaPlot$setState(list(
          loadings        = pca_res$loadings_df,
          variance_text   = vp$text,
          coord_flip      = coord_flip
        ))
      }, error = function(e) {
        stop(paste("Error in residual PCA:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # .runResidualPCA  -- inlined from easyRasch2::RMdimResidualPCA()
    # ---------------------------------------------------------------------
    .runResidualPCA = function(df, n_components) {

      data_mat      <- as.matrix(df)
      is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

      if (is_polytomous) {
        erm_fit      <- eRm::PCM(df)
        thresh_table <- eRm::thresholds(erm_fit)$threshtable[[1L]]
        if ("Location" %in% colnames(thresh_table)) {
          item_locations <- thresh_table[, "Location"]
        } else {
          item_locations <- rowMeans(thresh_table, na.rm = TRUE)
        }
        # rownames are already item names for PCM, but be defensive
        names(item_locations) <- sub("^beta\\s+", "", names(item_locations))
      } else {
        erm_fit        <- eRm::RM(df)
        item_locations <- stats::coef(erm_fit, "beta") * -1
        # `coef(fit, "beta")` names items as "beta I1", "beta I2", ...
        # Strip the prefix so `item_locations[item_name]` lookups work.
        names(item_locations) <- sub("^beta\\s+", "", names(item_locations))
      }

      pp        <- eRm::person.parameter(erm_fit)
      ifit      <- eRm::itemfit(pp)
      st_resids <- ifit$st.res

      if (anyNA(st_resids)) {
        keep_rows <- stats::complete.cases(st_resids)
        st_resids <- st_resids[keep_rows, , drop = FALSE]
      }

      # --- Variance partition (Linacre-style; CML item params, MLE thetas) ---
      thetas_all    <- pp$theta.table[["Person Parameter"]]
      finite_thetas <- is.finite(thetas_all)

      if (sum(finite_thetas) >= 2L) {
        if (is_polytomous) {
          thresh_only <- if ("Location" %in% colnames(thresh_table)) {
            thresh_table[, colnames(thresh_table) != "Location", drop = FALSE]
          } else {
            thresh_table
          }
          expected_mat <- pcm_expected_scores(
            thetas_all[finite_thetas],
            as.matrix(thresh_only)
          )
        } else {
          expected_mat <- outer(
            thetas_all[finite_thetas],
            as.numeric(item_locations),
            function(t, b) stats::plogis(t - b)
          )
        }
        data_finite     <- data_mat[finite_thetas, , drop = FALSE]
        var_total       <- sum(apply(data_finite,  2L, stats::var, na.rm = TRUE))
        var_explained   <- sum(apply(expected_mat, 2L, stats::var, na.rm = TRUE))
        var_unexplained <- max(var_total - var_explained, 0)
        pct_explained   <- if (var_total > 0) var_explained   / var_total else NA_real_
        pct_unexplained <- if (var_total > 0) var_unexplained / var_total else NA_real_
        n_partition     <- sum(finite_thetas)
        partition_avail <- TRUE
      } else {
        var_total <- var_explained <- var_unexplained <- NA_real_
        pct_explained <- pct_unexplained <- NA_real_
        n_partition <- 0L
        partition_avail <- FALSE
      }

      partition_text <- if (partition_avail) {
        paste0(
          "Total observed variance: ",
          round(pct_explained * 100, 1), "% explained by measures, ",
          round(pct_unexplained * 100, 1),
          "% unexplained\n(basis for PCA; n = ", n_partition,
          " non-extreme cases)."
        )
      } else {
        "Variance partition unavailable (too few persons with finite theta MLEs)."
      }

      # --- Run unrotated PCA ----------------------------------------------
      pca_fit  <- stats::prcomp(st_resids)
      eigvals  <- pca_fit$sdev^2
      total    <- sum(eigvals)
      prop_var <- eigvals / total

      k_show <- min(as.integer(n_components), length(eigvals))
      result_df <- data.frame(
        Component              = paste0("PC", seq_len(k_show)),
        Eigenvalue             = eigvals[seq_len(k_show)],
        Proportion_of_variance = prop_var[seq_len(k_show)],
        stringsAsFactors       = FALSE,
        row.names              = NULL
      )

      # --- Loadings data.frame (PC1 loading + item location) --------------
      loadings_df <- as.data.frame(pca_fit$rotation[, 1L, drop = FALSE])
      colnames(loadings_df) <- "PC1"
      loadings_df$Item     <- rownames(loadings_df)
      loadings_df$Location <- as.numeric(item_locations[loadings_df$Item])
      rownames(loadings_df) <- NULL

      list(
        result_df          = result_df,
        loadings_df        = loadings_df,
        variance_partition = list(
          available       = partition_avail,
          total           = var_total,
          explained       = var_explained,
          unexplained     = var_unexplained,
          pct_explained   = pct_explained,
          pct_unexplained = pct_unexplained,
          n_persons       = n_partition,
          text            = partition_text
        )
      )
    },

    # ---------------------------------------------------------------------
    # .runCutoffSim  -- inlined from easyRasch2::RMdimResidualPCACutoff()
    # ---------------------------------------------------------------------
    .runCutoffSim = function(df, iterations, seed) {
      if (!is.null(seed)) set.seed(seed)

      sim_seeds <- sample.int(.Machine$integer.max, iterations)
      data_mat  <- as.matrix(df)
      sample_n  <- nrow(data_mat)
      is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
      item_names_vec <- colnames(data_mat)

      if (is_polytomous) {
        pcm_fit     <- eRm::PCM(data_mat)
        pp          <- eRm::person.parameter(pcm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores  <- rowSums(data_mat, na.rm = TRUE)
        thetas      <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        thresh_mat  <- extract_item_thresholds(data_mat)
        deltaslist  <- lapply(seq_len(nrow(thresh_mat)), function(i) {
          as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
        })
        sim_data_list <- list(
          type       = "polytomous",
          thetas     = thetas,
          deltaslist = deltaslist,
          n_items    = ncol(data_mat),
          sample_n   = sample_n,
          item_names = item_names_vec
        )
      } else {
        rm_fit      <- eRm::RM(data_mat)
        pp          <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores  <- rowSums(data_mat, na.rm = TRUE)
        thetas      <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        sim_data_list <- list(
          type        = "dichotomous",
          thetas      = thetas,
          item_params = item_params,
          n_items     = ncol(data_mat),
          sample_n    = sample_n,
          item_names  = item_names_vec
        )
      }

      results_raw <- run_pca_sim_sequential(iterations, sim_seeds,
                                            sim_data_list, verbose = FALSE)

      ok         <- vapply(results_raw, is.numeric, logical(1L))
      successful <- results_raw[ok]

      # Guard against a degenerate cutoff: with very few successful
      # iterations the 99th percentile collapses onto a handful of
      # values. Require at least 20 successes and a 50% success rate;
      # otherwise report the dominant failure reason.
      n_ok <- length(successful)
      if (n_ok < 20L || n_ok < iterations / 2) {
        fail_msgs <- unlist(results_raw[!ok])
        top_reason <- if (length(fail_msgs) > 0L) {
          names(sort(table(fail_msgs), decreasing = TRUE))[1L]
        } else NULL
        stop(paste0(
          "Only ", n_ok, " of ", iterations, " simulation iterations ",
          "succeeded -- too few to estimate a reliable cutoff.",
          if (!is.null(top_reason))
            paste0(" Most common failure: ", top_reason, ".") else "",
          " Check your data: items must have sufficient response ",
          "variation, and the sample must be large enough for stable ",
          "Rasch model estimation."
        ), call. = FALSE)
      }

      actual_iterations <- length(successful)
      eig_vec <- as.numeric(unlist(successful))

      list(
        actual_iterations = actual_iterations,
        sample_n          = sample_n,
        p95               = stats::quantile(eig_vec, 0.95,  na.rm = TRUE),
        p99               = stats::quantile(eig_vec, 0.99,  na.rm = TRUE),
        p995              = stats::quantile(eig_vec, 0.995, na.rm = TRUE),
        p999              = stats::quantile(eig_vec, 0.999, na.rm = TRUE),
        max               = max(eig_vec, na.rm = TRUE),
        suggested_cutoff  = as.numeric(stats::quantile(eig_vec, 0.99,
                                                       na.rm = TRUE))
      )
    },

    # ---------------------------------------------------------------------
    # .pcaPlot -- loadings plot, optionally with coord_flip
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
