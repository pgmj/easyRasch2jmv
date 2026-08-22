#' @export
iteminfitmiClass <- R6::R6Class(
  "iteminfitmiClass",
  inherit = iteminfitmiBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met. With 2 items
      # the conditional infit is ~1 for both items by construction (no
      # degrees of freedom left for misfit once the total score is
      # conditioned on), so at least 3 items are required.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 3) {
        self$results$imputationNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. With only 2 ",
          "items the conditional infit equals 1 for both items by ",
          "construction and carries no information about item fit. ",
          "Select at least 3 items.</p>"
        ))
        return()
      }

      # 2. Required packages
      if (!requireNamespace("mice", quietly = TRUE))
        stop("Package 'mice' is required. Install it with: install.packages(\"mice\")")
      if (!requireNamespace("iarm", quietly = TRUE))
        stop("Package 'iarm' is required. Install it with: install.packages(\"iarm\")")

      # 3. Extract item data and convert to numeric
      data      <- self$data
      vars      <- self$options$vars
      aux_vars  <- self$options$auxVars
      if (is.null(aux_vars)) aux_vars <- character(0)
      aux_vars  <- setdiff(aux_vars, vars)

      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df_items <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df_items)
      if (!is.null(sparse_msg))
        self$results$infitTable$setNote("sparse", sparse_msg)

      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$infitTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df_items)
      if (!is.null(dup_msg))
        self$results$infitTable$setNote("duplicate", dup_msg)

      # 4. Detect missingness
      n_missing_total <- sum(is.na(df_items))
      if (n_missing_total == 0L) {
        stop(paste(
          "No missing values detected in the selected items.",
          "Use the standard 'Conditional Item Infit' analysis instead",
          "(no imputation needed)."
        ))
      }

      # 5. Read options
      method_choice  <- self$options$method
      m              <- self$options$m
      maxit          <- self$options$maxit
      seed           <- self$options$seed
      sort_by_infit  <- self$options$sortByInfit
      compute_cutoff <- isTRUE(self$options$computeCutoff)
      hdci_width     <- self$options$hdciWidth / 100
      sim_iterations <- self$options$iterations
      item_names     <- names(df_items)

      if (compute_cutoff && !requireNamespace("ggdist", quietly = TRUE)) {
        stop("Package 'ggdist' is required for HDCI cutoffs. Install with: install.packages(\"ggdist\")")
      }

      # 6. Build mice input (items + aux vars)
      mi_input <- df_items
      if (length(aux_vars) > 0L) {
        df_aux <- data[, aux_vars, drop = FALSE]
        mi_input <- cbind(mi_input, df_aux)
      }

      # For polr, items must be ordered factors with the observed integer levels
      if (method_choice == "polr") {
        for (col in item_names) {
          lvls <- sort(unique(stats::na.omit(mi_input[[col]])))
          mi_input[[col]] <- factor(mi_input[[col]],
                                    levels = as.character(lvls),
                                    ordered = TRUE)
        }
      }

      # 7. Per-column method vector
      methods_vec <- character(ncol(mi_input))
      names(methods_vec) <- names(mi_input)
      methods_vec[item_names] <- method_choice
      if (length(aux_vars) > 0L) {
        methods_vec[aux_vars] <- ""
      }

      # The imputation + pooling + cutoff-simulation pipeline below is the
      # expensive part, and jamovi reruns .run on every option change --
      # including the sort toggle, which does not affect any of it. The
      # hidden simCache element carries the pipeline outputs; jamovi
      # clears its state exactly when a pipeline-relevant option changes
      # (its clearWith list), and the signature + input-data check makes
      # reuse self-validating rather than relying on clearWith alone.
      sim_sig <- list(
        method         = method_choice,
        m              = m,
        maxit          = maxit,
        seed           = as.integer(seed),
        compute_cutoff = compute_cutoff,
        hdci_width     = hdci_width,
        sim_iterations = sim_iterations
      )
      cached <- self$results$simCache$state
      use_cache <- !is.null(cached) && !is.null(cached$results) &&
        identical(cached$sig, sim_sig) &&
        identical(cached$mi_input, mi_input)

      retry_note <- NULL
      mids_object <- NULL
      if (use_cache) {
        retry_note <- cached$retry_note
      } else {
        # 8. Run mice with retry-on-failure logic
        attempt <- private$.runMice(mi_input, methods_vec, m, maxit, seed)

        if (!attempt$ok && method_choice == "polr") {
          retry_maxit <- min(maxit * 2L, 100L)
          attempt <- private$.runMice(mi_input, methods_vec, m, retry_maxit, seed)
          if (attempt$ok) {
            retry_note <- paste0(
              "Initial imputation with maxit=", maxit, " failed to fully ",
              "impute the data. Succeeded after retrying with maxit=",
              retry_maxit, ". Consider increasing the iterations setting."
            )
          }
        }

        if (!attempt$ok) {
          msg <- paste0(
            "Imputation with method '", method_choice, "' failed",
            if (method_choice == "polr") " even after doubling maxit" else "",
            ". ",
            if (!is.null(attempt$message))
              paste0("Reason: ", attempt$message, ". ") else "",
            if (method_choice == "polr") {
              paste0("Try selecting method 'pmm' (recommended for difficult ",
                     "data) or 'cart'. Alternatively, check that all items ",
                     "have responses in all categories.")
            } else if (method_choice == "cart") {
              "Try selecting method 'pmm' or 'polr', or increase the iterations."
            } else {
              "Try increasing the iterations or selecting a different method."
            }
          )
          stop(msg)
        }

        mids_object <- attempt$imp
      }

      # 9. Pooled infit via easyRasch2. The imputation layer above (mice
      # with optional auxiliary variables) stays in the module -- the R
      # package takes a ready-made mids object instead of running mice
      # itself -- but the mids handed over must contain the item columns
      # only, so the auxiliary variables are stripped by round-tripping
      # through mice's long format. Results are numerically identical to
      # easyRasch2::RMitemInfitMI() / RMitemInfitCutoffMI() on the same
      # mids object and seed.
      tryCatch({
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        n_complete_first <- nrow(df_items)

        if (use_cache) {
          results      <- cached$results
          cutoff_res   <- cached$cutoff_res
          pooled_msq   <- cached$pooled_msq
          n_failed     <- cached$n_failed
          sim_fail_msg <- cached$sim_fail_msg
        } else {

        long_completed <- mice::complete(mids_object, action = "long",
                                         include = TRUE)
        mids_items <- mice::as.mids(
          long_completed[, c(".imp", ".id", item_names), drop = FALSE]
        )

        # 10. Optional: simulation-based cutoffs across imputations.
        # RMitemInfitCutoffMI() distributes the iterations over the imputed
        # datasets and stacks the simulated distributions, so the HDCI
        # bounds reflect both sampling and imputation uncertainty. If the
        # simulation cannot deliver reliable cutoffs, degrade gracefully:
        # pooled estimates are shown without the expected range, and the
        # note explains why.
        cutoff_res <- NULL
        sim_fail_msg <- NULL
        if (compute_cutoff) {
          cutoff_res <- tryCatch(
            suppressWarnings(suppressMessages(
              easyRasch2::RMitemInfitCutoffMI(
                mids_items,
                iterations = sim_iterations,
                parallel   = FALSE,
                seed       = as.integer(seed),
                hdci_width = hdci_width
              )
            )),
            error = function(e) {
              sim_fail_msg <<- e$message
              NULL
            }
          )
          # Guard against degenerate cutoffs: with very few successful
          # iterations in the stacked distribution the HDCI collapses.
          if (!is.null(cutoff_res) && cutoff_res$actual_iterations < 20L) {
            sim_fail_msg <- paste0(
              "Only ", cutoff_res$actual_iterations, " of ", sim_iterations,
              " simulation iterations succeeded across the imputed ",
              "datasets -- too few to estimate reliable cutoff intervals. ",
              "This typically happens when items have very low or very ",
              "high endorsement rates relative to the sample size."
            )
            cutoff_res <- NULL
          }
        }

        # Per-imputation CML fits + Rubin pooling. Failed imputations are
        # tolerated upstream (one warning per failure, at least 2 successes
        # required); the per-imputation warnings are counted here so the
        # note below can report them, matching the previous module
        # behaviour. Values are as reported by RMitemInfitMI(): pooled
        # infit and SE to 3 decimals, Rel. location to 2.
        n_failed <- 0L
        results <- withCallingHandlers(
          suppressMessages(
            easyRasch2::RMitemInfitMI(mids_items, cutoff = cutoff_res,
                                      output = "dataframe")
          ),
          warning = function(w) {
            if (grepl("^Model fitting failed for imputation",
                      conditionMessage(w)))
              n_failed <<- n_failed + 1L
            invokeRestart("muffleWarning")
          }
        )
        # Pooled observed infit per item, aligned to the item order, for
        # the plot overlay (extracted before any sorting below).
        pooled_msq <- results$Infit_MSQ[match(item_names, results$Item)]

        # Save the pipeline cache (results are pre-sort; the sort is
        # re-applied below on every run). sim_fail_msg is cached too so a
        # cached rerun replays the failure note instead of re-attempting
        # a simulation that would fail identically.
        self$results$simCache$setState(list(
          results      = results,
          cutoff_res   = cutoff_res,
          pooled_msq   = pooled_msq,
          n_failed     = n_failed,
          retry_note   = retry_note,
          sim_fail_msg = sim_fail_msg,
          sig          = sim_sig,
          mi_input     = mi_input
        ))

        } # end !use_cache

        m_ok <- m - n_failed

        # 11. Sort if requested
        if (isTRUE(sort_by_infit)) {
          results <- results[order(results$Infit_MSQ, decreasing = TRUE), ]
          rownames(results) <- NULL
        }

        # 12. Populate table
        table <- self$results$infitTable
        for (i in seq_len(nrow(results))) {
          vals <- list(
            item        = results$Item[i],
            infitMSQ    = results$Infit_MSQ[i],
            infitSE     = results$Infit_SE[i],
            relLocation = results$Relative_location[i]
          )
          if (!is.null(cutoff_res)) {
            vals$infitLow  <- results$Infit_low[i]
            vals$infitHigh <- results$Infit_high[i]
            vals$misfit    <- results$Flagged[i]
          }
          table$setRow(rowNo = i, values = vals)
        }

        # Footnotes: pooled-SE caveat, misfit rule, location definition
        table$setNote("se", paste0(
          "Pooled SE combines within-imputation variance (iarm's ",
          "asymptotic infit SE) and between-imputation variance via ",
          "Rubin's rules. Note that Mueller (2020) showed the asymptotic ",
          "infit SE can be unreliable -- base fit decisions on the ",
          "simulation-based expected range rather than the SE."
        ))
        table$setNote("loc", paste0(
          "Rel. location = pooled mean item (threshold) location ",
          "relative to the mean person location (weighted likelihood ",
          "estimates, WLE), in logits."
        ))
        if (!is.null(cutoff_res)) {
          table$setNote("misfit", paste0(
            "Flagged: pooled infit below the expected range = overfit ",
            "(item is more predictable than the model expects); above = ",
            "underfit (noisier than expected). Note the direction is ",
            "inverted relative to the item-restscore analyses."
          ))
        }

        # 13. Caption note
        aux_msg <- if (length(aux_vars) > 0L) {
          paste0(" Auxiliary predictor(s): ",
                 paste(aux_vars, collapse = ", "), ".")
        } else {
          ""
        }
        failed_msg <- if (n_failed > 0L) {
          paste0(" Note: ", n_failed, " of ", m,
                 " imputed datasets failed model fitting and were excluded.")
        } else {
          ""
        }
        retry_msg <- if (!is.null(retry_note)) {
          paste0(" <em>", retry_note, "</em>")
        } else {
          ""
        }
        cutoff_msg <- if (!is.null(cutoff_res)) {
          paste0(
            " Cutoff values based on ", cutoff_res$actual_iterations,
            " total simulation iterations across ",
            cutoff_res$n_imputations, " imputed datasets (",
            round(hdci_width * 100, 1), "% HDCI).",
            iteration_note(sim_iterations, 500L, corrected = TRUE),
            iteration_attrition_note(cutoff_res$actual_iterations,
                                     sim_iterations)
          )
        } else if (!is.null(sim_fail_msg)) {
          paste0(
            " <b>Simulation-based cutoffs unavailable:</b> ", sim_fail_msg
          )
        } else {
          ""
        }

        note_html <- paste0(
          "<p>Pooled MSQ values from ", m_ok,
          " successful imputation(s) (Rubin's rules), n = ",
          n_complete_first, " per imputed dataset. ",
          "Imputation method: '", method_choice, "', m = ", m,
          ", maxit = ", maxit, ", seed = ", seed, ".",
          aux_msg, failed_msg, retry_msg, cutoff_msg, "</p>"
        )
        self$results$imputationNote$setContent(note_html)

        # 14. Save plot state
        if (!is.null(cutoff_res)) {
          self$results$infitPlot$setState(list(
            results_df        = cutoff_res$results,
            item_names        = cutoff_res$item_names,
            actual_iterations = cutoff_res$actual_iterations,
            sample_n          = cutoff_res$sample_n,
            n_imputations     = cutoff_res$n_imputations,
            observed_infit    = pooled_msq,
            item_names_data   = item_names
          ))
        }

      }, error = function(e) {
        stop(paste("Error in pooled infit analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # .runMice — single mice attempt with diagnostic checks
    # ---------------------------------------------------------------------
    .runMice = function(mi_input, methods_vec, m, maxit, seed) {
      tryCatch({
        suppressWarnings({
          imp <- mice::mice(
            data             = mi_input,
            m                = m,
            maxit            = maxit,
            method           = methods_vec,
            seed             = seed,
            printFlag        = FALSE,
            remove.collinear = FALSE
          )
        })
        item_names <- names(methods_vec)[nzchar(methods_vec)]
        long <- mice::complete(imp, action = "long")
        if (any(is.na(long[, item_names, drop = FALSE]))) {
          return(list(
            ok = FALSE,
            imp = imp,
            message = "imputation left missing values in one or more items"
          ))
        }
        list(ok = TRUE, imp = imp, message = NULL)
      }, error = function(e) {
        list(ok = FALSE, imp = NULL, message = conditionMessage(e))
      })
    },

    # ---------------------------------------------------------------------
    # .infitPlot — same dot-plot style as iteminfit.b.R, with the
    # observed (= pooled MSQ) marker
    # ---------------------------------------------------------------------
    .infitPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)
      if (!requireNamespace("ggdist", quietly = TRUE)) return(FALSE)

      state <- image$state
      results_df       <- state$results_df
      item_names       <- state$item_names
      actual_iterations <- state$actual_iterations
      sample_n         <- state$sample_n
      n_imputations    <- state$n_imputations
      observed_infit   <- state$observed_infit
      item_names_data  <- state$item_names_data

      item_levels <- rev(item_names)

      lo_hi <- do.call(rbind, lapply(item_names, function(item) {
        sub <- results_df[results_df$Item == item, ]
        data.frame(
          Item            = item,
          min_infit_msq   = stats::quantile(sub$InfitMSQ, 0.001, na.rm = TRUE),
          max_infit_msq   = stats::quantile(sub$InfitMSQ, 0.999, na.rm = TRUE),
          p66lo_infit_msq = stats::quantile(sub$InfitMSQ, 0.167, na.rm = TRUE),
          p66hi_infit_msq = stats::quantile(sub$InfitMSQ, 0.833, na.rm = TRUE),
          median_infit    = stats::median(sub$InfitMSQ, na.rm = TRUE),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }))
      rownames(lo_hi) <- NULL

      observed_df <- data.frame(
        Item = item_names_data,
        observed_infit = observed_infit,
        stringsAsFactors = FALSE
      )

      infit_sim <- data.frame(
        Item  = results_df$Item,
        Value = results_df$InfitMSQ,
        stringsAsFactors = FALSE
      )
      infit_sim <- merge(infit_sim, observed_df[, c("Item", "observed_infit")],
                         by = "Item", sort = FALSE)
      infit_sim$Item <- factor(infit_sim$Item, levels = item_levels)
      lo_hi$Item_f   <- factor(lo_hi$Item, levels = item_levels)

      caption_text <- er2_caption(paste0(
        "Stacked results from ", actual_iterations,
        " simulated datasets across ", n_imputations,
        " imputations (n = ", sample_n, " per dataset).\n",
        "Orange diamonds indicate the pooled (Rubin's rules) observed infit.\n",
        "Black dots indicate median fit from simulations."
      ))

      p <- ggplot2::ggplot(infit_sim,
                           ggplot2::aes(x = .data$Value, y = .data$Item)) +
        ggdist::stat_dots(
          ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
          quantiles = actual_iterations,
          layout    = "weave",
          slab_color = NA,
          .width    = c(0.666, 0.999)
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$min_infit_msq, xend = .data$max_infit_msq,
                       y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 0.7
        ) +
        ggplot2::geom_segment(
          data = lo_hi,
          ggplot2::aes(x = .data$p66lo_infit_msq, xend = .data$p66hi_infit_msq,
                       y = .data$Item_f, yend = .data$Item_f),
          color = "black", linewidth = 1.2
        ) +
        ggplot2::geom_point(
          data = lo_hi,
          ggplot2::aes(x = .data$median_infit, y = .data$Item_f),
          size = 3.6
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = .data$observed_infit),
          color = "sienna2", shape = 18,
          position = ggplot2::position_nudge(y = -0.1),
          size = 7
        ) +
        ggplot2::labs(x = "Conditional Infit MSQ", y = "Item",
                      caption = caption_text) +
        ggplot2::scale_color_manual(
          values = scales::brewer_pal()(3)[-1],
          aesthetics = "slab_fill", guide = "none"
        ) +
        ggplot2::scale_x_continuous(minor_breaks = NULL) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
