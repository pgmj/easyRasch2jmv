#' @export
iteminfitmiClass <- R6::R6Class(
  "iteminfitmiClass",
  inherit = iteminfitmiBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2)
        return()

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

      df_items <- data[, vars, drop = FALSE]
      for (col in names(df_items)) {
        if (is.factor(df_items[[col]])) {
          df_items[[col]] <- as.numeric(as.character(df_items[[col]]))
        } else {
          df_items[[col]] <- as.numeric(df_items[[col]])
        }
      }

      # Identical-item check
      n_vars <- ncol(df_items)
      identical_pairs <- list()
      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          cor_val <- cor(df_items[[i]], df_items[[j]], use = "complete.obs")
          if (!is.na(cor_val) && cor_val == 1) {
            identical_pairs <- append(identical_pairs,
                                      list(c(names(df_items)[i], names(df_items)[j])))
          }
        }
      }
      if (length(identical_pairs) > 0) {
        pair_strings <- sapply(identical_pairs,
                               function(p) paste0("'", p[1], "' and '", p[2], "'"))
        pair_msg <- paste(pair_strings, collapse = "; ")
        if (ncol(df_items) == 2) {
          stop(paste("The two selected items are identical:", pair_msg,
                     "- please select different items."))
        } else {
          jmvcore::reject(
            "Warning: Some items appear to be identical ({pairs}). This may affect results.",
            pairs = pair_msg
          )
        }
      }

      all_na_cols <- sapply(df_items, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df_items)[all_na_cols]
        stop(paste("The following items contain no valid numeric data:",
                   paste(bad_vars, collapse = ", ")))
      }

      validate_response_data(df_items)

      for (col in names(df_items)) {
        unique_vals <- length(unique(stats::na.omit(df_items[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

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

      # 8. Run mice with retry-on-failure logic
      attempt <- private$.runMice(mi_input, methods_vec, m, maxit, seed)

      retry_note <- NULL
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

      # 9. Compute infit per imputation, pool with Rubin's rules
      tryCatch({
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        per_imp <- vector("list", m)
        n_failed <- 0L
        n_complete_first <- NULL

        for (i in seq_len(m)) {
          completed <- mice::complete(mids_object, action = i)
          completed <- completed[, item_names, drop = FALSE]
          for (col in item_names) {
            if (is.factor(completed[[col]])) {
              completed[[col]] <- as.numeric(as.character(completed[[col]]))
            } else {
              completed[[col]] <- as.numeric(completed[[col]])
            }
          }

          result_i <- tryCatch({
            data_mat <- as.matrix(completed)
            n_complete <- nrow(completed)

            if (max(data_mat, na.rm = TRUE) == 1L) {
              erm_out <- eRm::RM(completed)
              item_avg_locations <- stats::coef(erm_out, "beta") * -1
            } else {
              erm_out <- eRm::PCM(completed)
              thresh_table <- eRm::thresholds(erm_out)$threshtable[[1L]]
              if ("Location" %in% colnames(thresh_table)) {
                item_avg_locations <- thresh_table[, "Location"]
              } else {
                item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
              }
            }
            pp <- eRm::person.parameter(erm_out)
            person_avg_location <- mean(
              pp$theta.table[["Person Parameter"]], na.rm = TRUE
            )
            relative_locations <- item_avg_locations - person_avg_location

            cfit <- iarm::out_infit(erm_out)

            list(
              infit_msq         = cfit$Infit,
              infit_se          = cfit$Infit.se,
              relative_location = as.numeric(relative_locations),
              n_complete        = n_complete
            )
          }, error = function(e) {
            warning(sprintf("Model fitting failed for imputation %d: %s",
                            i, conditionMessage(e)), call. = FALSE)
            NULL
          })

          if (is.null(result_i)) {
            n_failed <- n_failed + 1L
            next
          }
          per_imp[[i]] <- result_i
          if (is.null(n_complete_first))
            n_complete_first <- result_i$n_complete
        }

        if (n_failed == m)
          stop("Model fitting failed for all ", m, " imputed datasets.")

        successful <- per_imp[!vapply(per_imp, is.null, logical(1L))]
        m_ok <- length(successful)
        if (m_ok < 2L)
          stop("Only ", m_ok, " imputation(s) succeeded. ",
               "At least 2 are required to estimate between-imputation variance.")

        # Pool with Rubin's rules
        n_items <- length(item_names)
        msq_mat <- matrix(NA_real_, nrow = n_items, ncol = m_ok)
        se_mat  <- matrix(NA_real_, nrow = n_items, ncol = m_ok)
        loc_mat <- matrix(NA_real_, nrow = n_items, ncol = m_ok)
        for (j in seq_len(m_ok)) {
          msq_mat[, j] <- successful[[j]]$infit_msq
          se_mat[, j]  <- successful[[j]]$infit_se
          loc_mat[, j] <- successful[[j]]$relative_location
        }
        pooled_msq      <- rowMeans(msq_mat)
        within_var      <- rowMeans(se_mat^2)
        between_var     <- apply(msq_mat, 1, stats::var)
        total_var       <- within_var + (1 + 1 / m_ok) * between_var
        pooled_se       <- sqrt(total_var)
        pooled_location <- rowMeans(loc_mat)

        results <- data.frame(
          Item              = item_names,
          Infit_MSQ         = round(pooled_msq, 3),
          Infit_SE          = round(pooled_se, 3),
          Relative_location = round(pooled_location, 2),
          stringsAsFactors  = FALSE,
          row.names         = NULL
        )

        # 10. Optional: simulation-based cutoffs across imputations
        cutoff_res <- NULL
        if (compute_cutoff) {
          cutoff_res <- private$.runCutoffSimMI(
            mids_object  = mids_object,
            item_names   = item_names,
            iterations   = sim_iterations,
            hdci_width   = hdci_width,
            seed         = seed
          )
          if (!is.null(cutoff_res)) {
            cutoff_df <- cutoff_res$item_cutoffs
            data_items <- results$Item
            cutoff_sub <- cutoff_df[, c("Item", "infit_low", "infit_high")]
            results <- merge(results, cutoff_sub, by = "Item", sort = FALSE)
            results <- results[match(data_items, results$Item), ]
            rownames(results) <- NULL
            results$Infit_low  <- round(results$infit_low,  3)
            results$Infit_high <- round(results$infit_high, 3)
            results$infit_low  <- NULL
            results$infit_high <- NULL
            results$Flagged <- results$Infit_MSQ < results$Infit_low |
                               results$Infit_MSQ > results$Infit_high
            results <- results[, c("Item", "Infit_MSQ", "Infit_SE",
                                   "Infit_low", "Infit_high", "Flagged",
                                   "Relative_location")]
          }
        }

        # 11. Sort if requested
        if (isTRUE(sort_by_infit)) {
          results <- results[order(results$Infit_MSQ, decreasing = TRUE), ]
          rownames(results) <- NULL
        }

        # 12. Populate table
        table <- self$results$infitTable
        for (i in seq_len(nrow(results))) {
          vals <- list(
            item     = results$Item[i],
            infitMSQ = results$Infit_MSQ[i],
            infitSE  = results$Infit_SE[i]
          )
          if (!is.null(cutoff_res)) {
            vals$infitLow  <- results$Infit_low[i]
            vals$infitHigh <- results$Infit_high[i]
            vals$flagged   <- ifelse(results$Flagged[i], "TRUE", "")
          }
          table$setRow(rowNo = i, values = vals)
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
            round(hdci_width * 100, 1), "% HDCI)."
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
    # .runCutoffSimMI — parametric bootstrap on each imputed dataset,
    # stack results, compute per-item HDCI bounds
    # ---------------------------------------------------------------------
    .runCutoffSimMI = function(mids_object, item_names, iterations,
                               hdci_width, seed) {

      m <- mids_object$m

      # Distribute iterations across imputations
      base_iter <- iterations %/% m
      remainder <- iterations %% m
      iters_per_imp <- rep(base_iter, m)
      if (remainder > 0L) {
        iters_per_imp[seq_len(remainder)] <-
          iters_per_imp[seq_len(remainder)] + 1L
      }

      # Per-imputation seeds
      set.seed(seed)
      imp_seeds <- sample.int(.Machine$integer.max, m)

      all_results          <- vector("list", m)
      actual_per_imp       <- integer(m)
      sample_n             <- NULL
      n_failed             <- 0L
      error_messages       <- character()

      for (i in seq_len(m)) {
        completed <- mice::complete(mids_object, action = i)
        completed <- completed[, item_names, drop = FALSE]
        for (col in item_names) {
          if (is.factor(completed[[col]])) {
            completed[[col]] <- as.numeric(as.character(completed[[col]]))
          } else {
            completed[[col]] <- as.numeric(completed[[col]])
          }
        }
        completed <- stats::na.omit(completed)
        if (nrow(completed) == 0L) {
          n_failed <- n_failed + 1L
          error_messages <- c(error_messages,
                              sprintf("imp %d: no complete cases", i))
          next
        }

        sim_res <- tryCatch({
          private$.runOneCutoffSim(
            data_complete = completed,
            iterations    = iters_per_imp[i],
            seed          = imp_seeds[i]
          )
        }, error = function(e) {
          error_messages <<- c(
            error_messages,
            sprintf("imp %d: %s", i, conditionMessage(e))
          )
          NULL
        })

        if (is.null(sim_res)) {
          n_failed <- n_failed + 1L
          next
        }

        sim_res$results_df$imputation <- i
        all_results[[i]]   <- sim_res$results_df
        actual_per_imp[i]  <- sim_res$actual_iterations
        if (is.null(sample_n)) sample_n <- sim_res$sample_n
      }

      if (n_failed == m) {
        # Surface the first few error messages so the failure isn't opaque
        sample_msgs <- utils::head(error_messages, 3L)
        stop(
          "Cutoff simulation failed for all ", m, " imputed datasets. ",
          "First error(s): ",
          paste(sample_msgs, collapse = " | ")
        )
      }

      stacked_df <- do.call(
        rbind,
        all_results[!vapply(all_results, is.null, logical(1L))]
      )
      rownames(stacked_df) <- NULL

      total_actual <- sum(actual_per_imp)
      n_imputations <- m - n_failed

      # Per-item HDCI cutoffs from stacked distribution
      unique_items <- unique(stacked_df$Item)
      item_cutoffs <- do.call(rbind, lapply(unique_items, function(item) {
        sub <- stacked_df[stacked_df$Item == item, ]
        infit_interval  <- ggdist::hdci(sub$InfitMSQ,  .width = hdci_width)
        outfit_interval <- ggdist::hdci(sub$OutfitMSQ, .width = hdci_width)
        data.frame(
          Item        = item,
          infit_low   = infit_interval[1L, 1L],
          infit_high  = infit_interval[1L, 2L],
          outfit_low  = outfit_interval[1L, 1L],
          outfit_high = outfit_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }))
      rownames(item_cutoffs) <- NULL

      list(
        results           = stacked_df,
        item_cutoffs      = item_cutoffs,
        actual_iterations = total_actual,
        sample_n          = sample_n,
        item_names        = item_names,
        n_imputations     = n_imputations,
        hdci_width        = hdci_width
      )
    },

    # ---------------------------------------------------------------------
    # .runOneCutoffSim — RMinfitcutoff()-style sequential simulation on a
    # single completed dataset, mirroring iteminfit.b.R$.runCutoffSim
    # ---------------------------------------------------------------------
    .runOneCutoffSim = function(data_complete, iterations, seed) {

      set.seed(seed)
      sim_seeds <- sample.int(.Machine$integer.max, iterations)

      data_mat       <- as.matrix(data_complete)
      sample_n       <- nrow(data_mat)
      is_polytomous  <- max(data_mat, na.rm = TRUE) > 1L
      item_names_vec <- colnames(data_mat)

      if (is_polytomous) {
        pcm_fit <- eRm::PCM(data_mat)
        pp <- eRm::person.parameter(pcm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        thresh_mat <- extract_item_thresholds(data_mat)
        deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
          as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
        })
        sim_data_list <- list(
          type = "polytomous", thetas = thetas, deltaslist = deltaslist,
          n_items = ncol(data_mat), sample_n = sample_n,
          item_names = item_names_vec
        )
      } else {
        rm_fit <- eRm::RM(data_mat)
        pp <- eRm::person.parameter(rm_fit)
        theta_table <- pp$theta.table[["Person Parameter"]]
        raw_scores <- rowSums(data_mat, na.rm = TRUE)
        thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
        item_params <- -rm_fit$betapar
        sim_data_list <- list(
          type = "dichotomous", thetas = thetas, item_params = item_params,
          n_items = ncol(data_mat), sample_n = sample_n,
          item_names = item_names_vec
        )
      }

      results_raw <- run_infit_sim_sequential(
        iterations, sim_seeds, sim_data_list, verbose = FALSE
      )

      ok <- vapply(results_raw, is.data.frame, logical(1L))
      successful <- results_raw[ok]
      if (length(successful) == 0L) {
        # Failed iterations come back as character strings (the error message)
        failed_msgs <- unlist(results_raw[!ok])
        sample_msg  <- if (length(failed_msgs) > 0L) {
          unique(failed_msgs)[1L]
        } else {
          "(no message captured)"
        }
        stop("all sim iterations failed; example: ", sample_msg)
      }

      actual_iterations <- length(successful)
      iter_dfs <- lapply(seq_along(successful), function(i) {
        d <- successful[[i]]
        d$iteration <- i
        d
      })
      results_df <- do.call(rbind, iter_dfs)
      rownames(results_df) <- NULL

      list(
        results_df        = results_df,
        actual_iterations = actual_iterations,
        sample_n          = sample_n
      )
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

      caption_text <- paste0(
        "Note: Stacked results from ", actual_iterations,
        " simulated datasets across ", n_imputations,
        " imputations (n = ", sample_n, " per dataset).\n",
        "Orange dots indicate the pooled (Rubin's rules) observed infit.\n",
        "Black dots indicate median fit from simulations."
      )

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
        ggplot2::scale_x_continuous(breaks = seq(0.5, 1.5, 0.1),
                                    minor_breaks = NULL) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"),
                       plot.caption  = ggplot2::element_text(size = 11))

      print(p)
      TRUE
    }
  )
)
