#' @export
personparamsClass <- R6::R6Class(
  "personparamsClass",
  inherit = personparamsBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- pre-create the fixed summary rows so the table renders its
    # structure immediately.
    # ---------------------------------------------------------------------
    .init = function() {
      st <- self$results$summaryTable
      rows <- list(
        n      = "Respondents used",
        mean   = "Mean theta",
        sd     = "SD theta",
        median = "Median theta",
        mad    = "MAD theta",
        iqr    = "IQR theta (25th-75th percentile)",
        min    = "Min theta",
        max    = "Max theta",
        sem    = "Mean SEM",
        exmin  = "Respondents at minimum score",
        exmax  = "Respondents at maximum score"
      )
      for (key in names(rows)) {
        st$addRow(rowKey = key,
                  values = list(statistic = rows[[key]], value = NA_real_))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {
      # Person locations are estimable from 2+ items (the CML item fit
      # needs at least 2); the residual-based analyses' 3-item floor does
      # not apply here.
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$personNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to estimate ",
          "the item parameters that person locations are based on. ",
          "Select at least 2 items.</p>"
        ))
        return()
      }

      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$summaryTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$summaryTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item cannot be
      # located; RMpersonParameters() drops them internally, so drop them
      # up front and keep the row mapping for the output variables.
      n_total <- nrow(df)
      keep    <- rowSums(!is.na(df)) > 0
      df_used <- df[keep, , drop = FALSE]
      n_used  <- nrow(df_used)
      if (n_used == 0)
        stop("No respondents with any responses on the selected items.")

      # Sparse response categories destabilise CML thresholds; switch to
      # MML with a note, mirroring the targeting analysis (and upstream
      # RMtargeting()).
      sparse_items <- sparse_category_items(df_used, min_n = 3L)
      use_mml      <- length(sparse_items) > 0L
      estimator    <- if (use_mml) "MML" else "CML"

      method    <- self$options$method
      theta_min <- self$options$thetaMin
      theta_max <- self$options$thetaMax
      if (theta_min >= theta_max)
        stop("Theta lower bound must be less than the upper bound.")

      tryCatch({
        # All computation is delegated to the easyRasch2 package, so
        # results are numerically identical to RMpersonParameters() with
        # the same method/estimator. Package messages (the all-NA drop
        # notice) are suppressed -- the module reports its own n clause.
        res <- suppressWarnings(suppressMessages(
          easyRasch2::RMpersonParameters(
            df_used,
            method      = method,
            estimator   = estimator,
            theta_range = c(theta_min, theta_max),
            output      = "dataframe"
          )
        ))
        if (nrow(res) != n_used)
          stop("Internal error: person-parameter rows do not match the data.")
        prior <- attr(res, "prior") # named numeric: mean, sd, estimated

        # --- Summary table --------------------------------------------------
        n_exmin <- sum(res$extreme & res$sum_score == 0)
        n_exmax <- sum(res$extreme & res$sum_score > 0)
        st <- self$results$summaryTable
        st$setRow(rowKey = "n",     values = list(value = n_used))
        st$setRow(rowKey = "mean",   values = list(value = mean(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "sd",     values = list(value = stats::sd(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "median", values = list(value = stats::median(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "mad",    values = list(value = stats::mad(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "iqr",    values = list(value = stats::IQR(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "min",    values = list(value = min(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "max",    values = list(value = max(res$theta, na.rm = TRUE)))
        st$setRow(rowKey = "sem",    values = list(value = mean(res$sem, na.rm = TRUE)))
        st$setRow(rowKey = "exmin",  values = list(value = n_exmin))
        st$setRow(rowKey = "exmax",  values = list(value = n_exmax))

        if (use_mml) {
          st$setNote("sparse", paste0(
            "Item parameters estimated by MML instead of CML because ",
            "item(s) ", paste(sparse_items, collapse = ", "), " have ",
            "response categories with fewer than 3 observations; MML is ",
            "more numerically stable under sparse categories."
          ))
        }

        # --- Output variables (written into the dataset) --------------------
        # Values are expanded back to the full row set: respondents with
        # no responses get empty cells. Row numbers come from the data's
        # row names so the columns stay aligned with the spreadsheet
        # (including under row filters).
        expand <- function(v, fill = NA_real_) {
          out <- rep(fill, n_total)
          out[keep] <- v
          out
        }
        row_nums <- rownames(df)
        set_output <- function(opt_name, values) {
          if (!isTRUE(self$options[[opt_name]])) return()
          out <- self$results[[opt_name]]
          out$setRowNums(row_nums)
          out$setValues(values)
        }
        set_output("outputTheta",     expand(res$theta))
        set_output("outputSem",       expand(res$sem))
        set_output("outputSumScore",  expand(res$sum_score))
        set_output("outputNAnswered", expand(res$n_answered))
        set_output("outputExtreme",   expand(as.integer(res$extreme),
                                             fill = NA_integer_))

        # --- Note -----------------------------------------------------------
        drop_clause <- if (n_used < n_total) {
          paste0(" (", n_total - n_used, " row(s) without any responses ",
                 "on the selected items excluded)")
        } else ""
        n_complete <- sum(stats::complete.cases(df_used))
        missing_clause <- if (n_complete < n_used) {
          paste0(" ", n_used - n_complete, " respondent(s) have partially ",
                 "missing responses; their locations are estimated from ",
                 "the items they answered (see the Items answered output ",
                 "variable), with correspondingly larger SEM.")
        } else ""
        method_clause <- if (method == "WLE") {
          paste0(
            "Person locations are Warm's weighted likelihood estimates ",
            "(WLE): bias-corrected, and finite for minimum and maximum ",
            "scores -- such extreme scores are extrapolated, carry large ",
            "SEM, and are marked in the extreme-score output variable."
          )
        } else {
          paste0(
            "Person locations are EAP estimates (posterior mean under a ",
            "normal prior; SEM = posterior SD), using a N(",
            round(prior[["mean"]], 2), ", ", round(prior[["sd"]], 2),
            " SD) prior",
            if (isTRUE(prior[["estimated"]] == 1))
              paste0(" with the SD estimated from the data by marginal ",
                     "maximum likelihood")
            else "",
            ". EAP shrinks extreme scores toward the prior mean."
          )
        }
        self$results$personNote$setContent(paste0(
          "<p>Person locations (theta, logits) for n = ", n_used,
          " respondents", drop_clause, ", based on ",
          if (use_mml) "MML" else "conditional maximum likelihood (CML)",
          " item parameters treated as fixed. ", method_clause,
          missing_clause,
          " Tick the checkboxes under <i>Save to dataset</i> to add the ",
          "estimates as variables in the spreadsheet for use in other ",
          "analyses or export. Results are identical to ",
          "easyRasch2::RMpersonParameters().</p>"
        ))

        # --- Plot state -------------------------------------------------
        self$results$thetaPlot$setState(list(
          theta   = res$theta,
          n_used  = n_used,
          n_exmin = n_exmin,
          n_exmax = n_exmax,
          method  = method
        ))

        # --- Optional sum-score-to-logit lookup (absorbed from the former
        # Sum Score to Logit Transformation analysis). One row per
        # possible raw sum score, a function of the item parameters only;
        # identical to easyRasch2::RMscoreSE() with the same method and
        # theta range.
        if (isTRUE(self$options$showScoreTable)) {
          score_table <- tryCatch({
            suppressWarnings(suppressMessages(
              easyRasch2::RMscoreSE(
                df_used,
                method      = method,
                output      = "dataframe",
                theta_range = c(theta_min, theta_max)
              )
            ))
          }, error = function(e) {
            hint <- if (grepl("degrees of freedom", e$message, fixed = TRUE)) {
              paste0(" With very few items the MML model cannot be ",
                     "estimated; use the WLE method or select more items.")
            } else ""
            stop(paste0("Error computing score-to-theta table: ",
                        e$message, hint))
          })

          # Rows added here rather than in .init() because the row count
          # (max sum score + 1) depends on the observed item maxima.
          sc <- self$results$scoreTable
          for (i in seq_len(nrow(score_table))) {
            sc$addRow(rowKey = i, values = list(
              rawScore   = score_table$raw_score[i],
              logitScore = score_table$logit_score[i],
              logitSE    = score_table$logit_se[i]
            ))
          }

          score_method_note <- if (method == "WLE") {
            paste0(
              "Person locations via Warm's WLE (CML item parameters via ",
              "psychotools), matching the person estimates above. The ",
              "standard error is the information-based 1 / sqrt(I(theta)) ",
              "evaluated at the estimate (as in catR / TAM); Warm's bias ",
              "correction yields finite estimates even at the lowest and ",
              "highest scores, where the SE is largest."
            )
          } else {
            paste0(
              "Score-to-theta lookup via EAPsum (MML item parameters from ",
              "mirt; standard normal prior; SEs are posterior SDs). Note ",
              "that this sum-score-based EAP uses a different engine and ",
              "prior than the pattern-based EAP person estimates above, ",
              "so values for the same sum score can differ slightly."
            )
          }
          sc$setNote("method", paste0(
            score_method_note, " Identical to easyRasch2::RMscoreSE()."
          ))
          if (use_mml && method == "WLE") {
            # RMscoreSE has no estimator argument: its WLE lookup always
            # rests on CML item parameters, while the person estimates
            # above switched to MML because of the sparse categories.
            sc$setNote("sparse", paste0(
              "Response categories with fewer than 3 observations were ",
              "detected (see the summary table note): this lookup's CML ",
              "item parameters may be unstable, and the person estimates ",
              "above use MML item parameters, so values for the same sum ",
              "score can differ."
            ))
          }

          if (isTRUE(self$options$showFigure)) {
            self$results$scorePlot$setState(list(
              score_table   = score_table,
              ci_multiplier = self$options$ciMultiplier,
              method        = method,
              n_used        = n_used
            ))
          }
        }

      }, error = function(e) {
        stop(paste("Error in person parameter estimation:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Histogram of the estimated person locations, with the number and
    # share of extreme (minimum/maximum) scores in the caption. Pure
    # presentation of the package-computed values in the plot state.
    # ---------------------------------------------------------------------
    .thetaPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      s <- image$state
      pct <- function(k) sprintf("%.1f", 100 * k / s$n_used)
      caption_text <- er2_caption(paste0(
        s$method, " estimates for n = ", s$n_used, " respondents. ",
        "Minimum score: n = ", s$n_exmin, " (", pct(s$n_exmin), "%); ",
        "maximum score: n = ", s$n_exmax, " (", pct(s$n_exmax), "%)."
      ))

      p <- ggplot2::ggplot(
        data.frame(theta = s$theta),
        ggplot2::aes(x = .data$theta)
      ) +
        ggplot2::geom_histogram(bins = 30, fill = "#78b0a0",
                                colour = "white") +
        ggplot2::labs(x = "Person location (logits)", y = "Respondents",
                      caption = caption_text) +
        ggplot2::theme_minimal(base_size = 15) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    },

    # ---------------------------------------------------------------------
    # .scorePlot -- sum-score-to-logit conversion: points with horizontal
    # CI bars (absorbed from the former Sum Score to Logit Transformation
    # analysis).
    # ---------------------------------------------------------------------
    .scorePlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state <- image$state
      d <- state$score_table
      d$lower <- d$logit_score - state$ci_multiplier * d$logit_se
      d$upper <- d$logit_score + state$ci_multiplier * d$logit_se

      # Drop boundary rows where SE is NA so the figure stays clean
      d <- d[is.finite(d$logit_score) & !is.na(d$logit_se), , drop = FALSE]

      caption_text <- er2_caption(paste0(
        if (state$method == "WLE") {
          "Warm's WLE (CML, psychotools)."
        } else {
          "EAPsum (MML, mirt)."
        },
        " Error bars: ±", state$ci_multiplier,
        " × logit SE. n = ", state$n_used, "."
      ))

      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(x = .data$logit_score, y = .data$raw_score)
      ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
          width = 0.5, colour = "darkgrey",
          orientation = "y"
        ) +
        ggplot2::geom_point(size = 3, shape = 18) +
        ggplot2::labs(
          x = "Logit interval score",
          y = "Ordinal sum score",
          caption = caption_text
        ) +
        ggplot2::theme_bw(base_size = 15) +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
