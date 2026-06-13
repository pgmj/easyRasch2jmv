#' @export
scoreseClass <- R6::R6Class(
  "scoreseClass",
  inherit = scoreseBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early / explain if requirements not met (eRm CML needs
      # at least 2 items)
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$scoreNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model. Select at least 2 items.</p>"
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
        self$results$scoreTable$setNote("sparse", sparse_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$scoreTable$setNote("duplicate", dup_msg)

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      # 3. Read options
      method        <- self$options$method
      theta_min     <- self$options$thetaMin
      theta_max     <- self$options$thetaMax
      ci_multiplier <- self$options$ciMultiplier

      if (theta_min >= theta_max)
        stop("Theta lower bound must be less than the upper bound.")

      # rgl workaround for any iarm/eRm dependency chains
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      # 4. Compute the score-to-theta lookup
      score_table <- tryCatch({
        if (method == "WLE") {
          private$.scoreSE_wle(df, c(theta_min, theta_max))
        } else {
          private$.scoreSE_eap(df)
        }
      }, error = function(e) {
        hint <- if (grepl("degrees of freedom", e$message, fixed = TRUE)) {
          paste0(" With very few items the MML model cannot be estimated; ",
                 "use the WLE method or select more items.")
        } else ""
        stop(paste0("Error computing score-to-theta table: ", e$message,
                    hint))
      })

      # 5. Populate table. Raw values (no pre-rounding) so the jamovi
      # frontend applies the user's "Number format" preferences.
      # NOTE: rows are added here rather than in .init() because the row
      # count (max sum score + 1) depends on the observed item maxima --
      # it cannot be derived from the options alone (defensible Level 3
      # case).
      table <- self$results$scoreTable
      for (i in seq_len(nrow(score_table))) {
        table$addRow(rowKey = i, values = list(
          rawScore   = score_table$raw_score[i],
          logitScore = score_table$logit_score[i],
          logitSE    = score_table$logit_se[i]
        ))
      }

      # Footnote: estimation method + sample basis. Both estimation
      # paths fit the model on all available responses (eRm CML / mirt
      # MML retain rows with partially missing responses).
      n_used <- sum(rowSums(!is.na(df)) > 0)
      method_text <- if (method == "WLE") {
        paste0("Person locations via Warm's WLE (CML item parameters from eRm). ",
               "Boundary scores searched within theta range [",
               theta_min, ", ", theta_max,
               "]; values outside that range are returned as Inf with NA SE.")
      } else {
        paste0("Person locations via EAPsum (MML item parameters from mirt). ",
               "Estimates are bounded by the standard normal prior; SEs are ",
               "posterior SDs.")
      }
      table$setNote(
        "method",
        paste0(method_text,
               " Item parameters fitted on N = ", n_used,
               " respondents (rows with partially missing responses are ",
               "retained by the estimation).")
      )

      # 6. Caption HTML (kept short; setNote already conveys most info)
      missing_msg <- if (n_used > n_complete) {
        paste0(" (", n_complete, " with complete responses)")
      } else {
        ""
      }
      self$results$scoreNote$setContent(
        paste0("<p>N = ", n_used, " respondents", missing_msg, ".</p>")
      )

      # 7. Save state for plot
      if (isTRUE(self$options$showFigure)) {
        self$results$scorePlot$setState(list(
          score_table   = score_table,
          ci_multiplier = ci_multiplier,
          method        = method,
          n_used        = n_used
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .scoreSE_wle  — eRm + patched iarm helpers (mirrors easyRasch2)
    # ---------------------------------------------------------------------
    .scoreSE_wle = function(df, theta_range) {
      data_mat <- as.matrix(df)
      if (max(data_mat, na.rm = TRUE) == 1L) {
        erm_out <- eRm::RM(df)
      } else {
        erm_out <- eRm::PCM(df)
      }
      score_list <- iarm_person_estimates(
        erm_out, properties = TRUE, sthetarange = theta_range
      )
      wle_mat <- score_list[[2L]]
      data.frame(
        raw_score   = as.integer(wle_mat[, "Raw Score"]),
        logit_score = as.numeric(wle_mat[, "WLE"]),
        logit_se    = as.numeric(wle_mat[, "SEM"]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    },

    # ---------------------------------------------------------------------
    # .scoreSE_eap  — mirt EAPsum
    # ---------------------------------------------------------------------
    .scoreSE_eap = function(df) {
      mirt_fit <- suppressMessages(
        mirt::mirt(
          data       = as.data.frame(df),
          model      = 1,
          itemtype   = "Rasch",
          verbose    = FALSE,
          accelerate = "squarem"
        )
      )
      sscores <- mirt::fscores(
        mirt_fit,
        method         = "EAPsum",
        full.scores    = FALSE,
        full.scores.SE = TRUE
      )
      sscores <- as.data.frame(sscores)

      raw_col   <- intersect(c("Sum.Scores", "Sum.Score"), names(sscores))[1L]
      theta_col <- intersect(c("F1", "Theta", "EAP"),       names(sscores))[1L]
      se_col    <- intersect(c("SE_F1", "SE", "SE_Theta"),  names(sscores))[1L]

      if (is.na(raw_col) || is.na(theta_col) || is.na(se_col)) {
        stop("Unexpected mirt::fscores() output structure; cannot locate ",
             "score / theta / SE columns. Columns returned: ",
             paste(names(sscores), collapse = ", "))
      }

      data.frame(
        raw_score   = as.integer(sscores[[raw_col]]),
        logit_score = as.numeric(sscores[[theta_col]]),
        logit_se    = as.numeric(sscores[[se_col]]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    },

    # ---------------------------------------------------------------------
    # .scorePlot — points with horizontal CI bars
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
          "Warm's WLE (CML, eRm)."
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
