#' @export
scoreseClass <- R6::R6Class(
  "scoreseClass",
  inherit = scoreseBase,
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

      # 2. Extract data and convert to numeric
      data <- self$data
      vars <- self$options$vars
      df   <- data[, vars, drop = FALSE]

      df <- to_numeric_responses_df(df)

      # All-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste("The following variables contain no valid numeric data:",
                   paste(bad_vars, collapse = ", ")))
      }

      # Sentinel-value sanity check
      max_obs <- max(as.matrix(df), na.rm = TRUE)
      if (is.finite(max_obs) && max_obs > 20) {
        bad_cols <- names(df)[
          vapply(df, function(x) {
            mx <- suppressWarnings(max(x, na.rm = TRUE))
            is.finite(mx) && mx > 20
          }, logical(1L))
        ]
        stop(paste0(
          "Item(s) ", paste(bad_cols, collapse = ", "),
          " contain values > 20, which look like missing-value codes ",
          "(e.g., 999, 8888) rather than ordinal responses. ",
          "Mark these codes as missing in the data editor, or recode your data."
        ))
      }

      validate_response_data(df)

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

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
        stop(paste("Error computing score-to-theta table:", e$message))
      })

      # 5. Populate table
      table <- self$results$scoreTable
      for (i in seq_len(nrow(score_table))) {
        table$addRow(rowKey = i, values = list(
          rawScore   = score_table$raw_score[i],
          logitScore = round(score_table$logit_score[i], 3),
          logitSE    = round(score_table$logit_se[i],    3)
        ))
      }

      # Footnote: estimation method + complete-case basis
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
               " Item parameters fitted on n = ", n_complete,
               " complete responses (rows with no missing values across the ",
               "selected items).")
      )

      # 6. Caption HTML (kept short; setNote already conveys most info)
      n_total    <- nrow(df)
      n_excluded <- n_total - n_complete
      excluded_msg <- if (n_excluded > 0L) {
        paste0(" (", n_excluded, " of ", n_total,
               " row(s) excluded due to missing responses)")
      } else {
        ""
      }
      self$results$scoreNote$setContent(
        paste0("<p>n = ", n_complete, " complete responses",
               excluded_msg, ".</p>")
      )

      # 7. Save state for plot
      if (isTRUE(self$options$showFigure)) {
        self$results$scorePlot$setState(list(
          score_table   = score_table,
          ci_multiplier = ci_multiplier,
          method        = method,
          n_complete    = n_complete
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

      caption_text <- paste0(
        if (state$method == "WLE") {
          "Warm's WLE (CML, eRm)."
        } else {
          "EAPsum (MML, mirt)."
        },
        " Error bars: ±", state$ci_multiplier,
        " × logit SE. n = ", state$n_complete,
        " complete responses."
      )

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
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(plot.caption = ggplot2::element_text(size = 10))

      print(p)
      TRUE
    }
  )
)
