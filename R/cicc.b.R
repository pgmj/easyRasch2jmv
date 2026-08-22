#' @export
ciccClass <- R6::R6Class(
  "ciccClass",
  inherit = ciccBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met. RMitemICCPlot()
      # needs at least 2 items (the class intervals are formed from the
      # restscore/theta, which is degenerate with a single item).
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$ciccNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model and form class intervals. Select at least 2 ",
          "items.</p>"
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

      sparse_msg <- sparse_note(df)
      dup_msg    <- duplicate_items_note(df)
      recode_msg <- recode_note(data, vars)

      n_total    <- nrow(df)
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data. The class-interval curves require complete responses.")

      # Optional DIF variable (aligned with the rows of df; complete-case
      # handling, including joint missingness with the DIF variable, is
      # done inside easyRasch2::RMitemICCPlot()).
      dif_var  <- self$options$difVar
      dif_vec  <- NULL
      dif_note <- ""
      if (!is.null(dif_var)) {
        dif_vec <- droplevels(as.factor(data[[dif_var]]))
        if (nlevels(dif_vec) < 2L) {
          stop(paste0(
            "The DIF variable '", dif_var, "' must have at least 2 ",
            "distinct non-missing levels."
          ))
        }
        dif_note <- paste0(
          " Observed class-interval averages are drawn separately per ",
          "level of <b>", dif_var, "</b> (",
          paste(levels(dif_vec), collapse = ", "),
          "); panels are annotated with the partial-gamma DIF statistic."
        )
      }

      # Save state for the plot: the figure is drawn by
      # easyRasch2::RMitemICCPlot() inside the render function.
      self$results$ciccPlot$setState(list(df = df, dif = dif_vec))

      # Explanatory note (jamovi users have no access to man pages).
      # All grouping happens on the total-score scale (Buchardt et al.,
      # 2023: empirical means per value of the total score or per value
      # of the grouped total score).
      # This note explains the grouping *rule*, which jamovi users cannot
      # look up in a man page. It deliberately does not state how many
      # groups were formed: the requested number is not always the realised
      # one, and the figure caption drawn by easyRasch2::RMitemICCPlot()
      # reports the realised grouping. Two statements that could disagree
      # would be worse than one.
      method_note <- switch(self$options$method,
        quantile = paste0(
          "Total-score grouping: quantile-based, aiming for ",
          self$options$classIntervals, " groups with approximately equal ",
          "numbers of respondents. Where total scores tie at a group ",
          "boundary the groups either side merge, so fewer may be formed."
        ),
        width = paste0(
          "Total-score grouping: ", self$options$classIntervals,
          " equal-width intervals over the observed total-score range. An ",
          "interval that no respondent falls into is still defined but ",
          "contributes no point."
        ),
        score = paste0(
          "Total-score grouping: each observed total score forms its own ",
          "group (the number-of-intervals setting does not apply)."
        )
      )
      if (self$options$method != "score") {
        method_note <- paste0(
          method_note,
          " The figure caption reports the grouping actually used."
        )
      }
      band_note <- if (isTRUE(self$options$errorBand)) paste0(
        " The shaded band around the expected curve is the model-implied ",
        self$options$confLevel, "% interval for the observed mean at each ",
        "total score, given how many respondents sit at that score: if ",
        "the Rasch model is true, the observed means should fall inside ",
        "it, and a point outside the band indicates localized misfit. It ",
        "complements the error bars (which show the sample-based ",
        "uncertainty of the observed group means)."
      ) else ""
      self$results$ciccNote$setContent(paste0(
        "<p>Model-expected item score curves (solid line; the exact ",
        "conditional expectation given the total score) with observed ",
        "mean item scores overlaid per total-score group (points",
        if (isTRUE(self$options$ci)) paste0(
          ", with ", self$options$confLevel, "% confidence intervals"
        ) else "",
        "). Observed averages that track the expected curve indicate ",
        "good graphical item fit; systematic deviations suggest over- or ",
        "under-discrimination. ", method_note,
        " Groups with fewer than ", self$options$minN,
        " observations are not drawn.", band_note, dif_note,
        " Analysis based on the ", n_complete, " of ", n_total,
        " respondents with complete responses. Model estimated with CML ",
        "(psychotools) via the easyRasch2 R package; identical to ",
        "easyRasch2::RMitemICCPlot().",
        if (!is.null(sparse_msg)) paste0(" ", sparse_msg) else "",
        if (!is.null(dup_msg)) paste0(" ", dup_msg) else "",
        if (!is.null(recode_msg)) paste0(" ", recode_msg) else "",
        "</p>"
      ))
    },

    # ---------------------------------------------------------------------
    # Conditional ICC panels — easyRasch2::RMitemICCPlot() (per-item
    # expected-score curves + observed class-interval means, patchwork).
    # ---------------------------------------------------------------------
    .ciccPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("patchwork", quietly = TRUE)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMitemICCPlot(
          image$state$df,
          dif_var         = image$state$dif,
          method          = self$options$method,
          class_intervals = self$options$classIntervals,
          ci              = isTRUE(self$options$ci),
          error_band      = isTRUE(self$options$errorBand),
          conf_level      = self$options$confLevel / 100,
          min_n           = self$options$minN,
          output          = "patchwork"
        )
      ))
      p <- er2_bump_text(p)

      print(p)
      TRUE
    }
  )
)
