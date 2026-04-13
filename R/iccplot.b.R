#' @export
iccplotClass <- R6::R6Class(
  "iccplotClass",
  inherit = iccplotBase,
  private = list(
    .run = function() {
      # Return early if no variables selected
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2)
        return()

      # Extract and validate data
      data <- self$data
      vars <- self$options$vars
      df <- data[, vars, drop = FALSE]

      # Convert factors to numeric
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          df[[col]] <- as.numeric(as.character(df[[col]]))
        } else {
          df[[col]] <- as.numeric(df[[col]])
        }
      }

      # Strip jmvcore S4 column wrappers
      col_names <- names(df)
      df <- as.data.frame(
        matrix(as.numeric(as.matrix(df)),
               nrow = nrow(df),
               ncol = ncol(df),
               dimnames = list(NULL, col_names))
      )

      # Check for all-NA columns
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste("The following variables contain no valid numeric data:",
                   paste(bad_vars, collapse = ", ")))
      }

      validate_response_data(df)

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # Get theta limits
      theta_min <- self$options$thetaMin
      theta_max <- self$options$thetaMax

      if (theta_min >= theta_max) {
        stop("Theta minimum must be less than Theta maximum.")
      }

      # Fit PCM model in .run() with clean data
      tryCatch({
        pcm_model <- eRm::PCM(df)

        self$results$iccPlot$setState(list(
          pcm_model = pcm_model,
          theta_min = theta_min,
          theta_max = theta_max,
          n_items = length(vars),
          show_legend = self$options$showLegend
        ))
      }, error = function(e) {
        stop(paste("Error fitting Partial Credit Model:", e$message))
      })
    },

    .iccPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      state <- image$state
      pcm_model <- state$pcm_model
      theta_min <- state$theta_min
      theta_max <- state$theta_max
      n_items <- state$n_items
      show_legend <- state$show_legend

      # Set up a grid layout for all items
      n_cols <- ceiling(sqrt(n_items))
      n_rows <- ceiling(n_items / n_cols)

      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))

      legpos <- if (isTRUE(show_legend)) "topleft" else FALSE

      eRm::plotICC(pcm_model, empICC = list("raw"),
                   xlim = c(theta_min, theta_max),
                   legpos = legpos, ylab = "Probability")

      TRUE
    }
  )
)
