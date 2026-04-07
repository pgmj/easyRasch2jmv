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
        stop("You need at least two variables to run an analysis.")

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

      # Save state for the plot render function
      self$results$iccPlot$setState(list(
        df = df,
        theta_min = theta_min,
        theta_max = theta_max
      ))
    },

    .iccPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      state <- image$state
      df <- state$df
      theta_min <- state$theta_min
      theta_max <- state$theta_max

      tryCatch({
        mirt_model <- mirt::mirt(
          df,
          model = 1,
          itemtype = "Rasch",
          verbose = FALSE
        )

        p <- plot(mirt_model, type = "trace", as.table = TRUE,
                  theta_lim = c(theta_min, theta_max))
        print(p)
        TRUE
      }, error = function(e) {
        stop(paste("Error fitting Rasch model:", e$message))
      })
    }
  )
)
