#' @export
iccplotClass <- R6::R6Class(
  "iccplotClass",
  inherit = iccplotBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met (the model needs
      # at least 2 items).
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$iccNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model. Select at least 2 items.</p>"
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

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      # Get theta limits
      theta_min <- self$options$thetaMin
      theta_max <- self$options$thetaMax

      if (theta_min >= theta_max) {
        stop("Theta minimum must be less than Theta maximum.")
      }

      tryCatch({
        max_score     <- max(as.matrix(df), na.rm = TRUE)
        is_polytomous <- max_score > 1L

        # Model-implied category probability curves via easyRasch2 (CML
        # item parameters via psychotools; a dichotomous item is a
        # 2-category PCM). The probabilities are numerically identical to
        # RMitemCatProb(); the dataframe feeds the module's joint-ICC
        # rendering for dichotomous data (a module-specific view -- the
        # package facets per item), while polytomous data are drawn by
        # the package plot directly in the render function.
        plot_df <- suppressWarnings(suppressMessages(
          easyRasch2::RMitemCatProb(
            df,
            theta_range = c(theta_min, theta_max),
            output      = "dataframe"
          )
        ))

        n_total <- sum(rowSums(!is.na(df)) > 0)

        self$results$iccPlot$setState(list(
          df          = df,
          plot_df     = plot_df,
          theta_range = c(theta_min, theta_max),
          max_score   = max_score,
          show_legend = self$options$showLegend,
          n_total     = n_total,
          model_label = if (is_polytomous) "partial credit model"
                        else "Rasch model"
        ))

        # Sample-size / estimation note (sparse warning folded in)
        self$results$iccNote$setContent(paste0(
          "<p>Model-implied ",
          if (is_polytomous) {
            "response category probabilities from a partial credit model"
          } else {
            "item characteristic curves from a Rasch model"
          },
          " estimated with CML (psychotools, via the easyRasch2 R ",
          "package) on N = ", n_total, " respondents (rows with ",
          "partially missing responses are retained). Probabilities are ",
          "identical to easyRasch2::RMitemCatProb().",
          if (!is.null(sparse_msg)) paste0(" ", sparse_msg) else "",
          if (!is.null(dup_msg)) paste0(" ", dup_msg) else "",
          if (!is.null(recode_msg)) paste0(" ", recode_msg) else "",
        if (!is.null(recode_msg)) paste0(" ", recode_msg) else "",
          "</p>"
        ))
      }, error = function(e) {
        stop(paste("Error fitting the Rasch model:", e$message))
      })
    },

    .iccPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state       <- image$state
      theta_range <- state$theta_range
      is_dicho    <- state$max_score == 1L

      if (!is_dicho) {
        # Polytomous: the faceted category-probability plot is drawn by
        # easyRasch2::RMitemCatProb() directly (its caption reports the
        # sample in the house style).
        p <- suppressWarnings(suppressMessages(
          easyRasch2::RMitemCatProb(
            state$df,
            theta_range  = theta_range,
            output       = "ggplot",
            label_curves = "legend"
          )
        ))
        p <- er2_bump_text(p)
        if (!isTRUE(state$show_legend)) {
          p <- p + ggplot2::theme(legend.position = "none")
        }
        print(p)
        return(TRUE)
      }

      # Dichotomous: joint ICC design (cf. eRm::plotjointICC) -- all items
      # in one panel, one curve per item showing P(X = 1). This view is
      # module-specific (the package facets per item; each panel would
      # show two mirror-image curves); the curves themselves come from
      # the package's probability grid.
      plot_df <- state$plot_df
      d1 <- plot_df[plot_df$Category == 1L, , drop = FALSE]

      caption_text <- er2_caption(paste0(
        "Model-implied item characteristic curves (",
        state$model_label, ", CML via psychotools). n = ",
        state$n_total, "."
      ))

      p <- ggplot2::ggplot(
        d1,
        ggplot2::aes(
          x     = .data$Theta,
          y     = .data$Probability,
          color = .data$Item
        )
      ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::scale_color_viridis_d(name = "Item", end = 0.95) +
        ggplot2::labs(
          x = expression(paste("Latent trait ", theta, " (logits)")),
          y = "P(response = 1)",
          caption = caption_text
        ) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = c(0, 0.25, 0.5, 0.75, 1)
        ) +
        ggplot2::scale_x_continuous(
          limits = theta_range,
          breaks = seq(ceiling(theta_range[1L]),
                       floor(theta_range[2L]), by = 1)
        ) +
        er2_axis_margins() +
        er2_plot_caption()

      if (!isTRUE(state$show_legend)) {
        p <- p + ggplot2::theme(legend.position = "none")
      }

      print(p)
      TRUE
    }
  )
)
