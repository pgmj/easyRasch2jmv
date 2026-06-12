#' @export
iccplotClass <- R6::R6Class(
  "iccplotClass",
  inherit = iccplotBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met (eRm CML needs at
      # least 2 items).
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

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")
      if (n_complete < 30)
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)

      # Get theta limits
      theta_min <- self$options$thetaMin
      theta_max <- self$options$thetaMax

      if (theta_min >= theta_max) {
        stop("Theta minimum must be less than Theta maximum.")
      }

      # Fit model and compute model-implied category probability curves.
      # Ports easyRasch2::RMitemCatProb(): RM for dichotomous data, PCM
      # for polytomous (the previous version always used PCM, which also
      # silently disabled eRm's empirical-ICC overlay).
      tryCatch({
        item_names    <- names(df)
        max_score     <- max(as.matrix(df), na.rm = TRUE)
        is_polytomous <- max_score > 1L

        if (is_polytomous) {
          fit <- eRm::PCM(df)
          thresh_mat <- eRm::thresholds(fit)$threshtable[[1L]]
          if ("Location" %in% colnames(thresh_mat)) {
            thresh_mat <- thresh_mat[, -1L, drop = FALSE]
          }
          item_thresholds <- lapply(seq_len(nrow(thresh_mat)), function(i) {
            v <- as.numeric(thresh_mat[i, ])
            v[!is.na(v)]
          })
          names(item_thresholds) <- rownames(thresh_mat)
        } else {
          fit <- eRm::RM(df)
          # Item difficulty delta_i = -beta_i (eRm parametrises beta as
          # easiness)
          deltas <- as.numeric(-fit$betapar)
          item_thresholds <- lapply(deltas, function(d) d)
          names(item_thresholds) <- item_names
        }
        item_thresholds <- item_thresholds[item_names]

        # Category probability curves on the theta grid (PCM formulation;
        # reduces to the Rasch ICC pair for dichotomous items)
        theta <- seq(theta_min, theta_max, length.out = 200L)
        pcm_probs <- function(thresholds, theta) {
          K <- length(thresholds)
          numerators <- matrix(0, nrow = length(theta), ncol = K + 1L)
          for (k in seq_len(K)) {
            numerators[, k + 1L] <- numerators[, k] + (theta - thresholds[k])
          }
          exp_num <- exp(numerators)
          exp_num / rowSums(exp_num)
        }
        per_item_dfs <- lapply(seq_along(item_thresholds), function(i) {
          thr   <- item_thresholds[[i]]
          probs <- pcm_probs(thr, theta)
          K_i   <- length(thr)
          data.frame(
            Item        = rep(item_names[i], times = length(theta) * (K_i + 1L)),
            Category    = rep(0L:K_i, each = length(theta)),
            Theta       = rep(theta, times = K_i + 1L),
            Probability = as.numeric(probs),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        })
        plot_df <- do.call(rbind, per_item_dfs)
        plot_df$Item <- factor(plot_df$Item, levels = item_names)

        n_total <- sum(rowSums(!is.na(df)) > 0)

        self$results$iccPlot$setState(list(
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
          " estimated with CML (eRm) on N = ", n_total, " respondents ",
          "(rows with partially missing responses are retained).",
          if (!is.null(sparse_msg)) paste0(" ", sparse_msg) else "",
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
      plot_df     <- state$plot_df
      theta_range <- state$theta_range
      max_score   <- state$max_score
      is_dicho    <- max_score == 1L

      caption_text <- er2_caption(paste0(
        "Model-implied ",
        if (is_dicho) "item characteristic curves"
        else "category probabilities",
        " (", state$model_label, ", CML). n = ", state$n_total, "."
      ))

      shared_scales <- list(
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = c(0, 0.25, 0.5, 0.75, 1)
        ),
        ggplot2::scale_x_continuous(
          limits = theta_range,
          breaks = seq(ceiling(theta_range[1L]),
                       floor(theta_range[2L]), by = 1)
        ),
        er2_axis_margins(),
        er2_plot_caption()
      )

      if (is_dicho) {
        # Joint ICC design (cf. eRm::plotjointICC): all items in one
        # panel, one curve per item showing P(X = 1). The faceted
        # category view is redundant for dichotomous items -- each panel
        # would show two mirror-image curves.
        d1 <- plot_df[plot_df$Category == 1L, , drop = FALSE]

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
          shared_scales
      } else {
        cat_breaks <- 0L:max_score

        p <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(
            x     = .data$Theta,
            y     = .data$Probability,
            group = .data$Category,
            color = .data$Category
          )
        ) +
          ggplot2::geom_line(linewidth = 0.9) +
          ggplot2::scale_color_viridis_c(
            name   = "Response\ncategory",
            breaks = cat_breaks,
            labels = as.character(cat_breaks),
            limits = c(0, max_score),
            end    = 0.95,
            guide  = ggplot2::guide_legend(
              override.aes = list(linewidth = 2)
            )
          ) +
          ggplot2::facet_wrap(~ Item) +
          ggplot2::labs(
            x = expression(paste("Latent trait ", theta, " (logits)")),
            y = "Category probability",
            caption = caption_text
          ) +
          ggplot2::theme_bw(base_size = 13) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "grey95",
                                                     colour = NA),
            strip.text       = ggplot2::element_text(face = "bold")
          ) +
          shared_scales
      }

      if (!isTRUE(state$show_legend)) {
        p <- p + ggplot2::theme(legend.position = "none")
      }

      print(p)
      TRUE
    }
  )
)
