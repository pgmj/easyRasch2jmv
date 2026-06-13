#' @export
targetingClass <- R6::R6Class(
  "targetingClass",
  inherit = targetingBase,
  private = list(
    .run = function() {
      # Return early / explain if requirements not met (eRm needs at
      # least 2 items to fit a model).
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$targetingNote$setContent(paste0(
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

      n_complete <- sum(complete.cases(df))
      if (n_complete == 0)
        stop("No complete cases found in the data.")

      tryCatch({
        data_mat  <- as.matrix(df)
        max_score <- max(data_mat, na.rm = TRUE)
        is_dicho  <- max_score == 1L
        item_names <- names(df)
        n_items    <- ncol(df)

        # Estimation method selection (mirrors easyRasch2::RMtargeting):
        # CML (eRm) when all response categories have at least 3
        # observations; MML (mirt) otherwise, which is more numerically
        # stable under sparse-category conditions. Person parameters
        # always come from the eRm fit.
        sparse_items <- sparse_category_items(df, min_n = 3L)
        use_mml      <- length(sparse_items) > 0L

        # Fit eRm model (thresholds in the CML branch; person parameters
        # in both branches)
        if (is_dicho) {
          erm_out <- eRm::RM(df)
        } else {
          erm_out <- eRm::PCM(df)
        }

        if (use_mml) {
          item_thresholds <- private$.thresholdsMML(df)
        } else if (is_dicho) {
          item_locs <- stats::coef(erm_out, "beta") * -1
          item_ses  <- erm_out$se.beta
          item_thresholds <- data.frame(
            Item      = item_names,
            Threshold = "T1",
            Location  = as.numeric(item_locs),
            SE        = as.numeric(item_ses),
            stringsAsFactors = FALSE
          )
        } else {
          thresh_obj   <- eRm::thresholds(erm_out)
          thresh_table <- thresh_obj$threshtable[[1]]

          # Remove Location column if present
          if ("Location" %in% colnames(thresh_table)) {
            thresh_table <- thresh_table[, colnames(thresh_table) != "Location",
                                         drop = FALSE]
          }

          # Center thresholds (subtract grand mean)
          grand_mean   <- mean(as.numeric(as.matrix(thresh_table)), na.rm = TRUE)
          thresh_table <- thresh_table - grand_mean

          # thresh_obj$se.thresh is a flat numeric vector matching
          # thresh_obj$threshpar (one SE per threshold, item-major)
          se_vec <- thresh_obj$se.thresh

          # Build long-format data.frame with per-threshold SE
          thresh_list <- vector("list", nrow(thresh_table))
          se_idx <- 1L
          for (i in seq_len(nrow(thresh_table))) {
            vals     <- thresh_table[i, ]
            non_na   <- !is.na(vals)
            n_thresh <- sum(non_na)
            thresh_list[[i]] <- data.frame(
              Item      = item_names[i],
              Threshold = paste0("T", seq_len(n_thresh)),
              Location  = as.numeric(vals[non_na]),
              SE        = as.numeric(se_vec[se_idx:(se_idx + n_thresh - 1L)]),
              stringsAsFactors = FALSE
            )
            se_idx <- se_idx + n_thresh
          }
          item_thresholds <- do.call(rbind, thresh_list)
          rownames(item_thresholds) <- NULL
        }

        # Compute CI bounds for thresholds (Location +/- z * SE).
        # `show_ci` is FALSE if the user turned it off or if any SE is
        # NA / non-finite (which can happen for degenerate items).
        ci_level <- self$options$ciLevel / 100
        show_ci  <- isTRUE(self$options$showCi) &&
                    all(is.finite(item_thresholds$SE))
        if (show_ci) {
          z_val <- stats::qnorm(1 - (1 - ci_level) / 2)
          item_thresholds$CI_low  <- item_thresholds$Location - z_val *
                                     item_thresholds$SE
          item_thresholds$CI_high <- item_thresholds$Location + z_val *
                                     item_thresholds$SE
        }

        # Person parameters
        pp <- eRm::person.parameter(erm_out)
        person_theta <- pp$theta.table[["Person Parameter"]]
        person_theta <- person_theta[!is.na(person_theta)]

        # Options
        bins       <- self$options$bins
        xlim       <- c(self$options$xlimLow, self$options$xlimHigh)
        robust     <- self$options$robust
        sort_items <- self$options$sortItems

        # Auto-expand xlim --- include CI bounds if shown so error bars
        # don't get clipped at the plot edge
        all_values <- c(person_theta, item_thresholds$Location)
        if (show_ci) {
          all_values <- c(all_values,
                          item_thresholds$CI_low,
                          item_thresholds$CI_high)
        }
        if (max(all_values, na.rm = TRUE) > xlim[2]) {
          xlim[2] <- ceiling(max(all_values, na.rm = TRUE))
        }
        if (min(all_values, na.rm = TRUE) < xlim[1]) {
          xlim[1] <- floor(min(all_values, na.rm = TRUE))
        }

        # Rows with no valid responses contribute nothing to estimation
        n_total <- sum(rowSums(!is.na(df)) > 0)

        # Save state for plot
        self$results$targetingPlot$setState(list(
          person_theta    = person_theta,
          item_thresholds = item_thresholds,
          bins            = bins,
          xlim            = xlim,
          robust          = robust,
          sort_items      = sort_items,
          item_names      = item_names,
          n_items         = n_items,
          is_dicho        = is_dicho,
          show_ci         = show_ci,
          ci_level        = ci_level,
          n_persons       = n_total
        ))

        # Populate threshold table.
        # NOTE: rows are added here rather than in .init() because the
        # number of rows equals the number of estimated thresholds, which
        # depends on the observed response categories per item -- it
        # cannot be derived from the options alone (defensible Level 3
        # case).
        table <- self$results$thresholdTable
        for (i in seq_len(nrow(item_thresholds))) {
          table$addRow(rowKey = i, values = list(
            item      = item_thresholds$Item[i],
            threshold = item_thresholds$Threshold[i],
            location  = item_thresholds$Location[i],
            se        = item_thresholds$SE[i]
          ))
        }

        # Estimation-method / sample-size note
        method_clause <- if (use_mml) {
          paste0(
            "Item thresholds were estimated with <b>MML (mirt)</b> ",
            "because item(s) ", paste(sparse_items, collapse = ", "),
            " have at least one response category with fewer than 3 ",
            "observations; MML is more numerically stable under sparse ",
            "categories. Person locations come from the eRm (CML) model. ",
            "This mirrors the estimation-method selection in ",
            "easyRasch2::RMtargeting()."
          )
        } else {
          paste0(
            "Item thresholds and person locations estimated with CML ",
            "(eRm); all response categories have at least 3 observations."
          )
        }
        self$results$targetingNote$setContent(paste0(
          "<p>Analysis based on N = ", n_total, " respondents (rows with ",
          "partially missing responses are retained by the estimation). ",
          method_clause, "</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in targeting plot analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Item thresholds via MML (mirt), used when response categories are
    # sparse. Ports easyRasch2:::.estimate_thresholds_mml(): IRT-pars
    # b-values + SEs from mirt, centered to grand mean zero so they are
    # comparable with the CML branch.
    # ---------------------------------------------------------------------
    .thresholdsMML = function(df) {
      mirt_out <- mirt::mirt(df, model = 1, itemtype = "Rasch",
                             SE = TRUE, verbose = FALSE)
      item_coefs <- mirt::coef(mirt_out, IRTpars = TRUE, printSE = TRUE,
                               simplify = FALSE)
      n_items <- ncol(df)
      thresh_list <- vector("list", n_items)
      for (i in seq_len(n_items)) {
        item_mat <- item_coefs[[i]]
        b_cols <- grep("^b\\d*$", colnames(item_mat), value = TRUE)
        if (length(b_cols) == 0L) b_cols <- "b"
        thresh_list[[i]] <- data.frame(
          Item      = names(df)[i],
          Threshold = paste0("T", seq_along(b_cols)),
          Location  = as.numeric(item_mat["par", b_cols]),
          SE        = as.numeric(item_mat["SE",  b_cols]),
          stringsAsFactors = FALSE
        )
      }
      thresh_df <- do.call(rbind, thresh_list)
      rownames(thresh_df) <- NULL

      # Center thresholds (SEs are invariant to location shift)
      thresh_df$Location <- thresh_df$Location -
        mean(thresh_df$Location, na.rm = TRUE)
      thresh_df
    },

    .targetingPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required but is not installed.")
      }

      state           <- image$state
      person_theta    <- state$person_theta
      item_thresholds <- state$item_thresholds
      bins            <- state$bins
      xlim            <- state$xlim
      robust          <- state$robust
      sort_items      <- state$sort_items
      item_names      <- state$item_names
      n_items         <- state$n_items
      is_dicho        <- state$is_dicho
      n_persons       <- state$n_persons
      show_ci         <- isTRUE(state$show_ci)
      ci_level        <- state$ci_level

      person_fill    <- "#0072B2"
      threshold_fill <- "#D55E00"

      # Summary statistics
      if (robust) {
        p_center <- stats::median(person_theta, na.rm = TRUE)
        p_spread <- stats::mad(person_theta, na.rm = TRUE)
        t_center <- stats::median(item_thresholds$Location, na.rm = TRUE)
        t_spread <- stats::mad(item_thresholds$Location, na.rm = TRUE)
        center_label <- "Median"
        spread_label <- "MAD"
      } else {
        p_center <- mean(person_theta, na.rm = TRUE)
        p_spread <- stats::sd(person_theta, na.rm = TRUE)
        t_center <- mean(item_thresholds$Location, na.rm = TRUE)
        t_spread <- stats::sd(item_thresholds$Location, na.rm = TRUE)
        center_label <- "Mean"
        spread_label <- "SD"
      }

      # Item ordering for bottom panel
      if (sort_items == "location") {
        item_means <- stats::aggregate(
          Location ~ Item, data = item_thresholds, FUN = mean, na.rm = TRUE
        )
        item_order <- rev(item_means$Item[order(item_means$Location)])
      } else {
        item_order <- rev(item_names)
      }

      # TOP PANEL: Person histogram
      person_df <- data.frame(theta = person_theta)

      p1 <- ggplot2::ggplot(person_df, ggplot2::aes(x = .data$theta)) +
        ggplot2::geom_histogram(
          bins = bins, fill = person_fill, colour = "white", alpha = 0.85
        ) +
        ggplot2::annotate(
          "rect",
          xmin = p_center - p_spread, xmax = p_center + p_spread,
          ymin = -Inf, ymax = Inf,
          fill = person_fill, alpha = 0.12
        ) +
        ggplot2::geom_vline(
          xintercept = p_center, linewidth = 0.8,
          linetype = "dashed", colour = "grey20"
        ) +
        ggplot2::annotate(
          "text",
          x = p_center, y = Inf, vjust = -0.5,
          label = paste0(
            center_label, " = ", round(p_center, 2),
            ", ", spread_label, " = ", round(p_spread, 2)
          ),
          size = 3.5, colour = "grey20"
        ) +
        ggplot2::scale_y_continuous(
          breaks = function(lim) {
            seq(0, floor(lim[2]), by = max(1, round(lim[2] / 6)))
          }
        ) +
        ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
        ggplot2::labs(x = NULL, y = "Persons") +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          plot.margin  = ggplot2::margin(5, 5, 0, 5)
        ) +
        er2_axis_margins()

      # MIDDLE PANEL: Inverted threshold histogram
      thresh_hist_df <- data.frame(location = item_thresholds$Location)

      p2 <- ggplot2::ggplot(thresh_hist_df, ggplot2::aes(x = .data$location)) +
        ggplot2::geom_histogram(
          bins = bins, fill = threshold_fill, colour = "white", alpha = 0.85
        ) +
        ggplot2::annotate(
          "rect",
          xmin = t_center - t_spread, xmax = t_center + t_spread,
          ymin = -Inf, ymax = Inf,
          fill = threshold_fill, alpha = 0.12
        ) +
        ggplot2::geom_vline(
          xintercept = t_center, linewidth = 0.8,
          linetype = "dashed", colour = "grey20"
        ) +
        ggplot2::annotate(
          "text",
          x = t_center, y = -Inf, vjust = 1.5,
          label = paste0(
            center_label, " = ", round(t_center, 2),
            ", ", spread_label, " = ", round(t_spread, 2)
          ),
          size = 4, colour = "grey20"
        ) +
        ggplot2::scale_y_reverse(
          minor_breaks = NULL,
          breaks = function(lim) {
            max_val <- abs(floor(lim[1]))
            seq(0, max_val, by = max(1, round(max_val / 4)))
          }
        ) +
        ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
        ggplot2::labs(x = NULL, y = "Thresholds") +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          plot.margin  = ggplot2::margin(0, 5, 0, 5)
        ) +
        er2_axis_margins()

      # BOTTOM PANEL: Item threshold dot plot
      item_thresholds$Item <- factor(item_thresholds$Item, levels = item_order)

      ci_caption <- if (show_ci) {
        paste0(" Horizontal error bars show ", round(ci_level * 100),
               "% CIs around item threshold locations.")
      } else {
        ""
      }
      caption_text <- er2_caption(paste0(
        "Person location ", tolower(center_label), ": ",
        round(p_center, 2), " (", spread_label, " ",
        round(p_spread, 2), "). Item threshold location ",
        tolower(center_label), ": ",
        round(t_center, 2), " (", spread_label, " ",
        round(t_spread, 2), "). n = ", n_persons, ".",
        ci_caption
      ))

      if (is_dicho) {
        p3 <- ggplot2::ggplot(
          item_thresholds,
          ggplot2::aes(x = .data$Location, y = .data$Item)
        )

        # CI error bars drawn first so the diamond points sit on top
        if (show_ci) {
          p3 <- p3 +
            ggplot2::geom_errorbar(
              ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
              width = 0.25, linewidth = 0.5, colour = threshold_fill,
              orientation = "y"
            )
        }

        p3 <- p3 +
          ggplot2::geom_point(size = 6, colour = threshold_fill, shape = 18) +
          ggplot2::coord_cartesian(xlim = xlim) +
          ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(8)) +
          ggplot2::labs(
            x       = "Location (logit scale)",
            y       = NULL,
            caption = caption_text
          ) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(
            plot.margin = ggplot2::margin(0, 5, 5, 5)
          ) +
          er2_axis_margins() +
          er2_plot_caption()
      } else {
        p3 <- ggplot2::ggplot(
          item_thresholds,
          ggplot2::aes(
            x      = .data$Location,
            y      = .data$Item,
            colour = .data$Threshold
          )
        )

        # CI error bars drawn first, dodged by Threshold so bars from
        # different thresholds for the same item don't overlap. Matches
        # the dodge width used by geom_point() below.
        if (show_ci) {
          p3 <- p3 +
            ggplot2::geom_errorbar(
              ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
              width = 0.25, linewidth = 0.5,
              position = ggplot2::position_dodge(width = 0.4),
              orientation = "y"
            )
        }

        p3 <- p3 +
          ggplot2::geom_point(
            size     = 6, shape = 18,
            position = ggplot2::position_dodge(width = 0.4)
          ) +
          ggplot2::scale_colour_viridis_d(end = 0.9) +
          ggplot2::coord_cartesian(xlim = xlim) +
          ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(8)) +
          ggplot2::labs(
            x       = "Location (logit scale)",
            y       = NULL,
            colour  = "Threshold",
            caption = caption_text
          ) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(
            legend.position = "bottom",
            plot.margin     = ggplot2::margin(0, 5, 5, 5)
          ) +
          er2_axis_margins() +
          er2_plot_caption()
      }

      print(p1 / p2 / p3 + patchwork::plot_layout(heights = c(3, 2, 5)))

      TRUE
    }
  )
)
