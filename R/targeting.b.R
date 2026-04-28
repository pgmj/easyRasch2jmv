#' @export
targetingClass <- R6::R6Class(
  "targetingClass",
  inherit = targetingBase,
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

      # Robust conversion: handles factors with text labels (SPSS),
      # haven_labelled vectors, and numerics.
      df <- to_numeric_responses_df(df)

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

      # Sentinel-value sanity check (e.g., 999, 8888 unmarked as missing)
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
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable.", n = n_complete)

      for (col in names(df)) {
        unique_vals <- length(unique(stats::na.omit(df[[col]])))
        if (unique_vals < 2)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      tryCatch({
        data_mat  <- as.matrix(df)
        max_score <- max(data_mat, na.rm = TRUE)
        is_dicho  <- max_score == 1L
        item_names <- names(df)
        n_items    <- ncol(df)

        # Fit eRm model
        if (is_dicho) {
          erm_out <- eRm::RM(df)
        } else {
          erm_out <- eRm::PCM(df)
        }

        # Extract item thresholds (CML approach, no SE needed for table)
        if (is_dicho) {
          item_locs <- stats::coef(erm_out, "beta") * -1
          item_thresholds <- data.frame(
            Item      = item_names,
            Threshold = "T1",
            Location  = as.numeric(item_locs),
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

          # Build long-format data.frame
          thresh_list <- vector("list", nrow(thresh_table))
          for (i in seq_len(nrow(thresh_table))) {
            vals     <- thresh_table[i, ]
            non_na   <- !is.na(vals)
            n_thresh <- sum(non_na)
            thresh_list[[i]] <- data.frame(
              Item      = item_names[i],
              Threshold = paste0("T", seq_len(n_thresh)),
              Location  = as.numeric(vals[non_na]),
              stringsAsFactors = FALSE
            )
          }
          item_thresholds <- do.call(rbind, thresh_list)
          rownames(item_thresholds) <- NULL
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

        # Auto-expand xlim
        all_values <- c(person_theta, item_thresholds$Location)
        if (max(all_values, na.rm = TRUE) > xlim[2]) {
          xlim[2] <- ceiling(max(all_values, na.rm = TRUE))
        }
        if (min(all_values, na.rm = TRUE) < xlim[1]) {
          xlim[1] <- floor(min(all_values, na.rm = TRUE))
        }

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
          n_persons       = nrow(df)
        ))

        # Populate threshold table
        table <- self$results$thresholdTable
        for (i in seq_len(nrow(item_thresholds))) {
          table$addRow(rowKey = i, values = list(
            item      = item_thresholds$Item[i],
            threshold = item_thresholds$Threshold[i],
            location  = item_thresholds$Location[i]
          ))
        }

      }, error = function(e) {
        stop(paste("Error in targeting plot analysis:", e$message))
      })
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
        )

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
        )

      # BOTTOM PANEL: Item threshold dot plot
      item_thresholds$Item <- factor(item_thresholds$Item, levels = item_order)

      caption_text <- paste0(
        "Person location ", tolower(center_label), ": ",
        round(p_center, 2), " (", spread_label, " ",
        round(p_spread, 2), "). Item threshold location ",
        tolower(center_label), ": ",
        round(t_center, 2), " (", spread_label, " ",
        round(t_spread, 2), "). n = ", n_persons, "."
      )

      if (is_dicho) {
        p3 <- ggplot2::ggplot(
          item_thresholds,
          ggplot2::aes(x = .data$Location, y = .data$Item)
        ) +
          ggplot2::geom_point(size = 9, colour = threshold_fill, shape = 18) +
          ggplot2::coord_cartesian(xlim = xlim) +
          ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(8)) +
          ggplot2::labs(
            x       = "Location (logit scale)",
            y       = NULL,
            caption = caption_text
          ) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(
            plot.caption = ggplot2::element_text(hjust = 0),
            plot.margin  = ggplot2::margin(0, 5, 5, 5)
          )
      } else {
        p3 <- ggplot2::ggplot(
          item_thresholds,
          ggplot2::aes(
            x      = .data$Location,
            y      = .data$Item,
            colour = .data$Threshold
          )
        ) +
          ggplot2::geom_point(
            size     = 6.5, shape = 18,
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
            plot.caption    = ggplot2::element_text(hjust = 0),
            plot.margin     = ggplot2::margin(0, 5, 5, 5)
          )
      }

      print(p1 / p2 / p3 + patchwork::plot_layout(heights = c(3, 2, 5)))

      TRUE
    }
  )
)
