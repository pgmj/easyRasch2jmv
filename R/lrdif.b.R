#' @export
lrdifClass <- R6::R6Class(
  "lrdifClass",
  inherit = lrdifBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run -- everything happens here (columns set up dynamically once
    # the DIF group structure is known from the data)
    # ---------------------------------------------------------------------
    .run = function() {

      # 1. Return early if requirements not met
      vars   <- self$options$vars
      difVar <- self$options$difVar
      if (is.null(vars) || length(vars) < 2 || is.null(difVar))
        return()

      # 2. Extract data + validate items
      data <- self$data
      df   <- data[, vars, drop = FALSE]
      df   <- to_numeric_responses_df(df)

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

      # 3. DIF variable + joint complete-case handling
      dif_raw    <- data[[difVar]]
      dif_factor <- droplevels(as.factor(dif_raw))

      # Drop rows where dif_var is NA, jointly with df
      na_mask <- is.na(dif_factor)
      if (any(na_mask)) {
        df         <- df[!na_mask, , drop = FALSE]
        dif_factor <- droplevels(dif_factor[!na_mask])
      }

      # Drop rows where any item is NA (eRm::LRtest needs complete cases)
      complete_mask <- stats::complete.cases(df)
      df         <- df[complete_mask, , drop = FALSE]
      dif_factor <- droplevels(dif_factor[complete_mask])

      groups <- levels(dif_factor)
      if (length(groups) < 2L)
        stop("DIF variable must have at least 2 distinct non-missing levels after dropping incomplete rows.")

      n_complete <- nrow(df)
      if (n_complete == 0L)
        stop("No complete cases remaining after dropping rows with NA in items or DIF variable.")
      if (n_complete < 30L)
        jmvcore::reject(
          "Warning: Only {n} complete cases found. Results may be unreliable.",
          n = n_complete
        )

      # Per-item variation check
      for (col in names(df)) {
        if (length(unique(df[[col]])) < 2L)
          stop(paste0("Item '", col, "' has no variation in responses."))
      }

      # Per-group variation check
      for (g in groups) {
        sub <- df[dif_factor == g, , drop = FALSE]
        bad <- vapply(sub, function(x) length(unique(x)) < 2L, logical(1L))
        if (any(bad)) {
          bad_items <- names(sub)[bad]
          stop(paste0(
            "Group '", g, "' has no variation on item(s): ",
            paste(bad_items, collapse = ", "),
            ". Inspect the response distribution per group before running ",
            "the LR test."
          ))
        }
      }

      # 4. Read options
      level         <- self$options$level
      cutoff_val    <- self$options$cutoff
      sort_by_max   <- isTRUE(self$options$sortByMaxDiff)
      show_figure   <- isTRUE(self$options$showFigure)
      conf_level    <- self$options$confLevel / 100

      # 5. Fit + LRtest
      tryCatch({
        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        data_mat      <- as.matrix(df)
        is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

        fit_full <- if (is_polytomous) eRm::PCM(df) else eRm::RM(df)

        lrt <- tryCatch(
          eRm::LRtest(fit_full, splitcr = dif_factor),
          error = function(e) {
            stop(paste0(
              "eRm::LRtest() failed: ", conditionMessage(e),
              ". This often indicates an empty response category in one ",
              "subgroup. Inspect the response distribution per group ",
              "before running the LR test."
            ))
          }
        )

        # 6. Long-form data: per-group + overall
        per_group <- private$.extractLrLocations(lrt, groups, names(df))
        overall   <- private$.oneFitLong(fit_full, "All", names(df))
        long_df   <- rbind(per_group, overall)
        long_df$Item <- factor(long_df$Item, levels = names(df))

        # 7. Aggregate to chosen level
        if (level == "item") {
          agg <- stats::aggregate(
            cbind(Location = long_df$Location, SE = long_df$SE),
            by  = list(Item = long_df$Item, DIFgroup = long_df$DIFgroup),
            FUN = function(x) mean(x, na.rm = TRUE)
          )
          agg$Item <- factor(agg$Item, levels = names(df))
          table_df <- private$.buildWideTable(agg, groups, level = "item")
          plot_df  <- agg
          plot_df$Threshold <- NA_character_
        } else {
          table_df <- private$.buildWideTable(long_df, groups, level = "threshold")
          plot_df  <- long_df
        }

        # 8. Flagged column
        if (cutoff_val > 0) {
          table_df$Flagged <- !is.na(table_df$MaxDiff) &
                              table_df$MaxDiff > cutoff_val
        } else {
          table_df$Flagged <- FALSE
        }

        # 9. Sort if requested
        if (sort_by_max) {
          table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
          rownames(table_df) <- NULL
        }

        # 10. Set up table columns now that the group structure is known.
        # Numeric columns use format = "zto" to match the formatting of
        # all other tables in the module (and avoids the per-column
        # sig-figs heuristic expanding decimals for small magnitudes).
        table <- self$results$lrtTable

        # When at threshold level, pack repeated Item values into a
        # single cell across the rows belonging to the same item.
        table$addColumn(name = "item", title = "Item", type = "text",
                        combineBelow = (level == "threshold"))
        if (level == "threshold") {
          table$addColumn(name = "threshold", title = "Threshold",
                          type = "text")
        }
        for (g in groups) {
          table$addColumn(
            name       = private$.locColName(g),
            title      = g,
            type       = "number",
            format     = "zto",
            superTitle = "Location"
          )
        }
        table$addColumn(name = "loc_overall", title = "All",
                        type = "number", format = "zto",
                        superTitle = "Location")
        table$addColumn(name = "maxDiff", title = "MaxDiff",
                        type = "number", format = "zto")
        if (cutoff_val > 0) {
          table$addColumn(name = "flagged", title = "Flagged", type = "text")
        }
        for (g in groups) {
          table$addColumn(
            name       = private$.seColName(g),
            title      = g,
            type       = "number",
            format     = "zto",
            superTitle = "SE"
          )
        }
        table$addColumn(name = "se_overall", title = "All",
                        type = "number", format = "zto",
                        superTitle = "SE")

        # 11. Populate rows. Pass raw numerics so the jamovi frontend
        # applies the user's "Number format" preferences (pt = decimal
        # places, sf = significant figures). Note: Location columns
        # may show one more decimal than SE for small-magnitude values
        # because jamovi's display heuristic adds digits to preserve
        # the chosen sig-figs count -- this is jamovi-wide behaviour
        # and matches every other analysis in the platform.
        for (i in seq_len(nrow(table_df))) {
          vals <- list(item = as.character(table_df$Item[i]))
          if (level == "threshold") {
            vals$threshold <- as.character(table_df$Threshold[i])
          }
          for (g in groups) {
            vals[[private$.locColName(g)]] <- table_df[[g]][i]
          }
          vals$loc_overall <- table_df$All[i]
          vals$maxDiff     <- table_df$MaxDiff[i]
          if (cutoff_val > 0) {
            vals$flagged <- if (isTRUE(table_df$Flagged[i])) "TRUE" else ""
          }
          for (g in groups) {
            vals[[private$.seColName(g)]] <- table_df[[paste0("SE_", g)]][i]
          }
          vals$se_overall <- table_df$SE_All[i]
          table$addRow(rowKey = i, values = vals)
        }

        # 12. LR test note (p-value rounded to 3 digits)
        model_name <- if (is_polytomous) "Partial Credit Model" else "Rasch Model"
        p_round <- round(lrt$pvalue, 3)
        p_str   <- if (p_round == 0) "&lt; 0.001" else format(p_round, nsmall = 3)
        lr_html <- paste0(
          "<p><b>Andersen LR test:</b> &chi;<sup>2</sup> = ",
          round(lrt$LR, 3),
          ", <i>df</i> = ", lrt$df,
          ", <i>p</i> = ", p_str,
          ". ", model_name, " split by `", difVar,
          "` (", length(groups), " groups: ",
          paste(groups, collapse = ", "), "). n = ", n_complete,
          " complete cases."
        )
        if (cutoff_val > 0) {
          lr_html <- paste0(
            lr_html,
            " Items flagged when MaxDiff &gt; ", cutoff_val, " logits."
          )
        }
        lr_html <- paste0(lr_html, "</p>")
        self$results$lrtNote$setContent(lr_html)

        # 13. Save state for the figure
        if (show_figure) {
          self$results$lrtPlot$setState(list(
            plot_df       = plot_df,
            groups        = groups,
            level         = level,
            conf          = conf_level,
            n_complete    = n_complete,
            is_polytomous = is_polytomous
          ))
        }

        # 14. Tileplot: per-item × category × DIF-group response counts
        if (isTRUE(self$options$showTileplot)) {
          tile_state <- private$.computeTileCounts(df, dif_factor)
          tile_state$cutoff  <- self$options$tileCutoff
          tile_state$percent <- isTRUE(self$options$tilePercent)
          self$results$tileplot$setState(tile_state)
        }
      }, error = function(e) {
        stop(paste("Error in LR-based DIF analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Internal helpers (inlined from easyRasch2::RMdifLR)
    # ---------------------------------------------------------------------

    .locColName = function(g) {
      paste0("loc_", private$.sanitizeName(g))
    },
    .seColName = function(g) {
      paste0("se_", private$.sanitizeName(g))
    },
    .sanitizeName = function(g) {
      out <- gsub("[^A-Za-z0-9]", "_", as.character(g))
      if (!nzchar(out) || grepl("^[0-9]", out)) out <- paste0("g", out)
      out
    },

    # Wrap long DIF-group labels onto multiple lines for the x-axis.
    # Base R only -- no stringr dependency. Width chosen so labels
    # like "Elementary school" wrap to two lines without crowding.
    .wrapLabels = function(x, width = 10L) {
      vapply(as.character(x), function(s) {
        if (is.na(s) || !nzchar(s)) return(s)
        paste(strwrap(s, width = width), collapse = "\n")
      }, character(1L))
    },

    # Pull per-fit thresholds + SEs into long-form
    .oneFitLong = function(fit, dif_group, item_names) {
      is_rm <- isTRUE(fit$model == "RM")
      if (is_rm) {
        # eRm::thresholds() is polytomous-only. For RM the "threshold"
        # is the item difficulty, available as -betapar with se.beta.
        beta <- fit$betapar
        se   <- fit$se.beta
        raw_names <- sub("^beta\\s+", "", names(beta))
        if (!all(raw_names %in% item_names)) {
          stop("Could not match RM item parameter names to data columns: ",
               paste(utils::head(setdiff(raw_names, item_names), 5L),
                     collapse = ", "))
        }
        return(data.frame(
          Item      = raw_names,
          Threshold = "1",
          DIFgroup  = dif_group,
          Location  = unname(-as.numeric(beta)),
          SE        = unname(as.numeric(se)),
          stringsAsFactors = FALSE
        ))
      }

      thr <- eRm::thresholds(fit)
      loc <- thr$threshpar
      se  <- thr$se.thresh
      if (is.list(loc)) loc <- unlist(loc, use.names = TRUE)
      if (is.list(se))  se  <- unlist(se,  use.names = TRUE)

      parsed <- private$.parseThresholdNames(names(loc), item_names)

      data.frame(
        Item      = parsed$Item,
        Threshold = parsed$Threshold,
        DIFgroup  = dif_group,
        Location  = unname(as.numeric(loc)),
        SE        = unname(as.numeric(se)),
        stringsAsFactors = FALSE
      )
    },

    # Parse "thresh beta I1.c1" / "beta I1.c1" / "I1.c1" / "I1" into
    # Item + Threshold parts.
    .parseThresholdNames = function(nm, item_names) {
      bare <- sub("^(thresh\\s+)?beta\\s+", "", nm)

      matched_item <- vapply(bare, function(b) {
        hits <- item_names[startsWith(b, item_names) |
                             startsWith(b, paste0(item_names, "."))]
        if (length(hits) == 0L) NA_character_
        else hits[which.max(nchar(hits))]
      }, character(1L))

      thr_part <- mapply(function(b, it) {
        if (is.na(it)) return(NA_character_)
        rest <- substr(b, nchar(it) + 1L, nchar(b))
        rest <- sub("^\\.", "", rest)
        if (!nzchar(rest)) "1" else rest
      }, bare, matched_item, USE.NAMES = FALSE)

      if (any(is.na(matched_item))) {
        bad <- nm[is.na(matched_item)]
        stop("Could not match threshold parameter name(s) to items: ",
             paste(utils::head(bad, 5L), collapse = ", "))
      }

      list(Item = unname(matched_item), Threshold = unname(thr_part))
    },

    # Walk LRtest's per-group fits
    .extractLrLocations = function(lrt, groups, item_names) {
      parts <- lapply(seq_along(groups), function(g) {
        fit_g <- lrt$fitobj[[g]]
        private$.oneFitLong(fit_g, dif_group = groups[g],
                            item_names = item_names)
      })
      out <- do.call(rbind, parts)
      rownames(out) <- NULL
      out
    },

    # Reshape long -> wide on DIFgroup; add MaxDiff
    .buildWideTable = function(long_df, groups, level) {
      id_cols <- if (level == "item") "Item" else c("Item", "Threshold")

      loc_wide <- stats::reshape(
        long_df[, c(id_cols, "DIFgroup", "Location"), drop = FALSE],
        idvar     = id_cols,
        timevar   = "DIFgroup",
        direction = "wide"
      )
      names(loc_wide) <- sub("^Location\\.", "", names(loc_wide))

      se_wide <- stats::reshape(
        long_df[, c(id_cols, "DIFgroup", "SE"), drop = FALSE],
        idvar     = id_cols,
        timevar   = "DIFgroup",
        direction = "wide"
      )
      names(se_wide) <- sub("^SE\\.", "SE_", names(se_wide))

      out <- merge(loc_wide, se_wide, by = id_cols, sort = FALSE)

      # MaxDiff across the per-group columns (excludes "All")
      group_mat <- as.matrix(out[, groups, drop = FALSE])
      out$MaxDiff <- apply(group_mat, 1L, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) NA_real_ else max(x) - min(x)
      })

      ord <- if (level == "threshold") {
        order(out$Item, out$Threshold)
      } else {
        order(out$Item)
      }
      out <- out[ord, , drop = FALSE]
      rownames(out) <- NULL
      out
    },

    # ---------------------------------------------------------------------
    # .lrtPlot -- facet-per-item ggplot
    # ---------------------------------------------------------------------
    .lrtPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state    <- image$state
      plot_df  <- state$plot_df
      groups   <- state$groups
      level    <- state$level
      conf     <- state$conf

      z <- stats::qnorm(1 - (1 - conf) / 2)

      # Show only the per-group locations (no overall/All diamond)
      plot_df    <- plot_df[plot_df$DIFgroup %in% groups, , drop = FALSE]
      plot_df$DIFgroup <- factor(plot_df$DIFgroup, levels = groups)

      caption <- er2_caption(paste0(
        "Error bars: ", round(conf * 100), "% CI. ",
        "n = ", state$n_complete, "."
      ))

      base_theme <- ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(
          legend.position = "none"
        ) +
        er2_axis_margins() +
        er2_plot_caption()

      if (level == "item") {
        p <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(x = .data$DIFgroup,
                       y = .data$Location,
                       group  = .data$Item,
                       colour = .data$Item)
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point() +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$Location - z * .data$SE,
                         ymax = .data$Location + z * .data$SE),
            width = 0.1
          ) +
          ggplot2::scale_x_discrete(labels = private$.wrapLabels) +
          ggplot2::facet_wrap(~ Item) +
          ggplot2::labs(
            title    = "DIF: item locations by group",
            subtitle = "Item locations are means of threshold locations",
            x        = "DIF group",
            y        = "Item location (logits)",
            caption  = caption
          ) +
          base_theme
      } else {
        p <- ggplot2::ggplot(
          plot_df,
          ggplot2::aes(x = .data$DIFgroup,
                       y = .data$Location,
                       group  = .data$Threshold,
                       colour = .data$Threshold)
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point(alpha = 0.9) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data$Location - z * .data$SE,
                         ymax = .data$Location + z * .data$SE),
            width = 0.1
          ) +
          ggplot2::scale_x_discrete(labels = private$.wrapLabels) +
          ggplot2::facet_wrap(~ Item) +
          ggplot2::labs(
            title   = "DIF: threshold locations by group",
            x       = "DIF group",
            y       = "Threshold location (logits)",
            caption = caption
          ) +
          base_theme
      }

      print(p)
      TRUE
    },

    # ------------------------------------------------------------------
    # Build per-(item x category x group) counts for the tileplot.
    # `df` is numeric, complete-cases item data; `dif_vec` is the
    # corresponding DIF variable values (length = nrow(df)).
    # Mirrors easyRasch2::RMplotTile logic.
    # ------------------------------------------------------------------
    .computeTileCounts = function(df, dif_vec) {
      item_names <- names(df)

      all_vals       <- unlist(df, use.names = FALSE)
      all_vals       <- all_vals[!is.na(all_vals)]
      min_val        <- min(all_vals)
      max_val        <- max(all_vals)
      all_categories <- seq(min_val, max_val)

      dif_factor <- if (is.factor(dif_vec)) {
        droplevels(dif_vec)
      } else {
        as.factor(dif_vec)
      }

      parts <- lapply(levels(dif_factor), function(g) {
        mask <- !is.na(dif_factor) & dif_factor == g
        sub  <- df[mask, , drop = FALSE]
        rows <- lapply(item_names, function(it) {
          vals <- sub[[it]]
          vals <- vals[!is.na(vals)]
          tab  <- table(factor(vals, levels = all_categories))
          data.frame(
            item     = it,
            category = as.integer(names(tab)),
            n        = as.integer(tab),
            group    = g,
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, rows)
      })
      count_df <- do.call(rbind, parts)
      rownames(count_df) <- NULL
      count_df$group <- factor(count_df$group, levels = levels(dif_factor))

      totals <- stats::aggregate(
        count_df$n,
        by  = count_df[, c("item", "group"), drop = FALSE],
        FUN = sum
      )
      colnames(totals)[ncol(totals)] <- "total"
      count_df <- merge(count_df, totals, by = c("item", "group"),
                       sort = FALSE)
      count_df$percentage <- round(count_df$n / count_df$total * 100, 1)

      # Item ordering: top of y-axis = first column of df
      count_df$item_label <- factor(count_df$item, levels = rev(item_names))

      list(
        count_df       = count_df,
        all_categories = all_categories,
        item_names     = item_names
      )
    },

    # ------------------------------------------------------------------
    # .tileplot -- faceted tile of item x category response counts,
    # faceted by the DIF grouping variable. Mirrors RMplotTile.
    # ------------------------------------------------------------------
    .tileplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      state          <- image$state
      count_df       <- state$count_df
      all_categories <- state$all_categories
      cutoff         <- state$cutoff
      use_pct        <- isTRUE(state$percent)

      p <- ggplot2::ggplot(
        count_df,
        ggplot2::aes(x = .data$category,
                     y = .data$item_label,
                     fill = .data$n)
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(
          expression(italic(n)),
          limits = c(0, NA)
        ) +
        ggplot2::scale_x_continuous(
          "Response category",
          expand = c(0, 0),
          breaks = all_categories
        )

      # Cell labels (counts or percentages) with low-count highlighting
      label_aes <- if (use_pct) {
        ggplot2::aes(label = paste0(.data$percentage, "%"),
                     color = ifelse(.data$n < cutoff, "red", "orange"))
      } else {
        ggplot2::aes(label = .data$n,
                     color = ifelse(.data$n < cutoff, "red", "orange"))
      }
      p <- p +
        ggplot2::geom_text(label_aes) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()

      # Wrap long DIF-group facet labels
      p <- p +
        ggplot2::facet_wrap(~ group,
                            labeller = ggplot2::labeller(
                              group = private$.wrapLabels
                            )) +
        ggplot2::labs(y = "Items") +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(
          axis.text.x   = ggplot2::element_text(size = 10),
          panel.grid    = ggplot2::element_blank(),
          panel.spacing = ggplot2::unit(0.7, "cm")
        ) +
        er2_axis_margins()

      print(p)
      TRUE
    }
  )
)
