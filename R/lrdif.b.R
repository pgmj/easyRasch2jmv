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

      # 1. Return early / explain if requirements not met (the model needs
      # at least 2 items; no DIF variable selected is the normal initial
      # state and stays silent)
      vars   <- self$options$vars
      difVar <- self$options$difVar
      if (is.null(vars) || length(vars) == 0 || is.null(difVar))
        return()
      if (length(vars) < 2) {
        self$results$lrtNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model in each DIF group. Select at least 2 items.</p>"
        ))
        return()
      }

      # 2. Extract data + validate items
      data <- self$data
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      # 3. DIF variable + joint complete-case handling
      n_total    <- nrow(df)
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

      # Sparse-category warning per DIF group -- the relevant split for
      # this analysis. Shown even when the tileplot is off, pointing
      # users to it for visual inspection.
      sparse_msg <- sparse_note_grouped(df, dif_factor)
      if (!is.null(sparse_msg)) {
        self$results$lrtTable$setNote("sparse", paste0(
          sparse_msg, " Enable 'Response distribution by DIF group' ",
          "to inspect the counts."
        ))
      }

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$lrtTable$setNote("duplicate", dup_msg)

      # 4. Read options
      level         <- self$options$level
      cutoff_val    <- self$options$cutoff
      sort_by_max   <- isTRUE(self$options$sortByMaxDiff)
      show_figure   <- isTRUE(self$options$showFigure)
      conf_level    <- self$options$confLevel / 100

      # 5. Run analysis via easyRasch2 (eRm::LRtest inside; RMdifLR()
      # deliberately remains eRm-based upstream as well, so the module and
      # the R package are identical by construction). Package messages are
      # suppressed -- the module pre-drops NA rows itself so the notes can
      # report the counts.
      tryCatch({
        # rgl workaround
        old_rgl <- getOption("rgl.useNULL")
        options(rgl.useNULL = TRUE)
        on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

        is_polytomous <- max(as.matrix(df), na.rm = TRUE) > 1L

        table_df <- tryCatch(
          suppressWarnings(suppressMessages(
            easyRasch2::RMdifLR(
              df,
              dif_var = dif_factor,
              level   = level,
              cutoff  = if (cutoff_val > 0) cutoff_val else NULL,
              output  = "dataframe"
            )
          )),
          error = function(e) {
            stop(paste0(
              conditionMessage(e),
              " This often indicates an empty response category in one ",
              "subgroup. Inspect the response distribution per group ",
              "before running the LR test."
            ))
          }
        )
        lr_summary <- attr(table_df, "lr_test")

        # 6. Sort if requested
        if (sort_by_max) {
          table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
          rownames(table_df) <- NULL
        }

        # 7. Set up table columns now that the group structure is known.
        # NOTE: these group columns depend on the levels of difVar, which
        # cannot be read until the data is loaded -- so they genuinely
        # cannot be moved to .init() (defensible Level 3 case). This causes
        # a one-time UI restructure when the analysis first opens: the table
        # appears blank, then gains its group columns once difVar resolves.
        # Numeric columns use format = "zto" to match the formatting of
        # all other tables in the module.
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

        # 8. Populate rows. Raw numerics so the jamovi frontend applies
        # the user's "Number format" preferences.
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

        # Footnotes: MaxDiff definition + flagging rule
        table$setNote("maxdiff", paste0(
          "MaxDiff = difference between the highest and lowest group ",
          "location (the All column is excluded)."
        ))
        if (cutoff_val > 0) {
          table$setNote("flag", paste0(
            "Flagged = TRUE when MaxDiff exceeds ", cutoff_val, " logits."
          ))
        }

        # 9. LR test note (p-value rounded to 3 digits)
        model_name <- if (is_polytomous) "Partial Credit Model" else "Rasch Model"
        p_round <- round(lr_summary$p_value, 3)
        p_str   <- if (p_round == 0) "&lt; 0.001" else format(p_round, nsmall = 3)
        n_excluded <- n_total - n_complete
        excluded_clause <- if (n_excluded > 0L) {
          paste0(" (", n_excluded, " of ", n_total, " row(s) excluded ",
                 "due to missing item responses or a missing DIF value; ",
                 "eRm::LRtest requires complete cases)")
        } else ""
        lr_html <- paste0(
          "<p><b>Andersen LR test:</b> χ<sup>2</sup> = ",
          round(lr_summary$LR, 3),
          ", <i>df</i> = ", lr_summary$df,
          ", <i>p</i> = ", p_str,
          ". ", model_name, " split by `", difVar,
          "` (", length(groups), " groups: ",
          paste(groups, collapse = ", "), "). n = ", n_complete,
          " complete cases", excluded_clause,
          ". Results are identical to easyRasch2::RMdifLR()."
        )
        if (cutoff_val > 0) {
          lr_html <- paste0(
            lr_html,
            " Items flagged when MaxDiff &gt; ", cutoff_val, " logits."
          )
        }
        lr_html <- paste0(lr_html, "</p>")
        self$results$lrtNote$setContent(lr_html)

        # 10. Save state for the figure (drawn by RMdifLR(output =
        # "ggplot") in the render function)
        if (show_figure) {
          self$results$lrtPlot$setState(list(df = df, dif = dif_factor))
        }

        # 11. Tileplot: per-item × category × DIF-group response counts,
        # drawn by easyRasch2::RMplotTile() in the render function.
        if (isTRUE(self$options$showTileplot)) {
          self$results$tileplot$setState(list(df = df, dif = dif_factor))
        }
      }, error = function(e) {
        stop(paste("Error in LR-based DIF analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # jamovi column names must be syntactically safe; group levels may
    # contain spaces or other characters.
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

    # ---------------------------------------------------------------------
    # Per-group locations figure — easyRasch2::RMdifLR(output = "ggplot")
    # ---------------------------------------------------------------------
    .lrtPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMdifLR(
          image$state$df,
          dif_var = image$state$dif,
          level   = self$options$level,
          cutoff  = if (self$options$cutoff > 0) self$options$cutoff else NULL,
          conf    = self$options$confLevel / 100,
          output  = "ggplot"
        )
      ))
      p <- er2_bump_text(p)

      print(p)
      TRUE
    },

    # ---------------------------------------------------------------------
    # Faceted tileplot of item × category response counts by DIF group,
    # drawn by easyRasch2::RMplotTile().
    # ---------------------------------------------------------------------
    .tileplot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      p <- suppressWarnings(suppressMessages(
        easyRasch2::RMplotTile(
          image$state$df,
          group   = image$state$dif,
          cutoff  = self$options$tileCutoff,
          percent = isTRUE(self$options$tilePercent)
        )
      ))
      p <- er2_bump_text(p)

      print(p)
      TRUE
    }
  )
)
