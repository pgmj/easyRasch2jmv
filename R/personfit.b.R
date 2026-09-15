#' @export
personfitClass <- R6::R6Class(
  "personfitClass",
  inherit = personfitBase,
  private = list(

    # Statistics requested via the checkboxes, in upstream order.
    .statistics = function() {
      stats <- character(0)
      if (isTRUE(self$options$statInfit))  stats <- c(stats, "infit")
      if (isTRUE(self$options$statOutfit)) stats <- c(stats, "outfit")
      if (isTRUE(self$options$statLz))     stats <- c(stats, "lz")
      stats
    },

    # ---------------------------------------------------------------------
    # .init -- pre-create the summary rows (fixed set + one flagged-count
    # row per requested statistic, known from the options).
    # ---------------------------------------------------------------------
    .init = function() {
      st <- self$results$summaryTable
      st$addRow(rowKey = "n",
                values = list(statistic = "Respondents assessed", value = NA_real_))
      st$addRow(rowKey = "extreme",
                values = list(statistic = "Extreme scorers (not assessed)", value = NA_real_))
      st$addRow(rowKey = "flagged",
                values = list(statistic = "Flagged (any statistic)", value = NA_real_))
      st$addRow(rowKey = "flaggedPct",
                values = list(statistic = "Flagged (%)", value = NA_real_))
      labels <- c(infit = "Flagged by infit", outfit = "Flagged by outfit",
                  lz = "Flagged by lz")
      for (s in private$.statistics()) {
        st$addRow(rowKey = paste0("by_", s),
                  values = list(statistic = labels[[s]], value = NA_real_))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {
      if (is.null(self$options$vars) || length(self$options$vars) == 0)
        return()
      if (length(self$options$vars) < 2) {
        self$results$personfitNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items</b> to fit a ",
          "Rasch model. Select at least 2 items.</p>"
        ))
        return()
      }
      stats <- private$.statistics()
      if (length(stats) == 0) {
        self$results$personfitNote$setContent(paste0(
          "<p>Select at least one person-fit statistic (infit, outfit, ",
          "or lz).</p>"
        ))
        return()
      }

      data <- self$data
      vars <- self$options$vars
      # Shared validation: conversion, all-NA / sentinel checks,
      # response validation, per-item variation, identical-items check
      df <- prepare_item_data(data, vars)

      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$summaryTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$summaryTable$setNote("duplicate", dup_msg)

      # Respondents with no responses on any selected item cannot be
      # assessed; drop them up front and keep the row mapping for the
      # output variables.
      n_total <- nrow(df)
      keep    <- rowSums(!is.na(df)) > 0
      df_used <- df[keep, , drop = FALSE]
      n_used  <- nrow(df_used)
      if (n_used == 0)
        stop("No respondents with any responses on the selected items.")

      # Sparse response categories destabilise CML thresholds; switch to
      # MML with a note, as in the person-parameters and targeting
      # analyses.
      sparse_items <- sparse_category_items(df_used, min_n = 3L)
      use_mml      <- length(sparse_items) > 0L
      estimator    <- if (use_mml) "MML" else "CML"

      flag_dir   <- self$options$flagDirection
      flag_alpha <- self$options$flagAlpha
      iterations <- self$options$iterations

      tryCatch({
        # The per-person resampling is the expensive part, and jamovi
        # reruns .run on every option change -- including the
        # output-variable toggles, which do not affect it. The hidden
        # simCache element carries the fit dataframe and the person-fit
        # maps (one RMpersonFit(output = "list") call); the signature
        # check makes reuse self-validating rather than relying on
        # clearWith alone. The figures read their plots from the cache.
        sim_sig <- list(
          statistics = stats,
          estimator  = estimator,
          flag       = flag_dir,
          flag_alpha = flag_alpha,
          iterations = iterations,
          seed       = as.integer(self$options$seed)
        )
        cached <- self$results$simCache$state
        if (!is.null(cached) && !is.null(cached$fit) &&
            identical(cached$sig, sim_sig) && identical(cached$df, df_used)) {
          res <- list(fit = cached$fit, grobs = cached$grobs)
        } else {
          # All computation is delegated to the easyRasch2 package, so
          # results are numerically identical to RMpersonFit() with the
          # same seed and iterations.
          res <- suppressWarnings(suppressMessages(
            easyRasch2::RMpersonFit(
              df_used,
              statistics   = stats,
              estimator    = estimator,
              theta_method = "WLE",
              iterations   = iterations,
              flag_alpha   = flag_alpha,
              flag         = flag_dir,
              parallel     = FALSE,
              seed         = as.integer(self$options$seed),
              output       = "list"
            )
          ))
          # The three maps are stored built. As ggplot objects they came to
          # 287, 287 and 752 KB compressed, over jmvcore's 500 KB warning
          # for one element and all of it written into the saved .omv; built,
          # they are 14, 14 and 12 KB.
          res$grobs <- lapply(res$plots, er2_plot_grob)
          res$plots <- NULL
          self$results$simCache$setState(list(
            fit = res$fit, grobs = res$grobs, sig = sim_sig, df = df_used
          ))
        }
        fit <- res$fit
        if (nrow(fit) != n_used)
          stop("Internal error: person-fit rows do not match the data.")

        # --- Summary table ----------------------------------------------
        # Extreme scorers get NA statistics and are not assessed.
        stat_cols <- c(infit = "infit_msq", outfit = "outfit_msq",
                       lz = "lz")[stats]
        p_cols    <- c(infit = "p_infit", outfit = "p_outfit",
                       lz = "p_lz")[stats]
        assessed  <- !is.na(fit[[stat_cols[1]]])
        n_assessed <- sum(assessed)
        n_flagged  <- sum(fit$flagged, na.rm = TRUE)

        st <- self$results$summaryTable
        st$setRow(rowKey = "n",       values = list(value = n_assessed))
        st$setRow(rowKey = "extreme", values = list(value = n_used - n_assessed))
        st$setRow(rowKey = "flagged", values = list(value = n_flagged))
        st$setRow(rowKey = "flaggedPct",
                  values = list(value = 100 * n_flagged / n_assessed))
        for (s in stats) {
          st$setRow(rowKey = paste0("by_", s), values = list(
            value = sum(fit[[p_cols[[s]]]] < flag_alpha, na.rm = TRUE)
          ))
        }
        st$setNote("flag", paste0(
          "Flagged: uncorrected resampled p-value < ", flag_alpha,
          " for any requested statistic (per-person screening; under ",
          "model fit, about ", round(100 * flag_alpha, 1), "% of ",
          "respondents are flagged by chance)."
        ))
        if (use_mml) {
          st$setNote("sparse", paste0(
            "Item parameters estimated by MML instead of CML because ",
            "item(s) ", paste(sparse_items, collapse = ", "), " have ",
            "response categories with fewer than 3 observations; MML is ",
            "more numerically stable under sparse categories."
          ))
        }

        # --- Output variables (written into the dataset) -----------------
        expand <- function(v, fill = NA_real_) {
          out <- rep(fill, n_total)
          out[keep] <- v
          out
        }
        row_nums <- rownames(df)
        set_output <- function(opt_name, values) {
          if (!isTRUE(self$options[[opt_name]])) return()
          out <- self$results[[opt_name]]
          out$setRowNums(row_nums)
          out$setValues(values)
        }
        get_col <- function(col) {
          if (col %in% names(fit)) fit[[col]] else rep(NA_real_, n_used)
        }
        set_output("outputFlagged",
                   expand(as.integer(fit$flagged), fill = NA_integer_))
        set_output("outputInfit",    expand(get_col("infit_msq")))
        set_output("outputOutfit",   expand(get_col("outfit_msq")))
        set_output("outputLz",       expand(get_col("lz")))
        set_output("outputPInfit",   expand(get_col("p_infit")))
        set_output("outputPOutfit",  expand(get_col("p_outfit")))
        set_output("outputPLz",      expand(get_col("p_lz")))

        # --- Note ---------------------------------------------------------
        drop_clause <- if (n_used < n_total) {
          paste0(" (", n_total - n_used, " row(s) without any responses ",
                 "on the selected items excluded)")
        } else ""
        dir_clause <- if (flag_dir == "underfit") {
          paste0(" MSQ p-values are one-sided (underfit only: noisy or ",
                 "erratic responding); benign overfit is ignored.")
        } else {
          paste0(" MSQ p-values are two-sided (underfit = noisy or ",
                 "erratic responding; overfit = overly deterministic, ",
                 "occasionally suspicious).")
        }
        self$results$personfitNote$setContent(paste0(
          "<p>Person-fit statistics for n = ", n_used, " respondents",
          drop_clause, ", based on ",
          if (use_mml) "MML" else "conditional maximum likelihood (CML)",
          " item parameters. The MSQ statistics are conditional on the ",
          "total score (no person estimate enters); lz uses WLE person ",
          "locations. Significance is assessed by ", iterations,
          " Monte-Carlo replications per person under the fitted model ",
          "(the asymptotic null distributions are unreliable).",
          dir_clause,
          " Extreme scorers cannot be assessed and have empty cells in ",
          "the saved variables -- note this when filtering on the ",
          "flag. Tick the checkboxes under <i>Save to dataset</i> to add ",
          "the flag, statistics, and p-values as variables, e.g. to ",
          "filter out aberrant respondents before rerunning other ",
          "analyses. Results are identical to easyRasch2::RMpersonFit().",
          "</p>"
        ))

      }, error = function(e) {
        stop(paste("Error in person fit analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Person-fit maps -- drawn by easyRasch2::RMpersonFit(output = "list")
    # in .run, restyled there and stored built in the simCache element
    # (single storage). See er2_plot_grob() for why they are stored built.
    # ---------------------------------------------------------------------
    .renderMap = function(stat) {
      state <- self$results$simCache$state
      if (is.null(state)) return(FALSE)
      er2_draw_grob(state$grobs[[stat]])
    },
    .infitMap  = function(image, ggtheme, theme, ...) private$.renderMap("infit"),
    .outfitMap = function(image, ggtheme, theme, ...) private$.renderMap("outfit"),
    .lzMap     = function(image, ggtheme, theme, ...) private$.renderMap("lz")
  )
)
