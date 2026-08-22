#' @export
diftreeClass <- R6::R6Class(
  "diftreeClass",
  inherit = diftreeBase,
  private = list(

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {
      vars <- self$options$vars
      covs <- self$options$covariates
      if ((is.null(vars) || length(vars) == 0) &&
          (is.null(covs) || length(covs) == 0))
        return()
      if (is.null(vars) || length(vars) < 3) {
        self$results$difNote$setContent(paste0(
          "<p>This analysis requires at least <b>3 items</b>. Select at ",
          "least 3 items and at least one covariate.</p>"
        ))
        return()
      }
      if (is.null(covs) || length(covs) == 0) {
        self$results$difNote$setContent(paste0(
          "<p>Assign at least one <b>covariate</b> (the variables the ",
          "tree may split on -- grouping factors and/or continuous ",
          "variables such as age).</p>"
        ))
        return()
      }

      data <- self$data
      # Shared validation for the items: conversion, all-NA / sentinel
      # checks, response validation, per-item variation, identical items
      df <- prepare_item_data(data, vars)

      sparse_msg <- sparse_note(df)
      if (!is.null(sparse_msg))
        self$results$difTable$setNote("sparse", sparse_msg)
      recode_msg <- recode_note(data, vars)
      if (!is.null(recode_msg))
        self$results$difTable$setNote("recode", recode_msg)

      dup_msg <- duplicate_items_note(df)
      if (!is.null(dup_msg))
        self$results$difTable$setNote("duplicate", dup_msg)

      # Covariates keep their measure types: factors stay factors,
      # continuous variables stay numeric (the tree machinery handles
      # both natively).
      cov_df <- as.data.frame(data[covs], stringsAsFactors = FALSE)
      for (nm in names(cov_df)) {
        if (is.factor(cov_df[[nm]])) cov_df[[nm]] <- droplevels(cov_df[[nm]])
      }
      n_total  <- nrow(df)
      n_cov_na <- sum(!stats::complete.cases(cov_df))

      # rgl workaround (iarm dependency chain)
      old_rgl <- getOption("rgl.useNULL")
      options(rgl.useNULL = TRUE)
      on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

      tryCatch({
        # All computation is delegated to the easyRasch2 package, so
        # results are identical to RMdifTree() (deterministic -- no
        # simulation, no seed). psychotree's per-node rescaling
        # diagnostic arrives as a message and is surfaced in the note.
        rescale_msgs <- character(0)
        res <- withCallingHandlers(
          suppressWarnings(
            easyRasch2::RMdifTree(
              df,
              covariates       = cov_df,
              purification     = if (isTRUE(self$options$purification))
                                   "iterative" else "none",
              p_adj            = self$options$pAdj,
              alpha            = self$options$alpha,
              prune_negligible = isTRUE(self$options$pruneNegligible),
              on_rescale       = "message",
              output           = "list"
            )
          ),
          message = function(m) {
            rescale_msgs <<- c(rescale_msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
          }
        )
        tbl <- res$table
        es  <- attr(tbl, "effect_size")
        is_mh <- identical(es, "MH")

        # --- Effect-size table ------------------------------------------
        table <- self$results$difTable
        table$getColumn("effectSize")$setTitle(
          if (is_mh) "MH (Delta scale)" else "Partial gamma"
        )
        for (i in seq_len(nrow(tbl))) {
          table$addRow(rowKey = i, values = list(
            node       = as.integer(tbl$NodeID[i]),
            split      = tbl$Split[i],
            item       = tbl$Item[i],
            effectSize = tbl$EffectSize[i],
            se         = tbl$SE[i],
            class      = tbl$Class[i],
            flagged    = if (isTRUE(tbl$Flagged[i])) "TRUE" else "",
            nLeft      = as.integer(tbl$n_left[i]),
            nRight     = as.integer(tbl$n_right[i])
          ))
        }

        alpha <- self$options$alpha
        if (nrow(tbl) > 0) {
          class_note <- if (is_mh) {
            paste0(
              "Mantel-Haenszel effect size on the ETS Delta scale. ",
              "Class A (negligible): |Delta| < 1 or not significantly ",
              "different from 0 at alpha = ", alpha, "; C (large): ",
              "|Delta| >= 1.5 and significantly > 1; B (moderate): ",
              "otherwise. Positive values = item more difficult for the ",
              "right-hand group of the split."
            )
          } else {
            paste0(
              "Partial gamma effect size (Bjorner et al., 1998). Class ",
              "A (negligible): |gamma| < 0.21 or not significantly ",
              "different from 0 at alpha = ", alpha, "; C (large): ",
              "|gamma| > 0.31 and significantly outside +/-0.21; B ",
              "(moderate): otherwise.",
              if (self$options$pAdj != "none") {
                paste0(" p-values adjusted (",
                       if (self$options$pAdj == "fdr")
                         "Benjamini-Hochberg FDR" else "Bonferroni",
                       ") before classification.")
              } else ""
            )
          }
          table$setNote("class", class_note)
          table$setNote("flag", "Flagged: Class B or C.")
          table$setNote("caveat", paste0(
            "The A/B/C boundaries are conventions from large-scale ",
            "educational testing, not values calibrated to this sample ",
            "and these items -- read the classification as a rough ",
            "magnitude guide. For a sample-calibrated DIF test with a ",
            "pre-specified grouping, use the Partial Gamma DIF analysis ",
            "with simulation-based cutoffs."
          ))
          if (is_mh && self$options$pAdj != "none") {
            table$setNote("padj", paste0(
              "The p-value adjustment option is ignored for dichotomous ",
              "data (the Mantel-Haenszel classification uses the ETS ",
              "test directly)."
            ))
          }
        }

        # --- Note -----------------------------------------------------------
        cov_na_clause <- if (n_cov_na > 0) {
          paste0(" ", n_cov_na, " row(s) with a missing covariate value ",
                 "are excluded by the tree machinery; missing item ",
                 "responses are handled by the models.")
        } else ""
        purif_clause <- if (isTRUE(self$options$purification)) {
          paste0(" Effect sizes use iterative purification: items ",
                 "already classified as DIF are excluded from the ",
                 "matching score and the effect sizes are recomputed.")
        } else ""
        rescale_clause <- if (length(rescale_msgs) > 0) {
          paste0(" <b>Note:</b> ", paste(rescale_msgs, collapse = " "),
                 " This concerns only the item parameters displayed in ",
                 "the tree figure; the effect sizes are computed on the ",
                 "raw responses and are not affected.")
        } else ""
        result_clause <- if (nrow(tbl) == 0) {
          paste0(
            " <b>The tree found no splits:</b> the parameter-instability ",
            "tests detected no covariate-related DIF, i.e. the item ",
            "parameters are stable across the supplied covariates. The ",
            "figure shows the single-group model."
          )
        } else {
          paste0(
            " The tree found ", length(unique(tbl$NodeID)), " split(s); ",
            "the table classifies every item's DIF effect size at each ",
            "split."
          )
        }
        self$results$difNote$setContent(paste0(
          "<p>Tree-based DIF analysis (",
          if (is_mh) {
            "Rasch model; Mantel-Haenszel effect sizes on the ETS Delta scale"
          } else {
            "partial credit model; partial gamma effect sizes"
          },
          ") for n = ", n_total, " respondents and ", length(covs),
          " covariate(s): the sample is recursively split wherever the ",
          "item parameters are unstable along a covariate, so groups ",
          "need not be pre-specified; continuous covariates are split ",
          "at data-driven cutpoints and interactions appear as nested ",
          "splits.", cov_na_clause, result_clause, purif_clause,
          rescale_clause,
          " The analysis is deterministic. Results are identical to ",
          "easyRasch2::RMdifTree().</p>"
        ))

        # --- Figure state -----------------------------------------------
        self$results$treePlot$setState(list(tree_ok = TRUE))
        private$.tree <- res$tree

      }, error = function(e) {
        stop(paste("Error in tree-based DIF analysis:", e$message))
      })
    },

    # The fitted tree is kept on the instance for the render pass of the
    # same run. partykit trees embed per-node models and data, which can
    # exceed jamovi's state-size budget, so the tree is deliberately NOT
    # stored in the image state; if the render happens without a
    # preceding .run in this process (e.g. reopening a saved file resizes
    # a figure), the tree is refitted -- deterministic, so the figure is
    # identical.
    .tree = NULL,

    # ---------------------------------------------------------------------
    # DIF tree figure -- partykit plot with item names on the
    # terminal-node axes (the upstream-documented tp_args).
    # ---------------------------------------------------------------------
    .treePlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)

      tree <- private$.tree
      if (is.null(tree)) {
        vars <- self$options$vars
        covs <- self$options$covariates
        df <- prepare_item_data(self$data, vars)
        cov_df <- as.data.frame(self$data[covs], stringsAsFactors = FALSE)
        for (nm in names(cov_df)) {
          if (is.factor(cov_df[[nm]]))
            cov_df[[nm]] <- droplevels(cov_df[[nm]])
        }
        tree <- suppressMessages(suppressWarnings(
          easyRasch2::RMdifTree(
            df,
            covariates       = cov_df,
            purification     = if (isTRUE(self$options$purification))
                                 "iterative" else "none",
            p_adj            = self$options$pAdj,
            alpha            = self$options$alpha,
            prune_negligible = isTRUE(self$options$pruneNegligible),
            on_rescale       = "message",
            output           = "tree"
          )
        ))
      }

      graphics::plot(tree, tp_args = list(names = TRUE))
      TRUE
    }
  )
)
