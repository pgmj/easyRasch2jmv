#' @export
personchangeClass <- R6::R6Class(
  "personchangeClass",
  inherit = personchangeBase,
  private = list(

    # ---------------------------------------------------------------------
    # .init -- the summary rows are fully determined by the options, so
    # they are built here and filled by .run().
    # ---------------------------------------------------------------------
    .init = function() {
      if (length(self$options$vars1) < 2 || length(self$options$vars2) < 2)
        return()

      st <- self$results$summaryTable
      rows <- list(
        n        = "Respondents tested",
        increase = "Increase in theta",
        none     = "No change detected",
        decrease = "Decrease in theta",
        critlo   = "Critical value (lower)",
        crithi   = "Critical value (upper)"
      )
      for (key in names(rows)) {
        st$addRow(rowKey = key, values = list(
          statistic = rows[[key]], value = NA_real_, notes = ""
        ))
      }
    },

    # ---------------------------------------------------------------------
    # .run
    # ---------------------------------------------------------------------
    .run = function() {
      v1 <- self$options$vars1
      v2 <- self$options$vars2
      if (is.null(v1) || is.null(v2) || length(v1) == 0 || length(v2) == 0)
        return()

      if (length(v1) < 2 || length(v2) < 2) {
        self$results$changeNote$setContent(paste0(
          "<p>This analysis requires at least <b>2 items per occasion</b> ",
          "to estimate the item parameters that person locations rest on. ",
          "Select at least 2 items in each box.</p>"
        ))
        return()
      }
      if (length(v1) != length(v2)) {
        self$results$changeNote$setContent(paste0(
          "<p><b>The two occasions must have the same number of items.</b> ",
          "Time 1 has ", length(v1), " and Time 2 has ", length(v2),
          ". Items are paired by their position in the two boxes, so the ",
          "lists have to be the same length and in the same order. Data ",
          "must be in wide format, with both occasions on the same row.</p>"
        ))
        return()
      }

      data <- self$data

      # --- Joint preparation ------------------------------------------------
      # Both occasions are prepared as ONE frame rather than separately.
      # prepare_item_data() decides the 1-based recode from the minimum of
      # every item it is given, so preparing the occasions apart could shift
      # one and not the other and silently destroy comparability. Names are
      # made unique so a variable used in both boxes does not collide.
      combined_names <- make.unique(c(v1, v2), sep = "_")
      raw <- data[, c(v1, v2), drop = FALSE]
      names(raw) <- combined_names
      df_all <- prepare_item_data(raw, combined_names)

      k <- length(v1)
      df1 <- df_all[, seq_len(k), drop = FALSE]
      df2 <- df_all[, k + seq_len(k), drop = FALSE]
      names(df1) <- v1
      names(df2) <- v2

      # Strip jmvcore S4 column wrappers, as elsewhere in the module
      as_plain <- function(d) {
        nm <- names(d)
        out <- as.data.frame(matrix(
          as.numeric(as.matrix(d)), nrow = nrow(d), ncol = ncol(d),
          dimnames = list(NULL, nm)
        ))
        out
      }
      df1 <- as_plain(df1)
      df2 <- as_plain(df2)

      # --- Pairing table (above the figure, deliberately) -------------------
      pt <- self$results$pairingTable
      cat_mismatch <- character(0)
      for (i in seq_len(k)) {
        m1 <- suppressWarnings(max(df1[[i]], na.rm = TRUE))
        m2 <- suppressWarnings(max(df2[[i]], na.rm = TRUE))
        cats <- if (!is.finite(m1) || !is.finite(m2)) {
          "—"
        } else if (m1 == m2) {
          paste0("0–", m1)
        } else {
          cat_mismatch <- c(cat_mismatch, paste0(v1[i], " / ", v2[i]))
          paste0("0–", m1, " vs 0–", m2)
        }
        pt$addRow(rowKey = i, values = list(
          pair = i, item1 = v1[i], item2 = v2[i], cats = cats
        ))
      }
      pt$setNote("pairing", paste0(
        "Items are matched by their position in the two variable boxes, ",
        "not by name. Check this table before reading the results: a ",
        "mispaired analysis produces a plausible-looking figure."
      ))
      if (length(cat_mismatch) > 0L) {
        pt$setNote("cats", paste0(
          "The highest observed response differs between occasions for ",
          paste(cat_mismatch, collapse = ", "),
          ". That is expected when a category simply went unused at one ",
          "occasion, but it also happens when the wrong items were paired."
        ))
      }
      recode_msg <- recode_note(raw, combined_names)
      if (!is.null(recode_msg)) pt$setNote("recode", recode_msg)
      dup_msg <- duplicate_items_note(df_all)
      if (!is.null(dup_msg)) pt$setNote("duplicate", dup_msg)

      # --- Who can be tested ------------------------------------------------
      # A respondent needs at least one response at EACH occasion. Rows that
      # are entirely missing at one occasion are dropped up front: the
      # bundled easyRasch2 release fails its CML fit on them ("subscript out
      # of bounds") rather than returning NA.
      n_total <- nrow(df1)
      has1 <- rowSums(!is.na(df1)) > 0
      has2 <- rowSums(!is.na(df2)) > 0
      n_any <- sum(has1 | has2)
      keep  <- has1 & has2
      n_used <- sum(keep)
      if (n_used == 0)
        stop("No respondents have responses at both occasions.")
      if (n_used < 10)
        stop(paste0(
          "Only ", n_used, " respondent(s) have responses at both ",
          "occasions. Item parameters cannot be calibrated from so few ",
          "people. A single-respondent analysis needs item parameters from ",
          "an external calibration, which this jamovi analysis cannot ",
          "accept; use easyRasch2::RMpersonChange() in R for that case."
        ))

      d1 <- df1[keep, , drop = FALSE]
      d2 <- df2[keep, , drop = FALSE]

      # RMpersonChange() requires the two occasions to carry the same item
      # names in the same order, which two sets of jamovi variables never
      # do. Both frames are given the Time 1 names: the pairing is by
      # position either way, and the pairing table above records which
      # Time 2 variable each name now stands for.
      names(d1) <- v1
      names(d2) <- v1

      ids <- if (!is.null(self$options$id)) {
        as.character(data[[self$options$id]])[keep]
      } else {
        as.character(seq_len(n_total))[keep]
      }

      sparse_msg <- sparse_note(df_all)
      if (!is.null(sparse_msg))
        self$results$summaryTable$setNote("sparse", sparse_msg)

      # --- Options ----------------------------------------------------------
      null_type <- self$options$nullType
      retest_sd <- if (null_type == "retest") self$options$retestSd else NULL
      theta_range <- c(self$options$thetaMin, self$options$thetaMax)
      if (theta_range[1L] >= theta_range[2L])
        stop("Theta lower bound must be less than the upper bound.")

      tryCatch({
        # --- Compute, or reuse the cache -----------------------------------
        # The enumeration is the expensive part and does not depend on the
        # output-variable toggles or the table switches, which jamovi
        # nevertheless reruns .run for. The figure is cached alongside the
        # result because RMpersonChange() offers no output mode returning
        # both, so a cold run costs two enumerations and a warm one none.
        sig <- list(
          anchor    = self$options$anchor,
          method    = self$options$method,
          null      = null_type,
          retest_sd = retest_sd,
          alpha     = self$options$alpha,
          direction = self$options$direction,
          cond      = isTRUE(self$options$conditionalCrit),
          theta     = theta_range,
          ids       = ids
        )
        cached <- self$results$changeCache$state
        if (!is.null(cached) && identical(cached$sig, sig) &&
            identical(cached$d1, d1) && identical(cached$d2, d2)) {
          res  <- cached$res
          plot <- cached$plot
        } else {
          call_args <- list(
            data_t1          = d1,
            data_t2          = d2,
            id               = ids,
            anchor           = self$options$anchor,
            method           = self$options$method,
            estimator        = "CML",
            null             = null_type,
            retest_sd        = retest_sd,
            critical         = "exact",
            alpha            = self$options$alpha,
            direction        = self$options$direction,
            conditional_crit = isTRUE(self$options$conditionalCrit),
            parallel         = FALSE,
            theta_range      = theta_range
          )
          res <- suppressWarnings(suppressMessages(
            do.call(easyRasch2::RMpersonChange,
                    c(call_args, list(output = "dataframe")))
          ))
          plot <- suppressWarnings(suppressMessages(
            do.call(easyRasch2::RMpersonChange,
                    c(call_args, list(output = "ggplot")))
          ))
          self$results$changeCache$setState(list(
            res = res, plot = plot, sig = sig, d1 = d1, d2 = d2
          ))
        }
        if (nrow(res) != n_used)
          stop("Internal error: change rows do not match the data.")

        # --- Summary table --------------------------------------------------
        cls <- as.character(res$change_class)
        n_inc  <- sum(cls == "increase", na.rm = TRUE)
        n_dec  <- sum(cls == "decrease", na.rm = TRUE)
        n_none <- sum(cls == "none detected", na.rm = TRUE)
        pct <- function(x) paste0(sprintf("%.1f", 100 * x / n_used), "% of ",
                                  n_used)

        st <- self$results$summaryTable
        st$setRow(rowKey = "n", values = list(
          value = n_used,
          notes = if (n_used < n_total)
            paste0("of ", n_total, " rows; ", n_total - n_used,
                   " lacked responses at one or both occasions")
          else ""
        ))
        st$setRow(rowKey = "increase", values = list(
          value = n_inc, notes = pct(n_inc)))
        st$setRow(rowKey = "none", values = list(
          value = n_none, notes = pct(n_none)))
        st$setRow(rowKey = "decrease", values = list(
          value = n_dec, notes = pct(n_dec)))

        cond <- isTRUE(self$options$conditionalCrit)
        crit_row <- function(v, label) {
          finite_v <- v[is.finite(v)]
          if (length(finite_v) == 0L)
            return(list(value = NA_real_,
                        notes = paste0("Not used: ", label, " tail not tested.")))
          if (cond) {
            list(value = stats::median(finite_v),
                 notes = paste0("Median of the per-respondent values, ",
                                "range ", sprintf("%.2f", min(finite_v)),
                                " to ", sprintf("%.2f", max(finite_v)),
                                ". See the figure below."))
          } else {
            list(value = finite_v[1L],
                 notes = "Pooled across respondents.")
          }
        }
        st$setRow(rowKey = "critlo", values = crit_row(res$crit_lower, "lower"))
        st$setRow(rowKey = "crithi", values = crit_row(res$crit_upper, "upper"))

        st$setNote("null", private$.nullText(null_type, retest_sd))
        st$setNote("units", paste0(
          "Change is in logits. The change index is unitless, and is a ",
          "decision statistic rather than a measure of how much someone ",
          "changed: it mixes the size of the move with how precisely each ",
          "of the two positions was pinned down, so respondents must not ",
          "be ranked by it."
        ))

        # --- Per-respondent table -------------------------------------------
        if (isTRUE(self$options$showTable)) {
          rows <- if (isTRUE(self$options$flaggedOnly)) {
            which(cls != "none detected")
          } else {
            seq_len(nrow(res))
          }
          ct <- self$results$changeTable
          fin <- function(x) if (is.finite(x)) x else NA_real_
          for (i in rows) {
            ct$addRow(rowKey = i, values = list(
              id          = as.character(res$id[i]),
              sumT1       = res$sum_t1[i],
              thetaT1     = res$theta_t1[i],
              sumT2       = res$sum_t2[i],
              thetaT2     = res$theta_t2[i],
              change      = res$change[i],
              seDiff      = res$se_diff[i],
              rci         = res$rci[i],
              pValue      = res$p_value[i],
              critLower   = fin(res$crit_lower[i]),
              critUpper   = fin(res$crit_upper[i]),
              changeClass = as.character(res$change_class[i]),
              retestTip   = res$retest_sd_tip[i]
            ))
          }
          if (isTRUE(self$options$flaggedOnly) && length(rows) == 0L) {
            ct$setNote("empty", paste0(
              "No respondent's change exceeded the critical value, so ",
              "there is nothing to list. Untick <i>Flagged respondents ",
              "only</i> to see every respondent."
            ))
          }
          ct$setNote("tip", paste0(
            "Retest SD tip: for a respondent whose change was detected, ",
            "the per-occasion retest SD that would bring their change ",
            "index back to the critical value. Empty for respondents ",
            "whose change was not detected, for whom the question does ",
            "not arise. A small value means the result would not survive ",
            "much occasion-to-occasion noise."
          ))
          if (self$options$direction != "two.sided") {
            ct$setNote("onesided", paste0(
              "One tail is not tested, so its critical value is empty."
            ))
          }
        }

        # --- Retest SD estimated from these data ----------------------------
        if (isTRUE(self$options$estimateRetestSd))
          private$.fillRetest(d1, d2, theta_range)

        # --- Output variables ------------------------------------------------
        expand <- function(v, fill = NA_real_) {
          out <- rep(fill, n_total)
          out[keep] <- v
          out
        }
        row_nums <- rownames(df1)
        set_output <- function(opt_name, values) {
          if (!isTRUE(self$options[[opt_name]])) return()
          out <- self$results[[opt_name]]
          out$setRowNums(row_nums)
          out$setValues(values)
        }
        set_output("outputThetaT1", expand(res$theta_t1))
        set_output("outputThetaT2", expand(res$theta_t2))
        set_output("outputChange",  expand(res$change))
        set_output("outputSeDiff",  expand(res$se_diff))
        set_output("outputRci",     expand(res$rci))
        set_output("outputPValue",  expand(res$p_value))
        if (isTRUE(self$options$outputClass)) {
          out <- self$results$outputClass
          out$setRowNums(row_nums)
          vals <- rep(NA_character_, n_total)
          vals[keep] <- cls
          out$setValues(vals)
        }

        # --- Note --------------------------------------------------------
        drop_clause <- if (n_used < n_total) {
          paste0(" ", n_total - n_used, " row(s) had no responses at one or ",
                 "both occasions and could not be tested",
                 if (n_any > n_used)
                   paste0(", ", n_any - n_used,
                          " of them having answered at one occasion only")
                 else "",
                 ".")
        } else ""
        anchor_clause <- switch(self$options$anchor,
          stack = paste0(
            "Item parameters are calibrated on both occasions stacked ",
            "together, which puts them on one metric by construction and ",
            "assumes the items behave the same way at both occasions. Test ",
            "that first with one of the DIF analyses, using occasion as the ",
            "grouping variable."
          ),
          t1 = paste0(
            "Item parameters are calibrated on occasion 1 only, so ",
            "occasion 2 is measured on the baseline metric and no ",
            "follow-up data enters the calibration."
          ),
          t2 = paste0(
            "Item parameters are calibrated on occasion 2 only, the mirror ",
            "of the baseline-only calibration."
          )
        )
        self$results$changeNote$setContent(paste0(
          "<p>Change tested for n = ", n_used, " respondents across ", k,
          " item pairs.", drop_clause, " ", anchor_clause, " ",
          private$.nullText(null_type, retest_sd),
          " Critical values are obtained by enumerating the null exactly, ",
          "so this analysis involves no simulation, no iteration count and ",
          "no random seed, and repeats identically. ",
          "A detected change is a necessary condition for real change, not ",
          "evidence of it. Read <i>change</i> for magnitude, in logits, and ",
          "the change index only for the decision. Results are identical to ",
          "easyRasch2::RMpersonChange().</p>"
        ))

        # --- Plot state ---------------------------------------------------
        self$results$changePlot$setState(list(ready = TRUE))
        if (cond) {
          self$results$critPlot$setState(list(
            res = res[, c("theta_t1", "se_t1", "theta_t2", "se_t2",
                          "crit_lower", "crit_upper")],
            key = paste(private$.answeredKey(d1), private$.answeredKey(d2),
                        sep = "|"),
            pooled = private$.pooledCrit(d1, d2, ids, theta_range, retest_sd)
          ))
        }

      }, error = function(e) {
        stop(paste("Error in person change analysis:", e$message))
      })
    },

    # ---------------------------------------------------------------------
    # Prose for the null in force. Used in both the table footnote and the
    # Html note, so it lives in one place.
    # ---------------------------------------------------------------------
    .nullText = function(null_type, retest_sd) {
      if (null_type == "measurement") {
        paste0(
          "Null: no change beyond measurement error, with the standard ",
          "error of the change the root of SE1² + SE2². This is ",
          "the narrower of the two nulls and makes no allowance for ",
          "occasion-to-occasion fluctuation that is not change in the ",
          "construct, so it detects more change than a retest-based null ",
          "would."
        )
      } else {
        paste0(
          "Null: no change beyond measurement error and occasion-to-",
          "occasion fluctuation, with the standard error of the change the ",
          "root of SE1² + SE2² + 2 × ",
          sprintf("%.3f", retest_sd), "². The retest SD is a ",
          "per-occasion quantity, which is why it enters twice."
        )
      }
    },

    # ---------------------------------------------------------------------
    # Answered-item pattern, one string per respondent. The critical value
    # is a deterministic function of the null location given which items
    # were answered, so respondents with different patterns lie on
    # different curves and must not be joined.
    # ---------------------------------------------------------------------
    .answeredKey = function(d) {
      apply(!is.na(as.matrix(d)), 1L,
            function(z) paste0(which(z), collapse = ","))
    },

    # ---------------------------------------------------------------------
    # The pooled critical value, for the reference line on the conditional
    # figure. A second enumeration, run only when the conditional figure is
    # shown and only for its reference line.
    # ---------------------------------------------------------------------
    .pooledCrit = function(d1, d2, ids, theta_range, retest_sd) {
      out <- try(suppressWarnings(suppressMessages(
        easyRasch2::RMpersonChange(
          d1, d2, id = ids,
          anchor           = self$options$anchor,
          method           = self$options$method,
          null             = self$options$nullType,
          retest_sd        = retest_sd,
          critical         = "exact",
          alpha            = self$options$alpha,
          direction        = self$options$direction,
          conditional_crit = FALSE,
          parallel         = FALSE,
          theta_range      = theta_range,
          output           = "dataframe"
        )
      )), silent = TRUE)
      if (inherits(out, "try-error")) return(NA_real_)
      v <- out$crit_upper[1L]
      if (is.finite(v)) v else NA_real_
    },

    # ---------------------------------------------------------------------
    # Retest SD estimated from the two occasions. Only defensible when the
    # data are a stability study; the UI says so and so does the footnote.
    # ---------------------------------------------------------------------
    .fillRetest = function(d1, d2, theta_range) {
      rt <- self$results$retestTable
      res <- try(suppressWarnings(suppressMessages(
        easyRasch2::RMretestSD(
          d1, d2,
          anchor      = self$options$anchor,
          method      = self$options$method,
          estimator   = "CML",
          sim_iter    = self$options$retestIter,
          parallel    = FALSE,
          boot        = isTRUE(self$options$retestBoot),
          boot_iter   = self$options$retestBootIter,
          conf_int    = self$options$confInt / 100,
          seed        = as.integer(self$options$seed),
          theta_range = theta_range,
          output      = "dataframe"
        )
      )), silent = TRUE)
      if (inherits(res, "try-error")) {
        rt$setNote("failed", paste0(
          "The retest SD could not be estimated: ",
          conditionMessage(attr(res, "condition"))
        ))
        return()
      }

      rows <- list(
        list("Variance of the change", res$var_change[1L], NA, NA),
        list("Variance from measurement error", res$var_error[1L], NA, NA),
        list("Occasion variance", res$variance[1L], NA, NA),
        list("Retest SD per occasion", res$sd[1L],
             res$lower[1L], res$upper[1L])
      )
      for (i in seq_along(rows)) {
        rt$addRow(rowKey = i, values = list(
          quantity = rows[[i]][[1L]],
          value    = rows[[i]][[2L]],
          lower    = rows[[i]][[3L]],
          upper    = rows[[i]][[4L]]
        ))
      }

      negative <- is.finite(res$variance[1L]) && res$variance[1L] < 0
      rt$setNote("use", paste0(
        "The occasion variance is the variance of the change minus the ",
        "part attributable to measurement error, halved because each ",
        "occasion carries one. To use this estimate, set <i>Null ",
        "hypothesis</i> to include occasion fluctuation and enter the ",
        "retest SD above. It is not applied automatically, so that the ",
        "value in force is always the one shown in the options."
      ))
      rt$setNote("valid", paste0(
        "Valid only if these two occasions are a stability study in which ",
        "no real change was expected. Estimated from data in which people ",
        "did change, this figure absorbs that change and the test then ",
        "detects almost nothing."
      ))
      if (negative) {
        rt$setNote("negative", paste0(
          "The occasion variance came out negative, which happens when ",
          "the observed change varies by less than measurement error ",
          "alone predicts. It is reported rather than set to zero. Treat ",
          "it as no detectable occasion fluctuation, estimated imprecisely."
        ))
      }
    },

    # ---------------------------------------------------------------------
    # Occasion 1 against occasion 2, with the no-change band, drawn by
    # easyRasch2::RMpersonChange(). Read from the cache rather than
    # recomputed: the enumeration behind it is the expensive part.
    # ---------------------------------------------------------------------
    .changePlot = function(image, ggtheme, theme, ...) {
      cached <- self$results$changeCache$state
      if (is.null(cached) || is.null(cached$plot)) return(FALSE)

      p <- er2_bump_text(cached$plot)
      print(p)
      TRUE
    },

    # ---------------------------------------------------------------------
    # Per-respondent critical values against the respondent's null
    # location, shown only when critical values are computed per
    # respondent. A curve rather than a histogram, because the critical
    # value is a deterministic function of that location given which items
    # were answered: a histogram would show the same curve marginalised
    # over wherever this sample happens to sit, with no way to tell the
    # items' behaviour from the sample's targeting.
    # ---------------------------------------------------------------------
    .critPlot = function(image, ggtheme, theme, ...) {
      if (is.null(image$state)) return(FALSE)
      if (!requireNamespace("ggplot2", quietly = TRUE)) return(FALSE)

      s <- image$state
      r <- s$res
      w1 <- 1 / r$se_t1^2
      w2 <- 1 / r$se_t2^2

      d <- data.frame(
        theta_null = (r$theta_t1 * w1 + r$theta_t2 * w2) / (w1 + w2),
        upper      = r$crit_upper,
        lower      = abs(r$crit_lower),
        pattern    = s$key,
        stringsAsFactors = FALSE
      )
      d <- d[is.finite(d$theta_null), , drop = FALSE]
      if (nrow(d) == 0L) return(FALSE)
      d <- d[order(d$theta_null), , drop = FALSE]

      keep_u <- is.finite(d$upper)
      keep_l <- is.finite(d$lower)
      long <- rbind(
        data.frame(d[keep_u, c("theta_null", "pattern")],
                   value = d$upper[keep_u], bound = "Upper"),
        data.frame(d[keep_l, c("theta_null", "pattern")],
                   value = d$lower[keep_l], bound = "Lower (absolute)")
      )
      if (nrow(long) == 0L) return(FALSE)
      long$bound <- factor(long$bound,
                           levels = c("Upper", "Lower (absolute)"))
      long$series <- paste(long$bound, long$pattern)

      n_pattern <- length(unique(d$pattern))
      symmetric <- sum(keep_u) == sum(keep_l) &&
        isTRUE(all.equal(d$upper[keep_u], d$lower[keep_l], tolerance = 1e-8))

      p <- ggplot2::ggplot(
        long, ggplot2::aes(x = .data$theta_null, y = .data$value)
      )

      if (is.finite(s$pooled)) {
        p <- p +
          ggplot2::geom_hline(yintercept = s$pooled, linetype = "dotted",
                              colour = "grey35", linewidth = 0.6) +
          ggplot2::annotate(
            "text", x = max(d$theta_null), y = s$pooled,
            label = paste0("Pooled  ", sprintf("%.2f", s$pooled)),
            hjust = 1, vjust = -0.7, size = 3.6, colour = "grey35"
          )
      }

      # Step lines only while there are few enough curves to tell apart.
      # Ragged missingness can produce nearly one pattern per respondent,
      # and a hundred overlaid step curves is worse than the points alone.
      if (n_pattern <= 3L) {
        p <- p +
          ggplot2::geom_step(
            ggplot2::aes(colour = .data$bound, linetype = .data$bound,
                         group = .data$series),
            linewidth = 0.8, alpha = 0.9
          ) +
          ggplot2::scale_linetype_manual(
            values = c(Upper = "solid", `Lower (absolute)` = "22")
          )
      }

      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(colour = .data$bound, shape = .data$bound),
          size = 1.7, alpha = 0.75
        ) +
        ggplot2::geom_rug(
          data = d, ggplot2::aes(x = .data$theta_null), inherit.aes = FALSE,
          sides = "b", alpha = 0.25, colour = "grey25"
        ) +
        ggplot2::scale_colour_manual(
          values = c(Upper = "#0072B2", `Lower (absolute)` = "#D55E00")
        ) +
        ggplot2::scale_shape_manual(
          values = c(Upper = 16, `Lower (absolute)` = 1)
        ) +
        ggplot2::labs(
          x = "Respondent's null location (logits)",
          y = "Critical value for the change index (unitless)",
          colour = NULL, linetype = NULL, shape = NULL,
          caption = er2_caption(paste0(
            "Each respondent's own critical value, from enumerating the ",
            "null at their location, the precision-weighted mean of their ",
            "two occasion estimates. Values pull in toward the ends of the ",
            "scale, where fewer scores remain attainable. ",
            if (n_pattern > 3L) {
              paste0(n_pattern, " distinct pairs of answered-item sets, ",
                     "too many to join into curves, so points only. ")
            } else if (n_pattern > 1L) {
              paste0("One step curve per pair of answered-item sets (",
                     n_pattern, " here). ")
            } else {
              ""
            },
            if (symmetric) {
              "The two bounds coincide while the occasions stay exchangeable. "
            } else {
              "The bounds part where the occasions differ in items answered. "
            },
            "Rug marks show where respondents sit. Dotted line: the single ",
            "value that would apply to everyone with per-respondent ",
            "critical values turned off."
          ))
        ) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(legend.position = "top") +
        er2_axis_margins() +
        er2_plot_caption()

      print(p)
      TRUE
    }
  )
)
