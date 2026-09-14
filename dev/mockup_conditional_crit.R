# Mockup: per-respondent critical values across the latent scale
#
# Proposed extra figure for the planned Person Change analysis in
# easyRasch2jmv, shown only when conditional_crit = TRUE.
#
# The point of the figure: under conditional_crit the critical value is a
# deterministic function of the respondent's null location (and of which items
# they answered at each occasion), so plotting it against that location shows
# the mechanism. A histogram of the same numbers would show the sample's
# targeting mixed in with the items' behaviour, with no way to separate them.
#
# Everything the figure needs is in the public RMpersonChange() output, so no
# package internals are touched. theta_null is reconstructed the way
# easyRasch2:::.pc_theta_null() builds it, as the precision-weighted mean of
# the two occasion estimates.
#
# Requires easyRasch2 >= 1.3.0.

if (!"easyRasch2" %in% loadedNamespaces()) {
  library(easyRasch2)
}
library(ggplot2)

if (utils::packageVersion("easyRasch2") < "1.3.0") {
  stop("This mockup needs easyRasch2 1.3.0 or later.")
}


# --- Module theme helpers (copied from R/utils-theme.R so this runs alone) ---

er2_axis_margins <- function() {
  theme(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12))
  )
}

er2_plot_caption <- function() {
  theme(plot.caption = element_text(hjust = 0, size = 10))
}

er2_caption <- function(text, width = 90L) {
  paste(strwrap(paste("Note.", text), width = width), collapse = "\n")
}


# --- A six-item, four-category scale measured twice --------------------------

sim_pcm <- function(theta, thr_list) {
  out <- matrix(NA_integer_, length(theta), length(thr_list))
  for (j in seq_along(thr_list)) {
    d <- thr_list[[j]]
    eta <- cbind(0, t(vapply(theta, function(th) cumsum(th - d), numeric(length(d)))))
    p <- exp(eta - apply(eta, 1, max))
    p <- p / rowSums(p)
    out[, j] <- apply(p, 1, function(pr) sample.int(length(pr), 1L, prob = pr)) - 1L
  }
  colnames(out) <- names(thr_list)
  as.data.frame(out)
}

thr_list <- list(
  i1 = c(-2.0, -0.8, 0.6),
  i2 = c(-1.5, -0.3, 1.0),
  i3 = c(-1.0, 0.1, 1.4),
  i4 = c(-0.6, 0.5, 1.8),
  i5 = c(-0.2, 0.9, 2.2),
  i6 = c(0.3, 1.4, 2.7)
)

set.seed(1234)
n <- 250
theta_1 <- rnorm(n, mean = 0, sd = 1.4)
theta_2 <- theta_1 + rnorm(n, mean = 0.35, sd = 0.55)

d1 <- sim_pcm(theta_1, thr_list)
d2 <- sim_pcm(theta_2, thr_list)


# --- The figure --------------------------------------------------------------

#' Answered-item pattern key, one string per respondent per occasion
#'
#' The critical value is a deterministic function of the null location *given*
#' which items were answered at each occasion. Respondents who answered
#' different item sets therefore sit on different curves, and joining across
#' them draws a zigzag that looks like noise. Grouping on the pair of keys
#' separates them.
answered_key <- function(data) {
  apply(!is.na(as.matrix(data)), 1L, function(z) paste0(which(z), collapse = ","))
}

#' Build the conditional-critical-value figure from an RMpersonChange result
#'
#' @param res_cond RMpersonChange(..., conditional_crit = TRUE,
#'   output = "dataframe")
#' @param data_t1,data_t2 The two occasions, for the answered-item grouping.
#' @param res_pool The same call with conditional_crit = FALSE, for the pooled
#'   reference line. Optional.
crit_curve_plot <- function(res_cond, data_t1, data_t2, res_pool = NULL,
                            base_size = 15) {

  w1 <- 1 / res_cond$se_t1^2
  w2 <- 1 / res_cond$se_t2^2
  theta_null <- (res_cond$theta_t1 * w1 + res_cond$theta_t2 * w2) / (w1 + w2)

  d <- data.frame(
    theta_null = theta_null,
    upper = res_cond$crit_upper,
    lower = abs(res_cond$crit_lower),
    pattern = paste(answered_key(data_t1), answered_key(data_t2), sep = "|")
  )
  d <- d[is.finite(d$theta_null) & is.finite(d$upper), ]
  d <- d[order(d$theta_null), ]

  long <- rbind(
    data.frame(d[c("theta_null", "pattern")], value = d$upper, bound = "Upper"),
    data.frame(d[c("theta_null", "pattern")], value = d$lower,
               bound = "Lower (absolute)")
  )
  long$bound <- factor(long$bound, levels = c("Upper", "Lower (absolute)"))
  long$series <- paste(long$bound, long$pattern)

  symmetric <- isTRUE(all.equal(d$upper, d$lower, tolerance = 1e-8))
  n_pattern <- length(unique(d$pattern))

  pooled <- if (!is.null(res_pool)) res_pool$crit_upper[1L] else NA_real_

  p <- ggplot(long, aes(x = theta_null, y = value))

  if (is.finite(pooled)) {
    p <- p +
      geom_hline(yintercept = pooled, linetype = "dotted",
                 colour = "grey35", linewidth = 0.6) +
      annotate("text", x = max(d$theta_null), y = pooled,
               label = paste0("Pooled  ", sprintf("%.2f", pooled)),
               hjust = 1, vjust = -0.7, size = base_size / 4.5,
               colour = "grey35")
  }

  # Step lines only while there are few enough curves to tell apart. Ragged
  # missingness can produce almost one pattern per respondent, and a hundred
  # overlaid step curves is worse than the points alone.
  if (n_pattern <= 3L) {
    p <- p +
      geom_step(aes(colour = bound, linetype = bound, group = series),
                linewidth = 0.8, alpha = 0.9) +
      scale_linetype_manual(values = c(Upper = "solid",
                                       `Lower (absolute)` = "22"))
  }

  p +
    geom_point(aes(colour = bound, shape = bound), size = 1.7, alpha = 0.75) +
    geom_rug(
      data = d, aes(x = theta_null), inherit.aes = FALSE,
      sides = "b", alpha = 0.25, length = unit(0.025, "npc"), colour = "grey25"
    ) +
    scale_colour_manual(values = c(Upper = "#0072B2",
                                   `Lower (absolute)` = "#D55E00")) +
    scale_shape_manual(values = c(Upper = 16, `Lower (absolute)` = 1)) +
    labs(
      x = "Respondent's null location (logits)",
      y = "Critical value for the change index (unitless)",
      colour = NULL, linetype = NULL, shape = NULL,
      caption = er2_caption(paste0(
        "Each respondent's own critical value, from enumerating the null at ",
        "their location, the precision-weighted mean of their two occasion ",
        "estimates. Values pull in toward the ends of the scale, where fewer ",
        "scores remain attainable. ",
        if (n_pattern > 3L) {
          paste0(n_pattern, " distinct pairs of answered-item sets, too many ",
                 "to join into curves, so points only. ")
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
        "Rug marks show where respondents sit. Dotted line: the single value ",
        "reported for everyone when per-respondent critical values are off."
      ), width = 100L)
    ) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "top") +
    er2_axis_margins() +
    er2_plot_caption()
}


# --- Case 1, complete data at both occasions ---------------------------------

res_cond <- suppressMessages(RMpersonChange(
  d1, d2, conditional_crit = TRUE, output = "dataframe"
))
res_pool <- suppressMessages(RMpersonChange(
  d1, d2, conditional_crit = FALSE, output = "dataframe"
))

p_complete <- crit_curve_plot(res_cond, d1, d2, res_pool)
print(p_complete)


# --- Case 2, occasion 2 missing two items for a third of respondents ---------
# Shows the two bounds separating when exchangeability fails.

d2_miss <- d2
drop <- sample(seq_len(n), size = floor(n / 3))
d2_miss[drop, c("i5", "i6")] <- NA

res_cond_m <- suppressMessages(RMpersonChange(
  d1, d2_miss, conditional_crit = TRUE, output = "dataframe"
))
res_pool_m <- suppressMessages(RMpersonChange(
  d1, d2_miss, conditional_crit = FALSE, output = "dataframe"
))

p_missing <- crit_curve_plot(res_cond_m, d1, d2_miss, res_pool_m)
print(p_missing)


# --- What the summary table would carry --------------------------------------

crit_summary <- function(res) {
  q <- stats::quantile(res$crit_upper, c(0.25, 0.5, 0.75), na.rm = TRUE)
  data.frame(
    metric = "Critical value (upper)",
    median = unname(q[2]),
    iqr_low = unname(q[1]),
    iqr_high = unname(q[3]),
    min = min(res$crit_upper, na.rm = TRUE),
    max = max(res$crit_upper, na.rm = TRUE)
  )
}

print(crit_summary(res_cond))
print(crit_summary(res_cond_m))


# --- Case 3, ragged missingness, one random item missing per occasion --------
# The worst case for this figure: nearly one answered-item pattern per
# respondent, so the step lines are dropped and only points remain.

d1_rag <- d1
d2_rag <- d2
for (i in seq_len(n)) {
  if (runif(1) < 0.5) d1_rag[i, sample(seq_len(6), 1L)] <- NA
  if (runif(1) < 0.5) d2_rag[i, sample(seq_len(6), 1L)] <- NA
}

res_cond_r <- suppressMessages(RMpersonChange(
  d1_rag, d2_rag, conditional_crit = TRUE, output = "dataframe"
))
res_pool_r <- suppressMessages(RMpersonChange(
  d1_rag, d2_rag, conditional_crit = FALSE, output = "dataframe"
))

p_ragged <- crit_curve_plot(res_cond_r, d1_rag, d2_rag, res_pool_r)
print(p_ragged)
