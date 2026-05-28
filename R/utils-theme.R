# Internal theme helpers shared across easyRasch2jmv plotting backends.
# Mirrors the helpers in the easyRasch2 R package so jamovi-rendered
# plots look the same as the R-package plots.

#' Standard easyRasch2jmv axis-title margins
#'
#' Adds a little breathing room around the x and y axis titles. Apply
#' to any ggplot via `p + er2_axis_margins()`.
#'
#' @return A `ggplot2::theme()` object.
#' @noRd
er2_axis_margins <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
  )
}

#' Italic plot.caption theme element
#'
#' Returns a `ggplot2::theme()` setting `plot.caption` to render
#' left-aligned at 9 pt. Pair with
#' \code{\link{er2_caption}} when building the caption text so it
#' starts with a "Note. " prefix and wraps at a reasonable line width.
#'
#' @return A `ggplot2::theme()` object.
#' @noRd
er2_plot_caption <- function() {
  ggplot2::theme(
    plot.caption = ggplot2::element_text(hjust = 0, size = 10)
  )
}

#' "Note." caption-text prefix with line wrapping
#'
#' Build a plot caption with the standard "Note. " prefix and wrap it
#' at `width` characters so long captions don't run off the right edge
#' of the plot. Returned with `\\n` line breaks, which
#' `ggplot2::element_text()` renders as a multi-line caption.
#'
#' @param text Character. The caption body text (everything after the
#'   "Note. " prefix).
#' @param width Integer. Maximum characters per line; passed to
#'   [strwrap()]. Default `90`.
#' @return A character string ready to pass to
#'   `ggplot2::labs(caption = ...)`.
#' @noRd
er2_caption <- function(text, width = 90L) {
  prefixed <- paste("Note.", text)
  paste(strwrap(prefixed, width = width), collapse = "\n")
}
