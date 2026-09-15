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

#' Left-aligned plot.caption theme element
#'
#' Returns a `ggplot2::theme()` setting `plot.caption` to render
#' left-aligned at `er2_caption_size()`. Pair with
#' \code{\link{er2_caption}} when building the caption text so it
#' starts with a "Note. " prefix and wraps at a reasonable line width.
#'
#' @return A `ggplot2::theme()` object.
#' @noRd
er2_plot_caption <- function() {
  ggplot2::theme(plot.caption = er2_caption_element())
}

#' Caption text size, in points
#'
#' One place, so module-drawn captions and the captions on package-drawn
#' figures come out the same size. Raised in 3.2.0 for readability on
#' jamovi's fixed-size canvases: module-drawn captions were 10 pt and
#' package-drawn ones 9 pt, and both are now 10.5.
#'
#' @return Numeric scalar.
#' @noRd
er2_caption_size <- function() 10.5

#' Caption theme element, matching the class easyRasch2 uses
#'
#' ggplot2 refuses to merge two theme elements of different classes, so a
#' caption element handed to a package-drawn figure has to be the same class
#' as the one already there. easyRasch2 uses `ggtext::element_markdown()` when
#' `ggtext` is installed, so that its italic "*Note.*" prefix can sit beside
#' roman body text, and a plain `element_text()` otherwise. This follows the
#' same test, which is why `ggtext` is declared in `Imports:`: without it in
#' the bundle the two would disagree and the size change would not apply.
#'
#' @return A ggplot2 theme element.
#' @noRd
er2_caption_element <- function() {
  if (requireNamespace("ggtext", quietly = TRUE)) {
    ggtext::element_markdown(hjust = 0, size = er2_caption_size())
  } else {
    ggplot2::element_text(hjust = 0, size = er2_caption_size())
  }
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

#' Enlarge text on a package-drawn plot to the module's base size
#'
#' Plots returned by easyRasch2 use the package's default text sizes
#' (ggplot2 base size 11-13), which are too small on jamovi's fixed-size
#' canvases; module-drawn plots use base size 15. This bumps the theme's
#' root `text` element additively, so every rel()-sized element (axis
#' text, titles, strips, legends) scales along while each plot's own
#' specific theme settings -- rotated axis labels, blanked grids, the
#' 10 pt caption from `er2_plot_caption()` -- survive untouched.
#' patchwork objects get the theme applied to every panel via `&`.
#'
#' The caption is the exception: package figures set `plot.caption` to their
#' own absolute size, which an additive bump to the root `text` element does
#' not reach, so it is restated here at `er2_caption_size()`. Without that a
#' package-drawn figure and a module-drawn one would disagree about caption
#' size on the same results panel.
#'
#' @param p A `ggplot` or `patchwork` object.
#' @param size Base text size in points. Default 15 (13 suits dense
#'   multi-facet grids, matching the module's previous faceted plots).
#' @return The re-themed plot object.
#' @noRd
er2_bump_text <- function(p, size = 15) {
  th <- ggplot2::theme(
    text = ggplot2::element_text(size = size),
    plot.caption = er2_caption_element()
  )
  if (inherits(p, "patchwork")) {
    # Panels that deliberately blank their caption cannot take a text
    # element, so the caption size is applied to the assembled plot and the
    # per-panel pass carries the text bump alone.
    out <- p & ggplot2::theme(text = ggplot2::element_text(size = size))
    return(out + ggplot2::theme(plot.caption = er2_caption_element()))
  }
  p + th
}

#' Wrap long axis / facet labels onto multiple lines
#'
#' Base R only -- no stringr dependency. Width chosen so labels like
#' "Elementary school" wrap to two lines without crowding. Used for
#' DIF-group labels in the LR-DIF and partial gamma DIF plots.
#'
#' @param x Character vector (or factor) of labels.
#' @param width Integer. Maximum characters per line; passed to
#'   [strwrap()].
#' @return Character vector with `\n` line breaks.
#' @noRd
er2_wrap_labels <- function(x, width = 10L) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1L))
}

#' Build a figure into a grob for storage in results state
#'
#' A `ggplot` object is far larger than the figure it describes. In ggplot2
#' 4.0.3 an empty one serialises to 102 KB compressed, and the module's real
#' figures run from 287 KB to 752 KB each, because the object carries the
#' whole ggproto scaffolding rather than the drawn result. jmvcore warns
#' above 500 KB of compressed element state, and everything stored goes into
#' the saved `.omv`.
#'
#' Building the plot first collapses that: the same figures come out at 12 to
#' 16 KB, a factor of about 25, and draw identically. Most analyses in the
#' module already store the data and draw in the render function, which is
#' cheaper still; this is for the three that cannot, because their figure
#' comes back from one expensive package call that also produces the numbers.
#'
#' The theme bump has to be applied before building, since a built grob can
#' no longer take a theme. That is fine: `er2_bump_text()` uses fixed sizes
#' and does not depend on the render context.
#'
#' @param p A `ggplot` or `patchwork` object.
#' @param size Text size passed to [er2_bump_text()].
#' @return A `gtable` / grob, drawn with [er2_draw_grob()].
#' @noRd
er2_plot_grob <- function(p, size = 15) {
  if (is.null(p)) return(NULL)
  p <- er2_bump_text(p, size = size)
  # Building a plot measures text, which needs a graphics device. .run()
  # has none, and without this ggplot2 opens the default one: on macOS an
  # empty Quartz window in front of the user, and a stray device left open
  # in the engine. A file-less pdf device stands in and is closed again.
  # Only when nothing is open, so a render function's own device is never
  # touched.
  if (is.null(grDevices::dev.list())) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  if (inherits(p, "patchwork")) patchwork::patchworkGrob(p)
  else ggplot2::ggplotGrob(p)
}

#' Draw a grob stored by er2_plot_grob()
#'
#' @param gt A grob, or NULL.
#' @return TRUE when something was drawn, FALSE otherwise, matching what a
#'   jamovi render function is expected to return.
#' @noRd
er2_draw_grob <- function(gt) {
  if (is.null(gt)) return(FALSE)
  grid::grid.newpage()
  grid::grid.draw(gt)
  TRUE
}
