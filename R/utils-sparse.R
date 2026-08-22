# ---------------------------------------------------------------------------
# Shared sparse-response-category checks.
#
# Mirrors the trigger used by easyRasch2::RMtargeting()
# (.has_sparse_categories): for each item, response counts are tabulated
# over 0..max and any category with fewer than `min_n` observations
# (including empty intermediate categories) marks the item as sparse.
# Used (a) by the targeting analysis to switch from CML to MML threshold
# estimation, and (b) by all item-based analyses to warn users via a
# table footnote.
# ---------------------------------------------------------------------------

#' Items with at least one sparse response category
#'
#' @param df data.frame of numeric item responses (0-indexed).
#' @param min_n Minimum per-category count.
#' @param max_cats Optional named integer vector giving the top category
#'   per item. Used by the grouped variant so that a subgroup that never
#'   reaches the top category is still flagged (its within-group max
#'   would otherwise hide the empty category).
#' @return Character vector of item names (possibly empty).
#' @noRd
sparse_category_items <- function(df, min_n = 3L, max_cats = NULL) {
  out <- character(0)
  for (j in seq_len(ncol(df))) {
    x <- df[[j]]
    x <- x[!is.na(x)]
    if (length(x) == 0L) next
    top <- if (!is.null(max_cats)) max_cats[[names(df)[j]]] else max(x)
    counts <- tabulate(x + 1L, nbins = top + 1L)
    if (any(counts < min_n)) out <- c(out, names(df)[j])
  }
  out
}

#' Footnote text for sparse categories, or NULL when none
#' @noRd
sparse_note <- function(df, min_n = 3L) {
  items <- sparse_category_items(df, min_n = min_n)
  if (length(items) == 0L) return(NULL)
  paste0(
    "Item(s) ", paste(items, collapse = ", "), " have at least one ",
    "response category with fewer than ", min_n, " observations. ",
    "Estimates may be unstable."
  )
}

#' Footnote text for the 1-based recode, or NULL when none was applied
#'
#' Rasch analysis needs items scored from 0, and jamovi users often bring
#' data coded from 1, so [prepare_item_data()] applies the shift rather than
#' refusing the data. It is silent data handling, so every analysis states it.
#'
#' Takes the *raw* input and re-derives the condition rather than reading a
#' flag off the prepared data.frame. An attribute would not survive: several
#' analyses rebuild or subset the frame immediately after preparing it (see
#' the `as.data.frame(matrix(...))` in `cicc`), which drops attributes and
#' would silently suppress the note in exactly the analyses that need it.
#'
#' @param data The raw analysis data.
#' @param vars Item variable names, as passed to [prepare_item_data()].
#' @return Character scalar, or NULL when no recode happened.
#' @noRd
recode_note <- function(data, vars) {
  converted <- try(
    to_numeric_responses_df(data[, vars, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(converted, "try-error") || !is_one_based(converted)) {
    return(NULL)
  }
  paste0(
    "Every item's lowest response was 1, so 1 was subtracted from all ",
    "items to give the 0-based scoring that Rasch analysis requires. ",
    "This does not change the results, but response categories are ",
    "reported from 0, so category 1 in your data appears as 0 here."
  )
}

#' Footnote text for sparse categories within DIF groups, or NULL
#'
#' Category ranges are fixed at the whole-sample per-item maximum, so a
#' group that never uses the top category is flagged as well.
#' @noRd
sparse_note_grouped <- function(df, group, min_n = 3L) {
  group <- droplevels(as.factor(group))
  max_cats <- vapply(df, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) 0L else as.integer(max(x))
  }, integer(1L))
  names(max_cats) <- names(df)

  out <- character(0)
  for (g in levels(group)) {
    sub <- df[!is.na(group) & group == g, , drop = FALSE]
    items <- sparse_category_items(sub, min_n = min_n, max_cats = max_cats)
    if (length(items) > 0L) {
      out <- c(out, paste0(items, " (", g, ")"))
    }
  }
  if (length(out) == 0L) return(NULL)
  paste0(
    "When split by the DIF variable, the following items have at least ",
    "one response category with fewer than ", min_n, " observations ",
    "(group in parentheses): ", paste(out, collapse = ", "),
    ". Estimates may be unstable."
  )
}
