# Internal helpers for input validation

#' Validate response data for Rasch analysis
#'
#' Checks that the supplied object is a data.frame or matrix, that all values
#' are non-negative integers (or `NA`), and that the minimum non-missing value
#' is 0 (i.e., items are scored starting at 0).
#'
#' @param data A data.frame or matrix of item responses.
#'
#' @return Invisibly returns `TRUE` if validation passes; otherwise stops with
#'   an informative error message.
#'
#' @noRd
validate_response_data <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }

  vals <- as.vector(as.matrix(data))
  vals_nonmissing <- vals[!is.na(vals)]

  if (length(vals_nonmissing) == 0L) {
    stop("`data` contains no non-missing values.", call. = FALSE)
  }

  if (any(vals_nonmissing < 0)) {
    stop("`data` contains negative values. Item responses must be non-negative integers.",
         call. = FALSE)
  }

  if (!all(vals_nonmissing == floor(vals_nonmissing))) {
    stop("`data` contains non-integer values. Item responses must be integers.",
         call. = FALSE)
  }

  min_val <- min(vals_nonmissing)
  if (min_val != 0) {
    stop(
      paste0(
        "The minimum value in `data` is ", min_val, ", but Rasch analysis ",
        "requires items scored starting at 0. Please recode your data."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Convert a single item vector to numeric responses
#'
#' Robustly converts ordinal/nominal data — including factors with text value
#' labels (e.g., from SPSS imports) and `haven_labelled` vectors — back to the
#' underlying numeric response codes that Rasch analysis expects.
#'
#' Handling:
#' * **`haven_labelled`** (or any vector with a `labels` attribute): the
#'   underlying numeric values are used. SPSS user-defined missing values
#'   recorded in `na_values` / `na_range` attributes are converted to `NA`.
#' * **Factor with numeric-string levels** (e.g., levels `c("0","1","2","3")`,
#'   or `c("0","1","2","3","8888","999")` when the SPSS file declared
#'   missing-value codes): levels are parsed as numbers and the original
#'   numeric values are recovered. Levels not present in the data don't
#'   affect the result.
#' * **Factor with text-label levels** (e.g., `c("None","Mild","Moderate",
#'   "Severe")`): levels can't be parsed numerically, so 0-indexed factor
#'   positions are returned. This preserves the ordinal structure as long as
#'   levels are in the correct ordinal order — which is the case for ordinal
#'   variables in jamovi.
#' * **Numeric / integer**: returned via `as.numeric()`.
#'
#' Note: jamovi-side, the user should still flag values like 8888 or 999 as
#' missing in the data editor (or upstream in their SPSS file). This helper
#' picks them up automatically for `haven_labelled` inputs but cannot infer
#' which codes are "missing" from a plain factor.
#'
#' @param x A vector (factor, `haven_labelled`, or numeric).
#'
#' @return A numeric vector of the same length as `x`.
#'
#' @noRd
to_numeric_responses <- function(x) {

  # haven_labelled (haven / SPSS) — use underlying numeric values, drop
  # SPSS-flagged missing codes
  if (inherits(x, "haven_labelled") || !is.null(attr(x, "labels"))) {
    vals     <- as.numeric(unclass(x))
    na_vals  <- attr(x, "na_values")
    na_range <- attr(x, "na_range")
    if (!is.null(na_vals)) {
      vals[vals %in% na_vals] <- NA
    }
    if (!is.null(na_range) && length(na_range) == 2L) {
      vals[!is.na(vals) & vals >= na_range[1L] & vals <= na_range[2L]] <- NA
    }
    return(vals)
  }

  if (is.factor(x)) {
    levs     <- levels(x)
    num_levs <- suppressWarnings(as.numeric(levs))
    if (length(levs) > 0L && all(!is.na(num_levs))) {
      # Numeric-string levels — recover the actual numeric values
      return(num_levs[as.integer(x)])
    }
    # Text-label levels — fall back to 0-indexed factor position
    return(as.integer(x) - 1L)
  }

  as.numeric(x)
}

#' Convert a data.frame of items to numeric responses
#'
#' Applies [to_numeric_responses()] column-wise. Returns a data.frame with the
#' same column names and order, but every column coerced to numeric per the
#' rules in [to_numeric_responses()].
#'
#' @param df A data.frame of items.
#'
#' @return A data.frame with all columns numeric.
#'
#' @noRd
to_numeric_responses_df <- function(df) {
  for (col in names(df)) {
    df[[col]] <- to_numeric_responses(df[[col]])
  }
  df
}

#' Standard item-data preparation and validation
#'
#' Consolidates the validation block previously copy-pasted across all
#' analyses: numeric conversion, all-NA column check, sentinel-value
#' check (> 20 looks like an unmarked missing-value code), response
#' validation (non-negative integers starting at 0), per-item variation
#' check, and the identical-items check (correlation == 1). Complete-case
#' handling stays in each analysis because it legitimately varies
#' (retained by CML/MML, dropped for lavaan/prcomp/LRtest, joint with a
#' DIF variable, ...).
#'
#' @param data The jamovi data frame (`self$data`).
#' @param vars Character vector of selected item variable names.
#' @return The converted numeric data.frame of item responses.
#' @noRd
prepare_item_data <- function(data, vars) {
  # Robust conversion: handles factors with text labels (SPSS),
  # haven_labelled vectors, and numerics.
  df <- to_numeric_responses_df(data[, vars, drop = FALSE])

  # All-NA columns
  all_na_cols <- vapply(df, function(x) all(is.na(x)), logical(1L))
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

  # Per-item variation
  for (col in names(df)) {
    if (length(unique(stats::na.omit(df[[col]]))) < 2L)
      stop(paste0("Item '", col, "' has no variation in responses. ",
                  "Each item needs at least two different response values."))
  }

  # Identical-items check (pairwise correlation == 1)
  n_vars <- ncol(df)
  if (n_vars >= 2L) {
    identical_pairs <- list()
    for (i in 1:(n_vars - 1)) {
      for (j in (i + 1):n_vars) {
        cor_val <- suppressWarnings(
          stats::cor(df[[i]], df[[j]], use = "complete.obs")
        )
        if (!is.na(cor_val) && cor_val == 1) {
          identical_pairs <- append(identical_pairs,
                                    list(c(names(df)[i], names(df)[j])))
        }
      }
    }
    if (length(identical_pairs) > 0) {
      pair_strings <- vapply(identical_pairs, function(p) {
        paste0("'", p[1], "' and '", p[2], "'")
      }, character(1L))
      pair_msg <- paste(pair_strings, collapse = "; ")
      if (n_vars == 2L) {
        stop(paste("The two selected items are identical:", pair_msg,
                   "- please select different items."))
      } else {
        jmvcore::reject(
          "Warning: Some items appear to be identical ({pairs}). This may affect results.",
          pairs = pair_msg
        )
      }
    }
  }

  df
}

#' Iteration-count advice for simulation notes
#'
#' Returns the "more iterations recommended" sentence when the user is
#' running at or below the analysis default (someone who lowered the
#' count needs the advice even more); returns "" when they raised it.
#' `default_iterations` must be kept in sync with the a.yaml default.
#' The infit variant adds the small-sample exception (detection power
#' can be better with ~100 iterations than with more; Johansson, 2025).
#'
#' @param iterations Current option value.
#' @param default_iterations The a.yaml default for this analysis.
#' @param infit Logical; use the conditional-infit variant.
#' @return Character scalar (possibly "") with a leading space.
#' @noRd
iteration_note <- function(iterations, default_iterations, infit = FALSE) {
  if (iterations > default_iterations) return("")
  if (infit) {
    paste0(" More iterations are generally recommended for ",
           "publication-ready results; however, for conditional infit ",
           "with small samples, around 100 iterations can yield better ",
           "detection power than a larger number (Johansson, 2025).")
  } else {
    paste0(" More iterations are generally recommended for ",
           "publication-ready results.")
  }
}
