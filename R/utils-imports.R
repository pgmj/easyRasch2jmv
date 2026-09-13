# Deliberate Imports: entries that carry no `pkg::` call of their own

#' Bundling pins for packages reached through easyRasch2
#'
#' `difR`, `mirt` and `partykit` are declared in `DESCRIPTION`'s `Imports:`
#' but never called from this package. They are not leftovers: jmvtools
#' bundles a module's `Imports:` and their recursive dependencies into the
#' `.jmo`, and ignores `Suggests:`. Every one of these is a hard runtime
#' requirement of an easyRasch2 function this module calls, so a jamovi user
#' with no personal R library must get them from the bundle.
#'
#' - `difR` and `partykit`: `easyRasch2::RMdifTree()` stops on entry without
#'   either (Mantel-Haenszel effect size and the tree itself). Both sit in
#'   easyRasch2's `Suggests:`, so `difR` reaches the bundle *only* through
#'   this declaration. Used by the Tree-Based DIF analysis.
#' - `mirt`: the MML fallback for sparse response categories in the Targeting
#'   plot, the RMU reliability estimate, and EAPsum person parameters.
#'
#' This function exists only so `R CMD check` sees the references and does not
#' report "All declared Imports should be used". It is never called.
#'
#' @noRd
ignore_unused_imports <- function() {
  difR::mantelHaenszel
  mirt::mirt
  partykit::ctree
}
