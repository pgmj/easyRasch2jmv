# Shared test helpers.
#
# Analysis wrapper functions and the internal validation/note helpers are
# reached through the package namespace rather than by bare name: the
# module's NAMESPACE only exports a subset of analyses, but every object
# is present in the loaded namespace under devtools::test() / R CMD check.
er2 <- asNamespace("easyRasch2jmv")

# Locate a bundled example dataset whether tests run from tests/testthat/
# (devtools::test) or the package root.
demo_path <- function(file) {
  p <- testthat::test_path("..", "..", "data", file)
  if (!file.exists(p)) p <- file.path("data", file)
  p
}

read_demo <- function(file) {
  p <- demo_path(file)
  if (!file.exists(p)) testthat::skip(paste("bundled dataset not found:", file))
  utils::read.csv(p)
}

# Convenience loaders for the three bundled datasets.
poly_data  <- function() read_demo("eRm_pcmdat2.csv")      # polytomous (PCM)
dich_data  <- function() read_demo("eRm_raschdat3.csv")    # dichotomous (RM)
dif_data   <- function() {
  d <- read_demo("eRm_pcmdat2dif.csv")
  d$dif <- factor(d$dif)
  d
}
dif_items  <- function() c("I1", "I2", "I3", "I4")

# Append a perfectly-correlated duplicate of the first item.
with_duplicate_item <- function(df) {
  df$Dup <- df[[1L]]
  df
}

# Inject NA into the first `cols` columns so MI analyses have something
# to impute (deterministic).
with_missing <- function(df, cols = 3L, n = 4L, seed = 6L) {
  set.seed(seed)
  for (j in seq_len(cols)) df[sample(nrow(df), n), j] <- NA
  df
}
