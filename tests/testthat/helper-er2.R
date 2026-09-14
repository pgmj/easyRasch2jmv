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

# Two-occasion data for the Person Change analysis. Simulated rather than
# taken from a bundled dataset, since none of them is longitudinal, and
# deterministic so the tests can assert on counts.
change_data <- function(n = 150, shift = 0.4, seed = 11L) {
  thr <- list(c(-1.5, -0.3, 1.0), c(-1, 0.1, 1.4), c(-0.6, 0.5, 1.8),
              c(-0.2, 0.9, 2.2), c(0.3, 1.4, 2.7), c(-2, -0.8, 0.6))
  sim <- function(theta) {
    out <- matrix(NA_integer_, length(theta), length(thr))
    for (j in seq_along(thr)) {
      d <- thr[[j]]
      eta <- cbind(0, t(vapply(theta, function(th) cumsum(th - d),
                               numeric(length(d)))))
      p <- exp(eta - apply(eta, 1, max))
      p <- p / rowSums(p)
      out[, j] <- apply(p, 1, function(pr)
        sample.int(length(pr), 1L, prob = pr)) - 1L
    }
    out
  }
  old <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  th <- stats::rnorm(n, 0, 1.3)
  df <- as.data.frame(cbind(sim(th), sim(th + shift)))
  if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
  names(df) <- c(paste0("T1_i", 1:6), paste0("T2_i", 1:6))
  df$pid <- paste0("P", seq_len(n))
  df
}

change_t1 <- function() paste0("T1_i", 1:6)
change_t2 <- function() paste0("T2_i", 1:6)

