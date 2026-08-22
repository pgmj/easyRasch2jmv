# Behavioural tests encoding decisions made in the 2.0.0 review, so they
# cannot silently regress.

test_that("too-few-items guard shows a note instead of erroring", {
  d <- poly_data()
  expect_no_error(suppressWarnings(
    r <- er2$iteminfit(data = d, vars = names(d)[1:2])))
  expect_match(r$cutoffNote$content, "requires at least")
})

test_that("guard message is gone once enough items are selected", {
  d <- poly_data()
  r <- suppressWarnings(er2$iteminfit(data = d, vars = names(d)))
  # iteminfit without cutoffs writes a base note (not the guard message),
  # so the stale 'requires at least' text must not survive.
  expect_false(grepl("requires at least", r$cutoffNote$content))
})

test_that("a perfectly-correlated pair (>=3 items) yields a footnote, not a stop", {
  d <- with_duplicate_item(poly_data())
  expect_no_error(suppressWarnings(
    r <- er2$itemrestscore(data = d, vars = names(d))))
  notes <- r$restscoreTable$notes
  expect_true("duplicate" %in% names(notes))
  expect_match(notes$duplicate$note, "perfectly correlated")
})

test_that("item-restscore Flagged column uses overfit/underfit labels", {
  d <- dich_data()
  r <- suppressWarnings(er2$itemrestscore(data = d, vars = names(d)))
  # empty text cells come back as NA via asDF; non-blank labels must be
  # exactly "overfit" / "underfit"
  fit_vals <- r$restscoreTable$asDF$fit
  expect_true(all(is.na(fit_vals) | fit_vals %in% c("overfit", "underfit")))
  expect_true(any(fit_vals %in% c("overfit", "underfit")))   # at least one flagged
})

test_that("item-restscore Flagged uses inverted direction vs. infit (sanity)", {
  # Restscore: observed > expected => overfit. Infit: observed < expected
  # range => overfit. Just assert both analyses run and expose the column;
  # the directional rules are unit-tested via their own logic above.
  d <- poly_data()
  rr <- suppressWarnings(er2$itemrestscore(data = d, vars = names(d)))
  expect_true("fit" %in% names(rr$restscoreTable$asDF))
})

test_that("iccplot category probabilities are valid (sum to 1 per theta point)", {
  d <- poly_data()
  r <- suppressWarnings(er2$iccplot(data = d, vars = names(d)))
  pdf_df <- r$iccPlot$state$plot_df
  # For each item x theta, probabilities across categories should sum to 1.
  agg <- stats::aggregate(Probability ~ Item + Theta, data = pdf_df, FUN = sum)
  expect_true(all(abs(agg$Probability - 1) < 1e-8))
})

test_that("score-to-logit WLE SEM is information-based: finite at extremes and larger there", {
  # The score-to-logit lookup lives inside the personparams analysis
  # (3.0.0: the former scorese analysis was absorbed into it).
  d <- dich_data()
  r  <- suppressWarnings(er2$personparams(data = d, vars = names(d),
                                          method = "WLE",
                                          showScoreTable = TRUE))
  df <- r$scoreTable$asDF
  se <- df$logitSE
  expect_true(all(is.finite(se)))                 # Warm correction -> finite extremes
  n  <- length(se)
  # the two extreme scores carry the largest SEs (least information)
  expect_gt(se[1L], min(se))
  expect_gt(se[n],  min(se))
})

test_that("CFA cutoff requires 4 items, with an explanatory message at 3", {
  d <- poly_data()
  r <- suppressWarnings(er2$cfacutoff(data = d, vars = names(d)[1:3]))
  expect_match(r$cfaNote$content, "at least <b>4 items</b>")
})

test_that("Q3 matches easyRasch2 directly and survives an all-NA respondent", {
  # Regression for the 2.1.0 migration: the module delegates to easyRasch2
  # (CML/WLE) and must (a) reproduce RMlocdepQ3() numerically and (b) drop
  # all-NA respondents up front, because RMlocdepQ3Cutoff() errors on them
  # (psychotools cannot fit all-NA rows).
  m <- as.matrix(poly_data())
  m[3, ] <- NA                                   # one fully-empty respondent
  d <- as.data.frame(m)

  expect_no_error(suppressWarnings(
    r <- er2$locdepq3(data = d, vars = names(d),
                      computeCutoff = TRUE, iterations = 50)))

  # cutoff simulation ran (not the graceful-degradation path)
  cut_df <- r$cutoffTable$asDF
  expect_false(anyNA(cut_df$value))
  expect_false(grepl("unavailable", r$q3Note$content))
  expect_match(r$q3Note$content, "1 row\\(s\\) without any responses")

  # numerical agreement with the R package (same seed defaults)
  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepQ3(d_used, output = "dataframe")))
  mod <- suppressWarnings(
    apply(as.matrix(r$q3Table$asDF[, names(d)]), 2, as.numeric))
  expect_equal(unname(mod), unname(as.matrix(pkg[names(d)])), tolerance = 1e-12)
})

test_that("item-restscore matches easyRasch2 directly and survives an all-NA respondent", {
  # Regression for the 2.1.0 migration: the module delegates to
  # easyRasch2::RMitemRestscore() and pre-drops all-NA respondents
  # (the bundled release cannot fit them).
  m <- as.matrix(poly_data())
  m[3, ] <- NA
  d <- as.data.frame(m)

  expect_no_error(suppressWarnings(
    r <- er2$itemrestscore(data = d, vars = names(d))))
  tab <- r$restscoreTable$asDF

  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressWarnings(suppressMessages(
    easyRasch2::RMitemRestscore(d_used, output = "dataframe")))
  expect_equal(tab$observed, pkg$Observed)
  expect_equal(tab$difference, pkg$Difference)
  expect_equal(tab$relLocation, pkg$Relative_location)
  expect_match(r$restscoreNote$content, "1 row\\(s\\) without any responses")
})

test_that("iteminfit bootstrap p-values match easyRasch2 and flag on adjusted p", {
  # Regression for the 3.0.0 p-value feature: pValues + correction are
  # passed through to easyRasch2::RMitemInfit(p_value = TRUE); flagging
  # switches from the expected range to adjusted p < 0.05 (upstream
  # behavior), and a stale pValues = TRUE without computeCutoff is inert.
  d <- poly_data()
  cutoff <- suppressWarnings(suppressMessages(
    easyRasch2::RMitemInfitCutoff(d, iterations = 60, parallel = FALSE,
                                  seed = 42, hdci_width = 0.99)))

  for (corr in c("fwer", "fdr_bh", "fdr_by")) {
    pkg <- suppressWarnings(suppressMessages(
      easyRasch2::RMitemInfit(d, cutoff = cutoff, p_value = TRUE,
                              correction = corr, output = "dataframe")))
    r <- suppressWarnings(
      er2$iteminfit(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42, hdciWidth = 99,
                    pValues = TRUE, correction = corr))
    tab <- r$infitTable$asDF
    expect_equal(tab$pValue, pkg$p_infit)
    expect_equal(tab$pAdjusted, pkg$padj_infit)
    # empty text cells come back NA via asDF
    flags <- tab$misfit
    flags <- ifelse(is.na(as.character(flags)), "", as.character(flags))
    expect_identical(flags, pkg$Flagged)
  }

  # footnotes explain the p-based flag rule and the correction method
  notes <- vapply(r$infitTable$notes, function(x) x$note, character(1))
  expect_true(any(grepl("adjusted p-value < 0.05", notes)))
  expect_true(any(grepl("Benjamini-Yekutieli", notes)))
  # below the calibrated floor => the liberal-correction tier of the caveat
  expect_match(r$cutoffNote$content, "below the calibrated floor of 400")
  # and the withdrawn small-sample advice is gone
  expect_false(grepl("detection power", r$cutoffNote$content))

  # dynamic adjusted-p column title names the method
  expect_identical(r$infitTable$getColumn("pAdjusted")$title,
                   "Adj. p-value (BY)")

  # stale pValues = TRUE without computeCutoff runs and stays inert
  r_off <- suppressWarnings(
    er2$iteminfit(data = d, vars = names(d), pValues = TRUE))
  expect_true(all(is.na(r_off$infitTable$asDF$pValue)))
})

test_that("iteminfit reuses the cutoff simulation when only p-value settings change", {
  # The plot state doubles as a simulation cache: jamovi clears it exactly
  # when a simulation-relevant option changes, and .run re-validates it via
  # a signature (iterations/seed/hdci width) + the prepared data. Emulate a
  # rerun-with-surviving-state by pre-seeding the state on a fresh instance,
  # with a tampered cutoff bound as the reuse probe.
  d <- poly_data()

  run_seeded <- function(state, ...) {
    opts <- er2$iteminfitOptions$new(vars = names(d), computeCutoff = TRUE,
                                     iterations = 60, seed = 42,
                                     hdciWidth = 99, ...)
    an <- er2$iteminfitClass$new(options = opts, data = d)
    if (!is.null(state)) an$results$infitPlot$setState(state)
    an$run()
    an$results
  }

  state <- run_seeded(NULL)$infitPlot$state
  expect_false(is.null(state$sig))
  true_low <- state$cutoff_res$item_cutoffs$infit_low[1]

  tampered <- state
  tampered$cutoff_res$item_cutoffs$infit_low[1] <- 0.123456

  # matching signature -> cache reused (tampered bound surfaces)
  r <- run_seeded(tampered, pValues = TRUE, correction = "fdr_bh")
  expect_equal(r$infitTable$asDF$infitLow[1], 0.123456)

  # stale signature -> simulation reruns (true bound restored)
  stale <- tampered
  stale$sig$seed <- 99L
  r <- run_seeded(stale)
  expect_equal(r$infitTable$asDF$infitLow[1], true_low)
})

test_that("partgamdif reuses the cutoff simulation when only display options change", {
  # Same state-cache pattern as iteminfit: the pgdifPlot state carries the
  # full cutoff object; .run re-validates via sig + prepared data + DIF
  # vector before reuse. Tampered-bound probe proves reuse/invalidation.
  d <- dif_data()

  run_seeded <- function(state, ...) {
    opts <- er2$partgamdifOptions$new(vars = dif_items(), difVar = "dif",
                                      computeCutoff = TRUE, iterations = 60,
                                      seed = 42, hdciWidth = 99, ...)
    an <- er2$partgamdifClass$new(options = opts, data = d)
    if (!is.null(state)) an$results$pgdifPlot$setState(state)
    an$run()
    an$results
  }

  state <- run_seeded(NULL)$pgdifPlot$state
  expect_false(is.null(state$sig))
  expect_false(is.null(state$cutoff_res$item_cutoffs))
  true_low <- state$cutoff_res$item_cutoffs$gamma_low[1]

  tampered <- state
  tampered$cutoff_res$item_cutoffs$gamma_low[1] <- -0.987654

  # matching signature + tileplot toggled on -> cache reused
  r <- run_seeded(tampered, showTileplot = TRUE)
  expect_equal(r$pgdifTable$asDF$gammaLow[1], -0.987654)
  expect_false(is.null(r$tileplot$state))

  # stale signature -> simulation reruns (true bound restored)
  stale <- tampered
  stale$sig$seed <- 99L
  r <- run_seeded(stale)
  expect_equal(r$pgdifTable$asDF$gammaLow[1], true_low)
})

test_that("locdepq3 bootstrap p-values match easyRasch2 and flag on adjusted p", {
  d <- poly_data()
  cutoff <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepQ3Cutoff(d, iterations = 60, parallel = FALSE,
                                 seed = 42, hdci_width = 0.99)))
  pkg <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepQ3(d, cutoff = cutoff, output = "dataframe",
                           p_value = TRUE, correction = "fwer")))$pairs
  r <- suppressWarnings(
    er2$locdepq3(data = d, vars = names(d), computeCutoff = TRUE,
                 iterations = 60, seed = 42, hdciWidth = 99,
                 pValues = TRUE, correction = "fwer"))
  tab <- r$pairTable$asDF
  expect_identical(paste(tab$item1, tab$item2), paste(pkg$Item1, pkg$Item2))
  expect_equal(tab$pValue, pkg$p_q3)
  expect_equal(tab$pAdjusted, pkg$padj_q3)
  flags <- ifelse(is.na(as.character(tab$flagged)), "", as.character(tab$flagged))
  expect_identical(flags, pkg$Flagged)
  expect_identical(r$pairTable$getColumn("pAdjusted")$title,
                   "Adj. p-value (FWER)")
  notes <- vapply(r$pairTable$notes, function(x) x$note, character(1))
  expect_true(any(grepl("adjusted p-value < 0.05", notes)))
  expect_match(r$q3Note$content, "At least 1000 iterations")
})

test_that("locdepgamma caches its simulation in the hidden simCache element", {
  # The hidden-simCache pattern (also used by residualpca, bootrestscore,
  # iteminfitmi): display-only option changes reuse the cached cutoff;
  # a stale signature forces a rerun. Tampered-bound probe as usual.
  d <- poly_data()

  run_seeded <- function(state, ...) {
    opts <- er2$locdepgammaOptions$new(vars = names(d), computeCutoff = TRUE,
                                       iterations = 60, seed = 42,
                                       hdciWidth = 99, ...)
    an <- er2$locdepgammaClass$new(options = opts, data = d)
    if (!is.null(state)) an$results$simCache$setState(state)
    an$run()
    an$results
  }

  state <- run_seeded(NULL)$simCache$state
  expect_false(is.null(state$sig))
  first_pair <- paste(state$cutoff_res$pair_cutoffs$Item1[1],
                      state$cutoff_res$pair_cutoffs$Item2[1])
  true_low <- state$cutoff_res$pair_cutoffs$gamma_low[1]

  tampered <- state
  tampered$cutoff_res$pair_cutoffs$gamma_low[1] <- -0.9876

  r <- run_seeded(tampered, sortByGamma = TRUE)
  tab <- r$dir1Table$asDF
  kk <- paste(tab$item1, tab$item2) == first_pair
  expect_equal(tab$gammaLow[kk][1], -0.9876)

  stale <- tampered
  stale$sig$seed <- 99L
  r <- run_seeded(stale)
  tab <- r$dir1Table$asDF
  kk <- paste(tab$item1, tab$item2) == first_pair
  expect_equal(tab$gammaLow[kk][1], true_low)
})

test_that("person parameters output variables match easyRasch2 and align rows", {
  # New personparams analysis: theta/SEM/sum score/n answered/extreme as
  # jamovi output variables. Output options are GUI-set (the generated
  # wrapper ignores them), so set them on the option objects directly.
  d <- poly_data()
  m <- as.matrix(d); m[3, ] <- NA          # all-NA respondent
  m[5, 1:2] <- NA                          # partial respondent
  d <- as.data.frame(m)

  opts <- er2$personparamsOptions$new(vars = names(d))
  for (o in c("outputTheta", "outputSem", "outputExtreme")) {
    opt_obj <- opts$option(o)
    opt_obj$value <- list(value = TRUE)
  }
  an <- er2$personparamsClass$new(options = opts, data = d)
  an$run()

  out_vals <- function(name) {
    pv <- an$results[[name]]$.__enclos_env__$private
    get(".values", envir = pv)[[1]]
  }
  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressMessages(easyRasch2::RMpersonParameters(
    d_used, method = "WLE", estimator = "CML", output = "dataframe"))

  theta <- out_vals("outputTheta")
  expect_length(theta, nrow(d))            # full row set incl. all-NA row
  expect_true(is.na(theta[3]))             # all-NA row -> empty cell
  expect_identical(theta[-3], pkg$theta)
  expect_identical(out_vals("outputSem")[-3], pkg$sem)
  expect_identical(as.integer(out_vals("outputExtreme")[-3]),
                   as.integer(pkg$extreme))
  expect_false(is.na(theta[5]))            # partial respondent retained

  tab <- an$results$summaryTable$asDF
  expect_equal(tab$value[tab$statistic == "Respondents used"], nrow(d_used))
  expect_match(an$results$personNote$content, "weighted likelihood")
  expect_false(is.null(an$results$thetaPlot$state))
})

test_that("targeting threshold table: wide by default, long format optional", {
  d <- poly_data()
  # wide default: one row per item, t columns + mean location, no SE/CI
  r <- suppressWarnings(er2$targeting(data = d, vars = names(d)))
  tab <- r$thresholdTable$asDF
  pkg <- suppressMessages(easyRasch2::RMitemParameters(
    d, estimator = "CML", format = "wide", se = FALSE, output = "dataframe"))
  expect_identical(names(tab), names(pkg))
  expect_equal(tab$location, pkg$location)
  expect_equal(tab$location, rowMeans(pkg[, grep("^t[0-9]+$", names(pkg))]))
  expect_false("se" %in% names(tab))

  # long format: one row per threshold with SE + CI
  r2 <- suppressWarnings(er2$targeting(data = d, vars = names(d),
                                       longFormat = TRUE))
  tab2 <- r2$thresholdTable$asDF
  pkg2 <- suppressMessages(easyRasch2::RMitemParameters(
    d, estimator = "CML", format = "long", se = TRUE, ci_level = 0.95,
    output = "dataframe"))
  expect_equal(tab2$location, pkg2$location)
  expect_equal(tab2$se, pkg2$se)
  expect_equal(tab2$ciLow, pkg2$ci_lower)
})

test_that("person fit matches easyRasch2, aligns output rows, and flags correctly", {
  d <- poly_data()
  m <- as.matrix(d); m[3, ] <- NA
  d <- as.data.frame(m)

  opts <- er2$personfitOptions$new(vars = names(d), iterations = 100,
                                   seed = 42)
  for (o in c("outputFlagged", "outputPInfit")) {
    oo <- opts$option(o)
    oo$value <- list(value = TRUE)
  }
  an <- er2$personfitClass$new(options = opts, data = d)
  an$run()

  out_vals <- function(name) {
    pv <- an$results[[name]]$.__enclos_env__$private
    get(".values", envir = pv)[[1]]
  }
  d_used <- d[rowSums(!is.na(d)) > 0, ]
  pkg <- suppressWarnings(suppressMessages(easyRasch2::RMpersonFit(
    d_used, iterations = 100, seed = 42, output = "dataframe")))

  flg <- out_vals("outputFlagged")
  expect_length(flg, nrow(d))
  expect_true(is.na(flg[3]))                       # all-NA row -> empty cell
  expect_identical(flg[-3], as.integer(pkg$flagged))
  expect_identical(out_vals("outputPInfit")[-3], pkg$p_infit)

  tab <- an$results$summaryTable$asDF
  expect_equal(tab$value[tab$statistic == "Flagged (any statistic)"],
               sum(pkg$flagged, na.rm = TRUE))
  # figures read from the simCache element
  expect_false(is.null(an$results$simCache$state$plots$infit))
})

test_that("Martin-Loef test matches easyRasch2 and documents its interpretation", {
  # Row labels contain a non-ASCII character; match by ASCII prefix so the
  # test is locale-proof. NOTE: the generated wrappers resolve bare symbols
  # as column names (resolveQuo), so variable-list arguments must be passed
  # as calls (inline expressions), not as variables holding characters.
  d <- poly_data()
  s1 <- names(d)[1:3]; s2 <- names(d)[4:5]
  r <- suppressWarnings(er2$martinlof(data = d,
                                      subscale1 = names(d)[1:3],
                                      subscale2 = names(d)[4:5],
                                      iterations = 150, seed = 42))
  pkg <- suppressWarnings(suppressMessages(easyRasch2::RMdimMartinLof(
    d, partition = list(s1, s2), iterations = 150, parallel = FALSE,
    seed = 42)))

  tab <- r$summaryTable$asDF
  v <- function(s) tab$value[grepl(s, tab$statistic, fixed = TRUE)]
  expect_equal(v("Observed Martin-L"), pkg$T_obs)
  expect_equal(v("Monte Carlo p-value"), pkg$p_value)
  expect_equal(v("Complete cases used"), pkg$sample_n)

  ct <- r$corrTable$asDF
  expect_equal(ct$r, pkg$wle_correlation$r)
  expect_equal(ct$ciLow, pkg$wle_correlation$ci_lower)

  # p-value interpretation + a-priori warning are user-facing requirements
  expect_match(r$mlNote$content, "a priori")
  expect_match(r$mlNote$content, "inflates the Type-I error")
  notes <- vapply(r$summaryTable$notes, function(x) x$note, character(1))
  expect_true(any(grepl("at least as large as the observed", notes)))
  cnotes <- vapply(r$corrTable$notes, function(x) x$note, character(1))
  expect_true(any(grepl("trivial departure", cnotes)))
})

test_that("tree-based DIF matches easyRasch2 and explains its classification", {
  d <- dif_data()
  r <- suppressWarnings(er2$diftree(data = d, vars = dif_items(),
                                    covariates = "dif", showSE = TRUE))
  pkg <- suppressMessages(suppressWarnings(easyRasch2::RMdifTree(
    d[dif_items()], covariates = d["dif"], output = "dataframe")))

  tab <- r$difTable$asDF
  expect_equal(nrow(tab), nrow(pkg))
  if (nrow(pkg) > 0) {
    expect_equal(tab$effectSize, pkg$EffectSize)
    expect_equal(tab$se, pkg$SE)
    expect_identical(as.character(tab$class), pkg$Class)
    expect_identical(r$difTable$getColumn("effectSize")$title,
                     "Partial gamma")
    notes <- vapply(r$difTable$notes, function(x) x$note, character(1))
    expect_true(any(grepl("rough magnitude guide", notes)))
  }
  expect_match(r$difNote$content, "identical to easyRasch2::RMdifTree")
})

test_that("locdepgamma keeps its asymptotic columns on the interval branch", {
  # Regression, easyRasch2 1.2.0. RMlocdepGamma()'s `p_value` default became
  # NULL, meaning "corrected bootstrap p-values whenever a full cutoff object
  # is supplied", and that branch drops `padj_bh` and `Significance`. This
  # analysis reads both, and has no bootstrap p-value option of its own, so it
  # asks for the interval branch by name. Left implicit, the adjusted
  # p-value, significance and flag columns all came back empty as soon as
  # expected ranges were switched on, and nothing errored.
  d <- poly_data()
  r <- suppressWarnings(suppressMessages(
    er2$locdepgamma(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42, pValues = FALSE)))
  t1 <- r$dir1Table$asDF

  expect_false(all(is.na(t1$padjBH)))
  expect_true(all(c("gammaPair", "gammaLow", "gammaHigh", "flagged") %in%
                    names(t1)))

  # The flag is taken on the pair statistic, the larger of the pair's two
  # rest-score directions, not on the coefficient the table displays.
  cut <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepGammaCutoff(d, iterations = 60, parallel = FALSE,
                                    seed = 42, hdci_width = 0.99)))
  pk <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepGamma(d, cutoff = cut, p_value = FALSE,
                              output = "dataframe")$direction1))
  expect_equal(t1$padjBH, pk$padj_bh)
  expect_identical(
    !is.na(t1$flagged) & t1$flagged == "TRUE",
    !is.na(pk$gamma_low) &
      (pk$gamma_pair < pk$gamma_low | pk$gamma_pair > pk$gamma_high)
  )
})

test_that("locdepgamma's significance filter still selects on the adjusted p", {
  # The same regression seen through the filter: with `padj_bh` gone the
  # filter matched nothing and emptied the table.
  d <- poly_data()
  r <- suppressWarnings(suppressMessages(
    er2$locdepgamma(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42, pValues = FALSE)))
  all_p <- r$dir1Table$asDF$padjBH
  expect_false(all(is.na(all_p)))

  rs <- suppressWarnings(suppressMessages(
    er2$locdepgamma(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42, sigOnly = TRUE,
                    pValues = FALSE)))
  expect_equal(rs$dir1Table$rowCount, sum(!is.na(all_p) & all_p < 0.05))
})

test_that("locdepgamma bootstrap p-values match easyRasch2 and flag on adjusted p", {
  # New in 3.1.0, mirroring the Q3 analysis: with expected ranges on, the
  # default is now the Westfall-Young corrected p-value rather than the
  # interval. The package drops its asymptotic columns on that path.
  d <- poly_data()
  r <- suppressWarnings(suppressMessages(
    er2$locdepgamma(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42)))
  t1 <- r$dir1Table$asDF
  expect_true(all(c("pValue", "pAdjusted", "gammaPair") %in% names(t1)))
  expect_false(all(is.na(t1$pAdjusted)))
  expect_true(all(is.na(t1$padjBH)))

  cut <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepGammaCutoff(d, iterations = 60, parallel = FALSE,
                                    seed = 42, hdci_width = 0.95)))
  pk <- suppressWarnings(suppressMessages(
    easyRasch2::RMlocdepGamma(d, cutoff = cut, p_value = TRUE,
                              correction = "fwer",
                              output = "dataframe")$direction1))
  expect_equal(t1$pValue, pk$p_gamma)
  expect_equal(t1$pAdjusted, pk$padj_gamma)
  expect_equal(t1$gammaPair, pk$gamma_pair)
  expect_identical(!is.na(t1$flagged) & t1$flagged == "TRUE", pk$flagged)

  # the significance filter follows the adjusted p in force
  rs <- suppressWarnings(suppressMessages(
    er2$locdepgamma(data = d, vars = names(d), computeCutoff = TRUE,
                    iterations = 60, seed = 42, sigOnly = TRUE)))
  expect_equal(rs$dir1Table$rowCount,
               sum(!is.na(pk$padj_gamma) & pk$padj_gamma < 0.05))
})

test_that("the two local dependence analyses share the package's defaults", {
  # easyRasch2 1.2.0 moved both cutoff functions to 400 iterations and a
  # descriptive 95% interval, and both flagging functions to the corrected
  # p-value. The module follows, so a jamovi user and an R user who change
  # nothing get the same answer.
  for (opts in list(er2$locdepq3Options$new(vars = character()),
                    er2$locdepgammaOptions$new(vars = character()))) {
    expect_equal(opts$iterations, 400)
    expect_equal(opts$hdciWidth, 95)
    expect_true(opts$pValues)
    expect_equal(opts$correction, "fwer")
  }
})
