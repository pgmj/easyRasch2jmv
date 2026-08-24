# easyRasch2jmv

<!-- badges: start -->
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<a href="https://doi.org/10.5281/zenodo.20136864" target="_blank"><img src="https://zenodo.org/badge/1199494308.svg" alt="DOI"></a>
<!-- badges: end -->

A [jamovi](https://www.jamovi.org/) module for Rasch Measurement Theory
analysis. All analyses delegate their computations to the
[easyRasch2](https://github.com/pgmj/easyRasch2) R package, so results
are numerically identical to easyRasch2 with the same seeds and
iterations.

A distinguishing feature, inherited from easyRasch2, is the use of
parametric-bootstrap critical values in place of rule-of-thumb cutoffs. In
the Conditional Item Infit and the two local dependence analyses, items and
item pairs are flagged on a multiplicity-corrected bootstrap *p*-value
(Westfall-Young family-wise, or FDR) once simulation-based cutoffs are
enabled, which controls the error rate across the whole set of tests rather
than leaving it to an interval width (Johansson, 2025, 2026).

## Analyses

### Item fit

- **Conditional Item Infit** — conditional infit MSQ via the
  `easyRasch2` package (`iarm::out_infit()` on a CML fit with WLE person
  estimates; Müller, 2020), with optional simulation-based cutoffs
  (Johansson, 2025) and a dot-plot of observed vs simulated
  distributions. With cutoffs enabled, items flag on the
  multiplicity-corrected bootstrap *p*-value by default, with the
  expected range shown alongside as a reference (Johansson, 2026).
  Numerically identical to `easyRasch2::RMitemInfit()` /
  `RMitemInfitCutoff()` with the same seed and iterations.
- **Conditional Item Infit (Multiple Imputation)** — pooled infit MSQ
  via Rubin's rules across `m` `mice` imputations (with optional
  auxiliary variables in the imputation model), pooled by
  `easyRasch2::RMitemInfitMI()`.
- **Bootstrap Item-Restscore** — non-parametric bootstrap of
  `iarm::item_restscore()` for use with large samples where the
  asymptotic test over-rejects. The bootstrap runs via
  `easyRasch2::RMitemRestscoreBoot()` (numerically identical draws and
  classifications with the same seed).
- **Item-Restscore Correlations** — observed vs model-expected
  item-restscore correlations via Goodman-Kruskal's gamma (Kreiner, 2011).
  Supports dichotomous and polytomous (Partial Credit Model) data.
  Delegates to `easyRasch2::RMitemRestscore()` (CML item locations via
  `psychotools`, WLE person-mean reference).

### Local dependence

- **Q3 Residual Correlation Matrix** — Yen's Q3 (Yen, 1984) via the
  `easyRasch2` package (CML item estimation with WLE person estimates),
  with optional simulation-based cutoffs (Christensen et al., 2017;
  [Johansson, 2024](https://pgmj.github.io/simcutoffs.html)), a Q3
  heatmap, and a per-pair distribution plot. With cutoffs enabled, item
  pairs flag on the multiplicity-corrected bootstrap *p*-value by
  default (Johansson, 2026). Numerically identical to
  `easyRasch2::RMlocdepQ3()` / `RMlocdepQ3Cutoff()` / `RMlocdepQ3Plot()`
  with the same seed and iterations.
- **Partial Gamma Local Dependence** — partial-gamma coefficients per
  item pair in both rest-score directions, with optional
  simulation-based per-pair expected ranges and a per-pair distribution
  plot; via `easyRasch2::RMlocdepGamma()` / `RMlocdepGammaCutoff()` /
  `RMlocdepGammaPlot()`. With cutoffs enabled, each pair is tested once
  on the larger of its two directions (the *Gamma pair (max)* column)
  and flagged on the multiplicity-corrected bootstrap *p*-value by
  default (Johansson, 2026).

### Dimensionality / unidimensionality

- **PCA of Standardized Residuals** — eigenvalues from PCA on the
  CML/WLE standardized Rasch residuals, with optional simulation-based
  cutoff for the first contrast (Chou & Wang, 2010) and a
  PC1-loading-vs-item-location plot; computed by
  `easyRasch2::RMdimResidualPCA()` / `RMdimResidualPCACutoff()`.
- **Dynamic/adaptive CFA Fit-Index Cutoffs** — observed one-factor
  categorical-CFA fit indices (CFI, RMSEA, SRMR) compared to a
  parametric-bootstrap null distribution simulated under the fitted
  PCM/RM, via `lavaan` WLSMV / ULSMV (avoids rule-of-thumb cutoffs),
  plus each item's standardized loading against its simulated expected
  range; computed by `easyRasch2::RMdimCFACutoff()` / `RMdimCFA()` with
  figures from `RMdimCFAPlot()`.
- **Martin-Löf Test** — likelihood-ratio test of unidimensionality
  against an a priori two-subscale partition of the items, generalised
  to polytomous models (Christensen, Bjorner, Kreiner & Petersen, 2002),
  with a Monte Carlo p-value under the unidimensional null (Christensen &
  Kreiner, 2007) and optional Besag-Clifford sequential stopping; via
  `easyRasch2::RMdimMartinLof()`.

### Differential item functioning

- **Andersen LR-test DIF** — Andersen's likelihood-ratio test via
  `easyRasch2::RMdifLR()` (which remains based on `eRm::LRtest()`), with
  per-group item or threshold locations, MaxDiff flagging, a faceted
  figure of group-by-item locations, and a response-distribution
  tileplot.
- **Partial Gamma DIF** — partial-gamma coefficients for categorical
  DIF variables, with optional simulation-based expected ranges and a
  response-distribution tileplot; via `easyRasch2::RMdifGamma()` /
  `RMdifGammaCutoff()` / `RMplotTile()`.
- **Tree-Based DIF** — model-based recursive partitioning for DIF via
  `easyRasch2::RMdifTree()` (psychotree Rasch/PCM trees; Strobl et al.,
  2015), splitting the sample wherever item parameters are unstable
  along one or more covariates, so groups need not be pre-specified.
  Each split's items get an A/B/C effect size (Mantel-Haenszel ETS Delta
  for dichotomous data, partial gamma for polytomous; Henninger et al.,
  2023, 2025), with optional purification and pruning.

### Reliability and targeting

- **Reliability** — Cronbach's α, the WLE-based PSI, marginal
  reliability (Green, 1984), and RMU from plausible values (Bignardi,
  Kievit & Bürkner, 2025), computed by `easyRasch2::RMreliability()`;
  optional bootstrap HDCIs for α, PSI, and Marginal.
- **Targeting Plot** — Wright-map style person-item targeting with
  back-to-back histograms of person and item threshold locations, plus
  a threshold-location dot plot, drawn by `easyRasch2::RMtargeting()`
  (CML item parameters via `psychotools`, WLE person locations; MML
  fallback under sparse categories). Threshold table with Wald CIs from
  `easyRasch2::RMitemParameters()`.

### Person statistics

- **Person Fit** — per-respondent conditional infit and outfit MSQ and
  the standardized log-likelihood lz, each with Monte-Carlo resampled
  p-values and a person-fit map; via `easyRasch2::RMpersonFit()`. The
  asymptotic person-fit nulls are unreliable, so significance comes from
  resampling under the fitted model (Müller, 2020; Sinharay, 2016).
- **Person Parameters** — per-respondent latent locations (theta,
  logits) with standard error of measurement, by WLE (Warm, 1989; the
  default) or EAP, via `easyRasch2::RMpersonParameters()`. Theta, SEM,
  sum score, items answered, and an extreme-score flag can be **saved as
  new variables** in the dataset. Includes the sum-score-to-logit lookup
  that was previously a separate analysis, identical to
  `easyRasch2::RMscoreSE()`.

### Visualization

- **Item Probability Curves** — model-implied category probability
  curves (polytomous, faceted per item) or joint item characteristic
  curves (dichotomous), computed by `easyRasch2::RMitemCatProb()`.
- **Item Characteristic Curves (CICC)** — expected item score curves
  with observed class-interval averages overlaid (Kreiner-style
  graphical item fit), via `easyRasch2::RMitemICCPlot()`; optional DIF
  variable with per-group observed averages and partial-gamma
  annotations.

## Requirements

All R-package dependencies are bundled as Imports and are installed
alongside the jamovi module: `easyRasch2`, `eRm`, `iarm`, `lavaan`,
`mirt`, `psychotools`, `mice`, `ggplot2`, `ggdist`, `patchwork`,
`scales`.

## Installation

In jamovi, click on **Modules** (far top right) and choose **"jamovi Library"**. Search for *easyRasch* and install from there.

For development version (see NEWS.md for changelog):

1. Download the latest `.jmo` file from the [Releases](https://github.com/pgmj/easyRasch2jmv/releases) page.
2. In jamovi, go to the **Modules** menu (⊞) → **Sideload**.
3. Select the downloaded `.jmo` file.

## Sample Data

Two sample datasets are bundled with the module, both from R package `eRm`:

- **pcmdat2** — Polytomous dataset (50 persons × 5 items, scored 0–3)
- **raschdat3** — Dichotomous dataset (50 persons × 8 items, scored 0–1)

These can be loaded from jamovi's **Open** → **Data Library** after installing
the module.


## How to cite

If you use easyRasch2jmv in published work, please cite it and jamovi itself. All analyses are computed by the easyRasch2 R package, so please cite that as well. If you use the Zotero reference manager, you can copy the bibtex below and use File -> Import from clipboard.

**easyRasch2jmv** (this jamovi module)

Johansson, M. (2026). easyRasch2jmv: A jamovi module based on easyRasch2
(Version 3.1.0) [Computer software].
<https://github.com/pgmj/easyRasch2jmv>

```bibtex
@Manual{easyRasch2jmv,
  title  = {{easyRasch2jmv}: A {jamovi} module based on {easyRasch2}},
  author = {Magnus Johansson},
  year   = {2026},
  note   = {jamovi module version 3.1.0},
  url    = {https://github.com/pgmj/easyRasch2jmv},
}
```

**jamovi** (intentionally with lower case 'j')

The jamovi project. (2026). jamovi [Computer Software] (Version 2.7) [Computer software]. <https://www.jamovi.org>

```bibtex
@Manual{thejamoviprojectJamoviComputerSoftware2026,
	title  = {jamovi [{Computer} {Software}]},
	author = {The jamovi project},
	year   = {2026},
	note   = {jamovi version 2.7},
	url    = {https://www.jamovi.org},
}
```

**easyRasch2** (the underlying R package)

Johansson, M. (2026). easyRasch2: Psychometric Analysis with Rasch Measurement
Theory (Version 1.2.0) [R].
<https://doi.org/10.32614/CRAN.package.easyRasch2>

```bibtex
@Manual{easyRasch2,
  title  = {{easyRasch2}: Psychometric Analysis with {Rasch} Measurement Theory},
  author = {Magnus Johansson},
  year   = {2026},
  note   = {R package version 1.2.0},
  doi    = {10.32614/CRAN.package.easyRasch2},
  url    = {https://doi.org/10.32614/CRAN.package.easyRasch2},
}
```


## References

- Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
  estimating reliability using Bayesian measurement uncertainty. PsyArXiv.
  <https://doi.org/10.31234/osf.io/h54k8>
- Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item-response
  models with principal component analysis on standardized residuals.
  *Educational and Psychological Measurement, 70*(5), 717–731.
  <https://doi.org/10.1177/0013164410379322>
- Christensen, K. B., Bjorner, J. B., Kreiner, S., & Petersen, J. H. (2002).
  Testing unidimensionality in polytomous Rasch models. *Psychometrika, 67*(4),
  563–574. <https://doi.org/10.1007/BF02295132>
- Christensen, K. B., & Kreiner, S. (2007). A Monte Carlo approach to
  unidimensionality testing in polytomous Rasch models. *Applied Psychological
  Measurement, 31*(1), 20–30. <https://doi.org/10.1177/0146621605286204>
- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for
  Yen's Q3: Identification of local dependence in the Rasch model using residual
  correlations. *Applied Psychological Measurement, 41*(3), 178–194.
  <https://doi.org/10.1177/0146621616677520>
- Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number of
  hypotheses simultaneously. *Journal of Statistical Theory and Practice,
  19*(6). <https://doi.org/10.1007/s42519-024-00412-4>
- Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D.
  (1984). Technical guidelines for assessing computerized adaptive tests.
  *Journal of Educational Measurement, 21*(4), 347–360.
  <https://doi.org/10.1111/j.1745-3984.1984.tb01039.x>
- Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping criterion for
  Rasch trees based on the Mantel-Haenszel effect size measure for differential
  item functioning. *Educational and Psychological Measurement, 83*(1), 181–212.
  <https://doi.org/10.1177/00131644221077135>
- Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial credit
  trees meet the partial gamma coefficient for quantifying DIF and DSF in
  polytomous items. *Behaviormetrika, 52*, 221–257.
  <https://doi.org/10.1007/s41237-024-00243-4>
- Johansson, M. (2024). Simulation-based cutoff values for Rasch item fit and
  residual correlations. <https://pgmj.github.io/simcutoffs.html>
- Johansson, M. (2025). Detecting item misfit in Rasch models. *Educational
  Methods & Psychometrics, 3*(18). <https://doi.org/10.61186/emp.2025.5>
- Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
  Rasch models: Iterations, multiplicity correction, and decision stability.
  *PsyArXiv*. <https://doi.org/10.31234/osf.io/7pqz4_v2>
- Kreiner, S. (2011). A note on item-restscore association in Rasch models.
  *Applied Psychological Measurement, 35*(7), 557–561.
  <https://doi.org/10.1177/0146621611410227>
- Mair, P., & Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for
  the application of IRT models in R. *Journal of Statistical Software, 20*(9).
  <https://doi.org/10.18637/jss.v020.i09>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust them?
  *Journal of Statistical Distributions and Applications, 7*(1), 5.
  <https://doi.org/10.1186/s40488-020-00108-7>
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
  *Journal of Statistical Software, 48*(2), 1–36.
  <https://doi.org/10.18637/jss.v048.i02>
- Sinharay, S. (2016). Assessment of person fit using resampling-based
  approaches. *Journal of Educational Measurement, 53*(1), 63–85.
  <https://doi.org/10.1111/jedm.12101>
- Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for
  detecting differential item functioning in the Rasch model. *Psychometrika,
  80*(2), 289–316. <https://doi.org/10.1007/s11336-013-9388-3>
- Warm, T. A. (1989). Weighted likelihood estimation of ability in item response
  theory. *Psychometrika, 54*(3), 427–450. <https://doi.org/10.1007/BF02294627>
- Westfall, P. H., & Young, S. S. (1993). *Resampling-based multiple testing:
  Examples and methods for p-value adjustment*. Wiley.
- Yen, W. M. (1984). Effects of local item dependence on the fit and equating
  performance of the three-parameter logistic model. *Applied Psychological
  Measurement, 8*(2), 125–145. <https://doi.org/10.1177/014662168400800201>


## Credits

This is largely based on my `easyRasch` package, and I am using Claude Opus 
to "transfer" functions to a more properly formatted package - `easyRasch2` - 
which is the foundation for all code in the Jamovi module. 

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed 
psychologist with a PhD in behavior analysis. He works as a research specialist 
at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), 
Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
